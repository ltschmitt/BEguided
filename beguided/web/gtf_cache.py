"""Process-wide cache for the pre-loaded GTF annotations.

The web server parses every configured GTF once at startup (one background
thread per genome) and stores the resulting ``dict[str, IntervalTree]`` here
so every pipeline run reuses it instead of re-parsing multi-hundred-MB
annotation files. Entries are keyed by the genome key from
:class:`beguided.web.config.GenomeOption`.

Each entry also records the CDS ``select_tag`` used to build it
(``"MANE_Select"`` for human Ensembl, ``"Ensembl_canonical"`` for the EBI
mouse annotation).

Concurrency model:
    * One daemon thread per key performs the actual
      ``load_gtf_to_interval_tree`` call.
    * A per-key ``threading.Condition`` lets worker threads block until that
      key's preload completes (or surfaces an error).
    * If the preload of the *default* genome's GTF fails the process exits
      (``os._exit(1)``); other keys fail soft and are marked ``error``.

State machine (per key):
    not_loaded -> loading -> ready
                          \\-> error
"""

from __future__ import annotations

import logging
import os
import re
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from beguided.bystander.annotation import load_gtf_to_interval_tree

logger = logging.getLogger("beguided.web.gtf_cache")

# Captures the assembly identifier from an Ensembl-style GTF filename.
# Examples:
#   Homo_sapiens.GRCh38.112.gtf                -> "GRCh38.112"
#   Homo_sapiens.GRCh38.chr_patch_hapl_scaff.gtf -> "GRCh38"
#   Mus_musculus.GRCm39.2026_04.gtf            -> "GRCm39.2026_04"
_GTF_ASSEMBLY_RE = re.compile(
    r"^(?:Homo_sapiens|Mus_musculus)\.([A-Za-z0-9._]+?)(?:\.chr_|\.gtf$|\.gff$|\.gff3$|$)"
)


def derive_gtf_label(gtf_path: str) -> str:
    """Best-effort extraction of the assembly label from an Ensembl-style
    GTF filename.

    Examples:
        ``Homo_sapiens.GRCh38.112.gtf`` -> ``"GRCh38.112"``
        ``Homo_sapiens.GRCh38.chr_patch_hapl_scaff.gtf`` -> ``"GRCh38"``
        ``/data/anything.gtf`` -> ``"GRCh38"`` (fallback)

    The fallback matches the on-disk fixtures shipped with BEguided so
    operators see a sensible label even for non-Ensembl filenames.
    """
    name = Path(gtf_path).name
    m = _GTF_ASSEMBLY_RE.match(name)
    if m:
        return m.group(1)
    return "GRCh38"


_lock = threading.Lock()


@dataclass
class _Entry:
    key: str
    status: str = "not_loaded"        # not_loaded | loading | ready | error
    path: Optional[str] = None
    label: str = ""
    chromosome_count: int = 0
    error: Optional[str] = None
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    select_tag: str = "MANE_Select"
    interval_trees: Optional[dict] = None
    cond: threading.Condition = field(default_factory=threading.Condition)
    thread: Optional[threading.Thread] = None


_entries: "dict[str, _Entry]" = {}


def _reset_for_test() -> None:
    """Reset module state. Only intended for unit tests."""
    with _lock:
        for entry in _entries.values():
            with entry.cond:
                entry.cond.notify_all()
        _entries.clear()


def _get_or_create(key: str) -> _Entry:
    with _lock:
        entry = _entries.get(key)
        if entry is None:
            entry = _Entry(key=key)
            _entries[key] = entry
        return entry


def _resolve_entry(key: Optional[str]) -> Optional[_Entry]:
    with _lock:
        if key is not None:
            return _entries.get(key)
        if len(_entries) == 1:
            return next(iter(_entries.values()))
        return None


def _snapshot(entry: _Entry) -> dict:
    return {
        "key": entry.key,
        "status": entry.status,
        "path": entry.path,
        "label": entry.label,
        "chromosome_count": entry.chromosome_count,
        "error": entry.error,
        "started_at": entry.started_at,
        "finished_at": entry.finished_at,
        "select_tag": entry.select_tag,
    }


def _empty_snapshot(key: Optional[str] = None) -> dict:
    return {
        "key": key,
        "status": "not_loaded",
        "path": None,
        "label": "",
        "chromosome_count": 0,
        "error": None,
        "started_at": None,
        "finished_at": None,
        "select_tag": "MANE_Select",
    }


def _aggregate_status(snapshots: "list[dict]") -> str:
    statuses = [s["status"] for s in snapshots]
    if not statuses:
        return "not_loaded"
    if all(s == "ready" for s in statuses):
        return "ready"
    if any(s == "error" for s in statuses):
        return "error" if all(s == "error" for s in statuses) else "partial"
    if any(s == "loading" for s in statuses):
        return "loading"
    if all(s == "not_loaded" for s in statuses):
        return "not_loaded"
    return "partial"


def get_status(key: Optional[str] = None) -> dict:
    """Return a snapshot of the preload state.

    With ``key``: that entry's snapshot (or a ``not_loaded`` snapshot when
    the key was never registered). Without ``key``: an aggregate dict with an
    overall ``status`` plus a ``gtfs`` list of per-key snapshots.
    """
    if key is not None:
        entry = _resolve_entry(key)
        if entry is None:
            return _empty_snapshot(key)
        with entry.cond:
            return _snapshot(entry)

    with _lock:
        entries = list(_entries.values())
    snapshots = []
    for entry in entries:
        with entry.cond:
            snapshots.append(_snapshot(entry))
    started = [s["started_at"] for s in snapshots if s["started_at"]]
    finished = [s["finished_at"] for s in snapshots if s["finished_at"]]
    return {
        "status": _aggregate_status(snapshots),
        "label": "",
        "path": None,
        "chromosome_count": sum(s["chromosome_count"] for s in snapshots),
        "error": next((s["error"] for s in snapshots if s["error"]), None),
        "started_at": min(started) if started else None,
        "finished_at": max(finished) if finished else None,
        "gtfs": snapshots,
    }


def is_ready(key: Optional[str] = None) -> bool:
    entry = _resolve_entry(key)
    if entry is None:
        return False
    with entry.cond:
        return entry.status == "ready"


def available_keys() -> "list[str]":
    with _lock:
        return list(_entries.keys())


def _load_blocking(entry: _Entry, fail_fast: bool) -> None:
    """Run inside the preload thread. Sets state and (optionally) exits."""
    logger.info(
        "Preloading GTF annotation %r from %s (tag=%s)",
        entry.key,
        entry.path,
        entry.select_tag,
    )
    try:
        interval_trees = load_gtf_to_interval_tree(
            entry.path, select_tag=entry.select_tag
        )
    except Exception as exc:
        logger.exception("GTF preload failed for %r: %s", entry.key, exc)
        with entry.cond:
            entry.status = "error"
            entry.error = f"{type(exc).__name__}: {exc}"
            entry.finished_at = time.time()
            entry.cond.notify_all()
        if fail_fast:
            logger.error(
                "Exiting process because GTF preload failed (fail-fast)."
            )
            os._exit(1)
        return

    # Exclude the reserved __tx_codons__ key from the chromosome count.
    from beguided.bystander.annotation import _TX_CODONS_KEY
    chrom_count = len(interval_trees) - (1 if _TX_CODONS_KEY in interval_trees else 0)
    with entry.cond:
        entry.interval_trees = interval_trees
        entry.status = "ready"
        entry.chromosome_count = chrom_count
        # Derive an assembly label (e.g. "GRCh38.112") from the filename now
        # that the load succeeded, mirroring genome_cache behaviour.
        entry.label = derive_gtf_label(entry.path)
        entry.finished_at = time.time()
        entry.cond.notify_all()
    logger.info(
        "GTF preload complete for %r: %d chromosomes (%s CDS, label=%s)",
        entry.key,
        chrom_count,
        entry.select_tag,
        entry.label,
    )


def preload_async(
    key: str,
    path: str,
    *,
    select_tag: str = "MANE_Select",
    fail_fast: bool = False,
) -> None:
    """Start background preload for ``key``. Idempotent for the same path.

    Args:
        key: Genome key (e.g. ``"human"``).
        path: Absolute path to the GTF annotation file.
        select_tag: CDS attribute tag marking the canonical transcript.
        fail_fast: When True, exit the process if the load fails.

    Raises:
        RuntimeError: If a preload for a different path is already in flight
            for the same key.
    """
    entry = _get_or_create(key)
    with entry.cond:
        if entry.status in ("loading", "ready") and entry.path == path:
            return
        if entry.status in ("loading", "ready"):
            raise RuntimeError(
                f"GTF {key!r} already loaded from {entry.path!r}; "
                f"cannot preload a different path {path!r}"
            )
        entry.status = "loading"
        entry.path = path
        entry.label = ""
        entry.chromosome_count = 0
        entry.error = None
        entry.started_at = time.time()
        entry.finished_at = None
        entry.select_tag = select_tag
        entry.interval_trees = None

    entry.thread = threading.Thread(
        target=_load_blocking,
        args=(entry, fail_fast),
        name=f"beguided-gtf-preload-{key}",
        daemon=True,
    )
    entry.thread.start()


def wait_until_ready(key: Optional[str] = None, timeout: Optional[float] = None) -> dict:
    """Block until the GTF for ``key`` is loaded. Returns the interval trees.

    When ``key`` is ``None`` and exactly one GTF is registered, that entry is
    used (backward-compatible single-genome call).

    Args:
        key: Genome key, or ``None`` for the only registered GTF.
        timeout: Max seconds to wait. ``None`` waits forever.

    Returns:
        The cached ``dict[str, IntervalTree]`` mapping chromosome id -> tree.

    Raises:
        RuntimeError: If preload errored, was never started, or the key is
            ambiguous/unknown.
        TimeoutError: If ``timeout`` elapses before preload completes.
    """
    entry = _resolve_entry(key)
    if entry is None:
        if key is None:
            with _lock:
                known = len(_entries)
            if known == 0:
                raise RuntimeError(
                    "GTF preload has not been started; call preload_async() first."
                )
            raise RuntimeError(
                "Multiple GTFs are registered; wait_until_ready() requires a key."
            )
        raise RuntimeError(
            f"GTF {key!r} has not been preloaded; call preload_async() first."
        )

    deadline = None if timeout is None else time.monotonic() + timeout
    with entry.cond:
        if entry.status == "ready":
            return entry.interval_trees  # type: ignore[return-value]
        if entry.status == "error":
            raise RuntimeError(
                f"GTF preload failed for {entry.key!r}: {entry.error}. "
                "Server cannot run the bystander step for this genome."
            )
        if entry.status != "loading":
            raise RuntimeError(
                "GTF preload has not been started; call preload_async() first."
            )
        while entry.status == "loading":
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"GTF preload did not complete within {timeout:.1f}s."
                    )
                entry.cond.wait(timeout=remaining)
            else:
                entry.cond.wait()
        if entry.status == "ready":
            return entry.interval_trees  # type: ignore[return-value]
        if entry.status == "error":
            raise RuntimeError(
                f"GTF preload failed for {entry.key!r}: {entry.error}. "
                "Server cannot run the bystander step for this genome."
            )
        raise RuntimeError(f"Unexpected preload state: {entry.status!r}")
