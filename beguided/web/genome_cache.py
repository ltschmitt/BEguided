"""Process-wide cache for the pre-loaded genome FASTAs.

The web server reads every configured genome FASTA once at startup (one
background thread per genome) and stores the resulting ``dict[str, str]``
here so every pipeline run reuses it instead of re-parsing multi-GB files
from disk. Entries are keyed by the genome key from
:class:`beguided.web.config.GenomeOption` (``"human"``, ``"mouse"``, ...).

Concurrency model:
    * One daemon thread per key performs the actual ``load_genome`` call.
    * A per-key ``threading.Condition`` lets worker threads block until that
      genome's preload completes (or surfaces an error).
    * If the preload of the *default* genome fails the process exits
      (``os._exit(1)``) so the server never advertises a half-functional
      pipeline. Non-default genomes fail soft: the key is marked ``error``
      and a job that selects it fails with a clear message.

State machine (per key):
    not_loaded -> loading -> ready
                          \\-> error
"""

from __future__ import annotations

import logging
import os
import threading
import time
from dataclasses import dataclass, field
from typing import Optional

from beguided.genome.io import derive_genome_label, load_genome

logger = logging.getLogger("beguided.web.genome_cache")

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
    genome: Optional[dict] = None
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
    }


def aggregate_status(snapshots: "list[dict]") -> str:
    """Collapse per-key statuses into one of
    ``not_loaded | loading | ready | partial | error``.

    ``partial`` means at least one key is ready and at least one failed.
    """
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
    overall ``status`` plus a ``genomes`` list of per-key snapshots.
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
        "status": aggregate_status(snapshots),
        "label": "",
        "path": None,
        "chromosome_count": sum(s["chromosome_count"] for s in snapshots),
        "error": next((s["error"] for s in snapshots if s["error"]), None),
        "started_at": min(started) if started else None,
        "finished_at": max(finished) if finished else None,
        "genomes": snapshots,
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
    logger.info("Preloading genome %r from %s", entry.key, entry.path)
    try:
        genome = load_genome(entry.path)
    except Exception as exc:
        logger.exception("Genome preload failed for %r: %s", entry.key, exc)
        with entry.cond:
            entry.status = "error"
            entry.error = f"{type(exc).__name__}: {exc}"
            entry.finished_at = time.time()
            entry.cond.notify_all()
        if fail_fast:
            logger.error(
                "Exiting process because genome preload failed (fail-fast)."
            )
            os._exit(1)
        return

    chrom_count = len(genome)
    with entry.cond:
        entry.genome = genome
        entry.status = "ready"
        entry.chromosome_count = chrom_count
        # If no explicit label was supplied, derive one from the filename now
        # that we know the load succeeded.
        if not entry.label:
            entry.label = derive_genome_label(entry.path)
        entry.finished_at = time.time()
        entry.cond.notify_all()
    logger.info(
        "Genome preload complete for %r: %d chromosomes (label=%s)",
        entry.key,
        chrom_count,
        entry.label,
    )


def preload_async(
    key: str,
    path: str,
    label: str = "",
    *,
    fail_fast: bool = False,
) -> None:
    """Start background preload for ``key``. Idempotent for the same path.

    Args:
        key: Genome key (e.g. ``"human"``).
        path: Absolute path to the genome FASTA.
        label: Display label for the genome (e.g. ``"GRCh38.112"``).
            Empty string means "derive from FASTA filename".
        fail_fast: When True, exit the process if the load fails. Used for
            the default genome so a broken deployment never serves jobs.

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
                f"Genome {key!r} already loaded from {entry.path!r}; "
                f"cannot preload a different path {path!r}"
            )
        entry.status = "loading"
        entry.path = path
        entry.label = label.strip() if label else ""
        entry.chromosome_count = 0
        entry.error = None
        entry.started_at = time.time()
        entry.finished_at = None
        entry.genome = None

    entry.thread = threading.Thread(
        target=_load_blocking,
        args=(entry, fail_fast),
        name=f"beguided-genome-preload-{key}",
        daemon=True,
    )
    entry.thread.start()


def wait_until_ready(key: Optional[str] = None, timeout: Optional[float] = None) -> dict:
    """Block until the genome for ``key`` is loaded. Returns its dict.

    When ``key`` is ``None`` and exactly one genome is registered, that
    genome is used (backward-compatible single-genome call).

    Args:
        key: Genome key, or ``None`` for the only registered genome.
        timeout: Max seconds to wait. ``None`` waits forever.

    Returns:
        The cached ``dict[str, str]`` mapping chromosome id -> sequence.

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
                    "Genome preload has not been started; call preload_async() first."
                )
            raise RuntimeError(
                "Multiple genomes are registered; wait_until_ready() requires a key."
            )
        raise RuntimeError(
            f"Genome {key!r} has not been preloaded; call preload_async() first."
        )

    deadline = None if timeout is None else time.monotonic() + timeout
    with entry.cond:
        if entry.status == "ready":
            return entry.genome  # type: ignore[return-value]
        if entry.status == "error":
            raise RuntimeError(
                f"Genome preload failed for {entry.key!r}: {entry.error}. "
                "Server cannot run the guidefinder step for this genome."
            )
        if entry.status != "loading":
            raise RuntimeError(
                "Genome preload has not been started; call preload_async() first."
            )
        while entry.status == "loading":
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"Genome preload did not complete within {timeout:.1f}s."
                    )
                entry.cond.wait(timeout=remaining)
            else:
                entry.cond.wait()
        if entry.status == "ready":
            return entry.genome  # type: ignore[return-value]
        if entry.status == "error":
            raise RuntimeError(
                f"Genome preload failed for {entry.key!r}: {entry.error}. "
                "Server cannot run the guidefinder step for this genome."
            )
        raise RuntimeError(f"Unexpected preload state: {entry.status!r}")
