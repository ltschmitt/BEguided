"""Process-wide cache for the pre-loaded GTF annotation.

The web server parses the configured GTF exactly once at startup (in a
background thread) and stores the resulting ``dict[str, IntervalTree]`` here
so every pipeline run reuses it instead of re-parsing a multi-hundred-MB
annotation file from disk.

Concurrency model:
    * One daemon thread performs the actual ``load_gtf_to_interval_tree`` call.
    * A ``threading.Condition`` lets worker threads block until preload
      completes (or surfaces a fatal error).
    * If preload fails the process exits (``os._exit(1)``) so the server never
      advertises a half-functional pipeline.

State machine:
    not_loaded -> loading -> ready
                          \\-> error
"""

from __future__ import annotations

import logging
import os
import re
import threading
import time
from pathlib import Path
from typing import Optional

from beguided.bystander.annotation import load_gtf_to_interval_tree

logger = logging.getLogger("beguided.web.gtf_cache")

# Captures the assembly identifier from an Ensembl-style GTF filename.
# Examples:
#   Homo_sapiens.GRCh38.112.gtf                -> "GRCh38.112"
#   Homo_sapiens.GRCh38.chr_patch_hapl_scaff.gtf -> "GRCh38"
_GTF_ASSEMBLY_RE = re.compile(
    r"^Homo_sapiens\.([A-Za-z0-9.]+?)(?:\.chr_|\.gtf$|\.gff$|\.gff3$|$)"
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


_state: dict = {
    "status": "not_loaded",        # not_loaded | loading | ready | error
    "path": None,
    "label": "",
    "chromosome_count": 0,
    "error": None,
    "started_at": None,
    "finished_at": None,
}
_lock = threading.Lock()
_cond = threading.Condition(_lock)
_interval_trees: Optional[dict] = None
_preload_thread: Optional[threading.Thread] = None


def _reset_for_test() -> None:
    """Reset module state. Only intended for unit tests."""
    global _interval_trees, _preload_thread
    with _cond:
        _state["status"] = "not_loaded"
        _state["path"] = None
        _state["label"] = ""
        _state["chromosome_count"] = 0
        _state["error"] = None
        _state["started_at"] = None
        _state["finished_at"] = None
        _interval_trees = None
        _preload_thread = None
        _cond.notify_all()


def get_status() -> dict:
    """Return a snapshot of the preload state (safe for /gtf/status)."""
    with _lock:
        return {
            "status": _state["status"],
            "path": _state["path"],
            "label": _state["label"],
            "chromosome_count": _state["chromosome_count"],
            "error": _state["error"],
            "started_at": _state["started_at"],
            "finished_at": _state["finished_at"],
        }


def is_ready() -> bool:
    with _lock:
        return _state["status"] == "ready"


def _load_blocking(path: str) -> None:
    """Run inside the preload thread. Sets state and exits on success/failure."""
    global _interval_trees

    logger.info("Preloading GTF annotation from %s", path)
    try:
        interval_trees = load_gtf_to_interval_tree(path)
    except Exception as exc:
        logger.exception("GTF preload failed: %s", exc)
        with _cond:
            _state["status"] = "error"
            _state["error"] = f"{type(exc).__name__}: {exc}"
            _state["finished_at"] = time.time()
            _cond.notify_all()
        logger.error("Exiting process because GTF preload failed (fail-fast).")
        os._exit(1)
        return

    # Exclude the reserved __tx_codons__ key from the chromosome count.
    from beguided.bystander.annotation import _TX_CODONS_KEY
    chrom_count = len(interval_trees) - (1 if _TX_CODONS_KEY in interval_trees else 0)
    with _cond:
        _interval_trees = interval_trees
        _state["status"] = "ready"
        _state["chromosome_count"] = chrom_count
        # Derive an assembly label (e.g. "GRCh38.112") from the filename now
        # that the load succeeded, mirroring genome_cache behaviour.
        _state["label"] = derive_gtf_label(path)
        _state["finished_at"] = time.time()
        _cond.notify_all()
    logger.info(
        "GTF preload complete: %d chromosomes with MANE_Select CDS (label=%s)",
        chrom_count,
        _state["label"],
    )


def preload_async(path: str) -> None:
    """Start background preload. Idempotent.

    Args:
        path: Absolute path to the GTF annotation file.

    Raises:
        RuntimeError: If a preload for a different path is already in flight.
    """
    global _preload_thread
    with _cond:
        if _state["status"] in ("loading", "ready") and _state["path"] == path:
            return
        if _state["status"] in ("loading", "ready"):
            raise RuntimeError(
                f"GTF already loaded from {_state['path']!r}; "
                f"cannot preload a different path {path!r}"
            )
        _state["status"] = "loading"
        _state["path"] = path
        _state["label"] = ""
        _state["chromosome_count"] = 0
        _state["error"] = None
        _state["started_at"] = time.time()
        _state["finished_at"] = None
        _interval_trees = None

    _preload_thread = threading.Thread(
        target=_load_blocking,
        args=(path,),
        name="beguided-gtf-preload",
        daemon=True,
    )
    _preload_thread.start()


def wait_until_ready(timeout: Optional[float] = None) -> dict:
    """Block until the GTF is loaded. Returns the cached interval trees.

    Args:
        timeout: Max seconds to wait. ``None`` waits forever.

    Returns:
        The cached ``dict[str, IntervalTree]`` mapping chromosome id -> tree.

    Raises:
        RuntimeError: If preload errored or no preload was started.
        TimeoutError: If ``timeout`` elapses before preload completes.
    """
    deadline = None if timeout is None else time.monotonic() + timeout
    with _cond:
        if _state["status"] == "ready":
            return _interval_trees  # type: ignore[return-value]
        if _state["status"] == "error":
            raise RuntimeError(
                f"GTF preload failed: {_state['error']}. "
                "Server cannot run the bystander step."
            )
        if _state["status"] != "loading":
            raise RuntimeError(
                "GTF preload has not been started; call preload_async() first."
            )
        while _state["status"] == "loading":
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"GTF preload did not complete within {timeout:.1f}s."
                    )
                _cond.wait(timeout=remaining)
            else:
                _cond.wait()
        if _state["status"] == "ready":
            return _interval_trees  # type: ignore[return-value]
        if _state["status"] == "error":
            raise RuntimeError(
                f"GTF preload failed: {_state['error']}. "
                "Server cannot run the bystander step."
            )
        raise RuntimeError(f"Unexpected preload state: {_state['status']!r}")
