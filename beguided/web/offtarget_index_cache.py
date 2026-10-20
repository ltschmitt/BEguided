"""Process-wide cache for pre-loaded BEguidedOFF FM-index handles.

The web server opens every configured fast-mode index (``.fmi``) once — at
startup when ``cfg.preload_offtarget_index`` is set, otherwise lazily on the
first screening job — and keeps the handle here so screening jobs search
through a warm handle instead of re-opening the multi-GB index per job.

Concurrency model:
    * Handles are read-only (``locate`` / ``extract``); concurrent finders
      may share one object, mirroring how the genome/GTF caches are shared
      between worker threads.
    * The index is mmap'd: the OS pages it in on demand, so the cost shows
      as page cache (~1.7 GB per genome for human) rather than a hard RSS
      commitment, and a warm handle removes the ~1 s open cost per job.
    * Failures are remembered per path so a broken index surfaces as a clear
      per-genome "unavailable" reason instead of failing every job.
"""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Optional

logger = logging.getLogger("beguided.web.offtarget_index_cache")

_lock = threading.Lock()
_handles: "dict[str, object]" = {}
_failed: "dict[str, str]" = {}


def _reset_for_test() -> None:
    """Reset module state. Only intended for unit tests."""
    with _lock:
        _handles.clear()
        _failed.clear()


def load_index(path: str) -> Optional[object]:
    """Load (once) and return the FM-index handle for *path*.

    Returns ``None`` when the index cannot be opened (the error is logged and
    remembered). A path already loaded — by the startup preload or a previous
    job — returns the existing handle immediately.
    """
    if not path:
        return None
    path_str = str(Path(path))
    with _lock:
        cached = _handles.get(path_str)
        if cached is not None:
            return cached
    try:
        from beguidedoff.fast.fmindex import SdslGenomeIndex

        handle = SdslGenomeIndex.load(path_str)
    except Exception as exc:  # pragma: no cover - exercised via tests with a bad path
        with _lock:
            _failed[path_str] = str(exc)
        logger.warning("Failed to load FM-index %s: %s", path_str, exc)
        return None
    with _lock:
        # Keep the first successful handle per path; a concurrent duplicate
        # load is simply dropped.
        return _handles.setdefault(path_str, handle)


def get_handle(path: str) -> Optional[object]:
    """Return the warm handle for *path*, loading it lazily if necessary."""
    return load_index(path)


def preload(cfg) -> "dict[str, str]":
    """Preload the configured fast-mode indexes (one per genome key).

    Called once at app startup when ``cfg.preload_offtarget_index`` is set.
    Returns a mapping of genome key to status (``"ready"`` or
    ``"error: <detail>"``); genomes without a usable fast-mode index are
    omitted. Never raises — a broken index only means that genome's
    screening button stays disabled with the load error as its reason.
    """
    statuses: "dict[str, str]" = {}
    for key, info in cfg.offtarget_availability().items():
        if info["mode"] != "fast" or not info["index"]:
            continue
        handle = load_index(info["index"])
        if handle is not None:
            statuses[key] = "ready"
        else:
            statuses[key] = f"error: {_failed.get(str(Path(info['index'])), 'load failed')}"
    if statuses:
        logger.info(
            "Off-target index preload: %s",
            ", ".join(f"{k}={v}" for k, v in sorted(statuses.items())),
        )
    return statuses
