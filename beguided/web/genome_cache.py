"""Process-wide cache for the pre-loaded genome FASTA.

The web server reads the configured genome FASTA exactly once at startup (in a
background thread) and stores the resulting ``dict[str, str]`` here so every
pipeline run reuses it instead of re-parsing a multi-GB file from disk.

Concurrency model:
    * One daemon thread performs the actual ``load_genome`` call.
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
import threading
import time
from typing import Optional

from beguided.genome.io import derive_genome_label, load_genome

logger = logging.getLogger("beguided.web.genome_cache")

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
_genome: Optional[dict] = None
_preload_thread: Optional[threading.Thread] = None


def _reset_for_test() -> None:
    """Reset module state. Only intended for unit tests."""
    global _genome, _preload_thread
    with _cond:
        _state["status"] = "not_loaded"
        _state["path"] = None
        _state["label"] = ""
        _state["chromosome_count"] = 0
        _state["error"] = None
        _state["started_at"] = None
        _state["finished_at"] = None
        _genome = None
        _preload_thread = None
        _cond.notify_all()


def get_status() -> dict:
    """Return a snapshot of the preload state (safe for /genome/status)."""
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
    global _genome

    logger.info("Preloading genome from %s", path)
    try:
        genome = load_genome(path)
    except Exception as exc:
        logger.exception("Genome preload failed: %s", exc)
        with _cond:
            _state["status"] = "error"
            _state["error"] = f"{type(exc).__name__}: {exc}"
            _state["finished_at"] = time.time()
            _cond.notify_all()
        logger.error("Exiting process because genome preload failed (fail-fast).")
        os._exit(1)
        return

    chrom_count = len(genome)
    with _cond:
        _genome = genome
        _state["status"] = "ready"
        _state["chromosome_count"] = chrom_count
        # If no explicit label was supplied, derive one from the filename now
        # that we know the load succeeded.
        if not _state["label"]:
            _state["label"] = derive_genome_label(path)
        _state["finished_at"] = time.time()
        _cond.notify_all()
    logger.info(
        "Genome preload complete: %d chromosomes (label=%s)",
        chrom_count,
        _state["label"],
    )


def preload_async(path: str, label: str = "") -> None:
    """Start background preload. Idempotent.

    Args:
        path: Absolute path to the genome FASTA.
        label: Display label for the genome (e.g. ``"GRCh38.112"``).
            Empty string means "derive from FASTA filename".

    Raises:
        RuntimeError: If a preload for a different path is already in flight.
    """
    global _preload_thread
    with _cond:
        if _state["status"] in ("loading", "ready") and _state["path"] == path:
            return
        if _state["status"] in ("loading", "ready"):
            raise RuntimeError(
                f"Genome already loaded from {_state['path']!r}; "
                f"cannot preload a different path {path!r}"
            )
        _state["status"] = "loading"
        _state["path"] = path
        _state["label"] = label.strip() if label else ""
        _state["chromosome_count"] = 0
        _state["error"] = None
        _state["started_at"] = time.time()
        _state["finished_at"] = None
        _genome = None

    _preload_thread = threading.Thread(
        target=_load_blocking,
        args=(path,),
        name="beguided-genome-preload",
        daemon=True,
    )
    _preload_thread.start()


def wait_until_ready(timeout: Optional[float] = None) -> dict:
    """Block until the genome is loaded. Returns the cached genome dict.

    Args:
        timeout: Max seconds to wait. ``None`` waits forever.

    Returns:
        The cached ``dict[str, str]`` mapping chromosome id -> sequence.

    Raises:
        RuntimeError: If preload errored or no preload was started.
        TimeoutError: If ``timeout`` elapses before preload completes.
    """
    deadline = None if timeout is None else time.monotonic() + timeout
    with _cond:
        if _state["status"] == "ready":
            return _genome  # type: ignore[return-value]
        if _state["status"] == "error":
            raise RuntimeError(
                f"Genome preload failed: {_state['error']}. "
                "Server cannot run the guidefinder step."
            )
        if _state["status"] != "loading":
            raise RuntimeError(
                "Genome preload has not been started; call preload_async() first."
            )
        while _state["status"] == "loading":
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"Genome preload did not complete within {timeout:.1f}s."
                    )
                _cond.wait(timeout=remaining)
            else:
                _cond.wait()
        if _state["status"] == "ready":
            return _genome  # type: ignore[return-value]
        if _state["status"] == "error":
            raise RuntimeError(
                f"Genome preload failed: {_state['error']}. "
                "Server cannot run the guidefinder step."
            )
        raise RuntimeError(f"Unexpected preload state: {_state['status']!r}")
