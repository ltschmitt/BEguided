"""Custom logging handler that writes to an asyncio queue for SSE streaming."""

import asyncio
import contextvars
import logging


# The scheduler worker sets this to the running job's id for the duration of
# that job. ``AsyncQueueHandler.emit`` uses it to route each log record to
# only that job's SSE queue, so parallel workers cannot cross-talk. Records
# emitted outside any job (startup, maintenance) fall back to broadcasting.
current_job_id: contextvars.ContextVar = contextvars.ContextVar(
    "beguided_current_job_id", default=None
)


# Terminal marker pushed into a job's SSE queue when the job finishes, so the
# ``/logs/stream`` generator can close the connection immediately instead of
# waiting for the next keepalive timeout. Identity-compared, never yielded.
SSE_DONE = object()


class AsyncQueueHandler(logging.Handler):
    """Thread-safe logging handler that puts records into per-job asyncio queues."""

    def __init__(self):
        super().__init__()
        self._queues = {}
        self._loop = None
        self._formatter = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        )

    def set_loop(self, loop):
        """Store a reference to the main event loop for thread-safe callbacks."""
        self._loop = loop

    def emit(self, record):
        """Route a formatted log record to the originating job's queue.

        Thread-safe. When ``current_job_id`` is set (the normal case inside a
        job) the record goes only to that job's queue; with several worker
        threads running, each browser therefore sees only its own job's logs.
        Records emitted with no active job fall back to broadcasting to every
        registered queue (startup / maintenance logs).
        """
        msg = self._formatter.format(record)
        loop = self._loop
        if loop is None or loop.is_closed():
            return
        job_id = current_job_id.get()
        if job_id is not None:
            q = self._queues.get(job_id)
            queues = [q] if q is not None else []
        else:
            queues = list(self._queues.values())
        for q in queues:
            try:
                loop.call_soon_threadsafe(q.put_nowait, msg)
            except Exception:
                pass

    def register_job(self, job_id):
        """Create or retrieve a queue for a new job.

        Must be called from the asyncio loop thread. Returns the queue
        regardless so callers can stash it for the worker thread.
        """
        if job_id in self._queues:
            return self._queues[job_id]
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError as exc:
            raise RuntimeError(
                "AsyncQueueHandler.register_job must be called from the "
                "asyncio loop thread."
            ) from exc
        self._loop = loop
        q = asyncio.Queue(maxsize=1000)
        self._queues[job_id] = q
        return q

    def get_queue(self, job_id):
        """Retrieve an existing queue for a job (returns None if unknown)."""
        return self._queues.get(job_id)

    def push(self, job_id, msg):
        """Push a log line to a single job's queue from any thread."""
        q = self._queues.get(job_id)
        if q is None:
            return False
        loop = self._loop
        if loop is None or loop.is_closed():
            return False
        try:
            loop.call_soon_threadsafe(q.put_nowait, msg)
        except Exception:
            return False
        return True

    def unregister_job(self, job_id):
        """Remove queue for a completed job."""
        if job_id in self._queues:
            del self._queues[job_id]

    def close_job(self, job_id):
        """Signal a finished job's SSE stream to end.

        Pushes the ``SSE_DONE`` sentinel into the job's queue, if one is
        registered (registration is lazy — the browser may never have
        connected). Safe to call from the worker thread; a no-op otherwise.
        """
        q = self._queues.get(job_id)
        loop = self._loop
        if q is None or loop is None or loop.is_closed():
            return False
        try:
            loop.call_soon_threadsafe(q.put_nowait, SSE_DONE)
        except Exception:
            return False
        return True


# Global handler instance
_sse_handler = AsyncQueueHandler()


def attach_sse_handler_to_beguided_logs():
    """Attach the global SSE handler to the root beguided logger.

    Uses propagation so child loggers (beguided.pipeline, etc.)
    automatically feed through the root handler.
    """
    try:
        loop = asyncio.get_running_loop()
        _sse_handler.set_loop(loop)
    except RuntimeError:
        pass

    root = logging.getLogger("beguided")
    if _sse_handler not in root.handlers:
        root.addHandler(_sse_handler)
    # Ensure propagation is enabled on child loggers that may
    # have previously been set up with propagate=False.
    for name in list(logging.root.manager.loggerDict):
        if name.startswith("beguided."):
            logger = logging.getLogger(name)
            logger.propagate = True


def get_sse_handler():
    """Return the global AsyncQueueHandler instance."""
    return _sse_handler
