"""Custom logging handler that writes to an asyncio queue for SSE streaming."""

import asyncio
import logging


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
        """Put formatted log record into every registered job queue (thread-safe)."""
        msg = self._formatter.format(record)
        loop = self._loop
        if loop is None or loop.is_closed():
            return
        for q in list(self._queues.values()):
            try:
                loop.call_soon_threadsafe(q.put_nowait, msg)
            except Exception:
                pass

    def register_job(self, job_id):
        """Create or retrieve a queue for a new job."""
        if job_id in self._queues:
            return self._queues[job_id]
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
        self._loop = loop
        q = asyncio.Queue(maxsize=1000)
        self._queues[job_id] = q
        return q

    def get_queue(self, job_id):
        """Retrieve an existing queue for a job."""
        return self._queues.get(job_id)

    def unregister_job(self, job_id):
        """Remove queue for a completed job."""
        if job_id in self._queues:
            del self._queues[job_id]


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
