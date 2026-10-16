"""Tests for per-job SSE log routing (parallel-worker correctness).

``AsyncQueueHandler.emit`` must send each record only to the queue of the job
that produced it (tracked via the ``current_job_id`` contextvar), otherwise
parallel workers would leak log lines into each other's browsers.
"""

from __future__ import annotations

import asyncio
import logging

from beguided.web.sse_logger import AsyncQueueHandler, SSE_DONE, current_job_id


def _record(msg: str) -> logging.LogRecord:
    return logging.LogRecord(
        "beguided.test", logging.INFO, __file__, 1, msg, None, None
    )


async def _drain():
    # Let call_soon_threadsafe callbacks scheduled by emit() run.
    await asyncio.sleep(0)
    await asyncio.sleep(0)


def test_emit_routes_only_to_current_job():
    async def scenario():
        handler = AsyncQueueHandler()
        handler.set_loop(asyncio.get_running_loop())
        q1 = handler.register_job("job1")
        q2 = handler.register_job("job2")

        token = current_job_id.set("job1")
        try:
            handler.emit(_record("hello from j1"))
        finally:
            current_job_id.reset(token)

        await _drain()
        assert not q1.empty(), "originating job should receive the record"
        assert q2.empty(), "other job must not receive the record"
        assert "hello from j1" in q1.get_nowait()

    asyncio.run(scenario())


def test_emit_without_job_broadcasts():
    async def scenario():
        handler = AsyncQueueHandler()
        handler.set_loop(asyncio.get_running_loop())
        q1 = handler.register_job("job1")
        q2 = handler.register_job("job2")

        # No active job: startup/maintenance records broadcast to everyone.
        handler.emit(_record("startup"))

        await _drain()
        assert "startup" in q1.get_nowait()
        assert "startup" in q2.get_nowait()

    asyncio.run(scenario())


def test_emit_for_unknown_job_is_dropped():
    async def scenario():
        handler = AsyncQueueHandler()
        handler.set_loop(asyncio.get_running_loop())
        q1 = handler.register_job("job1")

        token = current_job_id.set("ghost")
        try:
            handler.emit(_record("orphan"))
        finally:
            current_job_id.reset(token)

        await _drain()
        assert q1.empty()

    asyncio.run(scenario())


def test_close_job_pushes_done_sentinel():
    """close_job must enqueue the terminal marker so /logs/stream can end."""

    async def scenario():
        handler = AsyncQueueHandler()
        handler.set_loop(asyncio.get_running_loop())
        q = handler.register_job("job1")

        assert handler.close_job("job1") is True
        await _drain()
        assert q.get_nowait() is SSE_DONE

    asyncio.run(scenario())


def test_close_job_unknown_job_is_noop():
    async def scenario():
        handler = AsyncQueueHandler()
        handler.set_loop(asyncio.get_running_loop())
        assert handler.close_job("ghost") is False

    asyncio.run(scenario())


def test_close_job_without_loop_is_noop():
    handler = AsyncQueueHandler()
    # No event loop attached: nothing to schedule onto.
    assert handler.close_job("job1") is False
