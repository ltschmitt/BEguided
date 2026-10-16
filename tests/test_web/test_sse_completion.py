"""Integration tests for prompt SSE stream completion.

Regression: ``/logs/stream`` used to stay open for up to 30 s after a job
finished (it only noticed completion on the next keepalive timeout). The
worker now pushes an ``SSE_DONE`` sentinel, and the endpoint closes a stream
immediately when the job is already finished.
"""

from __future__ import annotations

import sqlite3
import threading
import time

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("uvicorn")


@pytest.fixture
def web_env(tmp_path):
    """Install a WebConfig + Scheduler backed by a temp DB."""
    from beguided.web.config import WebConfig, set_config
    from beguided.web.scheduler import Scheduler, install

    db = tmp_path / "jobs.db"
    set_config(
        WebConfig(
            genome_fasta="/nonexistent",
            gtf_path="/nonexistent",
            beguidedoff_genome="",
            db_path=str(db),
            queue_max=3,
            cooldown_seconds=0,
            max_upload_mb=5,
            usage_stats_enabled=False,
        )
    )
    sch = Scheduler(db_path=str(db), queue_max=3, cooldown_seconds=0)
    install(sch)
    sch.start()
    try:
        yield sch
    finally:
        sch.stop()


def _insert_job(sch, job_id: str, status: str) -> None:
    now = time.time()
    with sqlite3.connect(sch.db_path) as c:
        c.execute(
            "INSERT INTO jobs(job_id, revoke_token, status, payload_json, "
            "created_at, started_at, finished_at) VALUES(?,?,?,?,?,?,?)",
            (
                job_id,
                "tok",
                status,
                "{}",
                now,
                now,
                now if status == "complete" else None,
            ),
        )


def _client():
    from fastapi.testclient import TestClient

    from beguided.web.app import create_app

    return TestClient(create_app())


def test_logs_stream_closes_immediately_when_job_finished(web_env):
    """A late connect to a finished job must not wait for the keepalive."""
    _insert_job(web_env, "donejob", "complete")
    client = _client()

    t0 = time.time()
    resp = client.get("/logs/stream", params={"job_id": "donejob"})
    elapsed = time.time() - t0

    assert resp.status_code == 200
    assert resp.text == ""
    assert elapsed < 10, f"stream lingered {elapsed:.1f}s (keepalive is 30s)"


def test_logs_stream_closes_on_done_sentinel(web_env):
    """Once the worker signals completion, the stream ends promptly."""
    from beguided.web.sse_logger import get_sse_handler

    _insert_job(web_env, "runningjob", "running")

    handler = get_sse_handler()

    def _finish():
        # Give the endpoint time to register its queue, then emit a final
        # line and the completion sentinel from another thread.
        time.sleep(0.5)
        handler.push("runningjob", "final line")
        handler.close_job("runningjob")

    finisher = threading.Thread(target=_finish, daemon=True)
    finisher.start()

    client = _client()
    t0 = time.time()
    resp = client.get("/logs/stream", params={"job_id": "runningjob"})
    elapsed = time.time() - t0
    finisher.join(timeout=5)

    assert resp.status_code == 200
    assert "final line" in resp.text
    assert elapsed < 15, f"stream lingered {elapsed:.1f}s (keepalive is 30s)"
