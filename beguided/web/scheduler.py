"""Single-worker FIFO job scheduler for the BEguided web server.

Submits a queue.Queue (maxsize=cfg.queue_max) plus one persistent
``threading.Thread`` worker. Job state is persisted to SQLite (WAL mode) so
the queue survives a restart. Each running job owns a ``threading.Event``
that the worker checks at pipeline-step boundaries (cooperative cancel).

Concurrency invariants:
    * At most one pipeline runs at any time (one worker thread).
    * At most ``queue_max`` jobs can sit in the queue; further submits raise
      ``QueueFull`` (HTTP 503).
    * Per-IP cooldown ``cooldown_seconds`` blocks back-to-back submits (HTTP 429).

API:
    start()                          - call once at app startup
    submit(ip, payload) -> Job       - raises Cooldown or QueueFull
    cancel(job_id, token) -> bool    - token-gated
    get_status(job_id) -> str
    get_result(job_id) -> dict|None
"""

from __future__ import annotations

import json
import logging
import os
import queue
import secrets
import sqlite3
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from beguided.pipeline import PipelineCancelled
from beguided.web.tasks import execute_job

logger = logging.getLogger("beguided.web.scheduler")


DEFAULT_DB_PATH = "/var/lib/beguided/jobs.db"


class SchedulerError(Exception):
    """Base for scheduler-rejection errors."""


class QueueFull(SchedulerError):
    """Raised when the in-process queue is at capacity."""


class Cooldown(SchedulerError):
    """Raised when the submitter is still inside the per-IP cooldown window."""

    def __init__(self, retry_after: float):
        super().__init__(f"cooldown, retry in {retry_after:.0f}s")
        self.retry_after = retry_after


@dataclass
class Job:
    job_id: str
    revoke_token: str
    cancel: threading.Event = field(default_factory=threading.Event)
    payload: dict = field(default_factory=dict)


class Scheduler:
    """Single-worker FIFO queue backed by SQLite."""

    STATUSES = ("queued", "running", "complete", "error", "cancelled")

    def __init__(
        self,
        db_path: str = DEFAULT_DB_PATH,
        queue_max: int = 8,
        cooldown_seconds: float = 60.0,
    ):
        self.db_path = db_path
        self.queue_max = queue_max
        self.cooldown_seconds = cooldown_seconds

        self._q: "queue.Queue[Job]" = queue.Queue(maxsize=queue_max)
        self._last_submit: dict[str, float] = {}
        self._db_lock = threading.Lock()
        self._running: dict[str, Job] = {}
        self._worker: Optional[threading.Thread] = None

    # ----- lifecycle -------------------------------------------------------

    def start(self) -> None:
        """Idempotent: start the worker thread and ensure the DB schema."""
        self._init_db()
        if self._worker is None or not self._worker.is_alive():
            self._worker = threading.Thread(
                target=self._worker_loop, name="beguided-worker", daemon=True
            )
            self._worker.start()
            logger.info(
                "Scheduler started (queue_max=%d, cooldown=%ds, db=%s)",
                self.queue_max,
                int(self.cooldown_seconds),
                self.db_path,
            )

    # ----- public API ------------------------------------------------------

    def submit(self, client_ip: str, payload: dict) -> Job:
        """Enqueue a new job. Raises ``Cooldown`` or ``QueueFull``."""
        ip = client_ip or "unknown"
        now = time.monotonic()
        last = self._last_submit.get(ip, 0.0)
        remaining = self.cooldown_seconds - (now - last)
        if remaining > 0:
            raise Cooldown(retry_after=remaining)

        job = Job(
            job_id=str(uuid.uuid4()),
            revoke_token=secrets.token_urlsafe(32),
            payload=payload,
        )
        with self._db_lock:
            self._conn().execute(
                "INSERT INTO jobs(job_id, revoke_token, status, payload_json, "
                "result_json, created_at, started_at, finished_at) "
                "VALUES (?, ?, ?, ?, NULL, ?, NULL, NULL)",
                (
                    job.job_id,
                    job.revoke_token,
                    "queued",
                    json.dumps(_payload_to_jsonable(payload)),
                    time.time(),
                ),
            )
        try:
            self._q.put_nowait(job)
        except queue.Full:
            # Roll back the DB row so the count matches the queue.
            with self._db_lock:
                self._conn().execute(
                    "DELETE FROM jobs WHERE job_id = ? AND status = 'queued'",
                    (job.job_id,),
                )
            raise QueueFull(f"queue is full ({self.queue_max} pending jobs)")
        self._last_submit[ip] = now
        logger.info("Job queued: %s (ip=%s)", job.job_id, ip)
        return job

    def cancel(self, job_id: str, token: str) -> bool:
        """Cancel a queued or running job. Returns False on token/job mismatch."""
        with self._db_lock:
            row = self._conn().execute(
                "SELECT revoke_token, status FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
            if not row or row[0] != token:
                return False
            status = row[1]
            if status == "queued":
                self._conn().execute(
                    "UPDATE jobs SET status='cancelled', finished_at=? "
                    "WHERE job_id=? AND status='queued'",
                    (time.time(), job_id),
                )
                return True
            if status == "running":
                job = self._running.get(job_id)
                if job is not None:
                    job.cancel.set()
                    return True
        return False

    def get_status(self, job_id: str) -> str:
        with self._db_lock:
            row = self._conn().execute(
                "SELECT status FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
        return row[0] if row else "unknown"

    def get_result(self, job_id: str) -> Optional[dict]:
        with self._db_lock:
            row = self._conn().execute(
                "SELECT status, result_json FROM jobs WHERE job_id = ?",
                (job_id,),
            ).fetchone()
        if not row:
            return None
        status, raw = row[0], row[1]
        if raw is None:
            return {"status": status, "records": [], "message": None}
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            return {"status": status, "records": [], "message": None}

    def revoke_token_for(self, job_id: str) -> Optional[str]:
        """Return the revoke token for a job (used by tests; do not expose in HTTP)."""
        with self._db_lock:
            row = self._conn().execute(
                "SELECT revoke_token FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
        return row[0] if row else None

    def shutdown(self, drain: bool = False) -> None:
        """Stop the worker. Not normally needed; daemon thread dies with the process."""
        if not drain:
            return
        self._q.join()

    # ----- internals -------------------------------------------------------

    def _conn(self) -> sqlite3.Connection:
        conn = sqlite3.connect(
            self.db_path, check_same_thread=False, isolation_level=None, timeout=30.0
        )
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        return conn

    def _init_db(self) -> None:
        Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)
        with self._db_lock, self._conn() as c:
            c.executescript(
                """
                CREATE TABLE IF NOT EXISTS jobs(
                    job_id        TEXT PRIMARY KEY,
                    revoke_token  TEXT NOT NULL,
                    status        TEXT NOT NULL,
                    payload_json  TEXT NOT NULL,
                    result_json   TEXT,
                    created_at    REAL NOT NULL,
                    started_at    REAL,
                    finished_at   REAL
                );
                CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
                """
            )

    def _worker_loop(self) -> None:
        # Local imports: keep heavy optional deps (BEguidedOFF) out of the
        # import graph until a job actually starts.
        from beguided.web.sse_logger import get_sse_handler

        sse = get_sse_handler()
        while True:
            job: Job = self._q.get()
            try:
                # Atomic claim: skip if cancelled before pickup.
                with self._db_lock:
                    cur = self._conn().execute(
                        "UPDATE jobs SET status='running', started_at=? "
                        "WHERE job_id=? AND status='queued'",
                        (time.time(), job.job_id),
                    )
                if cur.rowcount == 0:
                    logger.info(
                        "Job %s was cancelled before pickup; skipping.",
                        job.job_id,
                    )
                    continue

                self._running[job.job_id] = job
                # The SSE queue for this job is registered lazily by
                # /logs/stream on first connect. Until then, log lines from
                # this worker simply have no consumer — that's fine.
                logger.info("Job running: %s", job.job_id)
                try:
                    if job.cancel.is_set():
                        raise PipelineCancelled("cancelled before start")
                    result = execute_job(
                        job_id=job.job_id,
                        targets_df=job.payload["targets_df"],
                        base_editor_names=job.payload["base_editor_names"],
                        filter_mode=job.payload["filter_mode"],
                        search_offtargets=job.payload.get("search_offtargets", True),
                        guide_score_threshold=job.payload.get(
                            "guide_score_threshold", 0.5
                        ),
                        passthrough_extras=job.payload.get(
                            "passthrough_extras", False
                        ),
                        mode=job.payload.get("mode", "spreadsheet"),
                        custom_genome_b64=job.payload.get("custom_genome_b64"),
                        cancel_event=job.cancel,
                    )
                    status = result.get("status", "error")
                except PipelineCancelled:
                    status = "cancelled"
                    result = {
                        "status": status,
                        "records": [],
                        "message": "Cancelled by user.",
                    }
                    logger.info("Job %s cancelled mid-run.", job.job_id)
                except Exception as exc:
                    status = "error"
                    result = {"status": status, "records": [], "message": str(exc)}
                    logger.exception("Job %s failed: %s", job.job_id, exc)
                finally:
                    self._running.pop(job.job_id, None)
                    with self._db_lock:
                        self._conn().execute(
                            "UPDATE jobs SET status=?, finished_at=?, result_json=? "
                            "WHERE job_id=?",
                            (status, time.time(), json.dumps(result), job.job_id),
                        )
            finally:
                # Always release the queue slot, even if the job was cancelled
                # before pickup; otherwise `q.join()` would block forever.
                self._q.task_done()


def _payload_to_jsonable(payload: dict) -> dict:
    """Replace any non-JSON-serializable objects (e.g. pandas DataFrame) with markers.

    The full payload is *also* kept in memory in ``Job.payload``; the SQLite
    copy is purely for audit/recovery, so dropping DataFrames is fine.
    """
    out: dict = {}
    for k, v in payload.items():
        if k == "targets_df":
            try:
                rows = len(v) if v is not None else 0
            except TypeError:
                rows = 0
            out[k] = {"__type__": "dataframe", "rows": rows}
        elif k == "base_editor_names":
            out[k] = list(v) if v is not None else []
        elif k == "custom_genome_b64":
            # Don't bloat the SQLite audit row with the full FASTA payload;
            # record just the byte length for diagnostics.
            out[k] = {"__type__": "fasta_b64", "bytes": len(v) if v else 0}
        else:
            try:
                json.dumps(v)
                out[k] = v
            except TypeError:
                out[k] = repr(v)
    return out


# Module-level singleton — replaced by ``install`` from the app factory.
_scheduler: Optional[Scheduler] = None


def install(scheduler: Scheduler) -> Scheduler:
    global _scheduler
    _scheduler = scheduler
    return _scheduler


def current() -> Scheduler:
    if _scheduler is None:
        raise RuntimeError("Scheduler has not been installed; call install() at startup.")
    return _scheduler
