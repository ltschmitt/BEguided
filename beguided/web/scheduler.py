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

import hashlib
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

# How often the retention sweeper wakes up and deletes finished jobs whose
# ``finished_at`` is older than the configured retention window. The value is
# deliberately shorter than the default 24-hour retention so the worst-case
# over-retention is bounded by this interval.
SWEEPER_INTERVAL_SECONDS = 15 * 60


def _anonymise_ip(ip: str, salt: str) -> str:
    """Return a short salted SHA-256 fingerprint of ``ip`` for log use.

    Plaintext IP addresses in long-lived log files are personal data under
    GDPR. With a salted hash the IP is not reversible from the log alone
    and the salt can be rotated independently of the log retention window.

    If ``salt`` is empty (e.g. unset in tests) the prefix ``ip:`` is
    returned so the field is still traceable in local development without
    writing the raw address to disk.
    """
    if not ip:
        return "ip:none"
    if not salt:
        return f"ip:raw:{ip[:8]}"
    digest = hashlib.sha256(f"{salt}:{ip}".encode("utf-8")).hexdigest()
    return f"ip:{digest[:16]}"


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
        job_retention_seconds: float = 24 * 3600,
        log_ip_salt: str = "",
        sweeper_interval_seconds: float = SWEEPER_INTERVAL_SECONDS,
    ):
        self.db_path = db_path
        self.queue_max = queue_max
        self.cooldown_seconds = cooldown_seconds
        self.job_retention_seconds = max(0.0, float(job_retention_seconds))
        self.log_ip_salt = log_ip_salt or ""
        self.sweeper_interval_seconds = max(1.0, float(sweeper_interval_seconds))

        self._q: "queue.Queue[Job]" = queue.Queue(maxsize=queue_max)
        self._last_submit: dict[str, float] = {}
        self._db_lock = threading.Lock()
        self._running: dict[str, Job] = {}
        self._worker: Optional[threading.Thread] = None
        self._sweeper: Optional[threading.Thread] = None
        self._sweeper_stop = threading.Event()

    # ----- lifecycle -------------------------------------------------------

    def start(self) -> None:
        """Idempotent: start the worker + sweeper threads and init the DB schema."""
        self._init_db()
        if self._worker is None or not self._worker.is_alive():
            self._worker = threading.Thread(
                target=self._worker_loop, name="beguided-worker", daemon=True
            )
            self._worker.start()
            logger.info(
                "Scheduler started (queue_max=%d, cooldown=%ds, retention=%ds, db=%s)",
                self.queue_max,
                int(self.cooldown_seconds),
                int(self.job_retention_seconds),
                self.db_path,
            )
        if self.job_retention_seconds > 0 and (
            self._sweeper is None or not self._sweeper.is_alive()
        ):
            self._sweeper_stop.clear()
            self._sweeper = threading.Thread(
                target=self._sweeper_loop,
                name="beguided-sweeper",
                daemon=True,
            )
            self._sweeper.start()
            logger.info(
                "Retention sweeper started (interval=%ds, retention=%ds)",
                int(self.sweeper_interval_seconds),
                int(self.job_retention_seconds),
            )

    def stop(self) -> None:
        """Signal the sweeper to exit. Worker thread is a daemon and dies with the process."""
        self._sweeper_stop.set()

    # ----- public API ------------------------------------------------------

    def submit(self, client_ip: str, payload: dict) -> Job:
        """Enqueue a new job. Raises ``Cooldown`` or ``QueueFull``."""
        ip = client_ip or "unknown"
        self._prune_last_submit(now=time.monotonic())
        last = self._last_submit.get(ip, 0.0)
        now = time.monotonic()
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
        logger.info(
            "Job queued: %s (%s)",
            job.job_id,
            _anonymise_ip(ip, self.log_ip_salt),
        )
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
                CREATE INDEX IF NOT EXISTS idx_jobs_finished_at ON jobs(finished_at);
                """
            )

    # ----- retention sweeper ----------------------------------------------

    def purge_old(self, older_than_seconds: Optional[float] = None) -> int:
        """Delete every finished job whose ``finished_at`` is older than the cutoff.

        Uses the scheduler's ``job_retention_seconds`` unless an explicit
        ``older_than_seconds`` is provided. Returns the number of deleted
        rows. Running or queued jobs are never deleted.
        """
        if older_than_seconds is None:
            older_than_seconds = self.job_retention_seconds
        if older_than_seconds <= 0:
            return 0
        cutoff = time.time() - older_than_seconds
        with self._db_lock:
            cur = self._conn().execute(
                "DELETE FROM jobs "
                "WHERE finished_at IS NOT NULL AND finished_at < ?",
                (cutoff,),
            )
        deleted = cur.rowcount
        if deleted:
            logger.info(
                "Retention sweep: deleted %d job(s) older than %.0fs.",
                deleted,
                older_than_seconds,
            )
        return deleted

    def purge_all_finished(self) -> int:
        """Delete every job whose status is not queued or running.

        Used by the admin ``/admin/purge-all`` endpoint to honour manual
        GDPR Art. 17 erasure requests and to wipe the table after a
        security incident. Returns the number of deleted rows.
        """
        with self._db_lock:
            cur = self._conn().execute(
                "DELETE FROM jobs "
                "WHERE status NOT IN ('queued', 'running')"
            )
        deleted = cur.rowcount
        return deleted

    def _prune_last_submit(self, now: float) -> None:
        """Drop expired entries from ``_last_submit``.

        The dict records the last submit time per client IP for cooldown
        enforcement. Without pruning it grows unbounded over the lifetime of
        the process. Entries are dropped once they are 4x the cooldown
        window old — far enough past the actual cooldown that no client
        could still be affected.
        """
        if not self._last_submit:
            return
        cutoff = now - max(4.0 * self.cooldown_seconds, 60.0)
        stale = [ip for ip, t in self._last_submit.items() if t < cutoff]
        for ip in stale:
            self._last_submit.pop(ip, None)

    def _sweeper_loop(self) -> None:
        """Wake every ``sweeper_interval_seconds`` and delete old jobs."""
        while not self._sweeper_stop.is_set():
            # Sleep first so the initial DB row inserts from a fresh start
            # are not deleted before they are queued. ``wait`` returns True
            # when the event is set, which we treat as "exit early".
            if self._sweeper_stop.wait(self.sweeper_interval_seconds):
                return
            try:
                self.purge_old()
            except Exception:
                logger.exception("Retention sweep failed; will retry next interval.")
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
                        keep_targets_without_guides=job.payload.get(
                            "keep_targets_without_guides", False
                        ),
                        mode=job.payload.get("mode", "spreadsheet"),
                        custom_genome_b64=job.payload.get("custom_genome_b64"),
                        use_custom_profiles=job.payload.get("use_custom_profiles", False),
                        custom_profiles_b64=job.payload.get("custom_profiles_b64"),
                        merged_profiles_b64=job.payload.get("merged_profiles_b64"),
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
        elif k == "custom_profiles_b64":
            out[k] = {"__type__": "tsv_b64", "bytes": len(v) if v else 0}
        elif k == "merged_profiles_b64":
            out[k] = {"__type__": "tsv_b64", "bytes": len(v) if v else 0}
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
