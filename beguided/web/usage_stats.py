"""Anonymous web-usage statistics.

The web server records one CSV row per finished job so the operator can
answer high-level questions later: is the service actually used, how
often, and do people submit a handful of targets or hundreds? Those
numbers inform whether further development and long-term maintenance are
worthwhile, and support resource requests (e.g. a GPU for off-target
search).

The recorded data is deliberately anonymous:

* No IP address, no user-agent, no session or job identifier.
* Date granularity only (the server-local calendar day).
* No base-editor *names* — only **how many** were selected. (Custom
  profile names are user-typed free text and could in principle contain
  identifying information, so the names are never written.)
* No guide counts, no output contents, no uploaded data.

Under the GDPR such a file contains no personal data (Recital 26), so it
needs no legal basis, no consent, no cookie banner, and may be retained
indefinitely. The one hard rule for future changes: never add an IP
address, user-agent, or sub-day timestamp to this file, or it becomes
personal data and needs DPO review.

Recording must never affect job processing: every write error is caught
and logged, never raised into the caller.
"""

from __future__ import annotations

import csv
import datetime as _dt
import logging
import os
import threading
from pathlib import Path
from typing import List

logger = logging.getLogger("beguided.web.usage_stats")

# Column order is part of the on-disk contract — keep it stable so old and
# new CSVs stay analysable together.
COLUMNS: List[str] = [
    "date",
    "mode",
    "n_targets",
    "n_base_editors",
    "duration_s",
    "status",
]


def _monthly_path(base_path: str, today: _dt.date) -> str:
    """Return ``base_path`` with a ``_YYYY-MM`` suffix before its extension.

    ``/var/lib/beguided/usage_stats.csv`` becomes
    ``/var/lib/beguided/usage_stats_2026-10.csv``.
    """
    p = Path(base_path)
    stamped = f"{p.stem}_{today.strftime('%Y-%m')}{p.suffix or '.csv'}"
    return str(p.with_name(stamped))


class UsageStatsLogger:
    """Append-only CSV writer for anonymous per-job usage statistics.

    Thread-safe (the scheduler uses a single worker thread, but the lock
    keeps the class safe under tests and future concurrency) and
    fail-open: :meth:`record` swallows every error so a full/read-only
    disk can never break a job.
    """

    def __init__(self, path: str, enabled: bool = True):
        self.path = path
        self.enabled = bool(enabled)
        self._lock = threading.Lock()

    def record(
        self,
        mode: str,
        n_targets: int,
        n_base_editors: int,
        duration_s: float,
        status: str,
    ) -> None:
        """Append one anonymous row for a finished web job. Never raises."""
        if not self.enabled:
            return
        try:
            self._write_row(
                {
                    "date": _dt.date.today().isoformat(),
                    "mode": str(mode),
                    "n_targets": int(n_targets),
                    "n_base_editors": int(n_base_editors),
                    "duration_s": round(float(duration_s), 1),
                    "status": str(status),
                }
            )
        except Exception:
            logger.warning(
                "Usage stats write failed; job processing unaffected.",
                exc_info=True,
            )

    def _write_row(self, row: dict) -> None:
        target = _monthly_path(self.path, _dt.date.today())
        parent = os.path.dirname(target)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with self._lock:
            needs_header = not os.path.exists(target) or os.path.getsize(target) == 0
            with open(target, "a", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=COLUMNS)
                if needs_header:
                    writer.writeheader()
                writer.writerow(row)
