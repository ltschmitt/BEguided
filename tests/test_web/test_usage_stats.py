"""Tests for the anonymous web-usage statistics recorder.

These cover the on-disk contract of ``beguided/web/usage_stats.py`` plus
the anonymity guarantee: the CSV must never contain anything beyond the
fixed, identifier-free column set.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from beguided.web.usage_stats import COLUMNS, UsageStatsLogger, _monthly_path


def _read_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _csvs(directory: Path) -> list[Path]:
    return sorted(directory.glob("usage_stats_*.csv"))


def test_first_write_creates_header_and_row(tmp_path):
    base = tmp_path / "usage_stats.csv"
    logger = UsageStatsLogger(str(base))
    logger.record("spreadsheet", 5, 2, 1.234, "complete")

    files = _csvs(tmp_path)
    assert len(files) == 1
    rows = _read_rows(files[0])
    assert len(rows) == 1
    row = rows[0]
    # Exact column set, in order, with no extras.
    assert list(row.keys()) == COLUMNS
    assert row["mode"] == "spreadsheet"
    assert row["n_targets"] == "5"
    assert row["n_base_editors"] == "2"
    assert row["duration_s"] == "1.2"
    assert row["status"] == "complete"
    # Date-only granularity: YYYY-MM-DD, never a time component.
    assert len(row["date"]) == 10
    assert "T" not in row["date"]
    assert ":" not in row["date"]


def test_appends_without_duplicating_header(tmp_path):
    base = tmp_path / "usage_stats.csv"
    logger = UsageStatsLogger(str(base))
    logger.record("spreadsheet", 1, 1, 0.1, "complete")
    logger.record("sequence", 3, 4, 0.2, "complete")

    files = _csvs(tmp_path)
    assert len(files) == 1
    text = files[0].read_text()
    assert text.count("date,mode,n_targets") == 1
    rows = _read_rows(files[0])
    assert [r["mode"] for r in rows] == ["spreadsheet", "sequence"]


def test_monthly_rotation_uses_suffix(tmp_path):
    import datetime as dt

    base = tmp_path / "usage_stats.csv"
    # Two different months must land in two different files.
    p1 = Path(_monthly_path(str(base), dt.date(2026, 9, 30)))
    p2 = Path(_monthly_path(str(base), dt.date(2026, 10, 1)))
    assert p1.name == "usage_stats_2026-09.csv"
    assert p2.name == "usage_stats_2026-10.csv"


def test_disabled_creates_no_file(tmp_path):
    base = tmp_path / "usage_stats.csv"
    logger = UsageStatsLogger(str(base), enabled=False)
    logger.record("spreadsheet", 1, 1, 0.1, "complete")
    assert _csvs(tmp_path) == []


def test_parent_directory_is_created(tmp_path):
    base = tmp_path / "nested" / "deeper" / "usage_stats.csv"
    logger = UsageStatsLogger(str(base))
    logger.record("spreadsheet", 1, 1, 0.1, "complete")
    assert len(_csvs(base.parent)) == 1


def test_write_error_is_swallowed(tmp_path, caplog):
    """A broken target path must never raise into job processing."""
    import logging

    # A directory in place of the file makes the open() fail.
    base = tmp_path / "as_dir"
    base.mkdir()
    logger = UsageStatsLogger(str(base))  # path is a dir → write fails
    logger.record("spreadsheet", 1, 1, 0.1, "complete")  # must not raise
    # (No assertion on the file; the point is that we got here.)


def test_anonymity_guard_no_identifiers_in_output(tmp_path):
    """Nothing resembling an identifier may reach the CSV.

    The recorder API has no parameter for an IP/user-agent, so we also
    assert the column set is exactly the anonymous contract. This is the
    tripwire for any future change that adds personal data.
    """
    assert COLUMNS == [
        "date",
        "mode",
        "n_targets",
        "n_base_editors",
        "duration_s",
        "status",
    ]

    base = tmp_path / "usage_stats.csv"
    logger = UsageStatsLogger(str(base))
    logger.record("spreadsheet", 42, 7, 12.5, "complete")

    text = _csvs(tmp_path)[0].read_text()
    for forbidden in ("ip", "user-agent", "user_agent", "cookie", "token", "job_id"):
        assert forbidden not in text.lower()
