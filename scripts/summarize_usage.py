#!/usr/bin/env python3
"""Summarize anonymous BEguided web usage statistics.

Reads the monthly CSV files written by ``beguided/web/usage_stats.py``
and prints per-month aggregates suitable for reports, maintenance
decisions, and resource requests (e.g. GPU for off-target search).

Columns consumed: ``date, mode, n_targets, n_base_editors,
duration_s, status``.

Usage:
    python scripts/summarize_usage.py
    python scripts/summarize_usage.py --glob '/var/lib/beguided/usage_stats_*.csv'
    python scripts/summarize_usage.py --last 6
"""

from __future__ import annotations

import argparse
import csv
import glob
import statistics
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_GLOB = "/var/lib/beguided/usage_stats_*.csv"

# Target-count buckets chosen to answer "do people run a couple of targets
# or full screens?". Labels are inclusive ranges; the last is open-ended.
BUCKETS = [
    ("1-2", lambda n: 1 <= n <= 2),
    ("3-10", lambda n: 3 <= n <= 10),
    ("11-50", lambda n: 11 <= n <= 50),
    ("51-100", lambda n: 51 <= n <= 100),
    (">100", lambda n: n > 100),
]


def _to_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _load(patterns: list[str]) -> list[dict]:
    rows: list[dict] = []
    for pattern in patterns:
        for path in sorted(glob.glob(pattern)):
            with open(path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    if not row.get("date"):
                        continue
                    rows.append(row)
    return rows


def _month(row: dict) -> str:
    return str(row.get("date", ""))[:7] or "unknown"


def _print_month(month: str, rows: list[dict]) -> None:
    total = len(rows)
    by_status: dict[str, int] = defaultdict(int)
    by_mode: dict[str, int] = defaultdict(int)
    bucket_counts: dict[str, int] = defaultdict(int)
    durations: list[float] = []
    editor_counts: list[int] = []

    for row in rows:
        by_status[str(row.get("status", "unknown"))] += 1
        by_mode[str(row.get("mode", "unknown"))] += 1
        n_targets = _to_int(row.get("n_targets", "0"))
        matched = False
        for label, pred in BUCKETS:
            if pred(n_targets):
                bucket_counts[label] += 1
                matched = True
                break
        if not matched:
            bucket_counts["0"] += 1
        durations.append(_to_float(row.get("duration_s", "0")))
        editor_counts.append(_to_int(row.get("n_base_editors", "0")))

    print(f"\n=== {month} ===")
    print(f"  jobs total : {total}")
    print(
        "  by status  : "
        + ", ".join(f"{k}={v}" for k, v in sorted(by_status.items()))
    )
    print(
        "  by mode    : "
        + ", ".join(f"{k}={v}" for k, v in sorted(by_mode.items()))
    )
    print("  target-count buckets:")
    for label, _ in BUCKETS:
        print(f"    {label:>6} : {bucket_counts.get(label, 0)}")
    if durations:
        print(
            f"  duration_s : mean={statistics.mean(durations):.1f} "
            f"median={statistics.median(durations):.1f} "
            f"max={max(durations):.1f}"
        )
    if editor_counts:
        print(f"  base editors / job: mean={statistics.mean(editor_counts):.1f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--glob",
        dest="patterns",
        action="append",
        default=None,
        help=f"Glob for usage CSVs (repeatable). Default: {DEFAULT_GLOB}",
    )
    parser.add_argument(
        "--last",
        type=int,
        default=0,
        help="Show only the most recent N months (0 = all).",
    )
    args = parser.parse_args()

    patterns = args.patterns or [DEFAULT_GLOB]
    rows = _load(patterns)
    if not rows:
        print(f"No usage rows found for: {', '.join(patterns)}", file=sys.stderr)
        return 1

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[_month(row)].append(row)

    months = sorted(grouped)
    if args.last > 0:
        months = months[-args.last:]

    total_jobs = 0
    for month in months:
        _print_month(month, grouped[month])
        total_jobs += len(grouped[month])

    print(f"\nTotal jobs across {len(months)} month(s): {total_jobs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
