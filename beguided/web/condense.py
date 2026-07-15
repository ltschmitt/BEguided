"""Condense rows of the BEguided results table that differ only in
``base_editor`` into a single row with newline-separated editor names.

The web UI exposes this as a "Condense results" toggle so a guide that
appears once per selected base editor collapses to one row listing every
editor that produces it. Both the in-browser renderer (a faithful JS port
in ``beguided/web/static/app.js``) and the ``/download/results`` endpoint
share this algorithm.

Algorithm:
    1. Drop ``base_editor`` from the grouping key (everything else).
    2. Group identical rows.
    3. Aggregate ``base_editor`` as ``"\\n".join(sorted(unique))``.
    4. Preserve the first occurrence's values for every other column.

If ``base_editor`` is missing the frame is returned unchanged.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

import pandas as pd

BASE_EDITOR_COLUMN = "base_editor"
CELL_SEPARATOR = "\n"


def condense_records(records: list[dict]) -> list[dict]:
    """Condense a list of result records by ``base_editor``.

    Mirrors ``condense_df`` for the JSON-record shape the web UI consumes.
    Records are grouped on every column except ``base_editor``; the merged
    row carries the first occurrence's values for those columns and the
    unique ``base_editor`` values joined with ``"\\n"`` in sorted order.

    Records that lack a ``base_editor`` key are passed through unchanged.
    """
    if not records:
        return list(records)
    if any(BASE_EDITOR_COLUMN not in r for r in records):
        return list(records)

    group_keys: list[tuple] = []
    buckets: dict[tuple, List[dict]] = {}
    for r in records:
        key = _row_key(r)
        if key not in buckets:
            buckets[key] = []
            group_keys.append(key)
        buckets[key].append(r)

    out: list[dict] = []
    for key in group_keys:
        members = buckets[key]
        merged = dict(members[0])
        editors = sorted({str(m[BASE_EDITOR_COLUMN]) for m in members})
        merged[BASE_EDITOR_COLUMN] = CELL_SEPARATOR.join(editors)
        out.append(merged)
    return out


def condense_df(df: pd.DataFrame) -> pd.DataFrame:
    """Condense a result DataFrame by ``base_editor``.

    Rows whose values match in every column except ``base_editor`` are
    collapsed into one row whose ``base_editor`` cell lists the merged
    editors separated by newlines in alphabetical order.

    Frames that lack a ``base_editor`` column are returned unchanged.
    """
    if df.empty or BASE_EDITOR_COLUMN not in df.columns:
        return df

    key_cols = [c for c in df.columns if c != BASE_EDITOR_COLUMN]

    def _join_editors(series: pd.Series) -> str:
        # pandas groupby preserves the group member order, so a sorted
        # unique is enough to get the documented alphabetical join.
        return CELL_SEPARATOR.join(sorted({str(v) for v in series}))

    grouped = df.groupby(key_cols, sort=False, dropna=False)
    agg = grouped.agg({BASE_EDITOR_COLUMN: _join_editors})
    # groupby.agg with a single column returns a Series with a possibly
    # different shape; restore the DataFrame layout.
    if isinstance(agg, pd.Series):
        agg = agg.to_frame(BASE_EDITOR_COLUMN)
    else:
        agg = agg.rename(columns={BASE_EDITOR_COLUMN: BASE_EDITOR_COLUMN})
    agg = agg.reset_index()

    # Reorder columns to match the original frame so downstream readers
    # (the openpyxl writer, TSV writer, etc.) see the same schema.
    return agg.reindex(columns=df.columns)


def _row_key(row: dict) -> tuple:
    """Build the grouping key for a record, ignoring ``base_editor``.

    Implemented explicitly so it stays in sync with ``condense_df``'s
    groupby columns and so None / missing values don't crash the tuple
    hash.
    """
    return tuple(
        (k, _hashable(row.get(k))) for k in row if k != BASE_EDITOR_COLUMN
    )


def _hashable(value):
    """Coerce a value into a hashable form for the grouping key.

    Pandas may emit ``numpy.float64`` / ``numpy.int64`` / ``Timestamp``
    scalars which don't compare equal to their Python counterparts in a
    tuple. Converting via ``item()`` / ``isoformat()`` keeps the keys
    stable and equality-comparable.
    """
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return tuple(_hashable(v) for v in value)
    if isinstance(value, dict):
        return tuple(sorted((k, _hashable(v)) for k, v in value.items()))
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return item()
        except (ValueError, TypeError):
            pass
    if hasattr(value, "isoformat"):
        try:
            return value.isoformat()
        except (TypeError, ValueError):
            pass
    return value


__all__ = [
    "BASE_EDITOR_COLUMN",
    "CELL_SEPARATOR",
    "condense_records",
    "condense_df",
]
