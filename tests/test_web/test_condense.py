"""Unit tests for ``beguided.web.condense``."""

from __future__ import annotations

import pandas as pd
import pytest


pytestmark = pytest.mark.skipif(
    __import__("importlib").util.find_spec("fastapi") is None,
    reason="web extras not installed",
)


def _records(*rows):
    """Tiny helper to build dict records without repeating keys."""
    keys = ("chromosome", "position", "guide_rna_sequence", "pam_sequence",
            "guide_score", "worst_offtarget_score", "base_editor")
    return [dict(zip(keys, r)) for r in rows]


def test_condense_records_merges_identical_rows():
    from beguided.web.condense import condense_records

    rows = _records(
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8e"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE9e"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8.20"),
        ("chr2", 200, "CCCCCCCCCCCCCCCCCCCC", "NGG", 0.4, 0.0, "ABE8e"),
    )
    out = condense_records(rows)
    assert len(out) == 2
    chr1, chr2 = out
    assert chr1["chromosome"] == "chr1"
    assert chr1["position"] == 100
    assert chr1["base_editor"].split("\n") == ["ABE8.20", "ABE8e", "ABE9e"]
    assert chr2["chromosome"] == "chr2"
    assert chr2["base_editor"] == "ABE8e"


def test_condense_records_keeps_rows_with_different_guides():
    from beguided.web.condense import condense_records

    rows = _records(
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8e"),
        ("chr1", 100, "TTTTTTTTTTTTTTTTTTTT", "NGG", 0.9, 1.2, "ABE8e"),
    )
    out = condense_records(rows)
    assert len(out) == 2
    assert all(r["base_editor"] == "ABE8e" for r in out)
    assert {r["guide_rna_sequence"] for r in out} == {
        "AAAAAAAAAAAAAAAAAAAA",
        "TTTTTTTTTTTTTTTTTTTT",
    }


def test_condense_records_preserves_first_occurrence_order():
    from beguided.web.condense import condense_records

    rows = _records(
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "Z-last"),
        ("chr2", 200, "CCCCCCCCCCCCCCCCCCCC", "NGG", 0.4, 0.0, "A-first"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "M-middle"),
    )
    out = condense_records(rows)
    assert [r["chromosome"] for r in out] == ["chr1", "chr2"]
    assert out[0]["base_editor"].split("\n") == ["M-middle", "Z-last"]


def test_condense_records_handles_empty_and_missing_column():
    from beguided.web.condense import condense_records

    assert condense_records([]) == []
    # When any record lacks `base_editor` the helper passes through unchanged
    # (mirrors condense_df's behaviour for frames missing the column).
    no_be = [{"chromosome": "chr1", "position": 100}]
    assert condense_records(no_be) is not no_be
    assert condense_records(no_be) == no_be


def test_condense_records_sorts_unique_editors_alphabetically():
    from beguided.web.condense import condense_records

    rows = _records(
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE9"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8.20"),
    )
    out = condense_records(rows)
    assert out[0]["base_editor"] == "ABE8\nABE8.20\nABE9"


def test_condense_df_matches_records_helper():
    """The DataFrame and record-list implementations should agree."""
    from beguided.web.condense import condense_df, condense_records

    df = pd.DataFrame(_records(
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE8e"),
        ("chr1", 100, "AAAAAAAAAAAAAAAAAAAA", "NGG", 0.9, 1.2, "ABE9e"),
        ("chr2", 200, "CCCCCCCCCCCCCCCCCCCC", "NGG", 0.4, 0.0, "ABE8e"),
    ))
    from_records = condense_records(df.to_dict(orient="records"))
    from_df = condense_df(df)
    assert from_df.shape[0] == len(from_records) == 2
    # base_editor is the only column that should differ between the two views
    df_dict = from_df.to_dict(orient="records")
    for left, right in zip(from_records, df_dict):
        assert left["base_editor"] == right["base_editor"]
        for k in left:
            if k != "base_editor":
                assert left[k] == right[k], k


def test_condense_df_returns_unchanged_when_no_base_editor_column():
    from beguided.web.condense import condense_df

    df = pd.DataFrame({"chromosome": ["chr1"], "position": [100]})
    out = condense_df(df)
    assert out.equals(df)


def test_condense_df_empty_input():
    from beguided.web.condense import condense_df

    df = pd.DataFrame(columns=["chromosome", "position", "base_editor"])
    assert condense_df(df).empty
