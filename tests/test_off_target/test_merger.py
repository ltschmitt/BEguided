"""Tests for off-target result aggregation with BEguidedOFF CCTOP output."""

import pytest
import pandas as pd

from beguided.off_target.merger import aggregate_beguidedoff_results, merge_with_input


class TestAggregateBEguidedOFF:
    """Tests for aggregate_beguidedoff_results."""

    def test_basic(self, synthetic_beguidedoff_output):
        summary = aggregate_beguidedoff_results(synthetic_beguidedoff_output)
        assert len(summary) == 3
        assert "guide_seq" in summary.columns
        assert "worst_offtarget_score" in summary.columns
        assert "offtarget_count" in summary.columns
        assert "worst_offtarget_strand" in summary.columns

    def test_picks_worst_lowest_score(self, synthetic_beguidedoff_output):
        summary = aggregate_beguidedoff_results(synthetic_beguidedoff_output)
        # T-guide has a perfect match (score 0.0) → worst off-target
        row = summary[summary["guide_seq"] == "TTTTTTTTTTTTTTTTTTTT"]
        assert len(row) == 1
        assert row.iloc[0]["worst_offtarget_score"] == 0.0
        assert row.iloc[0]["worst_offtarget_strand"] == "-"

    def test_empty(self):
        empty = pd.DataFrame()
        summary = aggregate_beguidedoff_results(empty)
        assert len(summary) == 0
        assert list(summary.columns) == [
            "guide_seq", "worst_offtarget_score",
            "worst_offtarget_mismatch_pos", "offtarget_count",
            "worst_offtarget_chr", "worst_offtarget_pos",
            "worst_offtarget_strand",
            "worst_offtarget_sequence",
        ]

    def test_offtarget_count_from_threshold(self, synthetic_beguidedoff_output):
        summary = aggregate_beguidedoff_results(synthetic_beguidedoff_output)
        row = summary[summary["guide_seq"] == "AAAAAAAAAAAAAAAAAAAA"]
        assert row.iloc[0]["offtarget_count"] == 5


class TestMergeWithInput:
    """Tests for merge_with_input with BEguidedOFF summary."""

    def test_basic(self, synthetic_bystander_output, synthetic_beguidedoff_output):
        summary = aggregate_beguidedoff_results(synthetic_beguidedoff_output)
        merged = merge_with_input(synthetic_bystander_output, summary)
        assert "worst_offtarget_score" in merged.columns
        assert "offtarget_count" in merged.columns
        assert "worst_offtarget_strand" in merged.columns

    def test_missing_column_raises(self):
        df = pd.DataFrame({"other": ["AAA"]})
        summary = pd.DataFrame({"guide_seq": ["AAA"]})
        with pytest.raises(ValueError):
            merge_with_input(df, summary)

    def test_guides_without_offtargets_get_zero(self, synthetic_bystander_output):
        import pandas as pd
        ot = pd.DataFrame({
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
            "score": [1.5],
            "n_offtargets_below_threshold": [2],
            "ot_chromosome": ["chr1"],
            "ot_position": [100],
            "ot_sequence": ["TTTTTTTTTTTTTTTTTTTT"],
        })
        summary = aggregate_beguidedoff_results(ot)
        merged = merge_with_input(synthetic_bystander_output, summary)

        # NNNNN... guide should get inf (no off-targets) since not in OT output
        nnn_row = merged[merged["guide_rna_sequence"] == "NNNNNNNNNNNNNNNNNNNN"]
        assert nnn_row.iloc[0]["worst_offtarget_score"] == float("inf")
        assert nnn_row.iloc[0]["offtarget_count"] == 0
        assert nnn_row.iloc[0]["worst_offtarget_strand"] == ""
