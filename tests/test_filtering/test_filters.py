"""Tests for filtering logic."""

import pytest
import pandas as pd

from beguided.filtering.filters import load_filter_criteria, filter_guides


class TestLoadFilterCriteria:
    """Tests for load_filter_criteria."""

    def test_defaults(self):
        criteria = load_filter_criteria()
        assert criteria["mode"] == "all"

    def test_override(self):
        criteria = load_filter_criteria(mode="silent")
        assert criteria["mode"] == "silent"


class TestFilterGuides:
    """Tests for filter_guides."""

    def test_no_bystanders_mode(self, synthetic_bystander_output):
        criteria = {"mode": "no_bystanders"}
        result = filter_guides(synthetic_bystander_output, criteria)
        assert len(result) == 3
        assert "TTTTTTTTTTTTTTTTTTTT" in result["guide_rna_sequence"].values

    def test_silent_mode(self, synthetic_bystander_output):
        criteria = {"mode": "silent"}
        result = filter_guides(synthetic_bystander_output, criteria)
        assert len(result) == 4  # 3 no_bystanders + 1 silent

    def test_no_bystanders_from_beguided(self, synthetic_beguided_output):
        criteria = load_filter_criteria(mode="no_bystanders")
        result = filter_guides(synthetic_beguided_output, criteria)
        assert len(result) == 2

    def test_removes_no_guide_possible(self, synthetic_bystander_output):
        df = synthetic_bystander_output.copy()
        no_guide = pd.DataFrame({
            "chromosome": ["chrX"],
            "position": [999],
            "target_base": ["A"],
            "desired_base": ["G"],
            "genome_base": ["G"],
            "strand_used": ["+"],
            "editing_position": [5],
            "bystander_positions": [""],
            "guide_rna_sequence": ["**No guide possible**"],
            "pam_sequence": ["XXX"],
            "base_editor": ["ABE8e"],
            "bystander_effect": ["no_bystanders"],
            "simulated_outcome": [""],
            "aa_wanted": ["X"],
            "aa_bystander": ["X"],
        })
        df = pd.concat([df, no_guide], ignore_index=True)
        criteria = {"mode": "all"}
        result = filter_guides(df, criteria)
        assert "**No guide possible**" not in result["guide_rna_sequence"].values

    def test_all_mode(self, synthetic_bystander_output):
        criteria = {"mode": "all"}
        result = filter_guides(synthetic_bystander_output, criteria)
        assert len(result) == 5

    def test_silent_mode_empty_df(self):
        empty_df = pd.DataFrame(columns=["chromosome", "position", "target_base"])
        criteria = {"mode": "silent"}
        result = filter_guides(empty_df, criteria)
        assert len(result) == 0

    def test_no_bystanders_mode_empty_df(self):
        empty_df = pd.DataFrame(columns=["chromosome", "position", "target_base"])
        criteria = {"mode": "no_bystanders"}
        result = filter_guides(empty_df, criteria)
        assert len(result) == 0

    def test_guide_score_threshold_removes_low_scores(self):
        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1", "chr1"],
            "guide_rna_sequence": ["A", "B", "C"],
            "guide_score": [0.8, 0.4, 0.5],
        })
        criteria = {"mode": "all", "guide_score_threshold": 0.5}
        result = filter_guides(df, criteria)
        assert len(result) == 2
        assert set(result["guide_rna_sequence"]) == {"A", "C"}

    def test_guide_score_threshold_keeps_nan(self):
        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1"],
            "guide_rna_sequence": ["A", "B"],
            "guide_score": [0.3, float("nan")],
        })
        criteria = {"mode": "all", "guide_score_threshold": 0.5}
        result = filter_guides(df, criteria)
        assert len(result) == 1
        assert result.iloc[0]["guide_rna_sequence"] == "B"

    def test_guide_score_threshold_skips_if_column_missing(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "guide_rna_sequence": ["A"],
        })
        criteria = {"mode": "all", "guide_score_threshold": 0.5}
        result = filter_guides(df, criteria)
        assert len(result) == 1

    def test_load_filter_criteria_reads_guide_score_threshold(self, tmp_path):
        import yaml
        cfg = {"guide_score_threshold": 0.7}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["guide_score_threshold"] == 0.7

    def test_load_filter_criteria_arg_overrides_file(self, tmp_path):
        import yaml
        cfg = {"guide_score_threshold": 0.7}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p), guide_score_threshold=0.9)
        assert criteria["guide_score_threshold"] == 0.9

    def test_skip_bystander_keeps_all_when_mode_is_no_bystanders(self):
        """Regression: when ``skip_bystander=True`` is set (sequence-mode
        jobs), ``mode='no_bystanders'`` must NOT drop guides whose
        ``bystander_positions`` is non-empty. This is the bug that
        previously caused every guide to be filtered out on A-rich
        sequences in sequence mode.
        """
        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1", "chr1"],
            "position": [10, 10, 10],
            "target_base": ["A", "A", "A"],
            "desired_base": ["G", "G", "G"],
            "guide_rna_sequence": ["AAAAA", "BBBBB", "CCCCC"],
            # Every guide has at least one bystander. Without the fix,
            # all three would be dropped under ``mode='no_bystanders'``.
            "bystander_positions": ["5", "3,7", "2,4,8"],
        })
        criteria = {"mode": "no_bystanders", "skip_bystander": True}
        result = filter_guides(df, criteria)
        assert len(result) == 3
        assert set(result["guide_rna_sequence"]) == {"AAAAA", "BBBBB", "CCCCC"}

    def test_skip_bystander_default_false(self):
        """Default behaviour: ``skip_bystander`` is False, so
        ``mode='no_bystanders'`` still drops guides with bystanders.
        """
        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1"],
            "position": [10, 10],
            "target_base": ["A", "A"],
            "desired_base": ["G", "G"],
            "guide_rna_sequence": ["AAAAA", "BBBBB"],
            "bystander_positions": ["", "5"],  # one clean, one with bystander
        })
        criteria = {"mode": "no_bystanders"}
        result = filter_guides(df, criteria)
        assert len(result) == 1
        assert result.iloc[0]["guide_rna_sequence"] == "AAAAA"

        """Regression: same guide sequence with different base editors must all survive filtering."""
        df = pd.DataFrame({
            "chromosome": ["1", "1", "1"],
            "position": [18700668, 18700668, 18700668],
            "target_base": ["T", "T", "T"],
            "desired_base": ["C", "C", "C"],
            "genome_base": ["C", "C", "C"],
            "strand_used": ["-", "-", "-"],
            "editing_position": [4, 4, 4],
            "bystander_positions": ["8", "8", "8,12"],
            "guide_rna_sequence": [
                "GGCAGTTACTGAACCAGACC",
                "GGCAGTTACTGAACCAGACC",
                "GGCAGTTACTGAACCAGACC",
            ],
            "pam_sequence": ["TGG", "TGG", "TGG"],
            "base_editor": ["ABE-NW1", "ABE-Umax", "ABE8e"],
            "bystander_effect": ["silent", "silent", "missense"],
            "simulated_outcome": ["", "", ""],
            "aa_wanted": ["X", "X", "X"],
            "aa_bystander": ["X", "X", "X"],
        })
        criteria = {"mode": "silent"}
        result = filter_guides(df, criteria)
        assert len(result) == 2
        assert set(result["base_editor"].values) == {"ABE-NW1", "ABE-Umax"}
