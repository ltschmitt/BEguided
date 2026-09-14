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


class TestGCFilter:
    """Tests for the gc_content_filter behaviour in filter_guides."""

    def _make_df(self, gc_values):
        return pd.DataFrame({
            "chromosome": [f"chr1"] * len(gc_values),
            "guide_rna_sequence": [f"g{i}" for i in range(len(gc_values))],
            "gc_content": gc_values,
        })

    def test_off_keeps_everything(self):
        df = self._make_df([0.0, 25.0, 55.0, 75.0, 100.0])
        result = filter_guides(df, {"mode": "all", "gc_content_filter": "off"})
        assert len(result) == 5

    def test_good_keeps_only_40_to_60(self):
        df = self._make_df([20.0, 39.9, 40.0, 50.0, 60.0, 60.1, 90.0])
        result = filter_guides(df, {"mode": "all", "gc_content_filter": "good"})
        assert sorted(result["gc_content"].tolist()) == [40.0, 50.0, 60.0]

    def test_acceptable_keeps_30_to_70(self):
        df = self._make_df([20.0, 29.9, 30.0, 55.0, 70.0, 70.1, 90.0])
        result = filter_guides(df, {"mode": "all", "gc_content_filter": "acceptable"})
        assert sorted(result["gc_content"].tolist()) == [30.0, 55.0, 70.0]

    def test_nan_is_kept(self):
        df = self._make_df([float("nan"), 50.0, 25.0])
        result = filter_guides(df, {"mode": "all", "gc_content_filter": "good"})
        assert len(result) == 2  # NaN + 50.0
        assert set(result["gc_content"].dropna().tolist()) == {50.0}

    def test_missing_column_is_skipped(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "guide_rna_sequence": ["g1"],
        })
        result = filter_guides(df, {"mode": "all", "gc_content_filter": "good"})
        assert len(result) == 1

    def test_invalid_mode_raises(self):
        df = self._make_df([50.0])
        with pytest.raises(ValueError):
            filter_guides(df, {"mode": "all", "gc_content_filter": "bogus"})


class TestHomopolymerFilter:
    """Tests for the filter_homopolymers behaviour in filter_guides."""

    def _make_df(self, hp_values):
        return pd.DataFrame({
            "chromosome": [f"chr1"] * len(hp_values),
            "guide_rna_sequence": [f"g{i}" for i in range(len(hp_values))],
            "homopolymer": hp_values,
        })

    def test_disabled_keeps_all(self):
        df = self._make_df([True, False, True])
        result = filter_guides(df, {"mode": "all", "filter_homopolymers": False})
        assert len(result) == 3

    def test_enabled_drops_true_rows(self):
        df = self._make_df([True, False, True])
        result = filter_guides(df, {"mode": "all", "filter_homopolymers": True})
        assert len(result) == 1
        assert result.iloc[0]["homopolymer"] == False  # noqa: E712

    def test_nan_is_kept(self):
        df = self._make_df([True, float("nan"), False])
        result = filter_guides(df, {"mode": "all", "filter_homopolymers": True})
        assert len(result) == 2

    def test_missing_column_is_skipped(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "guide_rna_sequence": ["g1"],
        })
        result = filter_guides(df, {"mode": "all", "filter_homopolymers": True})
        assert len(result) == 1


class TestLoadFilterCriteriaNewKeys:
    """Tests for the new keys in load_filter_criteria."""

    def test_defaults_are_off_and_false(self):
        criteria = load_filter_criteria()
        assert criteria["gc_content_filter"] == "off"
        assert criteria["filter_homopolymers"] is False

    def test_yaml_root_keys(self, tmp_path):
        import yaml
        cfg = {"gc_content_filter": "good", "filter_homopolymers": True}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["gc_content_filter"] == "good"
        assert criteria["filter_homopolymers"] is True

    def test_yaml_filtering_block(self, tmp_path):
        import yaml
        cfg = {"filtering": {"gc_content_filter": "acceptable"}}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["gc_content_filter"] == "acceptable"

    def test_invalid_gc_value_raises(self):
        with pytest.raises(ValueError):
            load_filter_criteria(gc_content_filter="bogus")

    def test_kwargs_override_defaults(self):
        criteria = load_filter_criteria(
            gc_content_filter="good", filter_homopolymers=True
        )
        assert criteria["gc_content_filter"] == "good"
        assert criteria["filter_homopolymers"] is True


class TestFilterCombination:
    """Combinations of GC / homopolymer / score / mode filters."""

    def test_all_filters_compose(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"] * 5,
            "guide_rna_sequence": ["g0", "g1", "g2", "g3", "g4"],
            "guide_score": [0.8, 0.4, 0.6, 0.9, 0.7],
            "gc_content": [50.0, 50.0, 75.0, 55.0, 50.0],
            "homopolymer": [False, False, False, True, False],
        })
        criteria = {
            "mode": "all",
            "guide_score_threshold": 0.5,
            "gc_content_filter": "good",
            "filter_homopolymers": True,
        }
        result = filter_guides(df, criteria)
        # Only g0 and g4 survive: score>=0.5, gc in [40,60], not homopolymer.
        assert set(result["guide_rna_sequence"]) == {"g0", "g4"}


class TestTrinucEditScoreFilter:
    """Tests for the trinuc_edit_score_min filter (ABE only; CBE = NaN kept)."""

    def _make_df(self, scores):
        return pd.DataFrame({
            "chromosome": ["chr1"] * len(scores),
            "guide_rna_sequence": [f"g{i}" for i in range(len(scores))],
            "trinuc_edit_score": scores,
        })

    def test_none_threshold_keeps_all(self):
        df = self._make_df([0.01, 0.5, 1.0, float("nan")])
        result = filter_guides(df, {"mode": "all", "trinuc_edit_score_min": None})
        assert len(result) == 4

    def test_threshold_drops_low_scores(self):
        df = self._make_df([0.01, 0.5, 0.99, 1.0])
        result = filter_guides(df, {"mode": "all", "trinuc_edit_score_min": 0.5})
        assert sorted(result["trinuc_edit_score"].tolist()) == [0.5, 0.99, 1.0]

    def test_nan_is_kept(self):
        """CBE rows (NaN) are never dropped by the trinuc filter."""
        df = self._make_df([0.01, float("nan"), 0.8])
        result = filter_guides(df, {"mode": "all", "trinuc_edit_score_min": 0.5})
        assert len(result) == 2  # NaN + 0.8
        assert set(result["guide_rna_sequence"]) == {"g1", "g2"}

    def test_missing_column_is_skipped(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "guide_rna_sequence": ["g1"],
        })
        result = filter_guides(df, {"mode": "all", "trinuc_edit_score_min": 0.5})
        assert len(result) == 1

    def test_boundary_value_kept(self):
        """A score exactly equal to the threshold is kept (>=)."""
        df = self._make_df([0.5, 0.49])
        result = filter_guides(df, {"mode": "all", "trinuc_edit_score_min": 0.5})
        assert len(result) == 1
        assert result.iloc[0]["guide_rna_sequence"] == "g0"

    def test_combines_with_other_filters(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"] * 4,
            "guide_rna_sequence": ["g0", "g1", "g2", "g3"],
            "guide_score": [0.8, 0.3, 0.8, 0.3],
            "trinuc_edit_score": [0.9, 0.9, 0.1, 0.1],
        })
        result = filter_guides(df, {
            "mode": "all",
            "guide_score_threshold": 0.5,
            "trinuc_edit_score_min": 0.5,
        })
        # Only g0 survives: score>=0.5 AND trinuc>=0.5.
        assert set(result["guide_rna_sequence"]) == {"g0"}


class TestLoadFilterCriteriaTrinucKey:
    """Tests for trinuc_edit_score_min in load_filter_criteria."""

    def test_default_is_none(self):
        criteria = load_filter_criteria()
        assert criteria["trinuc_edit_score_min"] is None

    def test_yaml_root_key(self, tmp_path):
        import yaml
        cfg = {"trinuc_edit_score_min": 0.5}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["trinuc_edit_score_min"] == 0.5

    def test_yaml_filtering_block(self, tmp_path):
        import yaml
        cfg = {"filtering": {"trinuc_edit_score_min": 0.7}}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["trinuc_edit_score_min"] == 0.7

    def test_yaml_null_disables(self, tmp_path):
        import yaml
        cfg = {"trinuc_edit_score_min": None}
        p = tmp_path / "cfg.yaml"
        p.write_text(yaml.dump(cfg))
        criteria = load_filter_criteria(config_file=str(p))
        assert criteria["trinuc_edit_score_min"] is None

    def test_kwarg_overrides_default(self):
        criteria = load_filter_criteria(trinuc_edit_score_min=0.8)
        assert criteria["trinuc_edit_score_min"] == 0.8

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            load_filter_criteria(trinuc_edit_score_min="not_a_number")

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            load_filter_criteria(trinuc_edit_score_min=1.5)
        with pytest.raises(ValueError):
            load_filter_criteria(trinuc_edit_score_min=-0.1)


class TestFillMissingTargets:
    """Tests for the ``fill_missing_targets`` option of filter_guides."""

    def test_off_by_default_no_placeholders(self, synthetic_bystander_output):
        result = filter_guides(synthetic_bystander_output, {"mode": "no_bystanders"})
        assert len(result) == 3
        assert "**Guides removed by filter**" not in set(result["guide_rna_sequence"])

    def test_adds_guides_removed_by_filter_rows(self, synthetic_bystander_output):
        result = filter_guides(
            synthetic_bystander_output,
            {"mode": "no_bystanders"},
            fill_missing_targets=True,
        )
        # 3 surviving guides + 2 placeholder rows (one per removed target).
        assert len(result) == 5
        placeholders = result[result["guide_rna_sequence"] == "**Guides removed by filter**"]
        assert len(placeholders) == 2
        assert set(zip(placeholders["chromosome"], placeholders["position"])) == {
            ("chr1", 100),
            ("chr2", 300),
        }
        # Identity/base columns preserved on the placeholder rows.
        assert set(placeholders["target_base"]) == {"A"}
        assert set(placeholders["desired_base"]) == {"G"}

    def test_preserves_pipeline_marker(self):
        """Rows the pipeline already flagged ("**No guide possible**") keep
        that marker; filter_guides re-adds them when a filter would drop them."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "position": [10, 20],
                "target_base": ["A", "A"],
                "desired_base": ["G", "G"],
                "guide_rna_sequence": [
                    "**No guide possible**",
                    "AAAAAAAAAAAAAAAAAAAA",
                ],
                "bystander_effect": [None, "silent"],
                "bystander_positions": ["", ""],
                "guide_score": [None, 0.9],
            }
        )
        result = filter_guides(df, {"mode": "all"}, fill_missing_targets=True)
        assert len(result) == 2
        markers = set(result["guide_rna_sequence"])
        assert "**No guide possible**" in markers
        assert "AAAAAAAAAAAAAAAAAAAA" in markers

    def test_empty_result_still_emits_rows(self):
        df = pd.DataFrame(
            {
                "chromosome": ["chr1"],
                "position": [10],
                "target_base": ["A"],
                "desired_base": ["G"],
                "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
                "guide_score": [0.1],
            }
        )
        result = filter_guides(
            df, {"mode": "all", "guide_score_threshold": 0.5},
            fill_missing_targets=True,
        )
        assert len(result) == 1
        assert result.iloc[0]["guide_rna_sequence"] == "**Guides removed by filter**"
        assert result.iloc[0]["chromosome"] == "chr1"
        assert result.iloc[0]["position"] == 10

    def test_fill_off_still_drops_everything(self):
        df = pd.DataFrame(
            {
                "chromosome": ["chr1"],
                "position": [10],
                "target_base": ["A"],
                "desired_base": ["G"],
                "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
                "guide_score": [0.1],
            }
        )
        result = filter_guides(
            df, {"mode": "all", "guide_score_threshold": 0.5}
        )
        assert len(result) == 0

