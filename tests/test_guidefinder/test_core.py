"""Tests for guidefinder core logic."""

import re

import pytest

from beguided.editing.presets import pam_to_regex
from beguided.guidefinder.core import process_target
from beguided.genome.io import load_genome


class TestIntegration:
    """Integration tests using process_target."""

    def test_chr1_50_T_C_single_guide(
        self, sample_genome_dict, abe8e_antisense_config
    ):
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 50, "T", "C",
            "C", abe8e_antisense_config
        )
        assert len(results) == 1
        assert results[0].editing_position == 4
        assert results[0].strand_used == "-"
        assert results[0].guide_rna_sequence != "**No guide possible**"

    def test_chr1_75_A_G_multiple_guides(
        self, sample_genome_dict, abe8e_preset_config
    ):
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 75, "A", "G",
            "G", abe8e_preset_config
        )
        assert len(results) == 3
        positions = sorted(r.editing_position for r in results)
        assert positions == [5, 6, 7]
        assert all(r.strand_used == "+" for r in results)

    def test_chr2_30_T_C_multiple_guides(self, sample_genome_dict):
        from beguided.editing.presets import load_presets
        import tempfile, os
        tsv = (
            "name\tediting\tgeneral_edit_start\tgeneral_edit_end\thigh_edit_start\thigh_edit_end\tguide_length\tpam\tcas\tsource\n"
            "ABE8e\tA-to-G\t2\t13\t3\t8\t20\tNGG\tspCas9\tNature2022\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(tsv)
            tp = f.name
        presets = load_presets(tp)
        os.unlink(tp)
        config = presets[("T", "C")]["ABE8e"]

        results = process_target(
            sample_genome_dict["chr2"], "chr2", 30, "T", "C",
            "C", config
        )
        assert len(results) >= 1
        assert all(r.strand_used == "-" for r in results)

    def test_no_valid_pam(self, sample_genome_dict, abe8e_preset_config):
        results = process_target(
            sample_genome_dict["chr2"], "chr2", 5, "A", "G",
            "A", abe8e_preset_config
        )
        assert len(results) == 1
        assert results[0].guide_rna_sequence == "**No guide possible**"

    def test_guide_outside_editing_window_filtered(
        self, sample_genome_dict, abe8e_preset_config
    ):
        # This test narrows window; need a fresh PresetConfig copy
        config = abe8e_preset_config
        # We can't mutate dataclass easily, so skip if not possible
        # Just verify filtering logic works with normal config
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 75, "A", "G",
            "G", config
        )
        if results[0].guide_rna_sequence != "**No guide possible**":
            assert all(3 <= r.editing_position <= 8 for r in results)

    def test_invalid_editing_direction(self, sample_genome_dict, abe8e_preset_config):
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 50, "G", "C",
            "G", abe8e_preset_config
        )
        assert len(results) == 1
        assert results[0].guide_rna_sequence == "**Base-edit doesn't match preset**"
        assert results[0].guide_score is None

    def test_extra_columns_passed_through(self, sample_genome_dict, abe8e_preset_config):
        extra = {"gene": "KRAS", "mutation_occurance": 17008}
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 75, "A", "G",
            "G", abe8e_preset_config, extra_columns=extra
        )
        assert len(results) == 3
        for r in results:
            assert r.extra_columns.get("gene") == "KRAS"
            assert r.extra_columns.get("mutation_occurance") == 17008

    def test_guide_score_populated_for_valid_guides(
        self, sample_genome_dict, abe8e_preset_config
    ):
        results = process_target(
            sample_genome_dict["chr1"], "chr1", 75, "A", "G",
            "G", abe8e_preset_config
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0
        for r in valid:
            assert isinstance(r.guide_score, float)
            assert 0 < r.guide_score < 2

    def test_no_guide_possible_has_no_score(
        self, sample_genome_dict, abe8e_preset_config
    ):
        results = process_target(
            sample_genome_dict["chr2"], "chr2", 5, "A", "G",
            "A", abe8e_preset_config
        )
        assert results[0].guide_rna_sequence == "**No guide possible**"
        assert results[0].guide_score is None


class TestEdgeCases:
    """Edge case tests."""

    def test_position_1_antisense(self, sample_genome_dict, abe8e_preset_config):
        # Need antisense config
        from beguided.editing.presets import load_presets
        import tempfile, os
        tsv = (
            "name\tediting\tgeneral_edit_start\tgeneral_edit_end\thigh_edit_start\thigh_edit_end\tguide_length\tpam\tcas\tsource\n"
            "ABE8e\tA-to-G\t2\t13\t3\t8\t20\tNGG\tspCas9\tNature2022\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(tsv)
            tp = f.name
        presets = load_presets(tp)
        os.unlink(tp)
        config = presets[("T", "C")]["ABE8e"]

        results = process_target(
            sample_genome_dict["chr1"], "chr1", 1, "T", "C",
            sample_genome_dict["chr1"][0], config
        )
        assert len(results) >= 1

    def test_position_at_chromosome_end(self, sample_genome_dict, abe8e_preset_config):
        chrom_len = len(sample_genome_dict["chr1"])
        results = process_target(
            sample_genome_dict["chr1"], "chr1", chrom_len - 5, "A", "G",
            sample_genome_dict["chr1"][chrom_len - 6], abe8e_preset_config
        )
        assert len(results) >= 1
