"""Tests for bystander analysis module."""

import pytest
import pandas as pd
from unittest.mock import patch

from beguided.bystander.analysis import analyze_guide, run_all


class TestAnalyzeGuideDetail:
    """Tests for the bystander_effect label returned by analyze_guide."""

    def _make_row(self, **overrides):
        defaults = {
            "chromosome": "chr1",
            "position": 100,
            "strand_used": "+",
            "editing_position": 5,
            "bystander_positions": "",
            "desired_base": "G",
            "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
            "pam_sequence": "AGG",
            "genome_base": "A",
            "genomic_context_upstream": "CCC",
            "genomic_context_downstream": "GGG",
        }
        defaults.update(overrides)
        return defaults

    @patch("beguided.bystander.analysis.get_cds_overlap")
    def test_no_cds_overlap_detail(self, mock_get_cds_overlap):
        mock_get_cds_overlap.return_value = ([], False)
        row = self._make_row()
        effect, *_ = analyze_guide(row, {})
        assert effect == "no_cds_overlap"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_cds_reconstruction_failed_detail(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        mock_get_cds_overlap.return_value = ([(90, 120, "+", "tx1")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = (None, None, None, None)
        row = self._make_row()
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "cds_reconstruction_failed"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_editing_position_non_coding_detail(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        # CDS is far downstream, so editing position falls outside it
        mock_get_cds_overlap.return_value = ([(200, 250, "+", "tx1")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAA", 0, 0, "+")
        row = self._make_row(editing_position=5)
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "editing_position_non_coding"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_bystander_outside_cds_detail(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        mock_get_cds_overlap.return_value = ([(90, 105, "+", "tx1")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        # editing_position=5 is inside CDS (90-105), but bystander at pos 20 is outside
        row = self._make_row(editing_position=5, bystander_positions="20")
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "bystander_outside_cds"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_success_effect_detail_matches(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        mock_get_cds_overlap.return_value = ([(90, 130, "+", "tx1")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ([], "M")
        mock_mark.return_value = "M"
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "no_bystanders"


class TestRunAll:
    """Tests for run_all DataFrame output."""

    def test_emits_single_bystander_effect_column(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [100],
            "strand_used": ["+"],
            "editing_position": [5],
            "bystander_positions": [""],
            "desired_base": ["G"],
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
            "pam_sequence": ["AGG"],
            "genome_base": ["A"],
            "genomic_context_upstream": ["CCC"],
            "genomic_context_downstream": ["GGG"],
        })

        with patch("beguided.bystander.analysis.get_cds_overlap") as mock_overlap:
            mock_overlap.return_value = ([], False)
            result = run_all(df, {})

        assert "bystander_effect" in result.columns
        assert "bystander_effect_detail" not in result.columns
        assert result.loc[0, "bystander_effect"] == "no_cds_overlap"
