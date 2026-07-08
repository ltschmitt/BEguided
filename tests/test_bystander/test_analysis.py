"""Tests for bystander analysis module."""

import numpy as np
import pytest
import pandas as pd
from unittest.mock import patch

from beguided.bystander.analysis import analyze_guide, run_all


def _make_trees_with_tx(tx_id="tx1", strand="+", exons=None):
    """Build a mock interval_trees dict with a __tx_codons__ entry.

    Args:
        tx_id: Transcript ID.
        strand: CDS strand.
        exons: List of (begin_0based, end_0based, exon_number_or_None).
               Defaults to a single exon at [90, 130).
    """
    if exons is None:
        exons = [(90, 130, 1)]
    lengths = [e[1] - e[0] for e in exons]
    cum_nt = np.cumsum(lengths, dtype=np.int64)
    return {
        "__tx_codons__": {
            tx_id: {
                "strand": strand,
                "exons": exons,
                "cum_nt": cum_nt,
            }
        }
    }


class TestAnalyzeGuideDetail:
    """Tests for the bystander_effect label returned by analyze_guide."""

    def _make_row(self, **overrides):
        defaults = {
            "chromosome": "chr1",
            "position": 100,
            "strand_used": "+",
            "editing_position": 5,
            "bystander_positions": "",
            "target_base": "A",
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
        mock_get_cds_overlap.return_value = ([(90, 120, "+", "tx1", "KRAS")], False)
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
        mock_get_cds_overlap.return_value = ([(200, 250, "+", "tx1", "KRAS")], False)
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
        mock_get_cds_overlap.return_value = ([(90, 105, "+", "tx1", "KRAS")], False)
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
        mock_get_cds_overlap.return_value = ([(90, 130, "+", "tx1", "KRAS")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ([], "M")
        mock_mark.return_value = "M"
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "no_bystanders"


class TestGeneEditDesired:
    """Tests for the gene_edit_desired field returned by analyze_guide."""

    def _make_row(self, **overrides):
        defaults = {
            "chromosome": "chr1",
            "position": 100,
            "strand_used": "+",
            "editing_position": 5,
            "bystander_positions": "",
            "target_base": "A",
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
    def test_no_cds_overlap_returns_empty(self, mock_get_cds_overlap):
        mock_get_cds_overlap.return_value = ([], False)
        row = self._make_row()
        result = analyze_guide(row, {})
        assert result[4] == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_plus_strand(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """+ strand, single exon at [90,130). position=100, editing_position=5.
        guide_start=95, genomic_editing_pos=99.
        pos_in_exon = 99-90 = 9, codon = 9//3+1 = 4.
        target_base='A', desired_base='G', coding_seq='AAAA...'
        mrna_codon='AAA', offset_in_codon=1
        target_codon='AAA'->K, desired_codon='AGA'->R.
        -> 'KRAS K4R'."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "KRAS K4R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_absolute_multi_exon(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Multi-exon transcript: exon1=[0,30) exon2=[90,130).
        Target at genomic 99 falls in exon2. Exon1 contributes 30 nt = 10 codons.
        pos_in_exon = 99-90 = 9, codon = (30+9)//3+1 = 14.
        target_codon='AAA'->K, desired_codon='AGA'->R. -> 'KRAS K14R'."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(0, 30, 1), (90, 130, 2)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "KRAS K14R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_empty_when_no_gene_name(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """CDS overlaps but gene_name is empty -> gene_edit_desired is ''."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_minus_strand_cds(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Minus-strand CDS at [90,130). Antisense guide (+ strand).
        abs_codon=11. coding_seq='AAAA...' (22 chars), mrna=rc='TTTT...'
        mrna_edit_pos=22-0-5=17, local_codon=6, offset=2
        mrna_codon='TTT', mrna_target=COMPLEMENT['A']='T', mrna_desired=COMPLEMENT['G']='C'
        target_codon='TTT'->F, desired_codon='TTC'->F.
        -> 'TP53 F11F'."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "-", "tx1", "TP53")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 2, "-")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "RRRRRRR")
        mock_mark.return_value = "R"
        trees = _make_trees_with_tx("tx1", "-", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="", strand_used="+")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "TP53 F11F"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_multi_transcript_uses_first(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Multi-transcript overlap: gene_edit_desired uses the first CDS entry."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS"), (90, 130, "+", "tx2", "NRAS")],
            True,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert effect == "multi-transcript"
        assert gene_edit == "KRAS K4R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_gene_edit_desired_empty_when_editing_pos_non_coding(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        mock_get_cds_overlap.return_value = (
            [(200, 250, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAA", 0, 0, "+")
        row = self._make_row(editing_position=5)
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, {"chr1": True})
        assert effect == "editing_position_non_coding"
        assert gene_edit == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    def test_gene_edit_desired_populated_when_bystander_outside_cds(
        self, mock_translate, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Target is in CDS but a bystander is outside -> gene_edit_desired
        is still populated, effect is 'bystander_outside_cds'."""
        # CDS at [90, 105). position=100, editing_position=5, strand='+'.
        # guide_start=95, genomic_editing_pos=99 -> in [90,105) ✓
        # Bystander at pos 20: genomic_pos=95+19=114 -> outside [90,105) ✓
        mock_get_cds_overlap.return_value = (
            [(90, 105, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_translate.return_value = ("AAAAAA", "MM")
        trees = _make_trees_with_tx("tx1", "+", [(90, 105, 1)])
        row = self._make_row(editing_position=5, bystander_positions="20")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert effect == "bystander_outside_cds"
        # abs_codon=4, target_codon='AAA'->K, desired_codon='AGA'->R
        assert gene_edit == "KRAS K4R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_empty_when_no_tx_index(
        self, mock_mark, mock_translate, mock_simulate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """When __tx_codons__ is absent, absolute_codon_number returns None
        and gene_edit_desired is ''."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+")
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        # No __tx_codons__ key — simulates an old-style interval_trees dict
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, {"chr1": True})
        assert gene_edit == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.simulate_editing")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_minus_strand_correct_aa(
        self, mock_mark, mock_simulate, mock_reconstruct, mock_get_cds_overlap
    ):
        """End-to-end - strand test with REAL extract_coding_sequence_from_region.

        Build a - strand CDS where the mRNA is known: 'ATGCAG' * 5 -> M,Q,M,Q,...
        Genomic = reverse_complement = 'CTGCAT' * 5 at [100, 130).
        Antisense guide (+ strand) targeting codon 10.
        target_base='A', desired_base='G', cds_strand='-'
        mrna_target=COMPLEMENT['A']='T', mrna_desired=COMPLEMENT['G']='C'
        mrna_codon='CAG', offset=0
        target_codon='TAG'->* (stop), desired_codon='CAG'->Q.
        -> 'TP53 *10Q'
        """
        from Bio.Seq import Seq

        mrna = "ATGCAG" * 5  # 30 nt
        genomic = str(Seq(mrna).reverse_complement())  # 'CTGCAT' * 5

        # CDS at [100, 130), strand '-'
        mock_get_cds_overlap.return_value = (
            [(100, 130, "-", "tx1", "TP53")], False,
        )
        # Region = genomic sequence starting at genomic 100
        mock_reconstruct.return_value = (genomic, 0)
        # Don't mock extract — use the real function
        mock_simulate.side_effect = lambda seq, offset, positions, base, strand: seq
        mock_mark.return_value = "X"
        trees = _make_trees_with_tx("tx1", "-", [(100, 130, 1)])

        # + strand guide (antisense), position=103, editing_position=3
        # guide_start = 103 - 3 = 100, guide_end = 120
        # genomic_editing_pos = 100 + 2 = 102
        # In tx_codons: pos_in_exon = (130-1) - 102 = 27, codon = 27//3 + 1 = 10
        # mRNA codon 10 = CAG (Q)
        # target_base='A' -> mrna_target='T', target_codon='TAG'->* (stop)
        # desired_base='G' -> mrna_desired='C', desired_codon='CAG'->Q
        row = self._make_row(
            position=103, editing_position=3, strand_used="+",
            bystander_positions="", desired_base="G",
        )
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "TP53 *10Q"


class TestRunAll:
    """Tests for run_all DataFrame output."""

    def test_emits_single_bystander_effect_column(self):
        df = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [100],
            "strand_used": ["+"],
            "editing_position": [5],
            "bystander_positions": [""],
            "target_base": ["A"],
            "desired_base": ["G"],
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
            "pam_sequence": ["AGG"],
            "genome_base": "A",
            "genomic_context_upstream": ["CCC"],
            "genomic_context_downstream": ["GGG"],
        })

        with patch("beguided.bystander.analysis.get_cds_overlap") as mock_overlap:
            mock_overlap.return_value = ([], False)
            result = run_all(df, {})

        assert "bystander_effect" in result.columns
        assert "bystander_effect_detail" not in result.columns
        assert "gene_edit_desired" in result.columns
        assert result.loc[0, "bystander_effect"] == "no_cds_overlap"
        assert result.loc[0, "gene_edit_desired"] == ""
