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
        mock_extract.return_value = (None, None, None, None, 0, 0)
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
        mock_extract.return_value = ("AAAA", 0, 0, "+", 0, 0)
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
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
        # editing_position=5 is inside CDS (90-105), but bystander at pos 20 is outside
        row = self._make_row(editing_position=5, bystander_positions="20")
        effect, *_ = analyze_guide(row, {"chr1": True})
        assert effect == "bystander_outside_cds"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_success_effect_detail_matches(
        self, mock_mark, mock_translate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        mock_get_cds_overlap.return_value = ([(90, 130, "+", "tx1", "KRAS")], False)
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
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
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_plus_strand(
        self, mock_mark, mock_translate,
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
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "KRAS K4R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_absolute_multi_exon(
        self, mock_mark, mock_translate,
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
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(0, 30, 1), (90, 130, 2)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "KRAS K14R"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_empty_when_no_gene_name(
        self, mock_mark, mock_translate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """CDS overlaps but gene_name is empty -> gene_edit_desired is ''."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        trees = _make_trees_with_tx("tx1", "+", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_minus_strand_cds(
        self, mock_mark, mock_translate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Minus-strand CDS at [90,130). Antisense guide (+ strand).
        abs_codon=11. coding_seq='AAAA...' (22 chars), mrna=rc='TTTT...'
        mep = anchor_mrna - mrna_pos(99) = 30 - 13 = 17, local_codon=6, offset=2
        mrna_codon='TTT', mrna_target=COMPLEMENT['A']='T', mrna_desired=COMPLEMENT['G']='C'
        target_codon='TTT'->F, desired_codon='TTC'->F.
        -> 'TP53 F11F'."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "-", "tx1", "TP53")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 2, "-", 0, 0)
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "RRRRRRR")
        mock_mark.return_value = "R"
        trees = _make_trees_with_tx("tx1", "-", [(90, 130, 1)])
        row = self._make_row(editing_position=5, bystander_positions="", strand_used="+")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)
        assert gene_edit == "TP53 F11F"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_multi_transcript_uses_first(
        self, mock_mark, mock_translate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """Multi-transcript overlap: gene_edit_desired uses the first CDS entry."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS"), (90, 130, "+", "tx2", "NRAS")],
            True,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
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
        mock_extract.return_value = ("AAAA", 0, 0, "+", 0, 0)
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
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
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
    @patch("beguided.bystander.analysis.translate_with_codons")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_empty_when_no_tx_index(
        self, mock_mark, mock_translate,
        mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """When __tx_codons__ is absent, absolute_codon_number returns None
        and gene_edit_desired is ''."""
        mock_get_cds_overlap.return_value = (
            [(90, 130, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("CCC" + "A" * 20 + "AGG" + "GGG", 3)
        mock_extract.return_value = ("AAAAAAAAAAAAAAAAAAAAAA", 0, 0, "+", 0, 0)
        mock_translate.return_value = ("AAAAAAAAAAAAAAA", "MMMMMMM")
        mock_mark.return_value = "M"
        # No __tx_codons__ key — simulates an old-style interval_trees dict
        row = self._make_row(editing_position=5, bystander_positions="")
        effect, _seq, _aa_w, _aa_b, gene_edit = analyze_guide(row, {"chr1": True})
        assert gene_edit == ""

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.mark_mutations")
    def test_gene_edit_desired_minus_strand_correct_aa(
        self, mock_mark, mock_reconstruct, mock_get_cds_overlap
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


class TestSplicingAwareCodonOffset:
    """Regression tests for the codon-offset bug fixed by switching from
    the genomic-distance ``codon_offset`` to ``cumulative_mrna_pos``.

    The bug: ``codon_offset_start`` / ``codon_offset_end`` are genomic
    distances (``guide_start - extract_start``, ``extract_end - guide_end``)
    that were used as indices into the *spliced* ``mrna_seq``. When the
    codon-aligned window spans an intron, the genomic distance is much
    larger than the corresponding ``mrna_seq`` index, so the computed
    edit position lands out of bounds, ``gene_edit_desired`` ends up
    ``"?"``, and the AA / simulated-outcome brackets are silently
    dropped by ``mark_mutations``.

    The fix: ``mep = cumulative_mrna_pos(edit) - cumulative_mrna_pos(anchor)``
    where ``anchor`` is ``extract_start`` for ``+`` CDS strand and
    ``extract_end - 1`` for ``-`` CDS strand. ``cumulative_mrna_pos``
    is splicing-aware (it skips introns), so ``mep`` is always a
    correct 0-based index into ``mrna_seq``.
    """

    def _make_row(self, **overrides):
        defaults = {
            "chromosome": "chr1",
            "position": 100,
            "strand_used": "-",
            "editing_position": 8,
            "bystander_positions": "",
            "target_base": "T",
            "desired_base": "C",
            "guide_rna_sequence": "X" * 20,
            "pam_sequence": "NNNNNN",
            "genome_base": "A",
            "genomic_context_upstream": "AAAA",
            "genomic_context_downstream": "TTTT",
        }
        defaults.update(overrides)
        return defaults

    def _two_exon_minus_strand_trees(self):
        """- strand transcript with two exons separated by an intron.

        Exon A (exon_number=1, mRNA 5' side):  [150, 200)  (50 nt)
        Intron:                                [100, 150)
        Exon B (exon_number=2, mRNA 3' side):  [ 20, 100)  (80 nt)
        Total mRNA = 130 nt.
        """
        return _make_trees_with_tx(
            "tx1", "-",
            [(150, 200, 1), (20, 100, 2)],
        )

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_minus_strand_codon_offset_window_spans_intron(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """- strand guide whose codon-aligned window crosses an intron.

        Setup:
        - Guide at - strand, position=100 (1-based, edit at + strand coord 99)
          editing_position=8.
        - Local CDS = exon B [20, 100).
        - Codon-aligned window in mRNA coords: mRNA pos 48..62 (5 codons).
        - Translated to + strand coords: window = [87, 152) = 65 bp.
        - Spliced coding_seq = 13 nt from exon B + 2 nt from exon A = 15 nt.
        - codon_offset_end (genomic distance from extract_end to guide_end)
          = 152 - 107 = 45.

        With the OLD ``codon_offset_end`` math:
          mep = 45 + 8 - 1 = 52, but len(mrna_seq) = 15 -> OOB.
          gene_edit_desired would be ``"TEST ?17?"`` and brackets missing.

        With ``cumulative_mrna_pos``:
          mrna_pos(99) = 50, anchor_mrna_pos(151) = 48,
          mep = 50 - 48 = 2 -> in bounds, gene_edit_desired is computed.
        """
        trees = self._two_exon_minus_strand_trees()
        mock_get_cds_overlap.return_value = (
            [(20, 100, "-", "tx1", "TEST")], False,
        )
        mock_reconstruct.return_value = ("placeholder_region", 0)
        # coding_seq, cos, coe, cds_strand, extract_start, extract_end.
        # cos=0 (guide at extract_start), coe=45 (extract_end far past
        # guide_end due to splicing). extract_start/extract_end match the
        # codon-aligned window in + strand coords (spans both exons).
        mock_extract.return_value = ("ATGCAGATGCAGATG", 0, 45, "-", 87, 152)

        row = self._make_row()
        _effect, sim, aa_w, aa_b, gene_edit = analyze_guide(row, trees)

        # Old code: "TEST ?17?" with no brackets.
        # New code: real AA letters at codon 17 (mRNA pos 50 / 3 + 1).
        assert "?" not in gene_edit, (
            f"Expected real AA letters in gene_edit_desired, got {gene_edit!r}"
        )
        assert gene_edit == "TEST Q17Q", (
            f"Expected 'TEST Q17Q', got {gene_edit!r}"
        )

        # mep = 2 in the spliced mrna_seq; the edit base (position 2 of
        # mrna_seq = position 12 of coding_seq) is mutated to mrna_desired_base
        # (= COMPLEMENT["C"] = "G") for cds_strand="-".
        #
        # mrna_seq_unmutated = reverse_complement("ATGCAGATGCAGATG")
        #                    = "CATCTGCATCTGCAT"
        # After mutation: "CAGCTGCATCTGCAT" (15 chars, no insertion).
        assert sim == "CA[G]CTGCATCTGCAT", (
            f"Expected bracketed mutation at mrna position 2, got {sim!r}"
        )

        # Codon 1 of the local AA sequence carries the bracket.
        # AAs: Q, L, H, L, H -> "[Q]LHLH"
        assert aa_w == "[Q]LHLH", f"Expected '[Q]LHLH', got {aa_w!r}"
        assert aa_b == "[Q]LHLH", f"Expected '[Q]LHLH', got {aa_b!r}"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_minus_strand_codon_offset_window_within_single_exon(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """- strand guide whose codon-aligned window stays in one exon.

        Control case: no splicing happens, codon_offset_end stays small,
        so the fix must produce the same result as the original logic.
        """
        # Single-exon - strand transcript at [100, 200), 100 nt.
        trees = _make_trees_with_tx("tx1", "-", [(100, 200, 1)])
        mock_get_cds_overlap.return_value = (
            [(100, 200, "-", "tx1", "TP53")], False,
        )
        mock_reconstruct.return_value = ("placeholder_region", 0)
        # coding_seq matches the actual 20-bp window: guide covers
        # + strand [118, 138), 20 bp, all within exon [100, 200).
        mock_extract.return_value = ("A" * 20, 0, 0, "-", 0, 0)

        # Edit at + strand coord 130 (1-based 131), ep=8.
        # guide_start = 118, guide_end = 138.
        row = self._make_row(position=131)
        _effect, _sim, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)

        # Edit at coord 130 = mRNA pos 69 = middle base of codon 24.
        # mrna_seq (reverse_complement of "A"*20) = "T"*20.
        # mep = 7 (mrna_seq[7]), local_codon = 3, offset = 1 (middle).
        # mrna_codon = mrna_seq[6:9] = "TTT".
        # mrna_target = COMPLEMENT["T"] = "A", mrna_desired = COMPLEMENT["C"] = "G".
        # target_codon = "T" + "A" + "T" = "TAT" -> Y
        # desired_codon = "T" + "G" + "T" = "TGT" -> C
        assert gene_edit == "TP53 Y24C"

    @patch("beguided.bystander.analysis.get_cds_overlap")
    @patch("beguided.bystander.analysis.reconstruct_genomic_region")
    @patch("beguided.bystander.analysis.extract_coding_sequence_from_region")
    def test_plus_strand_codon_offset_window_spans_intron(
        self, mock_extract, mock_reconstruct, mock_get_cds_overlap
    ):
        """+ strand guide whose codon-aligned window crosses an intron.

        Mirrors test_minus_strand_codon_offset_window_spans_intron for the
        + CDS strand to confirm the fix works for both strands.

        Setup:
        - + strand transcript: exon A [20, 100), intron [100, 150),
          exon B [150, 200). Total mRNA = 130 nt.
        - Guide extends upstream past the 5' end of exon B (into the intron)
          so the codon-aligned window reaches back into exon A.
        - Guide coords [145, 165), edit at coord 152 (1-based 153), ep=8.
        - Window in mRNA coords: 78..95 (18 nt = 6 codons).
        - Window in + strand coords: [98, 166) = 68 bp.
        - Spliced coding_seq = 2 (exon A) + 16 (exon B) = 18 nt.
        - codon_offset_start (genomic) = 47. With the OLD formula,
          mep = 47 + 7 = 54, far past len(mrna_seq) = 18 -> OOB.
        """
        trees = _make_trees_with_tx(
            "tx1", "+",
            [(20, 100, 1), (150, 200, 2)],
        )
        mock_get_cds_overlap.return_value = (
            [(150, 200, "+", "tx1", "KRAS")], False,
        )
        mock_reconstruct.return_value = ("placeholder_region", 0)
        # Spliced coding_seq matching the expected window shape:
        # 2 nt from exon A + 16 nt from exon B = 18 chars.
        mock_extract.return_value = ("AATGCAGATGCAGATGCA", 47, 1, "+", 98, 166)

        # + strand guide at coord 152 (1-based 153), ep=8.
        # guide_start = 145, guide_end = 165.
        row = self._make_row(
            position=153, strand_used="+", target_base="A", desired_base="G",
        )
        _effect, _sim, _aa_w, _aa_b, gene_edit = analyze_guide(row, trees)

        # Without fix: codon_offset_start=47 -> mep=54, OOB.
        # With fix: cumulative_mrna_pos correctly skips the intron.
        assert "?" not in gene_edit, (
            f"Expected real AA letters, got {gene_edit!r}"
        )
        # Format check; exact AA letters depend on the synthetic coding_seq.
        assert gene_edit.startswith("KRAS "), (
            f"Expected 'KRAS ...' format, got {gene_edit!r}"
        )
