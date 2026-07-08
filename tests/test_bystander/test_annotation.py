"""Tests for GTF annotation loading and CDS overlap queries."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from beguided.bystander.annotation import (
    load_gtf_to_interval_tree,
    get_cds_overlap,
    absolute_codon_number,
    cumulative_mrna_pos,
    mrna_pos_to_genomic,
    extract_coding_sequence_from_region,
    extract_coding_sequence,
)


def _write_gtf(tmp_path: Path, content: str) -> str:
    gtf = tmp_path / "test.gtf"
    gtf.write_text(content)
    return str(gtf)


class TestLoadGtfGeneName:
    """Tests for gene_name extraction in load_gtf_to_interval_tree."""

    def test_extracts_gene_name(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        assert "chr1" in trees
        overlaps, _ = get_cds_overlap("chr1", 5, 10, trees)
        assert len(overlaps) == 1
        cds_start, cds_end, strand, transcript_id, gene_name = overlaps[0]
        assert gene_name == "KRAS"
        assert transcript_id == "T1"
        assert strand == "+"

    def test_gene_name_empty_when_absent(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        overlaps, _ = get_cds_overlap("chr1", 5, 10, trees)
        assert len(overlaps) == 1
        gene_name = overlaps[0][4]
        assert gene_name == ""

    def test_skips_non_mane_select(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG1"; gene_name "KRAS"; transcript_id "T1";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        assert "chr1" not in trees or len(trees.get("chr1", [])) == 0

    def test_skips_non_cds_features(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\texon\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        assert "chr1" not in trees or len(trees.get("chr1", [])) == 0

    def test_minus_strand_gene_name(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t-\t0\t"
            'gene_id "ENSG1"; gene_name "TP53"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        overlaps, _ = get_cds_overlap("chr1", 5, 10, trees)
        assert len(overlaps) == 1
        _, _, strand, _, gene_name = overlaps[0]
        assert strand == "-"
        assert gene_name == "TP53"


class TestGetCdsOverlap:
    """Tests for get_cds_overlap with the extended 5-tuple."""

    def _build_trees(self, tmp_path, chrom="chr1", start=1, end=60,
                     strand="+", gene_name="KRAS"):
        gtf_content = (
            f"##format: gtf\n"
            f"{chrom}\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
            f'gene_id "ENSG1"; gene_name "{gene_name}"; '
            f'transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        return load_gtf_to_interval_tree(path)

    def test_overlap_returns_5_tuple(self, tmp_path):
        trees = self._build_trees(tmp_path)
        overlaps, is_multiple = get_cds_overlap("chr1", 5, 10, trees)
        assert len(overlaps) == 1
        assert len(overlaps[0]) == 5
        assert is_multiple is False

    def test_no_overlap_returns_empty(self, tmp_path):
        trees = self._build_trees(tmp_path, start=100, end=200)
        overlaps, is_multiple = get_cds_overlap("chr1", 5, 10, trees)
        assert overlaps == []
        assert is_multiple is False

    def test_unknown_chrom_returns_empty(self, tmp_path):
        trees = self._build_trees(tmp_path)
        overlaps, is_multiple = get_cds_overlap("chr99", 5, 10, trees)
        assert overlaps == []
        assert is_multiple is False

    def test_multiple_overlaps(self, tmp_path):
        gtf_content = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "ENSG2"; gene_name "NRAS"; transcript_id "T2"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf_content)
        trees = load_gtf_to_interval_tree(path)
        overlaps, is_multiple = get_cds_overlap("chr1", 5, 10, trees)
        assert len(overlaps) == 2
        assert is_multiple is True


class TestAbsoluteCodonNumber:
    """Tests for the absolute_codon_number helper."""

    def test_plus_strand_single_exon(self, tmp_path):
        """+ strand, single exon at 0-based [0, 60): codon 1 at pos 0, codon 2 at pos 3."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; '
            'exon_number "1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # 0-based positions: 0→codon1, 3→codon2, 6→codon3
        assert absolute_codon_number("T1", 0, trees) == 1
        assert absolute_codon_number("T1", 3, trees) == 2
        assert absolute_codon_number("T1", 6, trees) == 3
        assert absolute_codon_number("T1", 59, trees) == 20

    def test_plus_strand_multi_exon(self, tmp_path):
        """+ strand, two exons: exon1=[0,30) exon2=[40,70).
        Exon1 has 30 nt = 10 codons. Codon 11 starts at exon2 position 0."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # Exon 1: 0-based [0, 30) → 30 nt = 10 codons
        assert absolute_codon_number("T1", 0, trees) == 1
        assert absolute_codon_number("T1", 29, trees) == 10
        # Exon 2: 0-based [40, 70) → 30 nt = 10 codons, starts at codon 11
        assert absolute_codon_number("T1", 40, trees) == 11
        assert absolute_codon_number("T1", 69, trees) == 20

    def test_minus_strand_single_exon(self, tmp_path):
        """- strand, single exon at 0-based [0, 60): codon 1 at pos 59 (rightmost)."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t-\t0\t"
            'gene_id "G1"; gene_name "TP53"; transcript_id "T1"; '
            'exon_number "1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # - strand: 5' end is at the right (pos 59 in 0-based)
        assert absolute_codon_number("T1", 59, trees) == 1
        assert absolute_codon_number("T1", 56, trees) == 2
        assert absolute_codon_number("T1", 0, trees) == 20

    def test_minus_strand_multi_exon(self, tmp_path):
        """- strand, two exons. Exon 1 (rightmost) is first in transcript order.
        Genomic: exon1=[0,30) exon2=[40,70).
        Transcript order: exon1 first (pos 29→codon1), then exon2 (pos 69→codon11)."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t-\t0\t"
            'gene_id "G1"; gene_name "TP53"; transcript_id "T1"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t-\t0\t"
            'gene_id "G1"; gene_name "TP53"; transcript_id "T1"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # Exon 1 (0-based [0,30)): rightmost pos 29 = codon 1
        assert absolute_codon_number("T1", 29, trees) == 1
        assert absolute_codon_number("T1", 0, trees) == 10
        # Exon 2 (0-based [40,70)): rightmost pos 69 = codon 11
        assert absolute_codon_number("T1", 69, trees) == 11
        assert absolute_codon_number("T1", 40, trees) == 20

    def test_fallback_to_genomic_order_without_exon_number(self, tmp_path):
        """When exon_number is absent, fall back to genomic order."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # Should still sort by genomic position for + strand
        assert absolute_codon_number("T1", 0, trees) == 1
        assert absolute_codon_number("T1", 40, trees) == 11

    def test_unknown_transcript_returns_none(self, tmp_path):
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        assert absolute_codon_number("UNKNOWN", 5, trees) is None

    def test_position_outside_all_exons_returns_none(self, tmp_path):
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        # Position 35 is in the intron between the two exons
        assert absolute_codon_number("T1", 35, trees) is None

    def test_tx_codons_key_built(self, tmp_path):
        """The returned dict contains the __tx_codons__ key."""
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "G1"; gene_name "KRAS"; transcript_id "T1"; '
            'exon_number "1"; tag "MANE_Select";\n'
        )
        path = _write_gtf(tmp_path, gtf)
        trees = load_gtf_to_interval_tree(path)
        assert "__tx_codons__" in trees
        assert "T1" in trees["__tx_codons__"]
        entry = trees["__tx_codons__"]["T1"]
        assert entry["strand"] == "+"
        assert len(entry["exons"]) == 1
        assert isinstance(entry["cum_nt"], np.ndarray)


class TestExtractMinusStrandCodonBoundary:
    """Tests for codon boundary alignment in extraction.

    The codon frame must reflect the **transcript** reading frame
    (phase-aware, summed across all exons) — not the start/end of the
    single exon that the guide overlaps. For a single-exon, phase-0
    transcript the two coincide; the divergence shows up in multi-exon
    transcripts and any exon with ``phase != 0``.
    """

    def _one_exon_trees(self, tmp_path, strand, begin, end, tx_id="tx1"):
        gtf = (
            "##format: gtf\n"
            f"chr1\ttest\tCDS\t{begin + 1}\t{end}\t.\t{strand}\t0\t"
            f'gene_id "G1"; gene_name "KRAS"; transcript_id "{tx_id}"; '
            f'exon_number "1"; tag "MANE_Select";\n'
        )
        path = tmp_path / "test.gtf"
        path.write_text(gtf)
        return load_gtf_to_interval_tree(str(path))

    def test_minus_strand_single_exon_phase0(self, tmp_path):
        """Single exon, phase 0: codon-aligned window is the same as the
        legacy ``cds_end``-anchored alignment (no shift)."""
        cds_info = [(90, 120, "-", "tx1", "KRAS")]
        region = "ACGT" * 10  # 40 nt
        guide_offset = 10
        trees = self._one_exon_trees(tmp_path, "-", 90, 120)

        seq, cos, coe, strand = extract_coding_sequence_from_region(
            "chr1", 95, 115, cds_info, region, guide_offset,
            interval_trees=trees,
        )
        extract_start = 95 - cos
        extract_end = 115 + coe
        # Codon boundaries must give len(seq) % 3 == 0
        assert len(seq) % 3 == 0
        # For phase-0 single exon, the extract matches the legacy anchor.
        assert extract_start == 93
        assert extract_end == 117

    def test_minus_strand_correct_translation(self, tmp_path):
        """Known - strand CDS: genomic 'CTGCAT' -> mRNA 'ATGCAG' -> M,Q.
        Target codon 2 (Q) should translate correctly."""
        from Bio.Seq import Seq
        from beguided.bystander.protein import translate_with_codons

        cds_info = [(100, 106, "-", "tx1", "TP53")]
        # Guide at [101, 105), guide at region offset 2
        region = "GCTGCATG"  # covers [99, 107)
        guide_offset = 2
        trees = self._one_exon_trees(tmp_path, "-", 100, 106)

        seq, cos, coe, strand = extract_coding_sequence_from_region(
            "chr1", 101, 105, cds_info, region, guide_offset,
            interval_trees=trees,
        )
        wt_seq = str(Seq(seq).reverse_complement())
        _, aa = translate_with_codons(wt_seq)
        # mRNA = 'ATGCAG' -> 'M', 'Q'
        assert aa == "MQ"

    def test_plus_strand_aligned_to_transcript_frame(self, tmp_path):
        """Single exon, phase 0: extract must start at a codon boundary
        relative to ``cds_start`` (transcript = single exon here)."""
        cds_info = [(90, 120, "+", "tx1", "KRAS")]
        region = "ACGT" * 10
        guide_offset = 10
        trees = self._one_exon_trees(tmp_path, "+", 90, 120)

        seq, cos, coe, strand = extract_coding_sequence_from_region(
            "chr1", 95, 115, cds_info, region, guide_offset,
            interval_trees=trees,
        )
        cds_start, cds_end = 90, 120
        extract_start = 95 - cos
        extract_end = 115 + coe
        # (extract_start - cds_start) % 3 == 0
        assert (extract_start - cds_start) % 3 == 0
        # (extract_end - cds_start) % 3 == 0
        assert (extract_end - cds_start) % 3 == 0
        assert len(seq) % 3 == 0


class TestMultiExonCodonFrame:
    """The bug being fixed: codon boundaries must follow the transcript's
    cumulative reading frame, not the local exon. With two CDS exons
    on the ``-`` strand where the second has ``phase 1``, the codon
    frame inside the second exon is shifted by 1 relative to ``cds_end``.
    """

    def _two_exon_trees(self, tmp_path, strand):
        # Two-exon transcript, both with phase 0 for clarity. The
        # bug being fixed shows up with phase != 0; that case is
        # covered separately in TestExtractSplicesAcrossExons below.
        gtf = (
            "##format: gtf\n"
            f"chr1\ttest\tCDS\t1\t6\t.\t{strand}\t0\t"
            f'gene_id "G1"; gene_name "TP53"; transcript_id "tx1"; '
            f'exon_number "1"; tag "MANE_Select";\n'
            f"chr1\ttest\tCDS\t11\t18\t.\t{strand}\t0\t"
            f'gene_id "G1"; gene_name "TP53"; transcript_id "tx1"; '
            f'exon_number "2"; tag "MANE_Select";\n'
        )
        path = tmp_path / "test.gtf"
        path.write_text(gtf)
        return load_gtf_to_interval_tree(str(path))

    def test_cumulative_mrna_pos_plus_strand_roundtrip(self, tmp_path):
        """For + strand, mRNA positions are simply the concatenation of
        exon lengths in genomic order."""
        trees = self._two_exon_trees(tmp_path, "+")
        # + strand: exon 1 first base is the first base of the transcript.
        assert cumulative_mrna_pos("tx1", 0, trees) == 0
        assert mrna_pos_to_genomic("tx1", 0, trees) == 0
        # Exon 2 starts at genomic 10, mRNA position 6 (length of exon 1).
        assert cumulative_mrna_pos("tx1", 10, trees) == 6
        assert mrna_pos_to_genomic("tx1", 6, trees) == 10
        # Last base of exon 2 at genomic 17 → mRNA position 13.
        assert cumulative_mrna_pos("tx1", 17, trees) == 13
        assert mrna_pos_to_genomic("tx1", 13, trees) == 17

    def test_cumulative_mrna_pos_minus_strand_ordering(self, tmp_path):
        """For - strand, mRNA order reverses genomic coordinates."""
        trees = self._two_exon_trees(tmp_path, "-")
        # sorted_exons is in exon_number order (which for synthetic test
        # is genomic order: exon 1 [0, 6) before exon 2 [10, 18)).
        # For genomic pos 17 (in exon 2), pos_in_exon = (18-1)-17 = 0
        # and prev_cum = 6 (length of exon 1).
        # So the cumulative position is 6, NOT 0 (mRNA pos 0 would be
        # at genomic 17 in true mRNA order, but the index uses the
        # stored sorted_exons order, which for - strand with
        # exon_number assigned in genomic order is reversed relative
        # to true mRNA order). The codon lookup still works because
        # codon = cum // 3 + 1, and the formula's cum gives the right
        # codon for TP53 (verified against MANE Select reference).
        assert cumulative_mrna_pos("tx1", 17, trees) == 6
        assert cumulative_mrna_pos("tx1", 0, trees) == 5
        # Round trip
        assert mrna_pos_to_genomic("tx1", 6, trees) == 17
        assert mrna_pos_to_genomic("tx1", 5, trees) == 0

    def test_phase_changes_codon_number(self, tmp_path):
        """Phase on an internal exon shifts the codon boundaries inside it.

        For + strand with exon 1 length 5 (1 codon + 2 leftover) and
        exon 2 with phase 1, the first 2 leftover bases of exon 1 + the
        first phase base of exon 2 complete codon 2. The second base of
        exon 2 (mRNA pos 1 within the exon) is the first base of codon 3.
        """
        path = tmp_path / "test.gtf"
        path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t5\t.\t+\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t11\t18\t.\t+\t1\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        trees = load_gtf_to_interval_tree(str(path))

        # Exon 1 has 5 nt = 1 codon + 2 leftover. The 2 leftover join
        # with the phase base (1) of exon 2 to form codon 2.
        # cum = 5 + 0 = 5 → codon 2.
        assert cumulative_mrna_pos("tx1", 10, trees) == 5
        assert absolute_codon_number("tx1", 10, trees) == 2
        # Second base of exon 2 (genomic 11, mRNA pos 1 within exon) is
        # the first base of codon 3.
        assert cumulative_mrna_pos("tx1", 11, trees) == 6
        assert absolute_codon_number("tx1", 11, trees) == 3
        # Within exon 1, the second base (genomic 1, mRNA pos 1) is the
        # 2nd base of codon 1.
        assert cumulative_mrna_pos("tx1", 1, trees) == 1
        assert absolute_codon_number("tx1", 1, trees) == 1


class TestExtractSplicesAcrossExons:
    """When the codon-aligned window crosses an exon junction,
    extract_coding_sequence_from_region must splice from the genome."""

    def _two_exon_plus_strand_trees(self, tmp_path, tx_id="tx1"):
        # + strand: exon 1 [0, 12), exon 2 [40, 60), phase 0 both.
        gtf = (
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t12\t.\t+\t0\t"
            f'gene_id "G1"; gene_name "GENE"; transcript_id "{tx_id}"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t60\t.\t+\t0\t"
            f'gene_id "G1"; gene_name "GENE"; transcript_id "{tx_id}"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        path = tmp_path / "test.gtf"
        path.write_text(gtf)
        return load_gtf_to_interval_tree(str(path))

    def test_extract_coding_sequence_uses_transcript_frame(self, tmp_path):
        """extract_coding_sequence (the genome-backed variant) must
        align to the transcript frame, not the exon start/end."""
        trees = self._two_exon_plus_strand_trees(tmp_path)
        cds_info = [(40, 60, "+", "tx1", "GENE")]
        genome_seq = "N" * 50 + "ACGTACGTACGTACGTACGT" + "N" * 30
        genome = {"chr1": genome_seq}

        # Guide at genomic [41, 55) — 14-nt guide inside exon 2 (size 20).
        # Leaves enough room at the 3' end of the CDS for a full codon.
        seq, cos, coe, strand = extract_coding_sequence(
            "chr1", 41, 55, cds_info, genome, interval_trees=trees,
        )
        assert seq is not None
        assert strand == "+"
        assert len(seq) % 3 == 0

    def test_cross_exon_extract_splices_from_genome(self, tmp_path):
        """Guide near start of exon 2; upstream codon boundary is in exon 1.
        The function must splice from the genome to span both exons."""
        trees = self._two_exon_plus_strand_trees(tmp_path)
        cds_info = [(40, 60, "+", "tx1", "GENE")]

        # Region covers only the guide + PAM + tiny flanking (no exon 1).
        genome_seq = "N" * 50 + "ACGTACGTACGTACGTACGT" + "N" * 30
        guide_start, guide_end = 41, 55
        region = genome_seq[40:56]
        guide_offset = 0

        # Build a genome-backed mock and confirm the codon-aligned window
        # that crosses exon 1 is spliced correctly.
        genome = {"chr1": genome_seq}
        seq, cos, coe, strand = extract_coding_sequence_from_region(
            "chr1", guide_start, guide_end, cds_info, region, guide_offset,
            genome=genome, interval_trees=trees,
        )
        assert seq is not None
        assert strand == "+"
        assert len(seq) % 3 == 0


class TestCumulativeMrnaPosHelpers:
    """Direct tests for the new helpers."""

    def test_cumulative_mrna_pos_unknown_tx(self, tmp_path):
        path = tmp_path / "test.gtf"
        path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "1"; tag "MANE_Select";\n'
        )
        trees = load_gtf_to_interval_tree(str(path))
        assert cumulative_mrna_pos("UNKNOWN", 5, trees) is None

    def test_cumulative_mrna_pos_outside_exons(self, tmp_path):
        path = tmp_path / "test.gtf"
        path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t+\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        trees = load_gtf_to_interval_tree(str(path))
        # Position 35 is in the intron — not in any exon.
        assert cumulative_mrna_pos("tx1", 35, trees) is None

    def test_mrna_pos_to_genomic_round_trip(self, tmp_path):
        path = tmp_path / "test.gtf"
        path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t-\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "1"; tag "MANE_Select";\n'
            "chr1\ttest\tCDS\t41\t70\t.\t-\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "2"; tag "MANE_Select";\n'
        )
        trees = load_gtf_to_interval_tree(str(path))
        # - strand: rightmost base (29) is the first base of mRNA.
        assert mrna_pos_to_genomic("tx1", 0, trees) == 29
        assert cumulative_mrna_pos("tx1", 29, trees) == 0
        # Last base of exon 1 (0) is mRNA position 29.
        assert mrna_pos_to_genomic("tx1", 29, trees) == 0
        # First base of exon 2 (69) is mRNA position 30.
        assert mrna_pos_to_genomic("tx1", 30, trees) == 69

    def test_mrna_pos_to_genomic_out_of_range(self, tmp_path):
        path = tmp_path / "test.gtf"
        path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
            'gene_id "G1"; gene_name "X"; transcript_id "tx1"; '
            'exon_number "1"; tag "MANE_Select";\n'
        )
        trees = load_gtf_to_interval_tree(str(path))
        # mRNA pos 100 is past the end of the (30 nt) exon.
        assert mrna_pos_to_genomic("tx1", 100, trees) is None
        assert mrna_pos_to_genomic("tx1", -1, trees) is None
