"""Tests for bystander analysis with context columns."""

import pytest

from beguided.genome.context import reconstruct_genomic_region
from beguided.bystander.annotation import extract_coding_sequence_from_region


class TestGenomicRegionReconstruction:
    """Tests for reconstruct_genomic_region."""

    def test_sense_guide(self):
        row = {
            "strand_used": "+",
            "guide_rna_sequence": "AAAAGAAAAAAAAAAAAAAA",
            "pam_sequence": "AGG",
            "editing_position": 5,
            "genome_base": "G",
            "desired_base": "A",
            "genomic_context_upstream": "CCCC",
            "genomic_context_downstream": "TTT",
        }
        region, offset = reconstruct_genomic_region(row)
        assert region == "CCCCAAAAGAAAAAAAAAAAAAAAAGGTTT"
        assert offset == 4

    def test_sense_short_context(self):
        row = {
            "strand_used": "+",
            "guide_rna_sequence": "GAAAAAAAAAAAAAAAAAAA",
            "pam_sequence": "AGG",
            "editing_position": 1,
            "genome_base": "G",
            "desired_base": "A",
            "genomic_context_upstream": "A",
            "genomic_context_downstream": "CC",
        }
        region, offset = reconstruct_genomic_region(row)
        assert region == "AGAAAAAAAAAAAAAAAAAAAAGGCC"
        assert offset == 1

    def test_antisense_guide(self):
        """genome_guide restored from rc(guide_rna) with genome_base."""
        row = {
            "strand_used": "-",
            "guide_rna_sequence": "GAAAAAAAAAAAAAAAAAAA",
            "pam_sequence": "AGG",  # rc → CCT in genome
            "editing_position": 1,
            "genome_base": "C",
            "desired_base": "G",
            "genomic_context_upstream": "TT",
            "genomic_context_downstream": "AA",
        }
        region, offset = reconstruct_genomic_region(row)
        # rc("GAAAAAAAAAAAAAAAAAAA") = TTTTTTTTTTTTTTTTTTTC, replace idx15 w/ C -> same
        assert region == "AACCTTTTTTTTTTTTTTTTTTTTCTT"
        assert offset == 5

    def test_antisense_editing_position_4(self):
        row = {
            "strand_used": "-",
            "guide_rna_sequence": "TCTAATGATCGCTAGATCGA",
            "pam_sequence": "AGG",
            "editing_position": 4,
            "genome_base": "C",
            "desired_base": "C",
            "genomic_context_upstream": "GATA",
            "genomic_context_downstream": "GTA",
        }
        region, offset = reconstruct_genomic_region(row)
        # rc("TCTAATGATCGCTAGATCGA") = TCGATCTAGCGATCATTAGA
        # idx = 20-4 = 16; rc[16] = 'T' replaced with genome_base 'C'
        # genome_guide = TCGATCTAGCGATCATCAGA
        # genome_pam = CCT
        # region = GTA + CCT + TCGATCTAGCGATCATCAGA + GATA
        assert region == "GTACCTTCGATCTAGCGATCATCAGAGATA"
        assert offset == 6


class TestExtractCodingFromRegion:
    """Tests for extract_coding_sequence_from_region."""

    def test_full_region_covered(self):
        chrom = "chr1"
        guide_start = 4
        guide_end = 9
        cds_info = [(0, 20, "+", "tx1")]
        region = "ABCDEFGHIJKLMNO"
        guide_offset = 4
        seq, codon_start, codon_end, strand = extract_coding_sequence_from_region(
            chrom, guide_start, guide_end, cds_info, region, guide_offset
        )
        # offset=4→needed=1→extract_start=3; guide_end-extract_start=6→needed=3→extract_end=12
        # region_slice = [3:12]
        assert seq == "DEFGHIJKL"
        assert codon_start == 1
        assert codon_end == 3
        assert strand == "+"

    def test_region_too_short(self):
        chrom = "chr1"
        guide_start = 8
        guide_end = 12
        cds_info = [(0, 20, "+", "tx1")]
        region = "ABCDE"
        guide_offset = 3
        result = extract_coding_sequence_from_region(
            chrom, guide_start, guide_end, cds_info, region, guide_offset
        )
        # needed_upstream=2, needed_downstream=3 → extract_end=15,
        # region_extract_end=10 > len 5 → None
        assert result[0] is None

    def test_in_phase(self):
        chrom = "chr1"
        guide_start = 6
        guide_end = 10
        cds_info = [(0, 20, "+", "tx1")]
        region = "ABCDEFGHIJ"
        guide_offset = 2
        seq, codon_start, codon_end, strand = extract_coding_sequence_from_region(
            chrom, guide_start, guide_end, cds_info, region, guide_offset
        )
        # needed_upstream=0; guide_end-extract_start=4→needed=2→extract_end=12
        # region_slice = [2:8]
        assert seq == "CDEFGH"
        assert codon_start == 0
        assert codon_end == 2
