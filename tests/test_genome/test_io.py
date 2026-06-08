"""Tests for genome I/O and sequence utilities."""

import pytest

from beguided.genome.io import load_genome
from beguided.genome.sequence import reverse_complement, insert_mutation


class TestLoadGenome:
    """Tests for load_genome function."""

    def test_read_valid_fasta(self, sample_fasta_file):
        genome = load_genome(sample_fasta_file)
        assert "chr1" in genome
        assert "chr2" in genome
        assert len(genome["chr1"]) > 0
        assert all(c in "ACGTN" for c in genome["chr1"])

    def test_chromosome_id_extraction(self, tmp_path):
        fasta_content = (
            ">chr1 extra info here\nACGTACGT\n>chr2 more info\nTGCATGCA\n"
        )
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(fasta_content)
        genome = load_genome(str(fasta_file))

        assert "chr1" in genome
        assert "chr2" in genome
        assert "extra" not in genome


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_basic(self):
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAA") == "TTT"
        assert reverse_complement("ATGC") == "GCAT"

    def test_with_n(self):
        assert reverse_complement("ATNCG") == "CGNAT"


class TestInsertMutation:
    """Tests for insert_mutation function."""

    def test_sense(self):
        result = insert_mutation("AAAAAAAAAAAAAAAAAAAA", 5, "G", "+", 20)
        assert result == "AAAAGAAAAAAAAAAAAAAA"

    def test_antisense(self):
        result = insert_mutation("AAAAAAAAAAAAAAAAAAAA", 5, "G", "-", 20)
        # guide_length - editing_pos = 20 - 5 = 15
        expected = "AAAAAAAAAAAAAAAGAAAA"
        assert result == expected
