"""Tests for genome I/O and sequence utilities."""

import pytest

from beguided.genome.io import derive_genome_label, load_genome
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

    def test_skips_scaffold_records(self, tmp_path):
        """Records whose header contains 'dna:scaffold' must be dropped."""
        fasta = tmp_path / "mixed.fa"
        fasta.write_text(
            ">1 dna:chromosome chromosome:GRCh38:1:1:10:1 REF\nACGTACGTAC\n"
            ">KI270728.1 dna:scaffold scaffold:GRCh38:KI270728.1:1:8:1 REF\nGGGGGGGG\n"
            ">MT dna:chromosome chromosome:GRCh38:MT:1:4:1 REF\nTTTT\n"
            ">GL000225.1 dna:scaffold scaffold:GRCh38:GL000225.1:1:6:1 REF\nCCCCCC\n"
        )
        genome = load_genome(str(fasta))
        assert set(genome) == {"1", "MT"}
        assert genome["1"] == "ACGTACGTAC"
        assert genome["MT"] == "TTTT"

    def test_skips_mouse_unmasked_scaffolds(self, tmp_path):
        """The EBI mouse header uses 'unmasked:scaffold'; it must be dropped."""
        fasta = tmp_path / "mouse.fa"
        fasta.write_text(
            ">1 unmasked:chromosome GRCm39:1:1:10:1 REF\nACGTACGTAC\n"
            ">JH584299.1 unmasked:scaffold GRCm39:JH584299.1:1:8:1 REF\nGGGGGGGG\n"
            ">MT unmasked:chromosome GRCm39:MT:1:4:1 REF\nTTTT\n"
        )
        genome = load_genome(str(fasta))
        assert set(genome) == {"1", "MT"}
        assert genome["1"] == "ACGTACGTAC"

    def test_can_disable_scaffold_filter(self, tmp_path):
        """Passing skip_header_tokens=() restores the old all-records behaviour."""
        fasta = tmp_path / "mixed.fa"
        fasta.write_text(
            ">1 dna:chromosome chromosome:GRCh38:1:1:10:1 REF\nACGTACGTAC\n"
            ">KI270728.1 dna:scaffold scaffold:GRCh38:KI270728.1:1:8:1 REF\nGGGGGGGG\n"
        )
        genome = load_genome(str(fasta), skip_header_tokens=())
        assert "KI270728.1" in genome


class TestDeriveGenomeLabel:
    """Tests for the FASTA-filename label extractor."""

    def test_ensembl_no_version(self):
        assert (
            derive_genome_label("/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
            == "GRCh38"
        )

    def test_ensembl_with_version(self):
        assert (
            derive_genome_label("/data/Homo_sapiens.GRCh38.112.dna.primary_assembly.fa")
            == "GRCh38.112"
        )

    def test_unknown_fasta_falls_back_to_grch38(self):
        assert derive_genome_label("/tmp/random.fa") == "GRCh38"

    def test_mouse_release_date(self):
        assert (
            derive_genome_label("/data/Mus_musculus.GRCm39.2026_04.fa")
            == "GRCm39.2026_04"
        )


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
