"""Tests for the sequence property helpers used by the guidefinder."""

import pytest

from beguided.guidefinder.sequence_props import (
    gc_content,
    has_homopolymer,
    longest_homopolymer_run,
)


class TestGCContent:
    """Tests for ``gc_content``."""

    def test_all_at_returns_zero(self):
        assert gc_content("AAAA") == 0.0

    def test_all_gc_returns_hundred(self):
        assert gc_content("GGCC") == 100.0

    def test_mixed_string(self):
        # 4 GC out of 10 = 40.0
        assert gc_content("AAGGCCTTAA") == 40.0

    def test_rounded_to_one_decimal(self):
        # 4 / 13 = 30.769... -> 30.8
        result = gc_content("GCGCAAAAATTTT")
        assert result == 30.8

    def test_only_uppercase_cg_counted(self):
        # Document the behaviour: only upper-case CG are counted, matching
        # the convention in crisprater.py. Guidefinder input always arrives
        # as upper-case (FASTA), so this is not a problem in practice.
        assert gc_content("aaggcc") == 0.0
        # AAGGCC -> 4/6 GC = 66.7
        assert gc_content("AAGGCC") == 66.7

    def test_empty_string_returns_zero(self):
        assert gc_content("") == 0.0


class TestLongestHomopolymerRun:
    """Tests for ``longest_homopolymer_run``."""

    def test_no_runs_returns_one(self):
        assert longest_homopolymer_run("AGCTAGCTAG") == 1

    def test_three_mer_returns_three(self):
        assert longest_homopolymer_run("AAACCCGGGTTT") == 3

    def test_four_mer_at_start(self):
        assert longest_homopolymer_run("AAAAGCTAGCT") == 4

    def test_four_mer_in_middle(self):
        assert longest_homopolymer_run("GCTAAAAGCTA") == 4

    def test_five_mer_at_end(self):
        assert longest_homopolymer_run("AGCTAGCTAAAAA") == 5

    def test_uniform_string(self):
        assert longest_homopolymer_run("GGGGGGGG") == 8

    def test_empty_string_returns_zero(self):
        assert longest_homopolymer_run("") == 0


class TestHasHomopolymer:
    """Tests for ``has_homopolymer``."""

    def test_three_mer_is_false(self):
        assert has_homopolymer("AAATTTGGGCCC") is False

    def test_four_mer_is_true(self):
        assert has_homopolymer("AAAATTTGGGCC") is True

    def test_polyt_is_true(self):
        assert has_homopolymer("GCTTTTTTAGCT") is True

    def test_polyc_is_true(self):
        assert has_homopolymer("GCTCCCCCCAGCT") is True

    def test_custom_min_run(self):
        # 5-mer should be flagged at min_run=4 but not at min_run=6
        assert has_homopolymer("AAAAA") is True
        assert has_homopolymer("AAAAA", min_run=6) is False

    def test_default_min_run_is_four(self):
        # Sanity check that the default cutoff really is 4.
        assert has_homopolymer("AAA") is False
        assert has_homopolymer("AAAA") is True
