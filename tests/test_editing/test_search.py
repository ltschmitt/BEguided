"""Tests for PAM search operations."""

import re

import pytest

from beguided.editing.presets import pam_to_regex
from beguided.editing.search import find_pams_regex


class TestFindPamsRegex:
    """Tests for find_pams_regex function."""

    def test_simple_non_overlapping(self):
        seq = "AAGGTTCCGGAA"
        pam_regex = re.compile(pam_to_regex("NGG"))
        pams = find_pams_regex(seq, 0, len(seq), pam_regex, 3)
        assert len(pams) >= 1
        assert all(isinstance(pos, int) for pos in pams)

    def test_overlapping_pams(self):
        seq = "TAGGGG"
        pam_regex = re.compile(pam_to_regex("NGG"))
        pams = find_pams_regex(seq, 0, len(seq), pam_regex, 3)
        assert len(pams) == 3
        assert 1 in pams
        assert 2 in pams
        assert 3 in pams

    def test_no_matches(self):
        seq = "AAAAATTTTTCCCCC"
        pam_regex = re.compile(pam_to_regex("NGG"))
        pams = find_pams_regex(seq, 0, len(seq), pam_regex, 3)
        assert len(pams) == 0

    def test_search_range_subset(self):
        seq = "AAGGTTCCGGAA"
        pam_regex = re.compile(pam_to_regex("NGG"))
        pams = find_pams_regex(seq, 0, 6, pam_regex, 3)
        assert all(pos < 6 for pos in pams)

    def test_real_sample_chr1_75(self):
        seq = "ACGAAGACTACGAGTAGCCGTTGGATCGTACCTTCGATCTAGCGATCATCAGAGATATAGGGAACGGTAGCTAGGTAGCTAGCTAGCTAGGGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"  # noqa: E501
        pam_regex = re.compile(pam_to_regex("NGG"))
        pams = find_pams_regex(seq, 87, 93, pam_regex, 3)
        assert 88 in pams
        assert 89 in pams
        assert 90 in pams
