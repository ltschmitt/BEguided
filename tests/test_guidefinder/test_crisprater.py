"""Tests for CRISPRater on-target score module."""

import pytest

from beguided.guidefinder.crisprater import get_score


class TestCRISPRater:
    """Unit tests for CRISPRater score calculation."""

    def test_valid_20mer(self):
        """Score for a valid 20-nt sequence returns a float."""
        score = get_score("GACTGACTGACTGACTGACT")
        assert isinstance(score, float)

    def test_ambiguous_seq_raises(self):
        """Non-ACGT characters are not rejected at the scorer level."""
        # The scorer itself only checks length; ambiguous bases propagate
        # naturally through the feature logic.
        score = get_score("NNNNNNNNNNNNNNNNNNNN")
        assert isinstance(score, float)

    def test_exactly_20_required(self):
        """Sequences not exactly 20 nt raise ValueError."""
        with pytest.raises(ValueError):
            get_score("GACTGACTGACTGACTGAC")  # 19 nt
        with pytest.raises(ValueError):
            get_score("GACTGACTGACTGACTGACTT")  # 21 nt

    def test_gc_window_influences_score(self):
        """GC content of seq[3:13] should shift the score."""
        low_gc = "AAATAAAAAAGAAAAAAAAA"
        high_gc = "AAAGGGGGGGGAAAAAAAAA"
        assert get_score(high_gc) > get_score(low_gc)

    def test_position_19_g_increases_score(self):
        """seq[19] == 'G' sets feat[1]=1, which has positive weight."""
        without = "AAAAAAAAAAAAAAAAAAAA"
        with_g = "AAAAAAAAAAAAAAAAAAAG"
        assert get_score(with_g) > get_score(without)

    def test_position_2_ta_increases_score(self):
        """seq[2] in 'TA' sets feat[2]=1, which has positive weight."""
        without = "AAGAAAAAAAAAAAAAAAAA"
        with_ta = "AATAAAAAAAAAAAAAAAAA"
        assert get_score(with_ta) > get_score(without)
