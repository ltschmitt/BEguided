"""Tests for mutation and bystander logic."""

import pytest

from beguided.editing.mutations import (
    find_bystander_edits,
    simulate_editing,
    mark_mutations,
)


class TestFindBystanderEdits:
    """Tests for find_bystander_edits."""

    def test_no_bystanders(self):
        result = find_bystander_edits(
            "GGGGGGGGGGGGGGGGGGGG", "A", 5, 2, 13, "+"
        )
        assert result == ""

    def test_with_bystanders(self):
        result = find_bystander_edits(
            "AAATAAAAAATAAAAAAAAA", "A", 5, 2, 13, "+"
        )
        assert result == "2,3,6,7,8,9,10,12,13"

    def test_antisense(self):
        result = find_bystander_edits(
            "AATAATAATAATAATAATAA", "T", 5, 2, 13, "-"
        )
        # Antisense: complement of T is A, so look for A
        assert result != ""


class TestSimulateEditing:
    """Tests for simulate_editing."""

    def test_single_edit(self):
        result = simulate_editing("AAATAAAAAATAAAAAAAAA", 0, [5], "G", "+")
        assert result == "AAATGAAAAATAAAAAAAAA"

    def test_multiple_edits(self):
        result = simulate_editing(
            "AAATAAAAAATAAAAAAAAA", 0, [5, 11], "G", "+"
        )
        assert result == "AAATGAAAAAGAAAAAAAAA"


class TestMarkMutations:
    """Tests for mark_mutations."""

    def test_target_only(self):
        result = mark_mutations("ATCGATCG", 3, [], "+", "+", 0)
        assert result == "AT[C]GATCG"

    def test_with_bystander(self):
        result = mark_mutations("ATCGATCG", 3, [5], "+", "+", 0)
        assert result == "AT[C]G(A)TCG"
