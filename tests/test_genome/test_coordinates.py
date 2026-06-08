"""Tests for coordinate calculations."""

import pytest

from beguided.genome.coordinates import (
    calculate_search_range,
    calculate_editing_position,
    guide_coordinates,
)


class TestCalculateSearchRange:
    """Tests for calculate_search_range."""

    def test_antisense_upstream_basic(self):
        start, end = calculate_search_range(
            position=30, min_pos=3, max_pos=8,
            strand="-", guide_length=20, pam_length=3
        )
        assert start == 9
        assert end == 17

    def test_sense_downstream_basic(self):
        start, end = calculate_search_range(
            position=75, min_pos=3, max_pos=8,
            strand="+", guide_length=20, pam_length=3
        )
        assert start == 87
        assert end == 95

    def test_antisense_near_start(self):
        start, end = calculate_search_range(
            position=21, min_pos=3, max_pos=8,
            strand="-", guide_length=20, pam_length=3
        )
        assert start == 0
        assert end > start

    def test_sense_with_longer_pam(self):
        start, end = calculate_search_range(
            position=75, min_pos=5, max_pos=13,
            strand="+", guide_length=20, pam_length=6
        )
        assert start > 0
        assert end > start

    def test_narrow_editing_window(self):
        start, end = calculate_search_range(
            position=75, min_pos=5, max_pos=6,
            strand="+", guide_length=20, pam_length=3
        )
        assert end - start == 4


class TestCalculateEditingPosition:
    """Tests for calculate_editing_position."""

    def test_antisense_basic(self):
        editing_pos = calculate_editing_position(
            position=30, guide_start=14, guide_length=20, strand="-"
        )
        assert editing_pos == 5

    def test_sense_basic(self):
        editing_pos = calculate_editing_position(
            position=75, guide_start=68, guide_length=20, strand="+"
        )
        assert editing_pos == 7

    def test_antisense_edge_position_1(self):
        editing_pos = calculate_editing_position(
            position=23, guide_start=3, guide_length=20, strand="-"
        )
        assert editing_pos == 1

    def test_sense_edge_position_19(self):
        editing_pos = calculate_editing_position(
            position=87, guide_start=68, guide_length=20, strand="+"
        )
        assert editing_pos == 19

    def test_variable_guide_length(self):
        editing_pos = calculate_editing_position(
            position=30, guide_start=10, guide_length=21, strand="-"
        )
        assert 1 <= editing_pos <= 21


class TestGuideCoordinates:
    """Tests for guide_coordinates."""

    def test_sense(self):
        start, end = guide_coordinates(75, 7, "+", 20)
        assert start == 68
        assert end == 88

    def test_antisense(self):
        start, end = guide_coordinates(30, 5, "-", 20)
        assert start == 14
        assert end == 34
