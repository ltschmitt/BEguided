"""Tests for preset loading and PAM regex conversion."""

import re

import pytest

from beguided.editing.presets import (
    pam_to_regex,
    load_presets,
    get_available_preset_names,
    get_preset_names_by_type,
    get_supported_editor_keys,
    filter_compatible_targets,
)


class TestPAMToRegex:
    """Tests for pam_to_regex function."""

    @pytest.mark.parametrize(
        "pam_pattern,expected_substring",
        [("NGG", "[ACGT]GG"), ("NAG", "[ACGT]AG"), ("NGA", "[ACGT]GA")],
    )
    def test_single_n_wildcard(self, pam_pattern, expected_substring):
        result = pam_to_regex(pam_pattern)
        assert expected_substring in result
        assert result.startswith("(?=(")
        assert result.endswith("))")

    @pytest.mark.parametrize(
        "pam_pattern,expected",
        [("NNGRRT", "[ACGT][ACGT]G[AG][AG]T"), ("NNNRRT", "[ACGT][ACGT][ACGT][AG][AG]T")],
    )
    def test_multiple_wildcards(self, pam_pattern, expected):
        result = pam_to_regex(pam_pattern)
        assert expected in result

    @pytest.mark.parametrize(
        "pam_pattern,expected",
        [
            ("R", "[AG]"),
            ("Y", "[CT]"),
            ("W", "[AT]"),
            ("S", "[CG]"),
            ("M", "[AC]"),
            ("K", "[GT]"),
            ("B", "[CGT]"),
            ("D", "[AGT]"),
            ("H", "[ACT]"),
            ("V", "[ACG]"),
        ],
    )
    def test_all_iupac_wildcards(self, pam_pattern, expected):
        result = pam_to_regex(pam_pattern)
        assert expected in result

    def test_mixed_pattern(self):
        result = pam_to_regex("NGG")
        assert result == "(?=([ACGT]GG))"

    def test_lowercase_input(self):
        result = pam_to_regex("ngg")
        assert "[ACGT]GG" in result


class TestLoadPresets:
    """Tests for load_presets function."""

    def test_load_valid_presets(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        assert ("A", "G") in presets
        assert ("T", "C") in presets
        assert "ABE8e" in presets[("A", "G")]
        assert "ABE8e" in presets[("T", "C")]

    def test_sense_antisense_pair_generation(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        sense = presets[("A", "G")]["ABE8e"]
        assert sense.strand == "+"
        assert sense.pam_pattern == "NGG"
        assert not sense.reverse_complement

        antisense = presets[("T", "C")]["ABE8e"]
        assert antisense.strand == "-"
        assert antisense.pam_pattern == "CCN"
        assert antisense.reverse_complement

    def test_preset_parameters_loaded(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))
        abe8e = presets[("A", "G")]["ABE8e"]

        assert abe8e.high_edit_start == 3
        assert abe8e.high_edit_end == 8
        assert abe8e.guide_length == 20
        assert abe8e.pam_length == 3
        assert abe8e.name == "ABE8e"

    def test_duplicate_preset_names(self, tmp_path, duplicate_preset_tsv):
        preset_file = tmp_path / "duplicate_presets.tsv"
        preset_file.write_text(duplicate_preset_tsv)
        with pytest.raises(ValueError) as exc_info:
            load_presets(str(preset_file))
        assert "Duplicate preset names" in str(exc_info.value)
        assert "ABE8e" in str(exc_info.value)

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            load_presets("nonexistent_file.tsv")

    def test_multiple_presets(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        assert "ABE8e" in presets[("A", "G")]
        assert "saABE8e" in presets[("A", "G")]
        assert "ABE9" in presets[("A", "G")]


class TestSupportedEditorKeys:
    """Tests for get_supported_editor_keys."""

    def test_returns_sense_and_antisense_keys(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        keys = get_supported_editor_keys(presets, ["ABE8e"])
        assert ("A", "G") in keys
        assert ("T", "C") in keys
        assert ("C", "T") not in keys

    def test_subset_of_presets(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        keys = get_supported_editor_keys(presets, ["saABE8e"])
        assert ("A", "G") in keys
        assert ("T", "C") in keys
        assert ("C", "T") not in keys


class TestFilterCompatibleTargets:
    """Tests for filter_compatible_targets."""

    def test_filters_incompatible_targets(self, tmp_path, valid_preset_tsv):
        import pandas as pd
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        targets = pd.DataFrame({
            "chromosome": ["chr1", "chr1", "chr1"],
            "position": [10, 20, 30],
            "target_base": ["A", "T", "C"],
            "desired_base": ["G", "C", "T"],
        })

        compatible, incompatible = filter_compatible_targets(targets, presets, ["ABE8e"])
        assert len(compatible) == 2
        assert len(incompatible) == 1
        assert list(compatible["target_base"]) == ["A", "T"]
        assert list(incompatible["target_base"]) == ["C"]

    def test_all_targets_compatible(self, tmp_path, valid_preset_tsv):
        import pandas as pd
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        targets = pd.DataFrame({
            "chromosome": ["chr1", "chr1"],
            "position": [10, 20],
            "target_base": ["A", "T"],
            "desired_base": ["G", "C"],
        })

        compatible, incompatible = filter_compatible_targets(targets, presets, ["ABE8e"])
        assert len(compatible) == 2
        assert len(incompatible) == 0

    def test_empty_targets(self, tmp_path, valid_preset_tsv):
        import pandas as pd
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        targets = pd.DataFrame(columns=["chromosome", "position", "target_base", "desired_base"])
        compatible, incompatible = filter_compatible_targets(targets, presets, ["ABE8e"])
        assert len(compatible) == 0
        assert len(incompatible) == 0


class TestGetPresetNamesByType:
    """Tests for get_preset_names_by_type."""

    def test_classifies_abe_and_cbe(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        type_map = get_preset_names_by_type(presets)
        assert "abe" in type_map
        assert "cbe" in type_map
        assert "ABE8e" in type_map["abe"]
        assert "saABE8e" in type_map["abe"]
        assert "ABE9" in type_map["abe"]
        # valid_preset_tsv has no CBE presets, so cbe list should be empty
        assert type_map["cbe"] == []

    def test_no_duplicate_names(self, tmp_path, valid_preset_tsv):
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        type_map = get_preset_names_by_type(presets)
        # Antisense variants (T->C) should NOT appear in the ABE list
        abe_set = set(type_map["abe"])
        assert len(type_map["abe"]) == len(abe_set)


class TestResolvePresetNames:
    """Tests for _resolve_preset_names in pipeline (via presets on presets)."""

    def test_all_abes_expansion(self, tmp_path, valid_preset_tsv):
        from beguided.pipeline import _resolve_preset_names
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        names = _resolve_preset_names(presets, "all ABEs")
        assert "ABE8e" in names
        assert "saABE8e" in names
        assert "ABE9" in names

    def test_all_cbes_empty_when_no_cbe_presets(self, tmp_path, valid_preset_tsv):
        from beguided.pipeline import _resolve_preset_names
        preset_file = tmp_path / "valid_presets.tsv"
        preset_file.write_text(valid_preset_tsv)
        presets = load_presets(str(preset_file))

        names = _resolve_preset_names(presets, "all CBEs")
        assert names == []
