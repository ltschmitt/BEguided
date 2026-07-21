"""Tests for the trinucleotide edit-context lookup (ABE only)."""

from pathlib import Path

import pytest

from beguided.editing.context import (
    EditContextTable,
    load_edit_context_table,
    lookup,
    lookup_from_genome,
    is_abe,
)


# ---------------------------------------------------------------------------
# load_edit_context_table
# ---------------------------------------------------------------------------


class TestLoadEditContextTable:
    """Tests for load_edit_context_table (default bundled file + validation)."""

    def test_bundled_loads_16_keys(self):
        ect = load_edit_context_table()
        assert isinstance(ect, EditContextTable)
        assert len(ect.table) == 16
        # Max value is CAC = 1.00 per Marquart et al. 2021 (normalised).
        assert ect.table["CAC"] == 1.0
        assert ect.table["AAA"] == pytest.approx(0.02)

    def test_bundled_path_resolves(self):
        """The default path should be a real file under beguided/editing/data."""
        from beguided.editing.context import _default_table_path
        p = _default_table_path()
        assert p.exists()
        assert p.name == "edit_context_norm_proportion.tsv"

    def test_antisense_table_is_reverse_complement_of_sense(self):
        from beguided.genome.sequence import reverse_complement
        ect = load_edit_context_table()
        assert len(ect.antisense_table) == 16
        for k, v in ect.table.items():
            rc = reverse_complement(k)
            assert ect.antisense_table[rc] == v

    def test_source_label_documents_provenance(self):
        ect = load_edit_context_table()
        assert "Marquart" in ect.source_label
        assert "ABE8e" in ect.source_label
        assert ect.approximate_for_non_ABE8e_ABE is True

    def test_custom_path(self, tmp_path):
        p = tmp_path / "custom.tsv"
        p.write_text("trinuc\tproportion\n" + "".join(
            f"{k}\t{v}\n" for k, v in [
                ("AAA", 0.02), ("AAC", 0.11), ("AAG", 0.02), ("AAT", 0.07),
                ("CAA", 0.35), ("CAC", 1.00), ("CAG", 0.89), ("CAT", 0.70),
                ("GAA", 0.17), ("GAC", 0.84), ("GAG", 0.68), ("GAT", 0.60),
                ("TAA", 0.38), ("TAC", 0.91), ("TAG", 0.60), ("TAT", 0.79),
            ]
        ))
        ect = load_edit_context_table(str(p))
        assert ect.table["CAC"] == 1.0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_edit_context_table(str(tmp_path / "nope.tsv"))

    def test_wrong_columns_raise(self, tmp_path):
        p = tmp_path / "bad.tsv"
        p.write_text("foo\tbar\nAAA\t0.5\n")
        with pytest.raises(ValueError, match="trinuc"):
            load_edit_context_table(str(p))

    def test_invalid_key_pattern_raises(self, tmp_path):
        # Middle base must be A — swap one row to "ACA" and expect failure.
        p = tmp_path / "bad.tsv"
        p.write_text("trinuc\tproportion\n" + "".join(
            f"{k}\t{v}\n" for k, v in [
                ("ACA", 0.5),  # bad — middle is C
                ("AAC", 0.11), ("AAG", 0.02), ("AAT", 0.07),
                ("CAA", 0.35), ("CAC", 1.00), ("CAG", 0.89), ("CAT", 0.70),
                ("GAA", 0.17), ("GAC", 0.84), ("GAG", 0.68), ("GAT", 0.60),
                ("TAA", 0.38), ("TAC", 0.91), ("TAG", 0.60), ("TAT", 0.79),
            ]
        ))
        with pytest.raises(ValueError, match="invalid sense keys"):
            load_edit_context_table(str(p))

    def test_missing_key_raises(self, tmp_path):
        # Drop TAT to trigger the "missing sense keys" check.
        p = tmp_path / "incomplete.tsv"
        p.write_text("trinuc\tproportion\n" + "".join(
            f"{k}\t{v}\n" for k, v in [
                ("AAA", 0.02), ("AAC", 0.11), ("AAG", 0.02), ("AAT", 0.07),
                ("CAA", 0.35), ("CAC", 1.00), ("CAG", 0.89), ("CAT", 0.70),
                ("GAA", 0.17), ("GAC", 0.84), ("GAG", 0.68), ("GAT", 0.60),
                ("TAA", 0.38), ("TAC", 0.91), ("TAG", 0.60),
            ]
        ))
        with pytest.raises(ValueError, match="missing sense keys"):
            load_edit_context_table(str(p))


# ---------------------------------------------------------------------------
# is_abe
# ---------------------------------------------------------------------------


class TestIsABE:
    """Tests for is_abe — ABE detection from (from, to) bases."""

    @pytest.mark.parametrize(
        "fb,tb,expected",
        [
            ("A", "G", True),    # sense ABE
            ("T", "C", True),    # antisense ABE
            ("a", "g", True),    # case-insensitive
            ("t", "c", True),
            ("C", "T", False),   # sense CBE
            ("G", "A", False),   # antisense CBE
            ("A", "T", False),   # not a base edit
            ("C", "G", False),
        ],
    )
    def test_detection(self, fb, tb, expected):
        assert is_abe(fb, tb) is expected


# ---------------------------------------------------------------------------
# lookup
# ---------------------------------------------------------------------------


class TestLookup:
    """Tests for the direct sense-key lookup."""

    def test_known_keys(self):
        ect = load_edit_context_table()
        assert lookup(ect, "CAC") == 1.0
        assert lookup(ect, "AAA") == pytest.approx(0.02)
        assert lookup(ect, "TAC") == pytest.approx(0.91)

    def test_case_insensitive(self):
        ect = load_edit_context_table()
        assert lookup(ect, "cac") == 1.0
        assert lookup(ect, "Cac") == 1.0

    def test_unknown_key_returns_none(self):
        ect = load_edit_context_table()
        # The table only has XAY (middle = A) keys.
        assert lookup(ect, "ACA") is None
        assert lookup(ect, "NNN") is None

    def test_empty_returns_none(self):
        ect = load_edit_context_table()
        assert lookup(ect, "") is None


# ---------------------------------------------------------------------------
# lookup_from_genome
# ---------------------------------------------------------------------------


class TestLookupFromGenome:
    """Tests for the slice + reverse-complement + lookup helper."""

    def test_sense_strand_cac(self):
        ect = load_edit_context_table()
        # genome GTCACGT, target A at 1-based pos 4 -> slice[2:5] = CAC
        assert lookup_from_genome(ect, "GTCACGT", 4, "+", "A") == 1.0

    def test_sense_strand_aaa(self):
        ect = load_edit_context_table()
        # genome TAAAT, target A at 1-based pos 3 -> slice[1:4] = AAA
        assert lookup_from_genome(ect, "TAAAT", 3, "+", "A") == pytest.approx(0.02)

    def test_antistrand_returns_same_value_as_sense(self):
        """A `-`-strand guide editing T->C reads the reverse-complement of
        the genome slice, so an mRNA-sense CAC context scores the same on
        both strands."""
        ect = load_edit_context_table()
        # Sense: GTCACGT, pos 4, target A -> CAC -> 1.0
        sense = lookup_from_genome(ect, "GTCACGT", 4, "+", "A")
        # Antisense: reverse-complement of the same locus = ACGTGAC,
        # target T at 1-based pos 4 -> slice[2:5] = GTG -> revcomp = CAC -> 1.0
        antisense = lookup_from_genome(ect, "ACGTGAC", 4, "-", "T")
        assert sense == antisense == 1.0

    def test_target_base_used_as_middle(self):
        """When the genome base at the target position is the desired_base
        (e.g. already-edited 'G' at an A>G target), the lookup still uses
        target_base='A' as the middle of the trinucleotide."""
        ect = load_edit_context_table()
        # Genome has G at position 4 (already edited), but target_base='A'.
        # Neighbours at indices 2 and 4 are C and C, so context = CAC -> 1.0.
        assert lookup_from_genome(ect, "GTCGCGT", 4, "+", "A") == 1.0

    def test_oob_at_chromosome_start(self):
        ect = load_edit_context_table()
        # position=1 means we need slice[-1:2] -> OOB
        assert lookup_from_genome(ect, "AACCGGTT", 1, "+", "A") is None

    def test_oob_at_chromosome_end(self):
        ect = load_edit_context_table()
        # position=8 means slice[6:9] on an 8-char string -> OOB
        assert lookup_from_genome(ect, "AACCGGTT", 8, "+", "A") is None

    def test_non_acgt_neighbor_returns_none(self):
        ect = load_edit_context_table()
        # N as a neighbour (right of target at pos 2) -> non-ACGT -> None.
        # Genome "GANAA", position=2: left=G, middle=A(target_base), right=N.
        assert lookup_from_genome(ect, "GANAA", 2, "+", "A") is None

    def test_lowercase_input_is_normalised(self):
        ect = load_edit_context_table()
        assert lookup_from_genome(ect, "gtcacgt", 4, "+", "a") == 1.0


# ---------------------------------------------------------------------------
# EditContextTable dataclass
# ---------------------------------------------------------------------------


class TestEditContextTable:
    """Direct construction of EditContextTable bypasses file IO."""

    def test_construct_and_post_init_builds_antisense(self):
        from beguided.genome.sequence import reverse_complement
        ect = EditContextTable(table={"CAC": 1.0, "AAA": 0.02})
        assert ect.antisense_table[reverse_complement("CAC")] == 1.0
        assert ect.antisense_table[reverse_complement("AAA")] == 0.02

    def test_default_labels(self):
        ect = EditContextTable(table={"CAC": 1.0})
        assert "Marquart" in ect.source_label
        assert ect.approximate_for_non_ABE8e_ABE is True
