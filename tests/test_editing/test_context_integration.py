"""Integration tests for the trinuc_edit_score field on GuideResult.

These tests exercise the wiring between ``beguided.editing.context`` and
``beguided.guidefinder.core.process_target``, ensuring the score is
populated for ABE guides (sense and antisense), left at ``None`` for CBE
guides, and survives through ``run_genome_wide`` + ``results_to_dataframe``.
"""

import pytest

from beguided.editing.context import load_edit_context_table
from beguided.editing.presets import load_base_editor_profiles
from beguided.guidefinder.core import process_target, run_genome_wide, results_to_dataframe


# ---------------------------------------------------------------------------
# Synthetic genomes designed so the target sits in a known X-A-Y context.
# ---------------------------------------------------------------------------

# Sense ABE target: position 10 (1-based) is A, with C and C neighbours → CAC.
# A downstream NGG PAM at indices 25-27 ("AGG") so a + strand ABE8e guide
# with editing_position=5 picks up the CAC context.
#
# Layout (0-based indices):
#   0-7: filler "ACGTACGT"
#   8:   'C' (left neighbour of target)
#   9:   'A' (target at 1-based position 10)
#   10:  'C' (right neighbour of target)
#   11-24: filler "GTAGCCGTTGGACG" (14 chars)
#   25-27: "AGG" (PAM)
#   28-31: tail "TACG"
#
# Trinucleotide slice genome_seq[8:11] = "CAC" -> table[CAC] = 1.0.
# guide_start = 25 - 20 = 5; editing_pos = 9 - 5 + 1 = 5 (in [3, 8]). ✓
_SENSE_GENOME = "ACGTACGT" + "C" + "A" + "C" + "GTAGCCGTTGGACG" + "AGG" + "TACG"

# Antisense ABE target: position 10 (1-based) is T (genome base).
# Slice genome_seq[8:11] = "GTG"; reverse-complement = "CAC" -> 1.0.
# A CCN PAM (reverse complement of NGG) sits to the LEFT of the guide.
_ANTISENSE_GENOME = (
    "ACGTTGGCCG"     # positions 1-10 (index 0-9); index 8=G, 9=T(target), 10=G
    "AGTAGCCGTTGG"   # filler
    "CCN"            # not a real base; replaced below
)
# The above is just a sketch; we build the real antisense genome below so the
# PAM is correct. The key invariant is genome_seq[8:11] == "GTG".


def _make_abe_profiles(tmp_path, editing="ABE"):
    """Write a minimal ABE8e profiles TSV and return the loaded profiles dict."""
    tsv = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        f"ABE8e\t{editing}\t2\t13\t3\t8\tNGG\tSpCas9\ttest\n"
    )
    p = tmp_path / "profiles.tsv"
    p.write_text(tsv)
    return load_base_editor_profiles(str(p))


def _make_cbe_profiles(tmp_path):
    """Write a minimal CBE profiles TSV and return the loaded profiles dict."""
    tsv = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        f"CBE6d\tCBE\t2\t10\t6\t6\tNGG\tSpCas9\ttest\n"
    )
    p = tmp_path / "profiles.tsv"
    p.write_text(tsv)
    return load_base_editor_profiles(str(p))


# ---------------------------------------------------------------------------
# process_target integration
# ---------------------------------------------------------------------------


class TestProcessTargetTrinucScore:
    """Verify process_target populates trinuc_edit_score correctly."""

    def test_sense_abe_gets_cac_score(self, tmp_path):
        """A + strand ABE guide in a CAC context scores 1.0."""
        profiles = _make_abe_profiles(tmp_path)
        config = profiles[("A", "G")]["ABE8e"]
        ect = load_edit_context_table()

        results = process_target(
            _SENSE_GENOME, "chr1", 10, "A", "G",
            "A", config, edit_context_table=ect,
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0
        for r in valid:
            assert r.trinuc_edit_score == 1.0, (
                f"Expected CAC=1.0 for sense ABE guide at chr1:10, "
                f"got {r.trinuc_edit_score}"
            )

    def test_antistrand_abe_gets_reverse_complement_score(self, tmp_path):
        """A - strand ABE guide (T->C) on a GTG genome slice scores the same
        as the sense CAC context (1.0), because the slice is reverse-
        complemented before lookup."""
        profiles = _make_abe_profiles(tmp_path)
        config = profiles[("T", "C")]["ABE8e"]
        ect = load_edit_context_table()

        # Build a genome where position 10 (1-based) is T, with G neighbours
        # on each side (indices 8, 9, 10 = G, T, G). The reverse complement
        # of GTG is CAC -> 1.0. We need a CCN PAM (antisense of NGG) to the
        # left of the guide so a - strand guide can be generated.
        #
        # For position=10, editing_position=5 on the - strand:
        #   guide_end (0-based) = position + (editing_position - 1) - 1 = 10 + 4 - 1 = 13
        #   guide_start (0-based) = guide_end - guide_length = 13 - 20 = -7  (OOB)
        # So we need editing_position large enough that guide_start >= 0.
        # With editing_position=13: guide_end = 10 + 12 - 1 = 21, guide_start = 1.
        # But 13 is outside ABE8e's high window (3-8). Let's shift the target
        # further right instead.
        #
        # Use position=30 so there's room for a 20-nt guide to its left.
        # Trinucleotide at position 30: we want genome_seq[28:31] == "GTG".
        # CCN PAM at indices 5-7 (positions 6-8) so it's upstream of guide_start.
        # guide_end = 30 + (editing_pos - 1) - 1 (0-based). For editing_pos=5:
        #   guide_end = 30 + 4 - 1 = 33 (0-based), guide_start = 13 (0-based).
        # PAM for - strand is at [guide_start - pam_length : guide_start] = [10:13].
        # Hmm that overlaps the guide. Let me re-read the coordinates.
        #
        # Actually for - strand in calculate_search_range:
        #   base = position - 1 = 29
        #   offset = guide_length + pam_length = 23
        #   start = base - (offset - min_pos) = 29 - (23 - 3) = 9
        #   end = base - (offset - max_pos) + pam_length = 29 - (23 - 8) + 3 = 17
        # So the PAM search range is [9, 17] in 0-based coordinates.
        # For - strand: guide_start = pam_pos + pam_length
        # If pam_pos = 9, guide_start = 12, guide_end = 32.
        # editing_pos = guide_length - ((position - 1) - guide_start) = 20 - (29 - 12) = 3
        # That's in the high window [3, 8]. ✓
        #
        # So we need:
        #   - genome_seq[28:31] = "GTG" (position 30 = T target, neighbours G and G)
        #   - A "CCN" pattern (reverse complement of NGG) in the PAM search
        #     range [9, 17). CCN on the + strand = "CC[A/C/G/T]".
        genome = (
            "ACGAAGACTACGAGTAGCCG"   # indices 0-19 (positions 1-20)
            "TTGG"                   # indices 20-23
            "CC"                     # indices 24-25 — start of a CCN PAM
            "A"                      # index 26 — N of CCN
            "GG"                     # indices 27-28 — wait, need GTG at 28-30
        )
        # This is getting fiddly. Let me just construct it programmatically.
        genome = list("A" * 60)
        # Place GTG at indices 28, 29, 30 (positions 29, 30, 31)
        genome[28] = "G"
        genome[29] = "T"  # position 30 = target T
        genome[30] = "G"
        # Place a CCN PAM at indices 9, 10, 11 (positions 10, 11, 12)
        genome[9] = "C"
        genome[10] = "C"
        genome[11] = "A"  # N
        genome = "".join(genome)

        results = process_target(
            genome, "chr1", 30, "T", "C",
            "T", config, edit_context_table=ect,
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0, "Expected at least one valid - strand ABE guide"
        for r in valid:
            assert r.trinuc_edit_score == 1.0, (
                f"Expected reverse-complement CAC=1.0 for antisense ABE guide, "
                f"got {r.trinuc_edit_score}"
            )

    def test_cbe_guide_has_none_trinuc_score(self, tmp_path):
        """CBE guides must have trinuc_edit_score=None (table is ABE-only)."""
        profiles = _make_cbe_profiles(tmp_path)
        config = profiles[("C", "T")]["CBE6d"]
        ect = load_edit_context_table()

        # Use the sense genome but target C->T instead of A->G.
        # Position 11 (1-based) is C in _SENSE_GENOME.
        results = process_target(
            _SENSE_GENOME, "chr1", 11, "C", "T",
            "C", config, edit_context_table=ect,
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        if valid:
            for r in valid:
                assert r.trinuc_edit_score is None, (
                    f"CBE guide should have None trinuc_edit_score, "
                    f"got {r.trinuc_edit_score}"
                )

    def test_no_table_leaves_score_none(self, tmp_path):
        """When edit_context_table=None, every guide's score is None."""
        profiles = _make_abe_profiles(tmp_path)
        config = profiles[("A", "G")]["ABE8e"]

        results = process_target(
            _SENSE_GENOME, "chr1", 10, "A", "G",
            "A", config, edit_context_table=None,
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0
        for r in valid:
            assert r.trinuc_edit_score is None

    def test_no_guide_possible_has_none_score(self, tmp_path):
        """Invalid guides (**No guide possible**) never get a score."""
        profiles = _make_abe_profiles(tmp_path)
        config = profiles[("A", "G")]["ABE8e"]
        ect = load_edit_context_table()

        # Position 5 in chr2 of the sample genome has no valid PAM.
        from tests.conftest import sample_genome_dict  # noqa: F401
        # Use a tiny genome with no NGG PAM anywhere near the target.
        tiny_genome = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        results = process_target(
            tiny_genome, "chr1", 20, "A", "G",
            "A", config, edit_context_table=ect,
        )
        assert results[0].guide_rna_sequence == "**No guide possible**"
        assert results[0].trinuc_edit_score is None

    def test_score_survives_results_to_dataframe(self, tmp_path):
        """The trinuc_edit_score column appears in the DataFrame."""
        profiles = _make_abe_profiles(tmp_path)
        config = profiles[("A", "G")]["ABE8e"]
        ect = load_edit_context_table()

        results = process_target(
            _SENSE_GENOME, "chr1", 10, "A", "G",
            "A", config, edit_context_table=ect,
        )
        df = results_to_dataframe(results)
        assert "trinuc_edit_score" in df.columns
        valid = df[~df["guide_rna_sequence"].str.startswith("**")]
        assert len(valid) > 0
        assert valid["trinuc_edit_score"].notna().all()
        assert (valid["trinuc_edit_score"] == 1.0).all()


# ---------------------------------------------------------------------------
# run_genome_wide integration
# ---------------------------------------------------------------------------


class TestRunGenomeWideTrinucScore:
    """Verify run_genome_wide forwards edit_context_table to process_target."""

    def test_run_genome_wide_populates_score(self, tmp_path):
        profiles = _make_abe_profiles(tmp_path)
        ect = load_edit_context_table()

        import pandas as pd
        targets = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [10],
            "target_base": ["A"],
            "desired_base": ["G"],
        })
        results = run_genome_wide(
            {"chr1": _SENSE_GENOME}, targets, profiles, ["ABE8e"],
            edit_context_table=ect,
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0
        for r in valid:
            assert r.trinuc_edit_score == 1.0

    def test_run_genome_wide_without_table_leaves_none(self, tmp_path):
        profiles = _make_abe_profiles(tmp_path)

        import pandas as pd
        targets = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [10],
            "target_base": ["A"],
            "desired_base": ["G"],
        })
        results = run_genome_wide(
            {"chr1": _SENSE_GENOME}, targets, profiles, ["ABE8e"],
            # edit_context_table defaults to None
        )
        valid = [r for r in results if not r.guide_rna_sequence.startswith("**")]
        assert len(valid) > 0
        for r in valid:
            assert r.trinuc_edit_score is None
