"""Integration tests for the BEguidedOFF ``fast`` (FM-index) mode wiring.

Covers the BEguided-side seams only (the engine itself is tested in the
BEguidedOFF repo):

* pipeline config resolution (``offtargets.mode``, ``index_path``,
  ``max_mismatches``/``mismatch_tolerance`` alias),
* the fast-mode search-cutoff resolution and its soundness guard,
* an end-to-end ``_run_beguidedoff(mode="fast")`` search against a mini
  FM-index built from the sample-genome fixture.

These tests need the compiled FM-index binding; they skip when BEguidedOFF
(or the binding) is not installed in the active environment.
"""

import logging

import pandas as pd
import pytest

beguidedoff = pytest.importorskip(
    "beguidedoff", reason="BEguidedOFF is not installed"
)
pytest.importorskip(
    "beguidedoff._fmindex", reason="FM-index binding is not built"
)

from beguided.pipeline import (  # noqa: E402
    FAST_MAX_SEARCH_CUTOFF,
    PipelineError,
    _read_new_style_config,
    _resolve_fast_search_cutoff,
    _run_beguidedoff,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mini_index(tmp_path_factory):
    """Build a small FM-index from the sample-genome fixture once.

    The FASTA is copied into the tmp dir first: the index build writes a
    ``.fai`` sidecar next to the input, which must not pollute the checked-in
    fixtures directory.
    """
    import shutil
    from pathlib import Path

    from beguidedoff.fast.fmindex import build_index

    src = Path(__file__).parent.parent / "fixtures" / "sample_genome.fasta"
    work = tmp_path_factory.mktemp("fmindex")
    fasta = work / src.name
    shutil.copy(src, fasta)
    out = work / "sample_genome.fmi"
    build_index(str(fasta), str(out))
    return out


@pytest.fixture(scope="module")
def mini_index_guides(mini_index):
    """A guide that occurs in the sample genome (with an NGG PAM) plus one
    that occurs nowhere, as a full-format BEguidedOFF input TSV."""
    from pathlib import Path

    fasta = mini_index.with_suffix(".fasta")
    seqs = {}
    name = None
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip().upper()
    chr1 = seqs["chr1"]

    hit_guide = None
    for i in range(len(chr1) - 23):
        if chr1[i + 21 : i + 23] == "GG":
            guide = chr1[i : i + 20]
            # ABE (window 4-8) must see an editable base, otherwise the only
            # hit is non-editable and filtered out of the results.
            if any(base == "A" for base in guide[3:8]):
                hit_guide = guide
                break
    assert hit_guide, "sample genome has no editable NGG site?"

    miss_guide = "ACGTACGTACGTACGTACGT"
    assert miss_guide not in chr1

    tsv = mini_index.parent / "guides.tsv"
    tsv.write_text(
        "guide_seq\tpam\tedit_window_start\tedit_window_end\teditor_type\n"
        f"{hit_guide}\tNGG\t4\t8\tABE\n"
        f"{miss_guide}\tNGG\t4\t8\tABE\n"
    )
    return tsv


def _minimal_config(**offtargets):
    """A minimal new-style pipeline config with an offtargets block."""
    block = {
        "beguidedoff_genome": "/tmp/genome.npz",
        "count_threshold_score": 15.0,
        "remove_guides_score": 5.0,
        "remove_guides_n_offtargets": 100,
    }
    block.update(offtargets)
    return {
        "output_prefix": "/tmp/run",
        "input": {"targets": "/tmp/targets.tsv", "genome": "/tmp/genome.fa"},
        "offtargets": block,
    }


# ---------------------------------------------------------------------------
# config resolution
# ---------------------------------------------------------------------------

class TestOfftargetConfigResolution:
    def test_default_mode_is_extensive(self):
        cfg = _read_new_style_config(_minimal_config())
        assert cfg["beguidedoff_mode"] == "extensive"
        assert cfg["beguidedoff_index"] == ""
        assert cfg["max_mismatches"] == 5

    def test_fast_mode_and_index_path(self):
        cfg = _read_new_style_config(
            _minimal_config(mode="fast", index_path="/data/genome.fmi")
        )
        assert cfg["beguidedoff_mode"] == "fast"
        assert cfg["beguidedoff_index"] == "/data/genome.fmi"

    def test_mismatch_tolerance_alias(self):
        # The web UI historically sent the mismatch budget as
        # ``mismatch_tolerance`` (where it was silently ignored).
        cfg = _read_new_style_config(_minimal_config(mismatch_tolerance=3))
        assert cfg["max_mismatches"] == 3

    def test_max_mismatches_takes_precedence(self):
        cfg = _read_new_style_config(
            _minimal_config(max_mismatches=4, mismatch_tolerance=3)
        )
        assert cfg["max_mismatches"] == 4

    def test_invalid_mode_raises(self):
        with pytest.raises(PipelineError, match="offtargets.mode"):
            _read_new_style_config(_minimal_config(mode="turbo"))


# ---------------------------------------------------------------------------
# search-cutoff resolution (fast mode soundness guard)
# ---------------------------------------------------------------------------

class TestResolveFastSearchCutoff:
    @staticmethod
    def _resolved(**offtargets):
        return _read_new_style_config(_minimal_config(**offtargets))

    def test_removal_threshold_becomes_cutoff(self):
        assert _resolve_fast_search_cutoff(self._resolved()) == 5.0

    def test_disabled_removal_falls_back_to_count_threshold(self):
        cfg = self._resolved(remove_guides_score=None, count_threshold_score=3.0)
        assert _resolve_fast_search_cutoff(cfg) == 3.0

    def test_disabled_removal_clamps_high_count_threshold(self):
        cfg = self._resolved(remove_guides_score=None, count_threshold_score=15.0)
        assert _resolve_fast_search_cutoff(cfg) == pytest.approx(
            FAST_MAX_SEARCH_CUTOFF
        )

    def test_unsound_removal_threshold_raises(self):
        cfg = self._resolved(remove_guides_score=8.0)
        with pytest.raises(PipelineError, match="12 bp exact anchor"):
            _resolve_fast_search_cutoff(cfg)


# ---------------------------------------------------------------------------
# end-to-end fast-mode search through the pipeline seam
# ---------------------------------------------------------------------------

class TestRunBeguidedoffFastMode:
    def test_search_finds_only_real_matches(
        self, mini_index, mini_index_guides, tmp_path
    ):
        out = tmp_path / "ot.tsv"
        ot_df = _run_beguidedoff(
            mini_index_guides,
            "base_editor_profiles.tsv",  # unused by full-format input
            "/tmp/genome.npz",           # unused by fast mode
            max_mismatches=5,            # unused by fast mode
            count_threshold_score=15.0,
            output_file=out,
            mode="fast",
            index_path=str(mini_index),
            search_cutoff=5.0,
        )
        # The genome-occurring guide has exactly one perfect hit (the sampled
        # site itself); the absent guide has none and must not appear.
        assert len(ot_df) == 1
        assert ot_df.iloc[0]["score"] == 0.0
        assert ot_df.iloc[0]["n_offtargets_below_threshold"] == 1
        assert out.exists()

    def test_count_threshold_clamped_to_search_cutoff(
        self, mini_index, mini_index_guides, tmp_path, caplog
    ):
        with caplog.at_level(logging.WARNING, logger="pipeline"):
            ot_df = _run_beguidedoff(
                mini_index_guides,
                "base_editor_profiles.tsv",
                "/tmp/genome.npz",
                max_mismatches=5,
                count_threshold_score=15.0,
                output_file=None,
                mode="fast",
                index_path=str(mini_index),
                search_cutoff=5.0,
            )
        assert any("exceeds the fast-mode search cutoff" in r.message
                   for r in caplog.records)
        # Counts reflect the clamped threshold (score <= 5): the single
        # perfect match.
        assert ot_df.iloc[0]["n_offtargets_below_threshold"] == 1

    def test_missing_index_raises(self, mini_index_guides, tmp_path):
        with pytest.raises(PipelineError, match="FM-index"):
            _run_beguidedoff(
                mini_index_guides,
                "base_editor_profiles.tsv",
                "/tmp/genome.npz",
                max_mismatches=5,
                count_threshold_score=15.0,
                output_file=None,
                mode="fast",
                index_path="",
                search_cutoff=5.0,
            )

    def test_unknown_mode_raises(self, mini_index, mini_index_guides, tmp_path):
        with pytest.raises(PipelineError, match="Unknown BEguidedOFF mode"):
            _run_beguidedoff(
                mini_index_guides,
                "base_editor_profiles.tsv",
                "/tmp/genome.npz",
                max_mismatches=5,
                count_threshold_score=15.0,
                output_file=None,
                mode="turbo",
                index_path=str(mini_index),
                search_cutoff=5.0,
            )
