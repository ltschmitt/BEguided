"""Tests for ``beguided.web.tasks``.

Covers the TSV→records conversion used to feed the results UI and the
TSV/XLSX download, with a focus on the sequence-mode column-drop behaviour.
None of these tests require fastapi/uvicorn — they exercise pure helpers
and ``execute_job`` directly, so they run in the default pixi env.
"""

from __future__ import annotations

import base64
from pathlib import Path

import pandas as pd
import pytest

# Capture the genuine execute_job reference at module import time. Several
# tests in test_app.py patch ``tasks.execute_job`` (and ``scheduler.execute_job``)
# directly via raw ``setattr`` without using pytest's ``monkeypatch`` teardown,
# so by the time this file's tests run those modules' ``execute_job`` name
# may point to a stub. Capturing the real reference here, before any test
# function runs, gives us a stable handle we can re-attach before our tests.
from beguided.web.tasks import (  # noqa: E402
    SEQUENCE_MODE_DROP_COLUMNS,
    _df_to_records,
    _merge_target_desired,
    execute_job as _real_execute_job,
)


def _write_results_tsv(path: Path) -> Path:
    """Write a small canonical-looking results TSV with every column we
    care about, then return the path. Mirrors the schema produced by the
    pipeline's final step."""
    df = pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "position": [50, 75],
            "target_base": ["T", "A"],
            "desired_base": ["C", "G"],
            "genome_base": ["T", "A"],
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT"],
            "pam_sequence": ["AGG", "TGG"],
            "editing_position": [4, 5],
            "strand_used": ["-", "+"],
            "bystander_positions": ["", "3,7"],
            "bystander_effect": ["no_bystanders", "missense"],
            "base_editor": ["ABE8e", "ABE8e"],
            "guide_score": [0.0, 0.812345],
        }
    )
    df.to_csv(path, sep="\t", index=False)
    return path


class TestDfToRecords:
    def test_returns_all_columns_by_default(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        records = _df_to_records(path)
        assert len(records) == 2
        ref_cols = _write_results_tsv(tmp_path / "ref.tsv").read_text().splitlines()[
            0
        ].split("\t")
        assert set(records[0].keys()) == set(ref_cols)

    def test_missing_file_returns_empty_list(self, tmp_path):
        assert _df_to_records(tmp_path / "nope.tsv") == []

    def test_drop_columns_removes_keys(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        records = _df_to_records(
            path, drop_columns=("position", "genome_base", "strand_used")
        )
        assert records, "should still have rows"
        for record in records:
            assert "position" not in record
            assert "genome_base" not in record
            assert "strand_used" not in record
            # Untouched columns are still present.
            assert record["chromosome"] == "chr1"
            assert record["target_base"] in {"T", "A"}
            assert record["base_editor"] == "ABE8e"

    def test_drop_columns_ignores_missing_keys(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        # No KeyError even when asked to drop a column that isn't in the TSV.
        records = _df_to_records(path, drop_columns=("not_a_column", "position"))
        assert records
        for record in records:
            assert "position" not in record
            assert "not_a_column" not in record

    def test_drop_columns_empty_tuple_is_noop(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        all_records = _df_to_records(path)
        noop_records = _df_to_records(path, drop_columns=())
        assert all_records == noop_records

    def test_guide_score_rounded_to_three_decimals(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        records = _df_to_records(path)
        scores = [r["guide_score"] for r in records]
        # 0.0 stays 0.0; 0.812345 is rounded to 0.812.
        assert scores[0] == 0.0
        assert scores[1] == 0.812

    def test_does_not_touch_on_disk_tsv(self, tmp_path):
        path = _write_results_tsv(tmp_path / "results.tsv")
        before = path.read_text()
        _df_to_records(
            path, drop_columns=("position", "genome_base", "strand_used")
        )
        assert path.read_text() == before


class TestSequenceModeDropColumnsConstant:
    """Sanity-check the module-level drop list so that adding/removing a
    column here is a deliberate, visible change in the test suite."""

    def test_constant_lists_the_three_expected_columns(self):
        assert SEQUENCE_MODE_DROP_COLUMNS == ("position", "genome_base", "strand_used")

    def test_constant_is_tuple(self):
        assert isinstance(SEQUENCE_MODE_DROP_COLUMNS, tuple)


def _python_target_key(record: dict) -> str:
    """Mirror of the JS ``targetKey`` helper in ``app.js``.

    Used by the regression tests below to assert the frontend's count math
    stays in sync with the record shape that ``_df_to_records`` produces in
    each input mode.
    """
    chrom = record.get("chromosome") or record.get("sequence") or ""
    pos = record.get("position") or ""
    return f"{chrom}__{pos}"


class TestDistinctTargetsMath:
    """Regression coverage for the "Found guides for X/N targets" count.

    The frontend counts distinct (chromosome, position) pairs across the
    returned records. In spreadsheet mode that key uses ``chromosome`` +
    ``position`` directly. In sequence mode the backend renames
    ``chromosome`` -> ``sequence`` and drops ``position``, leaving only the
    ``sequence`` value to discriminate targets. Because the parser enforces
    unique sequence headers, each input target has a unique ``sequence``
    value and the count math must still return the correct total.

    Without the ``targetKey`` fallback in JS, every sequence-mode record
    collapses to the same key and the count is always 1.
    """

    def _write_two_distinct_targets(self, path: Path) -> Path:
        pd.DataFrame(
            {
                "chromosome": ["chrA", "chrB"],
                "position": [50, 50],
                "target_base": ["A", "A"],
                "desired_base": ["G", "G"],
                "genome_base": ["A", "A"],
                "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA"],
                "pam_sequence": ["AGG", "AGG"],
                "editing_position": [4, 4],
                "strand_used": ["+", "+"],
                "bystander_positions": ["", ""],
                "bystander_effect": ["no_bystanders", "no_bystanders"],
                "base_editor": ["ABE8e", "ABE8e"],
                "guide_score": [0.5, 0.5],
            }
        ).to_csv(path, sep="\t", index=False)
        return path

    def test_spreadsheet_mode_distinct_targets_count_is_two(self, tmp_path):
        path = self._write_two_distinct_targets(tmp_path / "results.tsv")
        records = _df_to_records(path)
        assert len({_python_target_key(r) for r in records}) == 2

    def test_sequence_mode_distinct_targets_count_is_two(self, tmp_path):
        path = self._write_two_distinct_targets(tmp_path / "results.tsv")
        records = _df_to_records(path, drop_columns=SEQUENCE_MODE_DROP_COLUMNS)
        for record in records:
            record["sequence"] = record.pop("chromosome")
        records = _merge_target_desired(records)
        assert len({_python_target_key(r) for r in records}) == 2

    def test_sequence_mode_records_lack_position_and_chromosome(self, tmp_path):
        path = self._write_two_distinct_targets(tmp_path / "results.tsv")
        records = _df_to_records(path, drop_columns=SEQUENCE_MODE_DROP_COLUMNS)
        for record in records:
            record["sequence"] = record.pop("chromosome")
        records = _merge_target_desired(records)
        for record in records:
            assert "position" not in record
            assert "chromosome" not in record
            assert "sequence" in record


class TestMergeTargetDesired:
    """Unit tests for the sequence-mode ``target_base``+``desired_base``
    merge that backs both the on-screen table and the downloaded files."""

    def test_basic_merge(self):
        records = [
            {"chromosome": "chr1", "position": 5, "target_base": "A", "desired_base": "G", "base_editor": "ABE8e"},
            {"chromosome": "chr1", "position": 7, "target_base": "T", "desired_base": "C", "base_editor": "ABE8e"},
        ]
        out = _merge_target_desired(records)
        assert len(out) == 2
        assert out[0]["desired_edit"] == "A>G"
        assert out[1]["desired_edit"] == "T>C"
        assert "target_base" not in out[0]
        assert "desired_base" not in out[0]

    def test_preserves_column_order(self):
        records = [{"chromosome": "chr1", "position": 5, "target_base": "A", "desired_base": "G", "base_editor": "ABE8e"}]
        out = _merge_target_desired(records)
        keys = list(out[0].keys())
        assert keys.index("desired_edit") == keys.index("position") + 1

    def test_missing_target_base_passes_through(self):
        records = [{"chromosome": "chr1", "desired_base": "G", "base_editor": "ABE8e"}]
        out = _merge_target_desired(records)
        assert "desired_edit" not in out[0]
        assert "desired_base" not in out[0]

    def test_missing_desired_base_passes_through(self):
        records = [{"chromosome": "chr1", "target_base": "A", "base_editor": "ABE8e"}]
        out = _merge_target_desired(records)
        assert "desired_edit" not in out[0]
        assert "target_base" not in out[0]

    def test_empty_string_treated_as_missing(self):
        records = [{"chromosome": "chr1", "target_base": "", "desired_base": "G", "base_editor": "ABE8e"}]
        out = _merge_target_desired(records)
        assert "desired_edit" not in out[0]

    def test_empty_list(self):
        assert _merge_target_desired([]) == []

    def test_none_row_preserved(self):
        assert _merge_target_desired([None]) == [None]


# ---------------------------------------------------------------------------
# Integration: execute_job end-to-end column-drop wiring
# ---------------------------------------------------------------------------


# 40-nt sequence with an 'A' at 1-based position 17 (matches the canonical
# example used across the web tests; valid A→G edit site).
_SEQUENCE = "ATAGAGTATAGACCCAA" + "CCTATCCAGTCTACCTAGGGAATC"


def _minimal_web_config(tmp_path: Path):
    """Build a WebConfig that is sufficient for a sequence-mode job.

    Genome/GTF paths point at non-existent files because sequence mode
    does not touch them.
    """
    from beguided.web.config import WebConfig, set_config

    presets = tmp_path / "presets.tsv"
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
    )
    cfg = WebConfig(
        base_editor_profiles_file=str(presets),
        genome_fasta="/nonexistent",
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        max_input_targets=100,
        db_path=str(tmp_path / "jobs.db"),
        queue_max=3,
        cooldown_seconds=0,
        max_upload_mb=5,
        preload_genome=False,
        preload_gtf=False,
    )
    set_config(cfg)
    return cfg


def _fake_pipeline_writer(*args, **kwargs) -> None:
    """Stand-in for ``run_pipeline`` used by the integration tests below.

    Writes the exact TSV the real pipeline would produce at
    ``<work_dir>/beguided_results.tsv`` so that ``execute_job`` can read it
    back via ``_df_to_records``. Accepts arbitrary positional/keyword args
    so it can stand in for the real ``run_pipeline`` signature.
    """
    # The first positional arg is the config dict; output_prefix is a key
    # inside it and that's where the TSV must land.
    config = args[0] if args else kwargs.get("config", {})
    work_dir = config.get("output_prefix") if isinstance(config, dict) else None
    if not work_dir:
        return
    out = Path(work_dir) / "beguided_results.tsv"
    pd.DataFrame(
        {
            "chromosome": ["custom"],
            "position": [17],
            "target_base": ["A"],
            "desired_base": ["G"],
            "genome_base": ["A"],
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA"],
            "pam_sequence": ["AGG"],
            "editing_position": [4],
            "strand_used": ["-"],
            "bystander_positions": [""],
            "bystander_effect": ["no_bystanders"],
            "base_editor": ["ABE8e"],
            "guide_score": [0.5],
        }
    ).to_csv(out, sep="\t", index=False)


@pytest.fixture(autouse=True)
def _restore_execute_job(monkeypatch):
    """Restore the real ``execute_job`` before each integration test.

    Several upstream tests in ``test_app.py`` patch ``tasks.execute_job``
    and ``scheduler.execute_job`` directly (without pytest's ``monkeypatch``
    teardown), which leaves a stub in place across tests. We re-attach
    both modules' ``execute_job`` to the genuine function captured at
    module import time so this test sees the real pipeline plumbing.
    """
    import beguided.web.scheduler as sched_mod
    import beguided.web.tasks as tasks_mod

    monkeypatch.setattr(tasks_mod, "execute_job", _real_execute_job)
    monkeypatch.setattr(sched_mod, "execute_job", _real_execute_job)
    yield


def _has_starlette() -> bool:
    try:
        import starlette  # noqa: F401
    except ImportError:
        return False
    return True


class TestExecuteJobColumnDrops:
    """End-to-end coverage that ``execute_job`` actually drops the
    sequence-mode columns from the records it returns.

    These tests require fastapi/starlette because ``execute_job`` lazily
    imports ``beguided.web.safety`` (which uses Starlette middleware).
    The pure ``_df_to_records`` unit tests above run in every env.
    """

    pytestmark = pytest.mark.skipif(
        not _has_starlette(),
        reason="starlette not installed (install with: pixi install --features web)",
    )

    def _run_execute_job(self, mode: str, monkeypatch, tmp_path) -> dict:
        from beguided.web import tasks as tasks_mod

        monkeypatch.setattr(tasks_mod, "run_pipeline", _fake_pipeline_writer)

        _minimal_web_config(tmp_path)

        fasta = f">custom\n{_SEQUENCE}\n"
        custom_genome_b64 = base64.b64encode(fasta.encode()).decode("ascii")

        targets_df = pd.DataFrame(
            {
                "chromosome": ["custom"],
                "position": [17],
                "target_base": ["A"],
                "desired_base": ["G"],
            }
        )

        return tasks_mod.execute_job(
            job_id="job-test",
            targets_df=targets_df,
            base_editor_names=["ABE8e"],
            filter_mode="all",
            search_offtargets=False,
            guide_score_threshold=0.0,
            mode=mode,
            custom_genome_b64=custom_genome_b64 if mode == "sequence" else None,
        )

    def test_sequence_mode_drops_three_columns(self, monkeypatch, tmp_path):
        result = self._run_execute_job("sequence", monkeypatch, tmp_path)
        assert result["status"] == "complete", result
        records = result["records"]
        assert len(records) == 1
        rec = records[0]
        for col in SEQUENCE_MODE_DROP_COLUMNS:
            assert col not in rec, f"{col} should be dropped in sequence mode"
        # target_base and desired_base are merged into desired_edit.
        assert "target_base" not in rec
        assert "desired_base" not in rec
        assert rec["desired_edit"] == "A>G"
        # chromosome is renamed to sequence; it must not be present under
        # its old key in sequence mode.
        assert "chromosome" not in rec
        assert rec["sequence"] == "custom"
        assert rec["base_editor"] == "ABE8e"

    def test_sequence_mode_target_count_survives_column_rename(
        self, monkeypatch, tmp_path
    ):
        """Regression: with N input targets the "Found guides for X/N"
        count must reach X=N when every target produces a guide. The
        frontend counts distinct ``(chromosome, position)`` keys; in
        sequence mode ``position`` is dropped and ``chromosome`` is
        renamed to ``sequence``, so the JS code falls back to
        ``sequence`` only. This test confirms the record shape supports
        that fallback."""
        from beguided.web import tasks as tasks_mod

        def _fake_pipeline_two_targets(config, **_kwargs):
            work_dir = config["output_prefix"] if isinstance(config, dict) else None
            if not work_dir:
                return
            out = Path(work_dir) / "beguided_results.tsv"
            pd.DataFrame(
                {
                    "chromosome": ["t1", "t2"],
                    "position": [17, 17],
                    "target_base": ["A", "A"],
                    "desired_base": ["G", "G"],
                    "genome_base": ["A", "A"],
                    "guide_rna_sequence": [
                        "AAAAAAAAAAAAAAAAAAAA",
                        "AAAAAAAAAAAAAAAAAAAA",
                    ],
                    "pam_sequence": ["AGG", "AGG"],
                    "editing_position": [4, 4],
                    "strand_used": ["+", "+"],
                    "bystander_positions": ["", ""],
                    "bystander_effect": ["no_bystanders", "no_bystanders"],
                    "base_editor": ["ABE8e", "ABE8e"],
                    "guide_score": [0.5, 0.5],
                }
            ).to_csv(out, sep="\t", index=False)

        monkeypatch.setattr(tasks_mod, "run_pipeline", _fake_pipeline_two_targets)
        _minimal_web_config(tmp_path)

        fasta = (
            f">t1\n{_SEQUENCE}\n"
            f">t2\n{_SEQUENCE}\n"
        )
        custom_genome_b64 = base64.b64encode(fasta.encode()).decode("ascii")
        targets_df = pd.DataFrame(
            {
                "chromosome": ["t1", "t2"],
                "position": [17, 17],
                "target_base": ["A", "A"],
                "desired_base": ["G", "G"],
            }
        )

        result = tasks_mod.execute_job(
            job_id="job-test",
            targets_df=targets_df,
            base_editor_names=["ABE8e"],
            filter_mode="all",
            search_offtargets=False,
            guide_score_threshold=0.0,
            mode="sequence",
            custom_genome_b64=custom_genome_b64,
        )
        assert result["status"] == "complete", result
        records = result["records"]
        assert len(records) == 2
        # The frontend's targetKey falls back to r.sequence when r.position
        # is missing — each record's sequence value must therefore differ.
        sequences = {r.get("sequence") for r in records}
        assert sequences == {"t1", "t2"}

    def test_spreadsheet_mode_keeps_all_columns(self, monkeypatch, tmp_path):
        result = self._run_execute_job("spreadsheet", monkeypatch, tmp_path)
        assert result["status"] == "complete", result
        records = result["records"]
        assert len(records) == 1
        rec = records[0]
        for col in SEQUENCE_MODE_DROP_COLUMNS:
            assert col in rec, f"{col} should be present in spreadsheet mode"
        assert rec["position"] == 17
        assert rec["genome_base"] == "A"
        assert rec["strand_used"] == "-"
