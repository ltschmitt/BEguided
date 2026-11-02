"""Tests for the single-target input parser and the /single/parse endpoint."""

from __future__ import annotations

import base64
import re

import pandas as pd
import pytest

pytest.importorskip("fastapi")

from beguided.web.single_input import (
    TARGETS_COLUMNS,
    SingleTargetError,
    parse_single_target,
    row_to_targets_dataframe,
)


@pytest.fixture(autouse=True)
def _reset_genome_cache():
    """Keep the genome-cache singleton clean around every test so the
    best-effort chromosome check behaves deterministically (no leftover
    ready entries from other test modules)."""
    from beguided.web import genome_cache

    genome_cache._reset_for_test()
    yield
    genome_cache._reset_for_test()


# ---------------------------------------------------------------------------
# parse_single_target: pure-parser tests
# ---------------------------------------------------------------------------


class TestParseSingleTarget:
    def test_valid_row(self):
        row = parse_single_target("chr1", "100", "A>G")
        assert row == {
            "chromosome": "chr1",
            "position": 100,
            "target_base": "A",
            "desired_base": "G",
        }

    def test_position_is_python_int(self):
        row = parse_single_target("1", "25245350", "C>T")
        assert row["position"] == 25245350
        assert isinstance(row["position"], int)

    def test_all_valid_edit_directions(self):
        for edit, (target, desired) in {
            "A>G": ("A", "G"),
            "C>T": ("C", "T"),
            "G>A": ("G", "A"),
            "T>C": ("T", "C"),
        }.items():
            row = parse_single_target("chr1", "10", edit)
            assert row["target_base"] == target
            assert row["desired_base"] == desired

    def test_whitespace_and_case_normalized(self):
        row = parse_single_target("  chr1  ", "  100  ", "  a > g  ")
        assert row == {
            "chromosome": "chr1",
            "position": 100,
            "target_base": "A",
            "desired_base": "G",
        }

    def test_x_y_chromosomes_accepted(self):
        for chrom in ("X", "Y", "chrX", "chrY"):
            row = parse_single_target(chrom, "5", "T>C")
            assert row["chromosome"] == chrom

    def test_empty_chromosome_rejected(self):
        with pytest.raises(SingleTargetError, match="Chromosome is required"):
            parse_single_target("", "100", "A>G")
        with pytest.raises(SingleTargetError, match="Chromosome is required"):
            parse_single_target("   ", "100", "A>G")

    def test_chromosome_with_internal_whitespace_rejected(self):
        with pytest.raises(SingleTargetError, match="whitespace"):
            parse_single_target("chr 1", "100", "A>G")

    def test_empty_position_rejected(self):
        with pytest.raises(SingleTargetError, match="Position is required"):
            parse_single_target("chr1", "", "A>G")

    def test_non_integer_positions_rejected(self):
        for bad in ("1.5", "abc", "-5", "1e3", "1,000", "1 2"):
            with pytest.raises(SingleTargetError, match="positive integer"):
                parse_single_target("chr1", bad, "A>G")

    def test_zero_position_rejected(self):
        with pytest.raises(SingleTargetError, match="positive integer"):
            parse_single_target("chr1", "0", "A>G")

    def test_malformed_edit_rejected(self):
        for bad in ("", "AG", "A>", ">G", "A", "A>G>C"):
            with pytest.raises(SingleTargetError, match="X>Y"):
                parse_single_target("chr1", "100", bad)

    def test_unsupported_edit_direction_rejected(self):
        with pytest.raises(SingleTargetError, match="not supported"):
            parse_single_target("chr1", "100", "A>T")

    def test_none_inputs_rejected(self):
        with pytest.raises(SingleTargetError):
            parse_single_target(None, "100", "A>G")


# ---------------------------------------------------------------------------
# row_to_targets_dataframe
# ---------------------------------------------------------------------------


class TestRowToTargetsDataframe:
    def test_columns_and_values(self):
        row = parse_single_target("chr1", "100", "A>G")
        df = row_to_targets_dataframe(row)
        assert df.columns.tolist() == TARGETS_COLUMNS
        assert len(df) == 1
        assert df.iloc[0]["chromosome"] == "chr1"
        assert df.iloc[0]["position"] == 100
        assert df.iloc[0]["target_base"] == "A"
        assert df.iloc[0]["desired_base"] == "G"

    def test_position_column_is_integer(self):
        # The pipeline expects an integer position column, exactly like the
        # spreadsheet path (pd.read_csv produces int64).
        df = row_to_targets_dataframe(parse_single_target("chr1", "100", "A>G"))
        assert pd.api.types.is_integer_dtype(df["position"])


# ---------------------------------------------------------------------------
# /single/parse + /run: FastAPI-level integration
# ---------------------------------------------------------------------------


@pytest.fixture
def cfg_and_db(tmp_path):
    """Minimal WebConfig + scheduler for FastAPI integration tests.

    Mirrors the fixture in test_sequence.py: ``genome_fasta`` is set (to a
    nonexistent path) so one ``human`` genome is registered and genome
    resolution succeeds, while no preload is started (the best-effort
    chromosome check is then skipped, same as in /upload tests).
    """
    from beguided.web import config as cfg_mod
    from beguided.web.config import WebConfig, set_config
    from beguided.web.scheduler import Scheduler, install

    db = tmp_path / "jobs.db"
    presets = tmp_path / "presets.tsv"
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
    )
    set_config(WebConfig(
        base_editor_profiles_file=str(presets),
        genome_fasta="/nonexistent",
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        max_input_targets=100,
        db_path=str(db),
        queue_max=3,
        cooldown_seconds=0,
        max_upload_mb=5,
        usage_stats_enabled=False,
    ))
    sch = Scheduler(db_path=str(db), queue_max=3, cooldown_seconds=0)
    install(sch)
    sch.start()
    yield cfg_mod.get_config(), sch


def _client():
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app
    return TestClient(create_app())


class TestSingleParseEndpoint:
    def test_valid_target_returns_upload_payload(self, cfg_and_db):
        res = _client().post(
            "/single/parse",
            data={"chromosome": "chr1", "position": "100", "edit": "A>G",
                  "genome": "human"},
        )
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is True
        assert body["row_count"] == 1
        assert body["columns"] == TARGETS_COLUMNS
        assert body["preview"] == [
            {"chromosome": "chr1", "position": 100,
             "target_base": "A", "desired_base": "G"}
        ]

        targets_tsv = base64.b64decode(body["raw_b64"]).decode()
        assert targets_tsv == (
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t100\tA\tG\n"
        )

    def test_default_genome_used_when_omitted(self, cfg_and_db):
        # genome_fasta is configured, so an empty genome field resolves to
        # the default (first) genome instead of failing.
        res = _client().post(
            "/single/parse",
            data={"chromosome": "chr1", "position": "100", "edit": "C>T"},
        )
        assert res.status_code == 200
        assert res.json()["ok"] is True

    def test_unknown_genome_rejected(self, cfg_and_db):
        # The fixture config registers "human" (and the repo-default "mouse"),
        # so an unregistered key must be rejected.
        res = _client().post(
            "/single/parse",
            data={"chromosome": "chr1", "position": "100", "edit": "A>G",
                  "genome": "zebrafish"},
        )
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is False
        assert "Unknown genome" in body["error"]

    def test_invalid_fields_reported(self, cfg_and_db):
        for data in (
            {"chromosome": "", "position": "100", "edit": "A>G"},
            {"chromosome": "chr1", "position": "abc", "edit": "A>G"},
            {"chromosome": "chr1", "position": "100", "edit": "A>T"},
            {"chromosome": "chr 1", "position": "100", "edit": "A>G"},
        ):
            res = _client().post("/single/parse", data=data)
            assert res.status_code == 200
            body = res.json()
            assert body["ok"] is False, data
            assert body["error"]

    def test_unknown_chromosome_rejected_when_genome_loaded(self, cfg_and_db, monkeypatch):
        """The best-effort chromosome check fires when the genome is ready."""
        import beguided.web.app as app_mod

        monkeypatch.setattr(
            app_mod, "_genome_chromosome_keys", lambda key: ["chr1", "1", "X", "Y"]
        )
        res = _client().post(
            "/single/parse",
            data={"chromosome": "chr2", "position": "100", "edit": "A>G",
                  "genome": "human"},
        )
        body = res.json()
        assert body["ok"] is False
        assert "Unknown chromosome" in body["error"]
        assert "chr1' vs '1'" in body["error"]

    def test_known_chromosome_passes_when_genome_loaded(self, cfg_and_db, monkeypatch):
        import beguided.web.app as app_mod

        monkeypatch.setattr(
            app_mod, "_genome_chromosome_keys", lambda key: ["chr1", "1", "X", "Y"]
        )
        res = _client().post(
            "/single/parse",
            data={"chromosome": "1", "position": "100", "edit": "T>C",
                  "genome": "human"},
        )
        assert res.json()["ok"] is True


class TestRunWithSingleMode:
    def test_run_in_single_mode_enqueues_job(self, cfg_and_db):
        """`/run` with mode=single must reach the worker like spreadsheet
        mode, carrying the parsed 1-row targets and the selected genome."""
        parse_res = _client().post(
            "/single/parse",
            data={"chromosome": "chr1", "position": "100", "edit": "A>G",
                  "genome": "human"},
        )
        body = parse_res.json()
        assert body["ok"] is True

        captured = {}

        def fake_execute_job(**kwargs):
            captured.update(kwargs)
            return {
                "job_id": kwargs["job_id"],
                "status": "complete",
                "message": "fake ok",
                "records": [{"chromosome": "chr1", "position": 100}],
            }

        import beguided.web.scheduler as sched_mod
        sched_mod.execute_job = fake_execute_job

        run_res = _client().post(
            "/run",
            data={
                "raw_b64": body["raw_b64"],
                "base_editor": "ABE8e",
                "mode": "single",
                "genome": "human",
            },
        )
        assert run_res.status_code == 200, run_res.text
        run_body = run_res.json()
        assert run_body["ok"] is True

        # Wait for the worker thread to invoke our fake.
        import time
        for _ in range(20):
            if captured:
                break
            time.sleep(0.05)
        assert captured.get("mode") == "single"
        assert captured.get("genome_key") == "human"
        assert captured.get("custom_genome_b64") is None
        df = captured.get("targets_df")
        assert df is not None and len(df) == 1
        assert df.iloc[0]["chromosome"] == "chr1"
        assert df.iloc[0]["position"] == 100
        assert df.iloc[0]["target_base"] == "A"
        assert df.iloc[0]["desired_base"] == "G"

    def test_run_single_mode_requires_known_genome(self, cfg_and_db):
        # An unregistered genome key must be rejected with a 400 instead of
        # silently running against another reference.
        raw = base64.b64encode(
            b"chromosome\tposition\ttarget_base\tdesired_base\nchr1\t100\tA\tG\n"
        ).decode()
        res = _client().post(
            "/run",
            data={"raw_b64": raw, "base_editor": "ABE8e",
                  "mode": "single", "genome": "zebrafish"},
        )
        assert res.status_code == 400
        assert res.json()["ok"] is False

    def test_check_accepts_single_mode_payload(self, cfg_and_db):
        """/check treats the /single/parse payload like any other upload."""
        parse_res = _client().post(
            "/single/parse",
            data={"chromosome": "chr1", "position": "100", "edit": "A>G",
                  "genome": "human"},
        )
        raw = parse_res.json()["raw_b64"]
        res = _client().post(
            "/check", data={"raw_b64": raw, "base_editor": "ABE8e"}
        )
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is True
        assert body["matched"] == 1
        assert body["can_run"] is True


class TestIndexPageRendersSinglePanel:
    def test_single_panel_elements_render(self, cfg_and_db):
        """The index page renders the Single mode button, panel, fields,
        and Load example button (inside the run group)."""
        html = _client().get("/").text

        # Mode switch button wired to the panel.
        assert 'id="modeSingle"' in html
        assert 'data-mode="single"' in html
        assert 'aria-controls="singlePanel"' in html

        # Extract the single panel block so assertions cannot be satisfied
        # by the spreadsheet panel's elements (it has its own genome select).
        panel = re.search(
            r'<div id="singlePanel".*?(?=<div id="sequencePanel")', html, re.S
        ).group(0)

        assert 'class="group single-field"' in panel
        assert 'id="singleChromosome"' in panel
        assert 'type="text"' in panel
        assert 'placeholder="1"' in panel
        # Position is a text input with a numeric keypad hint (no spinners).
        assert 'id="singlePosition"' in panel
        assert 'inputmode="numeric"' in panel
        assert 'type="number"' not in panel
        # Editing type dropdown with the four supported directions.
        assert 'id="singleEditType"' in panel
        for opt in ("A>G", "C>T", "G>A", "T>C"):
            assert f'value="{opt}"' in panel
        # Genome selector keeps the spreadsheet-panel styling hook.
        assert 'genome-select-group' in panel
        assert 'id="singleGenomeSelect"' in panel
        # Load example sits inside the run group, next to the shared Run.
        assert 'id="loadSingleExampleBtn"' in panel
        run_group = re.search(
            r'<div id="singleRunGroup".*?</div>', panel, re.S
        ).group(0)
        assert 'id="loadSingleExampleBtn"' in run_group

        # Guide-tab help section.
        assert "Single target input" in html
