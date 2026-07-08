"""Integration tests for the FastAPI app: HTTP layer around the scheduler."""

from __future__ import annotations

import asyncio
import base64
import io
import threading
import time
from pathlib import Path

import pandas as pd
import pytest


# Skip the whole module if web extras aren't installed.
pytest.importorskip("fastapi")
pytest.importorskip("uvicorn")


@pytest.fixture
def cfg_and_db(tmp_path, monkeypatch):
    """Set up a WebConfig + DB on a tmp path; restart the scheduler singleton."""
    from beguided.web import config as cfg_mod
    from beguided.web import scheduler as sched_mod
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
        genome_fasta="/nonexistent",  # web /run checks chromosomes; we won't trigger that
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        max_input_targets=1000,
        db_path=str(db),
        queue_max=3,
        cooldown_seconds=0,
        max_upload_mb=5,
    ))
    sch = Scheduler(db_path=str(db), queue_max=3, cooldown_seconds=0)
    install(sch)
    sch.start()
    yield cfg_mod.get_config(), sch


def _targets_tsv_b64(rows=("chr1\t50\tT\tC\n",)) -> str:
    content = "".join(rows).encode()
    return base64.b64encode(content).decode("ascii")


def _client(cfg_and_db):
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app

    return TestClient(create_app())


def test_run_returns_job_id_and_token(cfg_and_db):
    cfg, sch = cfg_and_db
    client = _client(cfg_and_db)

    # Patch tasks.execute_job to avoid touching the genome.
    from beguided.web import tasks as tasks_mod

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "fake ok",
            "records": [{"chromosome": "chr1", "position": 50}],
        }

    threading.Thread(
        target=lambda: setattr(tasks_mod, "execute_job", fake_execute_job),
        daemon=True,
    ).start()
    # Wait briefly for the monkey-patch above (TestClient uses no event loop here).
    import beguided.web.scheduler as sched_mod
    sched_mod.execute_job = fake_execute_job

    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert "job_id" in body
    assert "revoke_token" in body

    # Wait for completion
    for _ in range(40):
        s = client.get("/status", params={"job_id": body["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)
    assert s["status"] == "complete"

    r = client.get("/results", params={"job_id": body["job_id"]}).json()
    assert r["ok"] is True
    assert r["records"] == [{"chromosome": "chr1", "position": 50}]


def test_run_returns_429_on_cooldown(cfg_and_db):
    cfg, sch = cfg_and_db
    from beguided.web.scheduler import install, Scheduler
    # Re-install with a large cooldown so the second submit trips it.
    sch_big = Scheduler(
        db_path=cfg.db_path, queue_max=3, cooldown_seconds=60
    )
    install(sch_big)
    sch_big.start()
    try:
        # First submit succeeds.
        import beguided.web.scheduler as sched_mod
        sched_mod.execute_job = lambda **kw: {
            "job_id": kw["job_id"], "status": "complete",
            "message": "ok", "records": [],
        }

        client = _client(cfg_and_db)
        first = client.post(
            "/run",
            data={
                "raw_b64": _targets_tsv_b64(),
                "base_editor": "ABE8e",
            },
        )
        assert first.status_code == 200

        second = client.post(
            "/run",
            data={
                "raw_b64": _targets_tsv_b64(),
                "base_editor": "ABE8e",
            },
        )
        assert second.status_code == 429
        assert second.headers.get("Retry-After")
        assert second.json()["reason"] == "cooldown"
    finally:
        sch_big.shutdown(drain=True)


def test_run_rejects_oversized_upload(cfg_and_db):
    client = _client(cfg_and_db)
    big = b"a" * (6 * 1024 * 1024)  # 6 MB > 5 MB cap
    res = client.post(
        "/upload",
        files={"file": ("big.tsv", io.BytesIO(big), "text/tab-separated-values")},
    )
    # Either the middleware rejects it (413) or pandas parses it then we hit
    # a different cap. With our 5 MB middleware the middleware path triggers.
    assert res.status_code in (413, 400)


def test_run_rejects_unknown_base_editor(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "totally-fake-preset",
        },
    )
    assert res.status_code == 400
    assert "Unknown base editor" in res.json()["error"]


def test_cancel_requires_valid_token(cfg_and_db):
    import beguided.web.scheduler as sched_mod

    def slow_execute(**kw):
        import time as _t
        _t.sleep(0.3)  # give cancel a window
        return {
            "job_id": kw["job_id"], "status": "complete",
            "message": "ok", "records": [],
        }

    sched_mod.execute_job = slow_execute
    client = _client(cfg_and_db)
    submit = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
        },
    ).json()

    # Wait until the job is actually running (not queued-and-gone).
    import time as _t
    for _ in range(40):
        s = client.get("/status", params={"job_id": submit["job_id"]}).json()
        if s["status"] in ("running", "complete", "cancelled", "error"):
            break
        _t.sleep(0.05)

    bad = client.post(
        "/cancel",
        data={"job_id": submit["job_id"], "revoke_token": "wrong"},
    )
    assert bad.status_code == 403

    good = client.post(
        "/cancel",
        data={"job_id": submit["job_id"], "revoke_token": submit["revoke_token"]},
    )
    assert good.status_code == 200


def test_genome_status_endpoint_returns_metadata(cfg_and_db):
    """The /genome/status endpoint should always respond and reflect the
    preload config from WebConfig (preload_enabled flag)."""
    from beguided.web import genome_cache

    genome_cache._reset_for_test()
    client = _client(cfg_and_db)
    res = client.get("/genome/status")
    assert res.status_code == 200
    body = res.json()
    assert body["ok"] is True
    # cfg_and_db fixture doesn't toggle preload_genome explicitly; default is True.
    assert body["preload_enabled"] is True
    # With no preload started, status is "not_loaded".
    assert body["status"] == "not_loaded"
    # Label starts empty; it gets derived after a successful load.
    assert body["label"] == ""


def test_gtf_status_endpoint_returns_metadata(cfg_and_db):
    """The /gtf/status endpoint should always respond and reflect the
    preload config from WebConfig (preload_enabled flag)."""
    from beguided.web import gtf_cache

    gtf_cache._reset_for_test()
    client = _client(cfg_and_db)
    res = client.get("/gtf/status")
    assert res.status_code == 200
    body = res.json()
    assert body["ok"] is True
    # cfg_and_db fixture doesn't toggle preload_gtf explicitly; default is True.
    assert body["preload_enabled"] is True
    # With no preload started, status is "not_loaded".
    assert body["status"] == "not_loaded"


def test_config_endpoint_exposes_preload_genome(cfg_and_db):
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    body = client.get("/config").json()
    assert body["ok"] is True
    assert body["preload_genome"] == cfg.preload_genome


def test_config_endpoint_exposes_preload_gtf(cfg_and_db):
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    body = client.get("/config").json()
    assert body["ok"] is True
    assert body["preload_gtf"] == cfg.preload_gtf


def test_df_to_records_rounds_guide_score_to_3dp(tmp_path):
    """``guide_score`` should be rounded to 3 decimal places in web records
    so the UI and downloads show a compact ``0.73``-style value rather than
    the pipeline's full float precision."""
    from beguided.web.tasks import _df_to_records

    tsv = tmp_path / "results.tsv"
    tsv.write_text(
        "chromosome\tposition\tguide_score\n"
        "chr1\t10\t0.7345678\n"
        "chr2\t20\t0.5\n"
        "chr3\t30\t0.7355\n"
    )
    records = _df_to_records(tsv)
    scores = [r["guide_score"] for r in records]
    # Round half-to-even: 0.7345678 -> 0.735, 0.7355 -> 0.736 (banker's rounding
    # would give 0.736 anyway; either way must be ≤ 3 decimal places).
    assert scores[0] == round(0.7345678, 3)
    assert scores[1] == 0.5
    assert scores[2] == round(0.7355, 3)
    for s in scores:
        assert s is None or len(str(s).split(".")[-1]) <= 3


def test_df_to_records_preserves_non_score_columns(tmp_path):
    """Other columns must pass through untouched (including strings, ints,
    and missing values)."""
    from beguided.web.tasks import _df_to_records

    tsv = tmp_path / "results.tsv"
    tsv.write_text(
        "chromosome\tposition\tguide_score\tguide_rna_sequence\tbystander_effect\n"
        "chr1\t10\t0.123456789\tAAAAAAAAAAAAAAAAAAAA\tsilent\n"
    )
    records = _df_to_records(tsv)
    assert len(records) == 1
    r = records[0]
    assert r["chromosome"] == "chr1"
    assert r["position"] == 10
    assert r["guide_score"] == round(0.123456789, 3)
    assert r["guide_rna_sequence"] == "AAAAAAAAAAAAAAAAAAAA"
    assert r["bystander_effect"] == "silent"


def test_run_pipeline_invoked_with_unfiltered_defaults(cfg_and_db):
    """The web /run should hard-code filter_mode='all' and skip off-targets.

    The pipeline no longer filters during execution; the user filters
    dynamically in the Results tab. Verify the kwargs the scheduler sees.
    """
    cfg, _ = cfg_and_db
    captured = {}

    import beguided.web.scheduler as sched_mod

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "fake ok",
            "records": [{"chromosome": "chr1", "position": 50}],
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)

    res = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    )
    assert res.status_code == 200
    body = res.json()
    for _ in range(40):
        s = client.get("/status", params={"job_id": body["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)
    assert s["status"] == "complete"
    assert captured.get("filter_mode") == "all"
    assert captured.get("search_offtargets") is False
    assert captured.get("guide_score_threshold") == 0.0
    assert captured.get("passthrough_extras") is False


def test_run_ignores_legacy_form_fields(cfg_and_db):
    """The old filter_mode / score / search_offtargets form fields should be
    silently ignored (they're no longer part of the API)."""
    import beguided.web.scheduler as sched_mod

    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"], "status": "complete",
            "message": "ok", "records": [],
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)

    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
            "filter_mode": "silent",
            "guide_score_threshold": "0.9",
            "search_offtargets": "true",
        },
    )
    assert res.status_code == 200
    for _ in range(40):
        s = client.get("/status", params={"job_id": res.json()["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)
    # The web defaults should win, not the form fields.
    assert captured["filter_mode"] == "all"
    assert captured["search_offtargets"] is False
    assert captured["guide_score_threshold"] == 0.0


def test_results_returns_full_unfiltered_records(cfg_and_db):
    """With pipeline-level filtering disabled, /results must return every
    guide that the pipeline produced."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "silent",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "no_bystanders",
         "bystander_positions": "", "guide_score": None,
         "guide_rna_sequence": "GGGGGGGGGGGGGGGGGGGG"},
    ]

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"], "status": "complete",
            "message": "ok", "records": records,
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)

    submit = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    ).json()
    for _ in range(40):
        s = client.get("/status", params={"job_id": submit["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)

    r = client.get("/results", params={"job_id": submit["job_id"]}).json()
    assert r["ok"] is True
    assert len(r["records"]) == 3  # unfiltered


def test_download_results_applies_dynamic_filters(cfg_and_db):
    """The /download/results endpoint should respect the filter_mode and
    guide_score_threshold query params, mirroring the client-side filters."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "silent",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "no_bystanders",
         "bystander_positions": "", "guide_score": 0.6,
         "guide_rna_sequence": "GGGGGGGGGGGGGGGGGGGG"},
    ]

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"], "status": "complete",
            "message": "ok", "records": records,
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)
    submit = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    ).json()
    for _ in range(40):
        s = client.get("/status", params={"job_id": submit["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)

    # No filter -> all three rows.
    res_all = client.get(
        "/download/results",
        params={"job_id": submit["job_id"]},
    )
    assert res_all.status_code == 200
    lines = res_all.text.strip().split("\n")
    assert lines[0].startswith("chromosome")
    assert len(lines) == 1 + 3

    # silent: keep bystander_effect in {silent, no_bystanders}.
    res_silent = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "filter_mode": "silent"},
    )
    assert res_silent.status_code == 200
    rows = [l.split("\t") for l in res_silent.text.strip().split("\n")[1:]]
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50]

    # guide_score_threshold=0.5 -> drop the 0.2 row.
    res_score = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "guide_score_threshold": 0.5},
    )
    assert res_score.status_code == 200
    rows = [l.split("\t") for l in res_score.text.strip().split("\n")[1:]]
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50]

    # Combined: silent AND score >= 0.8 -> only the 0.9 silent row.
    res_combo = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "silent",
            "guide_score_threshold": 0.8,
        },
    )
    assert res_combo.status_code == 200
    rows = [l.split("\t") for l in res_combo.text.strip().split("\n")[1:]]
    assert len(rows) == 1
    assert rows[0][1] == "50"

    # Invalid filter_mode -> 400.
    bad = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "filter_mode": "bogus"},
    )
    assert bad.status_code == 400


def test_download_results_filename_reflects_filters(cfg_and_db):
    """The Content-Disposition filename should encode the active filters."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "silent",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
    ]

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"], "status": "complete",
            "message": "ok", "records": records,
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)
    submit = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    ).json()
    for _ in range(40):
        s = client.get("/status", params={"job_id": submit["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)

    base = client.get(
        "/download/results",
        params={"job_id": submit["job_id"]},
    )
    assert "beguided_results.tsv" in base.headers["content-disposition"]

    filtered = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "silent",
            "guide_score_threshold": 0.5,
        },
    )
    assert "beguided_results_silent_s0.5.tsv" in filtered.headers["content-disposition"]


def test_download_results_xlsx(cfg_and_db):
    """The /download/results endpoint with format=xlsx should serve an XLSX
    spreadsheet that round-trips back to the same DataFrame."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "silent",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "no_bystanders",
         "bystander_positions": "", "guide_score": 0.6,
         "guide_rna_sequence": "GGGGGGGGGGGGGGGGGGGG"},
    ]

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"], "status": "complete",
            "message": "ok", "records": records,
        }

    sched_mod.execute_job = fake_execute_job
    client = _client(cfg_and_db)
    submit = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    ).json()
    for _ in range(40):
        s = client.get("/status", params={"job_id": submit["job_id"]}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)

    # No filter -> all three rows.
    res_all = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "format": "xlsx"},
    )
    assert res_all.status_code == 200
    assert (
        res_all.headers["content-type"]
        == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
    assert "beguided_results.xlsx" in res_all.headers["content-disposition"]
    df = pd.read_excel(io.BytesIO(res_all.content), engine="openpyxl")
    assert len(df) == 3
    assert "chromosome" in df.columns
    assert "guide_score" in df.columns

    # silent filter applied -> 2 rows.
    res_silent = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "silent",
            "format": "xlsx",
        },
    )
    assert res_silent.status_code == 200
    df_silent = pd.read_excel(io.BytesIO(res_silent.content), engine="openpyxl")
    assert sorted(int(p) for p in df_silent["position"]) == [30, 50]

    # Filename reflects filters.
    filtered = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "silent",
            "guide_score_threshold": 0.5,
            "format": "xlsx",
        },
    )
    assert (
        "beguided_results_silent_s0.5.xlsx" in filtered.headers["content-disposition"]
    )

    # Default format (no param) is still TSV.
    res_default = client.get(
        "/download/results",
        params={"job_id": submit["job_id"]},
    )
    assert res_default.status_code == 200
    assert "text/tab-separated-values" in res_default.headers["content-type"]
    assert "beguided_results.tsv" in res_default.headers["content-disposition"]

    # Invalid format -> 400.
    bad = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "format": "csv"},
    )
    assert bad.status_code == 400
