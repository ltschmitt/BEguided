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
        usage_stats_enabled=False,
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


def test_run_rejects_too_many_targets(cfg_and_db):
    """The /run endpoint must enforce max_input_targets (direct API callers
    bypass /upload, which is where the UI-level cap is applied)."""
    import base64

    cfg, _ = cfg_and_db
    cfg.max_input_targets = 2
    try:
        client = _client(cfg_and_db)
        tsv = (
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "1\t10\tT\tC\n"
            "1\t20\tA\tG\n"
            "1\t30\tT\tC\n"
        )
        res = client.post(
            "/run",
            data={
                "raw_b64": base64.b64encode(tsv.encode()).decode(),
                "base_editor": "ABE8e",
                "mode": "spreadsheet",
            },
        )
        assert res.status_code == 400, res.text
        assert "Too many targets" in res.json()["error"]
    finally:
        cfg.max_input_targets = 1000


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


def test_upload_csv_normalizes_to_tsv(cfg_and_db):
    """Uploading a CSV must yield a TSV-normalized raw_b64 so /check and
    /run (which always parse raw_b64 as TSV) can read the columns."""
    client = _client(cfg_and_db)
    csv_bytes = (
        b"chromosome,position,target_base,desired_base\n"
        b"chr1,100,A,G\n"
        b"chr2,200,T,C\n"
    )
    res = client.post(
        "/upload",
        files={"file": ("targets.csv", io.BytesIO(csv_bytes), "text/csv")},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["row_count"] == 2

    decoded = base64.b64decode(body["raw_b64"]).decode()
    # Must be tab-separated with the 4 standard columns (not the original CSV).
    assert decoded.startswith("chromosome\tposition\ttarget_base\tdesired_base\n")
    assert "," not in decoded.splitlines()[1]
    # /decode_targets style round-trip must succeed.
    df = pd.read_csv(io.StringIO(decoded), sep="\t")
    assert list(df.columns) == ["chromosome", "position", "target_base", "desired_base"]
    assert df.iloc[0]["target_base"] == "A"
    assert df.iloc[1]["desired_base"] == "C"


def test_check_works_with_csv_upload(cfg_and_db):
    """Full flow: upload a CSV, then post the returned raw_b64 to /check.
    Regression test for the KeyError('target_base') bug where the original
    comma-separated bytes were forwarded instead of being re-encoded as TSV.
    """
    client = _client(cfg_and_db)
    csv_bytes = (
        b"chromosome,position,target_base,desired_base\n"
        b"chr1,100,A,G\n"
        b"chr2,200,A,G\n"
        b"chr3,300,T,C\n"
    )
    up = client.post(
        "/upload",
        files={"file": ("targets.csv", io.BytesIO(csv_bytes), "text/csv")},
    ).json()
    assert up["ok"] is True

    res = client.post(
        "/check",
        data={"raw_b64": up["raw_b64"], "base_editor": "ABE8e"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["total"] == 3
    # ABE8e supports A->G plus its antisense complement T->C.
    assert body["matched"] == 3
    assert body["mismatched"] == 0
    assert body["can_run"] is True


def test_run_accepts_csv_upload(cfg_and_db):
    """Full flow: upload a CSV and use the returned raw_b64 to start a job.
    Regression test for the same /upload -> /run KeyError.
    """
    import beguided.web.scheduler as sched_mod

    sched_mod.execute_job = lambda **kw: {
        "job_id": kw["job_id"], "status": "complete",
        "message": "ok", "records": [],
    }

    client = _client(cfg_and_db)
    csv_bytes = (
        b"chromosome,position,target_base,desired_base\n"
        b"chr1,100,A,G\n"
        b"chr2,200,A,G\n"
    )
    up = client.post(
        "/upload",
        files={"file": ("targets.csv", io.BytesIO(csv_bytes), "text/csv")},
    ).json()
    assert up["ok"] is True

    res = client.post(
        "/run",
        data={"raw_b64": up["raw_b64"], "base_editor": "ABE8e"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert "job_id" in body


def test_upload_tsv_returns_tsv_raw_b64(cfg_and_db):
    """Regression: a TSV upload should still produce a usable raw_b64 (TSV
    bytes, not the empty/garbled output of the bugged path)."""
    client = _client(cfg_and_db)
    tsv_bytes = b"chromosome\tposition\ttarget_base\tdesired_base\nchr1\t100\tA\tG\n"
    res = client.post(
        "/upload",
        files={"file": ("targets.tsv", io.BytesIO(tsv_bytes), "text/tab-separated-values")},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    decoded = base64.b64decode(body["raw_b64"]).decode()
    df = pd.read_csv(io.StringIO(decoded), sep="\t")
    assert list(df.columns) == ["chromosome", "position", "target_base", "desired_base"]
    assert df.iloc[0]["target_base"] == "A"


def test_upload_sara_csv_runs_check(cfg_and_db):
    """Exact regression test for the user's reported file: sara_targetlist.csv
    (3 A>G + 3 T>C). With the only configured profile being ABE8e, /check
    must return 200 and classify the rows by editing direction (not blow up
    with KeyError on the CSV-vs-TSV mismatch).
    """
    repo_root = Path(__file__).resolve().parents[2]
    csv_path = repo_root / "sara_targetlist.csv"
    if not csv_path.exists():
        pytest.skip(f"{csv_path} not present in this checkout")
    csv_bytes = csv_path.read_bytes()

    client = _client(cfg_and_db)
    up = client.post(
        "/upload",
        files={"file": ("sara_targetlist.csv", io.BytesIO(csv_bytes), "text/csv")},
    ).json()
    assert up["ok"] is True
    assert up["row_count"] == 6

    res = client.post(
        "/check",
        data={"raw_b64": up["raw_b64"], "base_editor": "ABE8e"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["total"] == 6
    # ABE8e supports A->G plus its antisense complement T->C, so all 6 rows
    # are compatible.
    assert body["matched"] == 6
    assert body["mismatched"] == 0


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


def test_config_endpoint_exposes_genomes(cfg_and_db):
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    body = client.get("/config").json()
    keys = {g["key"] for g in body["genomes"]}
    assert {"human", "mouse"} <= keys
    assert body["default_genome"] == cfg.default_genome
    human = next(g for g in body["genomes"] if g["key"] == "human")
    # The fixture points human at /nonexistent, so it is not available.
    assert human["available"] is False


def test_config_endpoint_exposes_scheduler_knobs(cfg_and_db):
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    body = client.get("/config").json()
    assert body["max_input_targets"] == cfg.max_input_targets
    assert body["workers"] == cfg.workers


def test_favicon_redirects_to_svg(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.get("/favicon.ico", follow_redirects=False)
    assert res.status_code == 307
    assert res.headers["location"] == "/static/favicon.svg"


def test_reference_status_endpoint_merges_genome_and_gtf(cfg_and_db):
    from beguided.web import genome_cache, gtf_cache

    genome_cache._reset_for_test()
    gtf_cache._reset_for_test()
    client = _client(cfg_and_db)
    res = client.get("/reference/status")
    assert res.status_code == 200
    body = res.json()
    assert body["ok"] is True
    assert body["status"] == "not_loaded"
    assert isinstance(body["genomes"], list)


def test_genome_status_endpoint_supports_key(cfg_and_db):
    from beguided.web import genome_cache

    genome_cache._reset_for_test()
    client = _client(cfg_and_db)
    body = client.get("/genome/status", params={"genome": "mouse"}).json()
    assert body["key"] == "mouse"
    assert body["status"] == "not_loaded"


def test_upload_accepts_genome_field(cfg_and_db):
    client = _client(cfg_and_db)
    tsv_bytes = b"chromosome\tposition\ttarget_base\tdesired_base\nchr1\t100\tA\tG\n"
    res = client.post(
        "/upload",
        files={"file": ("targets.tsv", io.BytesIO(tsv_bytes), "text/tab-separated-values")},
        data={"genome": "human"},
    )
    assert res.status_code == 200, res.text
    assert res.json()["ok"] is True


def test_run_rejects_unknown_genome(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
            "genome": "bogus",
        },
    )
    assert res.status_code == 400
    assert "Unknown genome" in res.json()["error"]


def test_run_forwards_genome_key_to_scheduler(cfg_and_db, monkeypatch):
    """The selected genome key must be carried in the job payload."""
    import beguided.web.scheduler as sched_mod

    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": [],
        }

    monkeypatch.setattr(sched_mod, "execute_job", fake_execute_job)
    client = _client(cfg_and_db)
    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
            "genome": "human",
        },
    )
    assert res.status_code == 200, res.text
    job_id = res.json()["job_id"]
    for _ in range(40):
        s = client.get("/status", params={"job_id": job_id}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            break
        time.sleep(0.1)
    assert captured.get("genome_key") == "human"


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
        "chr1\t10\t0.123456789\tAAAAAAAAAAAAAAAAAAAA\tsynonymous\n"
    )
    records = _df_to_records(tsv)
    assert len(records) == 1
    r = records[0]
    assert r["chromosome"] == "chr1"
    assert r["position"] == 10
    assert r["guide_score"] == round(0.123456789, 3)
    assert r["guide_rna_sequence"] == "AAAAAAAAAAAAAAAAAAAA"
    assert r["bystander_effect"] == "synonymous"


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
            "filter_mode": "synonymous",
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
        {"chromosome": "chr1", "position": 50, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "none",
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
        {"chromosome": "chr1", "position": 50, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "none",
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

    # synonymous: keep bystander_effect in {synonymous, none}. The "coding"
    # chr1:80 target loses its only guide, so the fill keeps one placeholder
    # row for it.
    res_syn = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "filter_mode": "synonymous"},
    )
    assert res_syn.status_code == 200
    rows = [l.split("\t") for l in res_syn.text.strip().split("\n")[1:]]
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50, 80]
    assert "**Guides removed by filter**" in res_syn.text

    # guide_score_threshold=0.5 -> drop the 0.2 row (placeheld for chr1:80).
    res_score = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "guide_score_threshold": 0.5},
    )
    assert res_score.status_code == 200
    rows = [l.split("\t") for l in res_score.text.strip().split("\n")[1:]]
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50, 80]

    # Combined: synonymous AND score >= 0.8 -> only the 0.9 synonymous chr1:50 row is
    # a real guide; chr2:30 and chr1:80 fall back to placeholder rows.
    res_combo = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "synonymous",
            "guide_score_threshold": 0.8,
        },
    )
    assert res_combo.status_code == 200
    rows = [l.split("\t") for l in res_combo.text.strip().split("\n")[1:]]
    assert len(rows) == 3
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50, 80]
    real = next(r for r in rows if int(r[1]) == 50)
    assert "Guides removed" not in "\t".join(real)
    assert "**No guide possible**" not in "\t".join(real)

    # Invalid filter_mode -> 400.
    bad = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "filter_mode": "bogus"},
    )
    assert bad.status_code == 400


def test_download_results_applies_trinuc_filter(cfg_and_db):
    """The /download/results endpoint should respect the
    trinuc_edit_score_min query param, keeping NaN (CBE) rows."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "trinuc_edit_score": 0.9},
        {"chromosome": "chr1", "position": 80, "guide_score": 0.9,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC",
         "trinuc_edit_score": 0.1},
        {"chromosome": "chr2", "position": 30, "guide_score": 0.9,
         "guide_rna_sequence": "GGGGGGGGGGGGGGGGGGGG",
         "trinuc_edit_score": None},  # CBE row — NaN, always kept
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
    assert len(lines) == 1 + 3

    # trinuc_edit_score_min=0.5 -> drop the 0.1 row; chr1:80 falls back to a
    # placeholder, keeping 0.9 and NaN (CBE).
    res_filtered = client.get(
        "/download/results",
        params={"job_id": submit["job_id"], "trinuc_edit_score_min": 0.5},
    )
    assert res_filtered.status_code == 200
    rows = [l.split("\t") for l in res_filtered.text.strip().split("\n")[1:]]
    positions = sorted(int(r[1]) for r in rows)
    assert positions == [30, 50, 80]
    assert "**Guides removed by filter**" in res_filtered.text


def test_download_results_filename_reflects_filters(cfg_and_db):
    """The Content-Disposition filename should encode the active filters."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "synonymous",
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
            "filter_mode": "synonymous",
            "guide_score_threshold": 0.5,
        },
    )
    assert "beguided_results_synonymous_s0.5.tsv" in filtered.headers["content-disposition"]


def test_download_results_xlsx(cfg_and_db):
    """The /download/results endpoint with format=xlsx should serve an XLSX
    spreadsheet that round-trips back to the same DataFrame."""
    import beguided.web.scheduler as sched_mod

    records = [
        {"chromosome": "chr1", "position": 50, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5,7", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC"},
        {"chromosome": "chr2", "position": 30, "bystander_effect": "none",
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

    # synonymous filter applied -> 2 guides + 1 placeholder for the "coding" target.
    res_syn = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "synonymous",
            "format": "xlsx",
        },
    )
    assert res_syn.status_code == 200
    df_syn = pd.read_excel(io.BytesIO(res_syn.content), engine="openpyxl")
    assert sorted(int(p) for p in df_syn["position"]) == [30, 50, 80]
    assert (df_syn["guide_rna_sequence"] == "**Guides removed by filter**").sum() == 1

    # Filename reflects filters.
    filtered = client.get(
        "/download/results",
        params={
            "job_id": submit["job_id"],
            "filter_mode": "synonymous",
            "guide_score_threshold": 0.5,
            "format": "xlsx",
        },
    )
    assert (
        "beguided_results_synonymous_s0.5.xlsx" in filtered.headers["content-disposition"]
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


def _submit_with_records(cfg_and_db, records):
    """Submit a /run job whose fake executor returns ``records`` and wait
    until the scheduler reports it complete. Returns the job_id."""
    import beguided.web.scheduler as sched_mod

    def fake_execute_job(**kwargs):
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": records,
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
    assert s["status"] == "complete"
    return submit["job_id"], client


def test_download_results_condense_merges_identical_rows(cfg_and_db):
    """condense=1 should collapse rows that differ only in base_editor
    into a single row whose base_editor cell is newline-separated, and
    the download filename should reflect the condense flag."""
    records = [
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE9e"},
        {"chromosome": "chr2", "position": 200, "bystander_effect": "none",
         "bystander_positions": "", "guide_score": 0.6,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    res = client.get(
        "/download/results",
        params={"job_id": job_id, "condense": "1"},
    )
    assert res.status_code == 200
    assert "beguided_results_condensed.tsv" in res.headers["content-disposition"]

    # Read with pandas so the embedded newline inside the quoted
    # base_editor cell is reconstructed cleanly.
    df = pd.read_csv(io.StringIO(res.text), sep="\t")
    assert len(df) == 2
    chr1 = df[df["chromosome"] == "chr1"].iloc[0]
    assert chr1["base_editor"].split("\n") == ["ABE8e", "ABE9e"]
    chr2 = df[df["chromosome"] == "chr2"].iloc[0]
    assert chr2["base_editor"] == "ABE8e"


def test_download_results_condense_keeps_different_guides_separate(cfg_and_db):
    """Rows whose guide_rna_sequence differs must NOT be merged, even
    when every other column matches."""
    records = [
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "TTTTTTTTTTTTTTTTTTTT",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    res = client.get(
        "/download/results",
        params={"job_id": job_id, "condense": "1"},
    )
    assert res.status_code == 200
    rows = [l.split("\t") for l in res.text.strip().split("\n")[1:]]
    assert len(rows) == 2  # not merged
    for r in rows:
        assert r[0] == "chr1"


def test_download_results_condense_off_keeps_every_row(cfg_and_db):
    """Without condense=1 the download must still return every row, even
    when several are identical except for base_editor."""
    records = [
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE9e"},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    res = client.get(
        "/download/results",
        params={"job_id": job_id},
    )
    assert res.status_code == 200
    assert "beguided_results.tsv" in res.headers["content-disposition"]
    rows = [l.split("\t") for l in res.text.strip().split("\n")[1:]]
    assert len(rows) == 2  # no merging without condense


def test_download_results_condense_combines_with_filter_and_score(cfg_and_db):
    """condense should compose with filter_mode and guide_score_threshold;
    the resulting filename should encode every active flag."""
    records = [
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE9e"},
        {"chromosome": "chr1", "position": 80, "bystander_effect": "coding",
         "bystander_positions": "5", "guide_score": 0.2,
         "guide_rna_sequence": "CCCCCCCCCCCCCCCCCCCC",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    res = client.get(
        "/download/results",
        params={
            "job_id": job_id,
            "filter_mode": "synonymous",
            "guide_score_threshold": 0.5,
            "condense": "1",
        },
    )
    assert res.status_code == 200
    cd = res.headers["content-disposition"]
    # Filename should mention both filters and condense (order is implementation
    # detail — check both tokens are present).
    assert "condensed" in cd
    assert "synonymous_s0.5" in cd

    df = pd.read_csv(io.StringIO(res.text), sep="\t")
    # After filter_mode=synonymous + score>=0.5 only the chr1:100 group remains
    # as a real guide (condensed into a single row); the "coding" chr1:80
    # target falls back to a "**Guides removed by filter**" placeholder row.
    assert len(df) == 2
    real = df[df["base_editor"].notna()].iloc[0]
    assert real["base_editor"].split("\n") == ["ABE8e", "ABE9e"]
    ph = df[df["guide_rna_sequence"] == "**Guides removed by filter**"].iloc[0]
    assert ph["position"] == 80
    assert pd.isna(ph["base_editor"])


def test_download_results_condense_xlsx(cfg_and_db):
    """condense=1 should also work for XLSX downloads."""
    records = [
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE8e"},
        {"chromosome": "chr1", "position": 100, "bystander_effect": "synonymous",
         "bystander_positions": "", "guide_score": 0.9,
         "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "pam_sequence": "NGG", "base_editor": "ABE9e"},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    res = client.get(
        "/download/results",
        params={"job_id": job_id, "condense": "1", "format": "xlsx"},
    )
    assert res.status_code == 200
    assert "beguided_results_condensed.xlsx" in res.headers["content-disposition"]
    df = pd.read_excel(io.BytesIO(res.content), engine="openpyxl")
    assert len(df) == 1
    assert df.iloc[0]["base_editor"].split("\n") == ["ABE8e", "ABE9e"]


def test_download_results_preserves_pipeline_placeholder_markers(cfg_and_db):
    """Pipeline-produced placeholder rows ("**No guide possible**" /
    "**Base change incompatible**") survive the dynamic download filters with
    their marker intact, so the web UI and the download agree on one row per
    target with no guide."""
    records = [
        # Target that could not produce a guide (compatible, no PAM).
        {"chromosome": "chr1", "position": 10, "target_base": "A",
         "desired_base": "G", "guide_rna_sequence": "**No guide possible**"},
        # Target whose mutation type does not fit the selected editors.
        {"chromosome": "chr2", "position": 20, "target_base": "C",
         "desired_base": "T", "guide_rna_sequence": "**Base change incompatible**"},
        # A real guide for a third target.
        {"chromosome": "chr3", "position": 30, "target_base": "A",
         "desired_base": "G", "guide_rna_sequence": "AAAAAAAAAAAAAAAAAAAA",
         "bystander_effect": "synonymous", "bystander_positions": "",
         "guide_score": 0.9},
    ]
    job_id, client = _submit_with_records(cfg_and_db, records)

    # Unfiltered download: all three rows present, markers intact.
    res = client.get(
        "/download/results",
        params={"job_id": job_id},
    )
    assert res.status_code == 200
    df = pd.read_csv(io.StringIO(res.text), sep="\t")
    assert len(df) == 3
    assert set(df["guide_rna_sequence"]) == {
        "**No guide possible**",
        "**Base change incompatible**",
        "AAAAAAAAAAAAAAAAAAAA",
    }

    # "synonymous" filter keeps the real guide, the no-guide placeholders survive.
    res_syn = client.get(
        "/download/results",
        params={"job_id": job_id, "filter_mode": "synonymous"},
    )
    assert res_syn.status_code == 200
    df_syn = pd.read_csv(io.StringIO(res_syn.text), sep="\t")
    assert len(df_syn) == 3
    markers = set(df_syn["guide_rna_sequence"])
    assert "**No guide possible**" in markers
    assert "**Base change incompatible**" in markers
    assert "AAAAAAAAAAAAAAAAAAAA" in markers

