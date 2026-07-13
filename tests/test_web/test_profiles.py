"""Tests for the custom base editor profiles upload + selection flow."""

from __future__ import annotations

import base64
import io

import pytest


pytest.importorskip("fastapi")
pytest.importorskip("uvicorn")


HEADER = (
    "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
    "high_edit_start\thigh_edit_end\tpam\n"
)
SAMPLE_PROFILES = (
    HEADER
    + "CustomABE1\tABE\t2\t13\t3\t8\tNGG\n"
    + "CustomCBE1\tCBE\t2\t10\t6\t6\tNGG\n"
)


@pytest.fixture
def cfg_and_db(tmp_path, monkeypatch):
    """Set up a WebConfig + DB on a tmp path; restart the scheduler singleton."""
    from beguided.web import config as cfg_mod
    from beguided.web import scheduler as sched_mod
    from beguided.web.config import WebConfig, set_config
    from beguided.web.profiles_store import get_store
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
        max_input_targets=1000,
        db_path=str(db),
        queue_max=3,
        cooldown_seconds=0,
        max_upload_mb=5,
    ))
    get_store().clear()
    sch = Scheduler(db_path=str(db), queue_max=3, cooldown_seconds=0)
    install(sch)
    sch.start()
    yield cfg_mod.get_config(), sch
    get_store().clear()


def _client(cfg_and_db):
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app

    return TestClient(create_app())


def _targets_tsv_b64(rows=("chr1\t50\tT\tC\n",)) -> str:
    content = "".join((
        "chromosome\tposition\ttarget_base\tdesired_base\n",
        *rows,
    )).encode()
    return base64.b64encode(content).decode("ascii")


def _wait_for_status(client, job_id, timeout=10.0):
    import time

    deadline = time.time() + timeout
    while time.time() < deadline:
        s = client.get("/status", params={"job_id": job_id}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            return s["status"]
        time.sleep(0.05)
    raise AssertionError("job did not finish in time")


def test_upload_valid_profiles(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("profiles.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["count"] == 2
    assert body["names"] == ["CustomABE1", "CustomCBE1"]
    assert body["key"] == "custom"
    assert body["label"] == "Custom (uploaded)"


def test_upload_missing_required_columns(cfg_and_db):
    client = _client(cfg_and_db)
    bad = "name\tediting\nABE8e\tABE\n"
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(bad.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 400
    assert "Missing required columns" in res.json()["error"]


def test_upload_duplicate_names_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    bad = HEADER + "Dup1\tABE\t2\t13\t3\t8\tNGG\n" + "Dup1\tCBE\t2\t10\t6\t6\tNGG\n"
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(bad.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 400
    assert "Invalid profiles" in res.json()["error"]


def test_upload_too_many_profiles(cfg_and_db):
    client = _client(cfg_and_db)
    bad = HEADER + "".join(
        f"BE{i}\tABE\t2\t13\t3\t8\tNGG\n" for i in range(11)
    )
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(bad.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 400
    assert "Too many profiles" in res.json()["error"]


def test_upload_pam_too_long(cfg_and_db):
    client = _client(cfg_and_db)
    bad = HEADER + "Custom1\tABE\t2\t13\t3\t8\tNNNNNNN\n"
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(bad.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 400
    assert "PAM too long" in res.json()["error"]


def test_upload_pam_exactly_six_bases_accepted(cfg_and_db):
    client = _client(cfg_and_db)
    ok = HEADER + "Custom1\tABE\t2\t13\t3\t8\tNNNNNN\n"
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(ok.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 200, res.text
    assert res.json()["ok"] is True


def test_upload_extra_columns_dropped(cfg_and_db):
    client = _client(cfg_and_db)
    body = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\textra\n"
        "Custom1\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\tdropped\n"
    )
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(body.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 200, res.text
    assert res.json()["ok"] is True

    table = client.get("/base-editor-table?source=custom").json()
    assert table["ok"] is True
    row_keys = set(table["rows"][0].keys())
    assert "extra" not in row_keys
    assert row_keys == {"name", "editing", "general", "high", "pam", "cas", "source"}


def test_max_upload_mb_is_one_by_default(cfg_and_db):
    cfg, _ = cfg_and_db
    assert cfg.max_profile_upload_mb == 1


def test_download_serves_profiles_template(cfg_and_db):
    """The /base-editor-profiles/download endpoint serves a template TSV
    containing a placeholder ``example`` row plus a real example row
    (ABE8e). The required schema columns are present, plus the optional
    ``cas`` column; the ``source`` column is intentionally omitted."""
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    res = client.get("/base-editor-profiles/download")
    assert res.status_code == 200
    assert res.headers["content-type"].startswith("text/tab-separated-values")
    assert "template" in res.headers["content-disposition"].lower()
    # Required header columns + the optional ``cas`` column.
    header_line = res.text.splitlines()[0]
    for col in (
        "name", "editing", "general_edit_start", "general_edit_end",
        "high_edit_start", "high_edit_end", "pam", "cas",
    ):
        assert col in header_line
    # The ``source`` column is intentionally not in the template.
    assert "source" not in header_line
    # Body rows: one ``example`` placeholder and one ABE8e example.
    assert "example" in res.text
    assert "ABE8e" in res.text


def test_upload_rejects_name_collision_with_default(cfg_and_db):
    """Uploading a custom profile whose name already exists in the default
    profiles file must be rejected with a clear collision error."""
    client = _client(cfg_and_db)
    # cfg_and_db's default file already contains ABE8e
    colliding = HEADER + "ABE8e\tABE\t2\t13\t3\t8\tNGG\n"
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(colliding.encode()), "text/tab-separated-values")},
    )
    assert res.status_code == 400
    assert "Name collision" in res.json()["error"]
    assert "ABE8e" in res.json()["error"]


def test_base_editor_table_combined_source(cfg_and_db):
    """source=combined returns default + custom profiles with a source tag."""
    client = _client(cfg_and_db)
    # Initial state: no custom profiles uploaded; combined == default.
    res = client.get("/base-editor-table?source=combined").json()
    assert res["ok"] is True
    assert all(row["source"] == "default" for row in res["rows"])
    default_names = [row["name"] for row in res["rows"]]

    # Upload custom profiles; combined now contains both.
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    res = client.get("/base-editor-table?source=combined").json()
    assert res["ok"] is True
    sources = {row["name"]: row["source"] for row in res["rows"]}
    # Custom profiles are listed first.
    assert sources["CustomABE1"] == "custom"
    assert sources["CustomCBE1"] == "custom"
    assert sources["ABE8e"] == "default"
    assert len(res["rows"]) == len(default_names) + 2


def test_status_initial_loaded_false(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.get("/base-editor-profiles/status")
    assert res.status_code == 200
    body = res.json()
    assert body["ok"] is True
    assert body["loaded"] is False
    assert body["names"] == []


def test_status_reflects_upload(cfg_and_db):
    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    res = client.get("/base-editor-profiles/status")
    body = res.json()
    assert body["loaded"] is True
    assert body["names"] == ["CustomABE1", "CustomCBE1"]


def test_check_uses_custom_profiles(cfg_and_db):
    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    res = client.post(
        "/check",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "custom"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["matched"] == 1
    assert body["can_run"] is True


def test_check_unknown_base_editor_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/check",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "totally-fake"},
    )
    assert res.status_code == 200
    assert res.json()["ok"] is False


def test_check_custom_without_upload_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/check",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "custom"},
    )
    assert res.status_code == 200
    assert res.json()["ok"] is False
    assert "No custom profiles" in res.json()["error"]


def test_run_payload_with_custom_profiles(cfg_and_db, monkeypatch):
    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": [],
        }

    import beguided.web.scheduler as sched_mod
    sched_mod.execute_job = fake_execute_job

    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )

    res = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "custom"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert "job_id" in body

    status = _wait_for_status(client, body["job_id"])
    assert status == "complete"

    assert captured["use_custom_profiles"] is True
    assert captured["custom_profiles_b64"] is not None
    decoded = base64.b64decode(captured["custom_profiles_b64"])
    assert b"CustomABE1" in decoded
    assert captured["base_editor_names"] == ["CustomABE1", "CustomCBE1"]


def test_run_default_after_upload_does_not_use_custom(cfg_and_db, monkeypatch):
    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": [],
        }

    import beguided.web.scheduler as sched_mod
    sched_mod.execute_job = fake_execute_job

    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )

    res = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "ABE8e"},
    )
    assert res.status_code == 200, res.text
    body = res.json()
    status = _wait_for_status(client, body["job_id"])
    assert status == "complete"

    assert captured["use_custom_profiles"] is False
    assert captured["custom_profiles_b64"] is None
    assert captured["base_editor_names"] == ["ABE8e"]


def test_run_custom_without_upload_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/run",
        data={"raw_b64": _targets_tsv_b64(), "base_editor": "custom"},
    )
    assert res.status_code == 400
    assert "No custom profiles" in res.json()["error"]


def test_store_lives_in_memory_only(cfg_and_db, tmp_path):
    from beguided.web.profiles_store import get_store

    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    assert get_store().is_loaded() is True
    assert get_store().get_bytes() == SAMPLE_PROFILES.encode()

    profile_paths = [p for p in tmp_path.rglob("*profiles*")]
    assert profile_paths == [], f"upload wrote unexpected files: {profile_paths}"

    get_store().clear()
    assert get_store().is_loaded() is False
    assert client.get("/base-editor-profiles/status").json()["loaded"] is False


def test_build_config_writes_custom_profiles_to_workdir(cfg_and_db, tmp_path):
    import base64 as _b64

    cfg, _ = cfg_and_db
    from beguided.web.tasks import _build_config
    import pandas as pd

    work_dir = tmp_path / "custom_workdir"
    work_dir.mkdir()

    targets = pd.DataFrame({
        "chromosome": ["chr1"],
        "position": [50],
        "target_base": ["T"],
        "desired_base": ["C"],
    })
    custom_b64 = _b64.b64encode(SAMPLE_PROFILES.encode()).decode("ascii")
    config = _build_config(
        targets,
        base_editor_names=["CustomABE1", "CustomCBE1"],
        filter_mode="all",
        work_dir=str(work_dir),
        use_custom_profiles=True,
        custom_profiles_b64=custom_b64,
    )
    expected = str(work_dir / "custom_profiles.tsv")
    assert config["input"]["base_editor_profiles"] == expected
    assert config["use_base_editors"] == "all"

    from pathlib import Path
    written = Path(expected)
    assert written.exists()
    assert "CustomABE1" in written.read_text()


def test_check_accepts_comma_separated_editor_names(cfg_and_db):
    """Multi-select UI sends comma-separated individual editor names."""
    from pathlib import Path
    from beguided.web.profiles_store import get_store

    cfg, _ = cfg_and_db
    # Inject a second preset so the resolver has more than one editor to pick.
    presets = Path(cfg.base_editor_profiles_file)
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
        "CBE6d\tCBE\t2\t10\t6\t6\tNGG\tSpCas9\tNature2023\n"
    )
    get_store().clear()

    client = _client(cfg_and_db)
    res = client.post(
        "/check",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e,CBE6d",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["can_run"] is True


def test_check_comma_separated_rejects_unknown_name(cfg_and_db):
    """Unknown editor names in a comma-separated list are rejected."""
    from beguided.web.profiles_store import get_store
    get_store().clear()

    client = _client(cfg_and_db)
    res = client.post(
        "/check",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e,totally-fake",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is False
    assert "totally-fake" in body["error"]


def test_run_resolves_comma_separated_editor_names(cfg_and_db, monkeypatch):
    """Comma-separated names flow through to the pipeline as a list."""
    from pathlib import Path
    from beguided.web.profiles_store import get_store

    cfg, _ = cfg_and_db
    presets = Path(cfg.base_editor_profiles_file)
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
        "CBE6d\tCBE\t2\t10\t6\t6\tNGG\tSpCas9\tNature2023\n"
    )
    get_store().clear()

    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": [],
        }

    import beguided.web.scheduler as sched_mod
    sched_mod.execute_job = fake_execute_job

    client = _client(cfg_and_db)
    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "CBE6d,ABE8e",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    status = _wait_for_status(client, body["job_id"])
    assert status == "complete"

    assert captured["use_custom_profiles"] is False
    assert sorted(captured["base_editor_names"]) == ["ABE8e", "CBE6d"]


def test_check_single_name_still_works_as_backward_compat(cfg_and_db):
    """A single editor name (no commas) is still accepted as a 1-element list."""
    from beguided.web.profiles_store import get_store
    get_store().clear()

    client = _client(cfg_and_db)
    res = client.post(
        "/check",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["can_run"] is True


def test_check_mixed_custom_and_default_names(cfg_and_db):
    """Comma-separated names spanning custom + default sources resolve
    correctly via a merged profiles dict."""
    from beguided.web.profiles_store import get_store
    get_store().clear()

    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    # Default fixture has ABE8e (ABE); SAMPLE_PROFILES has CustomABE1 + CustomCBE1.
    # Sending a mix should succeed and match both editing directions.
    targets = (
        "chromosome\tposition\ttarget_base\tdesired_base\n"
        "chr1\t50\tA\tG\n"      # matched by ABE8e (sense A→G)
        "chr1\t100\tC\tT\n"     # matched by CustomCBE1 (sense C→T)
    )
    raw_b64 = base64.b64encode(targets.encode()).decode("ascii")
    res = client.post(
        "/check",
        data={
            "raw_b64": raw_b64,
            "base_editor": "ABE8e,CustomCBE1",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["can_run"] is True
    assert body["matched"] == 2


def test_run_mixed_custom_and_default_includes_merged_profiles_b64(
    cfg_and_db, monkeypatch
):
    """A run with mixed sources carries a merged_profiles_b64 payload."""
    from beguided.web.profiles_store import get_store
    get_store().clear()

    captured = {}

    def fake_execute_job(**kwargs):
        captured.update(kwargs)
        return {
            "job_id": kwargs["job_id"],
            "status": "complete",
            "message": "ok",
            "records": [],
        }

    import beguided.web.scheduler as sched_mod
    sched_mod.execute_job = fake_execute_job

    client = _client(cfg_and_db)
    client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(SAMPLE_PROFILES.encode()), "text/tab-separated-values")},
    )
    res = client.post(
        "/run",
        data={
            "raw_b64": _targets_tsv_b64(),
            "base_editor": "ABE8e,CustomABE1",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    status = _wait_for_status(client, body["job_id"])
    assert status == "complete"

    assert captured["use_custom_profiles"] is False
    assert captured["merged_profiles_b64"] is not None
    merged_bytes = base64.b64decode(captured["merged_profiles_b64"])
    assert b"ABE8e" in merged_bytes
    assert b"CustomABE1" in merged_bytes
    assert sorted(captured["base_editor_names"]) == ["ABE8e", "CustomABE1"]
