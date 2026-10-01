"""Tests for the custom base editor profiles upload + selection flow.

The web upload accepts XLSX (preferred, matches the downloaded template)
and TSV for backwards compatibility. The server normalises both to TSV
in memory before storing, so all tests below use an XLSX fixture to
exercise the new code path while still asserting the on-disk format
the pipeline reads stays TSV.
"""

from __future__ import annotations

import base64
import io

import pandas as pd
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
SAMPLE_PROFILES_DF = pd.DataFrame(
    [
        {"name": "CustomABE1", "editing": "ABE",
         "general_edit_start": 2, "general_edit_end": 13,
         "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG"},
        {"name": "CustomCBE1", "editing": "CBE",
         "general_edit_start": 2, "general_edit_end": 10,
         "high_edit_start": 6, "high_edit_end": 6, "pam": "NGG"},
    ]
)


def _xlsx_bytes(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_excel(buf, index=False, engine="openpyxl")
    return buf.getvalue()


def _csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


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
        usage_stats_enabled=False,
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


def _post_profiles(client, name, payload, mime, df=None, suffix=None):
    """Upload helper: pass either raw ``payload`` bytes or a DataFrame."""
    if df is not None and suffix is not None:
        if suffix == ".xlsx":
            body = _xlsx_bytes(df)
        elif suffix == ".csv":
            body = _csv_bytes(df)
        else:
            body = df.to_csv(sep="\t", index=False).encode("utf-8")
        return client.post(
            "/base-editor-profiles/upload",
            files={"file": (name, io.BytesIO(body), mime)},
        )
    return client.post(
        "/base-editor-profiles/upload",
        files={"file": (name, io.BytesIO(payload), mime)},
    )


XLSX_MIME = (
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)
TSV_MIME = "text/tab-separated-values"
CSV_MIME = "text/csv"


def _wait_for_status(client, job_id, timeout=10.0):
    import time

    deadline = time.time() + timeout
    while time.time() < deadline:
        s = client.get("/status", params={"job_id": job_id}).json()
        if s["status"] in ("complete", "error", "cancelled"):
            return s["status"]
        time.sleep(0.05)
    raise AssertionError("job did not finish in time")


# ---------------------------------------------------------------------------
# Download endpoint — XLSX template
# ---------------------------------------------------------------------------

def test_download_serves_xlsx_template(cfg_and_db):
    """The /base-editor-profiles/download endpoint serves an XLSX
    template (matching the new preferred upload format). It must contain
    a placeholder ``example`` row and a real ``ABE8e`` example, plus the
    required schema columns plus the optional ``cas`` column. The
    ``source`` column is intentionally omitted."""
    cfg, _ = cfg_and_db
    client = _client(cfg_and_db)
    res = client.get("/base-editor-profiles/download")
    assert res.status_code == 200
    assert "spreadsheetml.sheet" in res.headers["content-type"]
    assert res.headers["content-disposition"].endswith("template.xlsx")
    df = pd.read_excel(io.BytesIO(res.content), engine="openpyxl")
    assert list(df.columns) == [
        "name", "editing", "general_edit_start", "general_edit_end",
        "high_edit_start", "high_edit_end", "pam", "cas",
    ]
    assert "source" not in df.columns
    assert set(df["name"]) == {"example", "ABE8e"}


# ---------------------------------------------------------------------------
# Upload — happy path XLSX
# ---------------------------------------------------------------------------

def test_upload_valid_xlsx(cfg_and_db):
    client = _client(cfg_and_db)
    res = _post_profiles(
        client, "profiles.xlsx",
        payload=None, mime=XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["ok"] is True
    assert body["count"] == 2
    assert body["names"] == ["CustomABE1", "CustomCBE1"]
    assert body["key"] == "custom"
    assert body["label"] == "Custom (uploaded)"


def test_upload_valid_tsv_still_accepted(cfg_and_db):
    """Backward compatibility: plain TSV uploads keep working."""
    client = _client(cfg_and_db)
    res = client.post(
        "/base-editor-profiles/upload",
        files={
            "file": (
                "profiles.tsv",
                io.BytesIO(SAMPLE_PROFILES.encode()),
                TSV_MIME,
            )
        },
    )
    assert res.status_code == 200, res.text
    assert res.json()["names"] == ["CustomABE1", "CustomCBE1"]


def test_upload_valid_csv_accepted(cfg_and_db):
    """CSV uploads are accepted defensively (mirrors /upload for targets)."""
    client = _client(cfg_and_db)
    res = _post_profiles(
        client, "profiles.csv",
        payload=None, mime=CSV_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".csv",
    )
    assert res.status_code == 200, res.text
    assert res.json()["names"] == ["CustomABE1", "CustomCBE1"]


def test_upload_xls_rejected_with_clear_message(cfg_and_db):
    """Legacy .xls files cannot be parsed by openpyxl and lack xlrd in
    the dependency set; the route must reject with a re-save hint."""
    client = _client(cfg_and_db)
    body = SAMPLE_PROFILES.encode()  # pretend TSV bytes
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("profiles.xls", io.BytesIO(body), "application/vnd.ms-excel")},
    )
    assert res.status_code == 400
    assert ".xlsx" in res.json()["error"]


def test_upload_empty_file_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("profiles.xlsx", io.BytesIO(b""), XLSX_MIME)},
    )
    assert res.status_code == 400
    assert "empty" in res.json()["error"].lower()


def test_upload_unknown_extension_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    body = SAMPLE_PROFILES.encode()
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("profiles.xlsx.bin", io.BytesIO(body), "application/octet-stream")},
    )
    assert res.status_code == 400
    assert "Unsupported file type" in res.json()["error"]


def test_upload_extension_mismatch_rejected(cfg_and_db):
    """An XLSX-magic-byte payload posted as ``.tsv`` must be rejected with
    a clear magic-byte mismatch message — same protection as /upload."""
    client = _client(cfg_and_db)
    xlsx = _xlsx_bytes(SAMPLE_PROFILES_DF)
    res = client.post(
        "/base-editor-profiles/upload",
        files={"file": ("p.tsv", io.BytesIO(xlsx), TSV_MIME)},
    )
    assert res.status_code == 400
    assert "do not match" in res.json()["error"]


# ---------------------------------------------------------------------------
# Upload — validation rules (originally TSV; now exercise them via XLSX)
# ---------------------------------------------------------------------------

def test_upload_missing_required_columns(cfg_and_db):
    client = _client(cfg_and_db)
    bad = pd.DataFrame([
        {"name": "ABE8e", "editing": "ABE"},
    ])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=bad, suffix=".xlsx",
    )
    assert res.status_code == 400
    assert "Missing required columns" in res.json()["error"]


def test_upload_duplicate_names_rejected(cfg_and_db):
    client = _client(cfg_and_db)
    bad = pd.DataFrame([
        {"name": "Dup1", "editing": "ABE",
         "general_edit_start": 2, "general_edit_end": 13,
         "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG"},
        {"name": "Dup1", "editing": "CBE",
         "general_edit_start": 2, "general_edit_end": 10,
         "high_edit_start": 6, "high_edit_end": 6, "pam": "NGG"},
    ])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=bad, suffix=".xlsx",
    )
    assert res.status_code == 400
    assert "Invalid profiles" in res.json()["error"]


def test_upload_too_many_profiles(cfg_and_db):
    client = _client(cfg_and_db)
    bad = pd.DataFrame([
        {"name": f"BE{i}", "editing": "ABE",
         "general_edit_start": 2, "general_edit_end": 13,
         "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG"}
        for i in range(11)
    ])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=bad, suffix=".xlsx",
    )
    assert res.status_code == 400
    assert "Too many profiles" in res.json()["error"]


def test_upload_pam_too_long(cfg_and_db):
    client = _client(cfg_and_db)
    bad = pd.DataFrame([{
        "name": "Custom1", "editing": "ABE",
        "general_edit_start": 2, "general_edit_end": 13,
        "high_edit_start": 3, "high_edit_end": 8, "pam": "NNNNNNN",
    }])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=bad, suffix=".xlsx",
    )
    assert res.status_code == 400
    assert "PAM too long" in res.json()["error"]


def test_upload_pam_exactly_six_bases_accepted(cfg_and_db):
    client = _client(cfg_and_db)
    ok = pd.DataFrame([{
        "name": "Custom1", "editing": "ABE",
        "general_edit_start": 2, "general_edit_end": 13,
        "high_edit_start": 3, "high_edit_end": 8, "pam": "NNNNNN",
    }])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=ok, suffix=".xlsx",
    )
    assert res.status_code == 200, res.text
    assert res.json()["ok"] is True


def test_upload_extra_columns_dropped(cfg_and_db):
    client = _client(cfg_and_db)
    df = pd.DataFrame([{
        "name": "Custom1", "editing": "ABE",
        "general_edit_start": 2, "general_edit_end": 13,
        "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG",
        "cas": "SpCas9", "source": "Nature2022", "extra": "dropped",
    }])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=df, suffix=".xlsx",
    )
    assert res.status_code == 200, res.text
    assert res.json()["ok"] is True

    table = client.get("/base-editor-table?source=custom").json()
    assert table["ok"] is True
    row_keys = set(table["rows"][0].keys())
    assert "extra" not in row_keys
    assert row_keys == {"name", "editing", "general", "high", "pam", "cas", "origin", "source"}


def test_max_upload_mb_is_one_by_default(cfg_and_db):
    cfg, _ = cfg_and_db
    assert cfg.max_profile_upload_mb == 1


def test_upload_rejects_name_collision_with_default(cfg_and_db):
    """Uploading a custom profile whose name already exists in the default
    profiles file must be rejected with a clear collision error."""
    client = _client(cfg_and_db)
    # cfg_and_db's default file already contains ABE8e
    df = pd.DataFrame([{
        "name": "ABE8e", "editing": "ABE",
        "general_edit_start": 2, "general_edit_end": 13,
        "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG",
    }])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=df, suffix=".xlsx",
    )
    assert res.status_code == 400
    assert "Name collision" in res.json()["error"]
    assert "ABE8e" in res.json()["error"]


# ---------------------------------------------------------------------------
# Base-editor table endpoint
# ---------------------------------------------------------------------------

def test_base_editor_table_combined_source(cfg_and_db):
    """source=combined returns default + custom profiles with an origin tag."""
    client = _client(cfg_and_db)
    # Initial state: no custom profiles uploaded; combined == default.
    res = client.get("/base-editor-table?source=combined").json()
    assert res["ok"] is True
    assert all(row["origin"] == "default" for row in res["rows"])
    default_names = [row["name"] for row in res["rows"]]

    # Upload custom profiles; combined now contains both.
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    res = client.get("/base-editor-table?source=combined").json()
    assert res["ok"] is True
    origins = {row["name"]: row["origin"] for row in res["rows"]}
    # Custom profiles are listed first.
    assert origins["CustomABE1"] == "custom"
    assert origins["CustomCBE1"] == "custom"
    assert origins["ABE8e"] == "default"
    assert len(res["rows"]) == len(default_names) + 2


def test_base_editor_table_parses_markdown_link(tmp_path, cfg_and_db):
    """Default profiles with a ``[text](url)`` source cell must surface a
    parsed ``{text, url}`` payload on every row."""
    from pathlib import Path
    from beguided.web.config import get_config

    cfg = get_config()
    presets = Path(cfg.base_editor_profiles_file)
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\t"
        "[Richter et al, 2020](https://doi.org/10.1038/s41587-020-0453-z)\n"
    )

    client = _client(cfg_and_db)
    res = client.get("/base-editor-table?source=default").json()
    assert res["ok"] is True
    rows = {row["name"]: row for row in res["rows"]}
    assert rows["ABE8e"]["source"] == {
        "text": "Richter et al, 2020",
        "url": "https://doi.org/10.1038/s41587-020-0453-z",
    }
    assert rows["ABE8e"]["origin"] == "default"


def test_base_editor_table_plain_text_source(tmp_path, cfg_and_db):
    """A source cell with no ``(...)`` part is returned as plain text with
    ``url=None`` so the frontend renders it as plain text."""
    from pathlib import Path
    from beguided.web.config import get_config

    cfg = get_config()
    presets = Path(cfg.base_editor_profiles_file)
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
    )

    client = _client(cfg_and_db)
    res = client.get("/base-editor-table?source=default").json()
    assert res["ok"] is True
    assert res["rows"][0]["source"] == {"text": "Nature2022", "url": None}


def test_base_editor_table_rejects_non_http_url(tmp_path, cfg_and_db):
    """Non-http(s) URLs (e.g. ``javascript:...``) must collapse to a plain
    text cell; the entire raw cell is kept as ``text`` so nothing is lost."""
    from pathlib import Path
    from beguided.web.config import get_config

    cfg = get_config()
    presets = Path(cfg.base_editor_profiles_file)
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\t"
        "[click](javascript:alert(1))\n"
    )

    client = _client(cfg_and_db)
    res = client.get("/base-editor-table?source=default").json()
    assert res["ok"] is True
    src = res["rows"][0]["source"]
    assert src["url"] is None
    assert "javascript:alert(1)" in src["text"]


def test_base_editor_table_custom_profiles_have_empty_source(cfg_and_db):
    """Custom uploads drop the ``source`` column, so the rendered row must
    carry an empty ``{text, url}`` payload."""
    client = _client(cfg_and_db)
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    res = client.get("/base-editor-table?source=custom").json()
    assert res["ok"] is True
    assert res["rows"][0]["source"] == {"text": "", "url": None}
    assert res["rows"][0]["origin"] == "custom"


# ---------------------------------------------------------------------------
# Status endpoint
# ---------------------------------------------------------------------------

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
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    res = client.get("/base-editor-profiles/status")
    body = res.json()
    assert body["loaded"] is True
    assert body["names"] == ["CustomABE1", "CustomCBE1"]


# ---------------------------------------------------------------------------
# /check and /run with custom profiles
# ---------------------------------------------------------------------------

def test_check_uses_custom_profiles(cfg_and_db):
    client = _client(cfg_and_db)
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
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
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
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
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
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


# ---------------------------------------------------------------------------
# Store-level invariants
# ---------------------------------------------------------------------------

def test_store_holds_normalised_tsv_after_xlsx_upload(cfg_and_db):
    """Regardless of upload format, ``get_store().get_bytes()`` returns
    canonical TSV bytes so downstream readers (workdir writer, scheduler
    audit, /base-editor-table) keep working unchanged."""
    from beguided.web.profiles_store import get_store

    client = _client(cfg_and_db)
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    stored = get_store().get_bytes()
    assert stored is not None
    df = pd.read_csv(io.BytesIO(stored), sep="\t")
    assert sorted(df["name"].tolist()) == ["CustomABE1", "CustomCBE1"]


def test_store_lives_in_memory_only(cfg_and_db, tmp_path):
    from beguided.web.profiles_store import get_store

    client = _client(cfg_and_db)
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    assert get_store().is_loaded() is True
    assert get_store().get_bytes() is not None

    profile_paths = [p for p in tmp_path.rglob("*profiles*")]
    assert profile_paths == [], f"upload wrote unexpected files: {profile_paths}"

    get_store().clear()
    assert get_store().is_loaded() is False
    assert client.get("/base-editor-profiles/status").json()["loaded"] is False


# ---------------------------------------------------------------------------
# Re-upload behaviour
# ---------------------------------------------------------------------------

def test_reupload_same_xlsx_is_idempotent(cfg_and_db):
    """The browser-side localStorage cache re-POSTs the same bytes on every
    page load. The server must accept the re-upload without complaint and
    end up holding the same profiles."""
    client = _client(cfg_and_db)
    for _ in range(2):
        res = _post_profiles(
            client, "p.xlsx", None, XLSX_MIME,
            df=SAMPLE_PROFILES_DF, suffix=".xlsx",
        )
        assert res.status_code == 200, res.text
        assert res.json()["ok"] is True
        assert res.json()["names"] == ["CustomABE1", "CustomCBE1"]

    status = client.get("/base-editor-profiles/status").json()
    assert status["loaded"] is True
    assert status["names"] == ["CustomABE1", "CustomCBE1"]


def test_reupload_different_xlsx_replaces(cfg_and_db):
    """A new upload supersedes whatever was previously stored. The
    localStorage cache mirrors this: every successful upload overwrites the
    stored bytes, so the next page load re-POSTs the latest file only."""
    client = _client(cfg_and_db)
    first = pd.DataFrame([{
        "name": "FirstBE", "editing": "ABE",
        "general_edit_start": 2, "general_edit_end": 13,
        "high_edit_start": 3, "high_edit_end": 8, "pam": "NGG",
    }])
    second = pd.DataFrame([{
        "name": "SecondBE", "editing": "CBE",
        "general_edit_start": 2, "general_edit_end": 10,
        "high_edit_start": 6, "high_edit_end": 6, "pam": "NGG",
    }])

    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=first, suffix=".xlsx",
    )
    assert res.status_code == 200
    assert res.json()["names"] == ["FirstBE"]

    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=second, suffix=".xlsx",
    )
    assert res.status_code == 200
    assert res.json()["names"] == ["SecondBE"]

    status = client.get("/base-editor-profiles/status").json()
    assert status["names"] == ["SecondBE"]
    table = client.get("/base-editor-table?source=custom").json()
    custom_names = [row["name"] for row in table["rows"]]
    assert "FirstBE" not in custom_names
    assert "SecondBE" in custom_names


def test_invalid_reupload_does_not_clear_existing(cfg_and_db):
    """If the browser's localStorage holds a stale payload that fails
    server validation, the existing server-side store must remain
    intact. The JS hydrate path drops the local entry on validation
    failure, so this guarantees a failed bootstrap never erases another
    user's upload on the same server."""
    client = _client(cfg_and_db)
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )

    bad = pd.DataFrame([
        {"name": "NoPam", "editing": "ABE"},
    ])
    res = _post_profiles(
        client, "p.xlsx", None, XLSX_MIME, df=bad, suffix=".xlsx",
    )
    assert res.status_code == 400

    status = client.get("/base-editor-profiles/status").json()
    assert status["loaded"] is True
    assert status["names"] == ["CustomABE1", "CustomCBE1"]


def test_build_config_writes_custom_profiles_to_workdir(cfg_and_db, tmp_path):
    cfg, _ = cfg_and_db
    from beguided.web.tasks import _build_config

    work_dir = tmp_path / "custom_workdir"
    work_dir.mkdir()

    targets = pd.DataFrame({
        "chromosome": ["chr1"],
        "position": [50],
        "target_base": ["T"],
        "desired_base": ["C"],
    })
    # Use the TSV-form bytes the store hands out (the canonical wire format
    # for the scheduler/workdir writer).
    from beguided.web.profiles_store import get_store
    get_store().clear()
    _post_profiles(
        client_tmp := _client(cfg_and_db), "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
    )
    custom_b64 = base64.b64encode(get_store().get_bytes()).decode("ascii")
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


# ---------------------------------------------------------------------------
# /check — comma-separated editor names
# ---------------------------------------------------------------------------

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
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
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
    _post_profiles(
        client, "p.xlsx", None, XLSX_MIME,
        df=SAMPLE_PROFILES_DF, suffix=".xlsx",
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
