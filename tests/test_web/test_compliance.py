"""Compliance & safety tests for the BEguided web server.

These tests cover the security and GDPR-related work tracked in
``compliance_safety.md``: response headers, XSS regressions, legal
pages, retention sweeper, IP anonymisation in logs, admin-token
gating, MIME sniffing, and the static-asset SRI hash.
"""

from __future__ import annotations

import base64
import io
import logging
import time
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("uvicorn")


WEB_DIR = Path(__file__).resolve().parents[2] / "beguided" / "web"
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"
APP_JS = STATIC_DIR / "app.js"


# ---------------------------------------------------------------------------
# Static-asset hygiene
# ---------------------------------------------------------------------------


def test_chartjs_is_self_hosted_in_index_html():
    """Chart.js must load from /static/vendor/ with an SRI hash, never from a CDN."""
    html = (TEMPLATE_DIR / "index.html").read_text()
    assert "cdn.jsdelivr.net" not in html, "index.html must not reference any CDN"
    assert "googleapis.com" not in html
    assert "unpkg.com" not in html
    assert "/static/vendor/chart.umd.min.js" in html
    # Subresource Integrity + crossorigin are mandatory for any
    # locally-loaded third-party bundle so a compromised file swap would
    # be detected by the browser.
    assert 'integrity="sha384-' in html
    assert 'crossorigin="anonymous"' in html


def test_chartjs_bundle_exists_and_matches_sri_hash():
    """The vendored Chart.js bundle must exist and its SRI hash must match."""
    bundle = STATIC_DIR / "vendor" / "chart.umd.min.js"
    assert bundle.exists(), "vendored Chart.js bundle is missing"
    assert bundle.stat().st_size > 50_000, "Chart.js bundle looks too small"

    import hashlib

    digest = hashlib.sha384(bundle.read_bytes()).digest()
    import base64 as _b64

    expected = "sha384-" + _b64.b64encode(digest).decode("ascii")
    html = (TEMPLATE_DIR / "index.html").read_text()
    assert expected in html, (
        f"SRI hash in index.html ({expected}) does not match the vendored bundle"
    )


def test_appjs_escapes_user_supplied_error_text():
    """Regression: data.error was interpolated raw into innerHTML in app.js.

    Search the source for any line that BOTH assigns to an ``innerHTML``
    and interpolates ``${data.error}`` without wrapping the value in
    ``escapeHtml(...)``. Lines that pass ``data.error`` to a
    ``textContent`` sink (e.g. ``showToast``) are not in scope — those
    are auto-escaped by the DOM.
    """
    src = APP_JS.read_text()
    bad = []
    for i, line in enumerate(src.splitlines(), start=1):
        if "innerHTML" not in line:
            continue
        if "${data.error" not in line:
            continue
        if "escapeHtml(data.error" in line:
            continue
        bad.append((i, line.strip()))
    assert not bad, f"unescaped data.error in innerHTML: {bad}"


# ---------------------------------------------------------------------------
# Response headers (security middleware)
# ---------------------------------------------------------------------------


@pytest.fixture
def cfg_and_client(tmp_path, monkeypatch):
    """Build a TestClient with the real middleware stack wired up."""
    from beguided.web import config as cfg_mod
    from beguided.web import scheduler as sched_mod
    from beguided.web.config import WebConfig, set_config
    from beguided.web.scheduler import Scheduler, install

    presets = tmp_path / "presets.tsv"
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\n"
    )
    db = tmp_path / "jobs.db"
    set_config(WebConfig(
        base_editor_profiles_file=str(presets),
        genome_fasta="/nonexistent",
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        max_input_targets=1000,
        db_path=str(db),
        queue_max=2,
        cooldown_seconds=0,
        max_upload_mb=5,
        trusted_proxy_ips=["127.0.0.1"],
        admin_token="test-admin-token",
        log_ip_salt="test-salt",
        usage_stats_enabled=False,
    ))
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=0,
                    log_ip_salt="test-salt")
    install(sch)
    sch.start()
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app

    client = TestClient(create_app())
    yield cfg_mod.get_config(), sch, client
    sch.stop()


def test_security_headers_on_root(cfg_and_client):
    _, _, client = cfg_and_client
    res = client.get("/")
    assert res.status_code == 200
    headers = res.headers
    assert "default-src 'self'" in headers.get("content-security-policy", "")
    assert headers.get("x-content-type-options") == "nosniff"
    assert headers.get("referrer-policy") == "strict-origin-when-cross-origin"
    assert "camera=()" in headers.get("permissions-policy", "")
    assert headers.get("cross-origin-opener-policy") == "same-origin"
    assert headers.get("server") == "beguided"
    # HSTS is gated on HTTPS; over HTTP it must not be emitted.
    assert "strict-transport-security" not in {k.lower() for k in headers}


def test_security_headers_on_json(cfg_and_client):
    _, _, client = cfg_and_client
    res = client.get("/config")
    assert res.status_code == 200
    assert "content-security-policy" in {k.lower() for k in res.headers}
    assert res.headers.get("x-content-type-options") == "nosniff"


def test_proxy_headers_set_real_client_ip(tmp_path, monkeypatch):
    """ProxyHeadersMiddleware must surface the X-Forwarded-For client IP."""
    from beguided.web import config as cfg_mod
    from beguided.web import scheduler as sched_mod
    from beguided.web.config import WebConfig, set_config
    from beguided.web.scheduler import Scheduler, install

    presets = tmp_path / "presets.tsv"
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\n"
    )
    db = tmp_path / "jobs.db"
    set_config(WebConfig(
        base_editor_profiles_file=str(presets),
        genome_fasta="/nonexistent",
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        db_path=str(db),
        queue_max=2,
        cooldown_seconds=0,
        max_upload_mb=5,
        # ``testclient`` is the hostname Starlette's TestClient uses for
        # the simulated upstream peer. Including it here exercises the
        # real middleware path the way Caddy will use it in production.
        trusted_proxy_ips=["127.0.0.1", "::1", "testclient"],
        usage_stats_enabled=False,
    ))
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=0)
    install(sch)
    sch.start()
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app

    client = TestClient(create_app())

    captured = {}
    real_submit = sch.submit

    def spy_submit(ip, payload):
        captured["ip"] = ip
        return real_submit(ip, payload)

    sch.submit = spy_submit

    targets_b64 = base64.b64encode(
        b"chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
    ).decode("ascii")
    res = client.post(
        "/run",
        data={"raw_b64": targets_b64, "base_editor": "ABE8e"},
        headers={"X-Forwarded-For": "203.0.113.42"},
    )
    assert res.status_code in (200, 503)
    assert captured.get("ip") == "203.0.113.42"


# ---------------------------------------------------------------------------
# Legal pages
# ---------------------------------------------------------------------------


def test_legal_routes_render_with_empty_config(cfg_and_client):
    """All legal pages must render with HTTP 200 even when ``legal.yaml``
    is missing or empty (the default for local dev and CI). No
    ``TODO_LEGAL`` tripwires should remain in the templates.
    """
    _, _, client = cfg_and_client
    for path in ("/impressum", "/datenschutz", "/agb", "/accessibility"):
        res = client.get(path)
        assert res.status_code == 200, f"{path} returned {res.status_code}"
        assert "text/html" in res.headers.get("content-type", "")
        body = res.text
        assert "<html" in body.lower()
        # No tripwire placeholders must remain in any legal page.
        assert "TODO_LEGAL" not in body, (
            f"{path} still contains TODO_LEGAL placeholders"
        )


def test_load_legal_config_missing_returns_empty(tmp_path):
    from beguided.web.legal_config import load_legal_config

    cfg = load_legal_config(tmp_path / "does-not-exist.yaml")
    assert cfg.institution_name == ""
    assert cfg.address_lines == []
    assert cfg.contact_email == ""
    assert cfg.dpo_block == []
    assert cfg.log_retention_days == 0


def test_load_legal_config_directory_returns_empty(tmp_path, caplog):
    """A directory at the config path (the docker bind-mount footgun when
    the host file is missing) must be treated like a missing file."""
    import logging

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "legal.yaml"
    p.mkdir()
    with caplog.at_level(logging.WARNING, logger="beguided.web.legal_config"):
        cfg = load_legal_config(p)
    assert cfg.institution_name == ""
    assert cfg.address_lines == []
    assert cfg.contact_email == ""
    assert cfg.dpo_block == []
    assert cfg.log_retention_days == 0
    assert any("not found or is not a file" in rec.message for rec in caplog.records)


def test_load_legal_config_none_returns_empty():
    from beguided.web.legal_config import load_legal_config

    cfg = load_legal_config(None)
    assert cfg.institution_name == ""
    assert cfg.address_lines == []


def test_load_legal_config_empty_file_returns_empty(tmp_path):
    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "empty.yaml"
    p.write_text("")
    cfg = load_legal_config(p)
    assert cfg.institution_name == ""


def test_load_legal_config_malformed_raises(tmp_path):
    import pytest

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "bad.yaml"
    p.write_text("institution_name: ACME\n  contact_email: oops\n")  # bad indent
    with pytest.raises(ValueError, match="Failed to parse legal config"):
        load_legal_config(p)


def test_load_legal_config_wrong_top_level_raises(tmp_path):
    import pytest

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "list.yaml"
    p.write_text("- one\n- two\n")
    with pytest.raises(ValueError, match="must contain a YAML mapping"):
        load_legal_config(p)


def test_load_legal_config_bad_log_retention_raises(tmp_path):
    import pytest

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "bad_int.yaml"
    p.write_text("log_retention_days: seven\n")
    with pytest.raises(ValueError, match="log_retention_days"):
        load_legal_config(p)


def test_load_legal_config_string_field_must_be_string(tmp_path):
    import pytest

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "bad_field.yaml"
    p.write_text("contact_email:\n  - a\n  - b\n")
    with pytest.raises(ValueError, match="contact_email"):
        load_legal_config(p)


def test_load_legal_config_block_from_multiline_string(tmp_path):
    """A block field given as a multi-line string must normalize to the
    same lines as a list of strings."""
    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "block_string.yaml"
    p.write_text(
        "dpo_block: |\n"
        "  Line one\n"
        "  \n"
        "    Line two\n"
    )
    cfg = load_legal_config(p)
    assert cfg.dpo_block == ["Line one", "Line two"]

    p2 = tmp_path / "block_list.yaml"
    p2.write_text("dpo_block: ['Line one', '', 'Line two']\n")
    cfg = load_legal_config(p2)
    assert cfg.dpo_block == ["Line one", "Line two"]


def test_load_legal_config_jinja_reuse_forward_reference(tmp_path):
    """A block may reference a variable that is defined later in the
    file (resolution iterates to a fixed point)."""
    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "forward.yaml"
    p.write_text(
        "controller_block:\n"
        "  - '{{ institution_name }}'\n"
        "institution_name: 'ACME University'\n"
    )
    cfg = load_legal_config(p)
    assert cfg.controller_block == ["ACME University"]


def test_load_legal_config_jinja_loop_over_list_field(tmp_path):
    """A block can loop over another list field to reuse its lines."""
    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "loop.yaml"
    p.write_text(
        "controller_block:\n"
        "  - '{{ institution_name }}'\n"
        "  - \"{% for line in address_lines %}{{ line }}\\n{% endfor %}\"\n"
        "institution_name: 'ACME University'\n"
        "address_lines: ['Street 1', '12345 Town']\n"
    )
    cfg = load_legal_config(p)
    assert cfg.controller_block == [
        "ACME University",
        "Street 1",
        "12345 Town",
    ]


def test_load_legal_config_undefined_variable_raises(tmp_path):
    """Referencing a variable that is not defined anywhere in the file
    must fail startup (a typo must not render a blank paragraph)."""
    import pytest

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "undefined.yaml"
    p.write_text("controller_block:\n  - '{{ instution_name }}'\n")
    with pytest.raises(ValueError, match="instution_name"):
        load_legal_config(p)


def test_load_legal_config_legacy_fields_warn(tmp_path, caplog):
    """Legacy per-facet fields (dpo_name, competent_dpa_*, ...) are no
    longer part of the block schema; a file still using them must log a
    warning naming the fields so the operator notices."""
    import logging

    from beguided.web.legal_config import load_legal_config

    p = tmp_path / "legacy.yaml"
    p.write_text(
        "institution_name: 'ACME University'\n"
        "dpo_name: 'old dpo field'\n"
        "competent_dpa: 'old dpa field'\n"
        "hosting_legal_postal: ['old postal']\n"
    )
    with caplog.at_level(logging.WARNING, logger="beguided.web.legal_config"):
        cfg = load_legal_config(p)
    assert cfg.institution_name == "ACME University"
    assert cfg.dpo_block == []
    messages = " ".join(rec.getMessage() for rec in caplog.records)
    assert "dpo_name" in messages
    assert "competent_dpa" in messages
    assert "hosting_legal_postal" in messages


def test_linkify_escapes_html_and_links_tokens():
    """Block lines are operator data but still untrusted input: the
    line must be HTML-escaped, and only mailto:/http(s):// tokens are
    turned into links (mailto labels drop the scheme prefix)."""
    from beguided.web.app import _linkify

    assert (
        _linkify("<script>alert(1)</script>")
        == "&lt;script&gt;alert(1)&lt;/script&gt;"
    )
    assert (
        _linkify("Email: mailto:dpa@example.de")
        == 'Email: <a href="mailto:dpa@example.de">dpa@example.de</a>'
    )
    assert (
        _linkify("Web: https://example.de/a?b=1&c=2")
        == 'Web: <a href="https://example.de/a?b=1&amp;c=2">'
        "https://example.de/a?b=1&amp;c=2</a>"
    )


def test_legal_routes_linked_from_base_template():
    """Every legal page must be reachable from the site footer on every page."""
    base = (TEMPLATE_DIR / "base.html").read_text()
    for link in ("/impressum", "/datenschutz", "/agb"):
        assert link in base, f"base.html footer is missing link to {link}"
    # Skip-link required for BFSG/EAA keyboard accessibility.
    assert 'class="skip-link"' in base
    assert 'href="#main"' in base


# ---------------------------------------------------------------------------
# Scheduler: retention sweeper + IP anonymisation
# ---------------------------------------------------------------------------


def test_purge_old_removes_finished_jobs(tmp_path):
    from beguided.web.scheduler import Scheduler

    db = tmp_path / "jobs.db"
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=0,
                    job_retention_seconds=60)
    sch.start()
    try:
        sch.submit("1.2.3.4", {"dummy": True})
        # Wait for the job to be queued (we don't actually run anything
        # here; manually mark it complete + old in the DB).
        import sqlite3
        with sqlite3.connect(str(db)) as c:
            old = time.time() - 3600
            c.execute(
                "UPDATE jobs SET status='complete', finished_at=?", (old,)
            )
        deleted = sch.purge_old()
        assert deleted == 1
    finally:
        sch.stop()


def test_purge_all_finished_keeps_running(tmp_path):
    from beguided.web.scheduler import Scheduler

    db = tmp_path / "jobs.db"
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=0)
    sch.start()
    try:
        sch.submit("1.2.3.4", {"dummy": True})
        # The job is in the DB with status=queued. purge_all_finished
        # must not touch it.
        deleted = sch.purge_all_finished()
        assert deleted == 0
        # A complete job should be deleted.
        import sqlite3
        with sqlite3.connect(str(db)) as c:
            c.execute(
                "UPDATE jobs SET status='complete', finished_at=?",
                (time.time(),),
            )
        deleted = sch.purge_all_finished()
        assert deleted == 1
    finally:
        sch.stop()


def test_ip_is_anonymised_in_log(tmp_path, caplog):
    from beguided.web.scheduler import Scheduler

    db = tmp_path / "jobs.db"
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=0,
                    log_ip_salt="very-secret-salt")
    sch.start()
    try:
        with caplog.at_level(logging.INFO, logger="beguided.web.scheduler"):
            sch.submit("203.0.113.42", {"dummy": True})
        text = "\n".join(record.getMessage() for record in caplog.records)
        assert "203.0.113.42" not in text, "plaintext IP leaked into log"
        assert "ip:" in text, "anonymised IP marker missing from log"
    finally:
        sch.stop()


def test_last_submit_is_ttl_pruned(tmp_path):
    """The cooldown map must not grow without bound."""
    from beguided.web.scheduler import Scheduler

    db = tmp_path / "jobs.db"
    sch = Scheduler(db_path=str(db), queue_max=2, cooldown_seconds=1)
    sch.start()
    try:
        # First submit from IP "x" goes through and seeds the map. The map
        # is nested: cooldown kind -> IP -> last submit time.
        sch.submit("10.0.0.1", {"dummy": True})
        assert "10.0.0.1" in sch._last_submit["pipeline"]
        # Roll the monotonic clock back far enough for the prune to
        # consider the entry stale.
        sch._last_submit["pipeline"]["10.0.0.1"] = 0.0
        sch._prune_last_submit(now=time.monotonic())
        assert "10.0.0.1" not in sch._last_submit.get("pipeline", {})
    finally:
        sch.stop()


def test_job_retention_hours_env_knob(monkeypatch):
    """``BEGUIDED_JOB_RETENTION_HOURS`` must override the 24 h default.

    Regression: the field was a plain attribute that ignored the env var
    (``--job-retention-hours`` was the only working override), so the knob
    advertised in the CLI help and docs was a silent no-op.
    """
    from beguided.web.config import WebConfig

    monkeypatch.setenv("BEGUIDED_JOB_RETENTION_HOURS", "1")
    assert WebConfig().job_retention_hours == 1

    monkeypatch.delenv("BEGUIDED_JOB_RETENTION_HOURS")
    assert WebConfig().job_retention_hours == 24


# ---------------------------------------------------------------------------
# Admin purge endpoint
# ---------------------------------------------------------------------------


def test_admin_purge_all_requires_token(cfg_and_client):
    _, sch, client = cfg_and_client

    # First seed a completed job.
    sch.submit("1.2.3.4", {"dummy": True})
    import sqlite3
    with sqlite3.connect(sch.db_path) as c:
        c.execute(
            "UPDATE jobs SET status='complete', finished_at=?", (time.time(),)
        )

    # No token → 403.
    res = client.post("/admin/purge-all")
    assert res.status_code == 403

    # Wrong token → 403.
    res = client.post(
        "/admin/purge-all", headers={"X-Admin-Token": "wrong"}
    )
    assert res.status_code == 403

    # Correct token → 200 + reports a deletion.
    res = client.post(
        "/admin/purge-all", headers={"X-Admin-Token": "test-admin-token"}
    )
    assert res.status_code == 200
    assert res.json()["ok"] is True
    assert res.json()["deleted"] >= 1


def test_admin_purge_all_disabled_without_token(cfg_and_client):
    """When admin_token is empty, the endpoint refuses every request."""
    from beguided.web import config as cfg_mod

    cfg, _, client = cfg_and_client
    cfg_mod.set_config(
        cfg_mod.WebConfig(
            base_editor_profiles_file=cfg.base_editor_profiles_file,
            genome_fasta=cfg.genome_fasta,
            gtf_path=cfg.gtf_path,
            beguidedoff_genome=cfg.beguidedoff_genome,
            db_path=cfg.db_path,
            queue_max=cfg.queue_max,
            cooldown_seconds=cfg.cooldown_seconds,
            max_upload_mb=cfg.max_upload_mb,
            admin_token="",
        )
    )
    res = client.post(
        "/admin/purge-all", headers={"X-Admin-Token": "anything"}
    )
    assert res.status_code == 403


# ---------------------------------------------------------------------------
# MIME sniffing on /upload
# ---------------------------------------------------------------------------


def test_upload_rejects_binary_with_tsv_extension(cfg_and_client):
    _, _, client = cfg_and_client
    # Pretend to upload a ZIP file masquerading as a .tsv. The magic
    # bytes PK\x03\x04 must trigger a 400.
    fake = b"PK\x03\x04" + b"\x00" * 100
    res = client.post(
        "/upload",
        files={"file": ("targets.tsv", io.BytesIO(fake),
                         "application/octet-stream")},
    )
    assert res.status_code == 200
    body = res.json()
    assert body["ok"] is False
    assert "do not match" in body["error"].lower()


def test_upload_accepts_real_tsv(cfg_and_client):
    _, _, client = cfg_and_client
    payload = (
        b"chromosome\tposition\ttarget_base\tdesired_base\n"
        b"chr1\t50\tT\tC\n"
    )
    res = client.post(
        "/upload",
        files={"file": ("targets.tsv", io.BytesIO(payload),
                         "text/tab-separated-values")},
    )
    assert res.status_code == 200
    assert res.json()["ok"] is True


# ---------------------------------------------------------------------------
# Index page accessibility hooks
# ---------------------------------------------------------------------------


def test_index_html_has_no_input_disclaimer():
    """The spreadsheet panel must not surface a disclaimer or terms-of-use
    notice inline — the footer link to /agb is the single source of
    acceptable-use policy. No consent checkbox either.
    """
    html = (TEMPLATE_DIR / "index.html").read_text()
    # No input-mode disclaimer under the spreadsheet Run button.
    assert "input-mode-disclaimer" not in html
    assert "auto-deleted" not in html
    assert "Acceptable Use" not in html
    # No consent checkbox must be present.
    assert "consentCheckbox" not in html
    assert "consentRow" not in html


def test_base_template_has_skip_link():
    """BFSG/EAA: every page must provide a skip-to-content link."""
    base = (TEMPLATE_DIR / "base.html").read_text()
    assert "skip-link" in base
    assert "#main" in base
    assert 'id="main"' in base
    # lang attribute set on the root html element (split on the <html tag,
    # not on the doctype's closing bracket).
    html_tag = base[base.index("<html"): base.index(">", base.index("<html"))]
    assert 'lang="' in html_tag


def test_style_css_defines_reduced_motion_and_focus_visible():
    css = (STATIC_DIR / "style.css").read_text()
    assert "prefers-reduced-motion: reduce" in css
    assert ":focus-visible" in css
    assert "skip-link" in css
    assert ".site-footer" in css
    # The removed input-mode disclaimer and consent checkbox styles must
    # both be gone.
    assert ".input-mode-disclaimer" not in css
    assert ".consent-row" not in css


# ---------------------------------------------------------------------------
# App.js integrity: no unescaped ${...} inside innerHTML strings
# ---------------------------------------------------------------------------


def test_appjs_innerHTML_templates_escape_user_input():
    """Every ``innerHTML = `` assignment in app.js that interpolates a
    non-literal expression must either escape that value via
    ``escapeHtml(...)`` or use only literals / local-constant ternaries.

    The exception is ``escapeHtml`` itself which sets
    ``div.innerHTML`` on a freshly created div — the value there is
    already-escaped HTML by construction.
    """
    import re

    src = APP_JS.read_text()
    pattern = re.compile(r"\b([a-zA-Z_][\w]*)\.innerHTML\s*=\s*(`[^`]*`)")
    # Patterns that indicate the expression touches user data:
    #   - member access: ``.foo`` or ``?.foo`` (not ``.`` inside a string)
    #   - function call: ``foo(`` (not ``(`` inside a string)
    member_access = re.compile(r"\?\s*\.[A-Za-z_]|\.[A-Za-z_]")
    function_call = re.compile(r"[A-Za-z_]\w*\s*\(")
    violations = []
    for match in pattern.finditer(src):
        line_no = src.count("\n", 0, match.start()) + 1
        target, body = match.group(1), match.group(2)
        # Skip the escapeHtml helper itself (target is "div" inside that fn).
        if target == "div" and "div.textContent" in src[max(0, match.start() - 80): match.start()]:
            continue
        for inner in re.findall(r"\$\{([^}]+)\}", body):
            expr = inner.strip()
            # Safe: already wrapped in escapeHtml(...).
            if expr.startswith("escapeHtml("):
                continue
            # Safe: literal / local-constant expression. Detect literals
            # by stripping all string contents and looking for any
            # remaining identifier-followed-by-paren or member access.
            no_strings = re.sub(r'"(?:\\.|[^"\\])*"', '""', expr)
            no_strings = re.sub(r"'(?:\\.|[^'\\])*'", "''", no_strings)
            if member_access.search(no_strings) or function_call.search(no_strings):
                violations.append((line_no, expr))
    assert not violations, (
        "unescaped interpolations in innerHTML: " + repr(violations)
    )
