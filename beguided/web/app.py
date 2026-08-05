"""FastAPI application for the BEguided web interface."""

from __future__ import annotations

import asyncio
import base64
import io
import logging
import re
import secrets
from pathlib import Path
from typing import Annotated, Optional
from urllib.parse import urlparse

import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, Query, Request, Response, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from jinja2 import Environment, FileSystemLoader

from beguided.filtering.filters import filter_guides, VALID_GC_FILTERS
from beguided.logging_config import setup_logging
from beguided.web.condense import condense_df
from beguided.web.config import (
    get_config,
    get_base_editor_display_names,
    load_base_editor_groups,
)
from beguided.editing.presets import (
    get_available_base_editor_names,
    get_supported_editor_keys,
    load_base_editor_profiles,
)
from beguided.web import genome_cache, gtf_cache
from beguided.web.profiles_store import (
    ProfilesUploadError,
    SUPPORTED_PROFILE_SUFFIXES,
    get_store as get_profiles_store,
    validate_uploaded_profiles,
)
from beguided.web.sequence_input import (
    SequenceParseError,
    parse_sequence_text,
    records_to_fasta,
    records_to_targets_dataframe,
)
from beguided.web.scheduler import (
    Cooldown,
    Job,
    QueueFull,
    Scheduler,
    current as current_scheduler,
    install as install_scheduler,
)
from beguided.web.safety import (
    MaxBodySizeMiddleware,
    ProxyHeadersMiddleware,
    SecurityHeadersMiddleware,
    validate_chromosomes,
    validate_operator_paths,
)
from beguided.web.legal_config import LegalConfig, load_legal_config
from beguided.web.sse_logger import attach_sse_handler_to_beguided_logs, get_sse_handler

logger = logging.getLogger("beguided.web.app")

REQUIRED_COLUMNS = {"chromosome", "position", "target_base", "desired_base"}
VALID_FILTER_MODES = ("all", "no_bystanders", "silent")
VALID_INPUT_MODES = ("spreadsheet", "sequence")
CUSTOM_EDITOR_KEY = "custom"
CUSTOM_EDITOR_LABEL = "Custom (uploaded)"

# Module-level singleton for the legal-page content. Loaded once at
# ``create_app`` startup from the YAML file referenced by
# ``WebConfig.legal_config_path`` (env ``BEGUIDED_LEGAL_CONFIG``). The
# route handlers read this rather than re-parsing the file per request.
_legal_config: LegalConfig = LegalConfig()


def install_legal_config(cfg: LegalConfig) -> LegalConfig:
    """Override the module-level legal-config singleton."""
    global _legal_config
    _legal_config = cfg
    return _legal_config


def current_legal_config() -> LegalConfig:
    """Return the legal config installed at app startup."""
    return _legal_config

# Magic bytes used by ``_sniff_upload_kind`` to reject uploads whose
# extension lies about their contents. The XLSX / XLS formats are ZIP /
# OLE compound documents; TSV / CSV / XLSX / XLS are the only extensions
# accepted by ``/upload``. ASCII-printable bytes with at least one tab or
# newline are accepted as TSV/CSV.
_XLSX_MAGIC = b"PK\x03\x04"
_XLS_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"

# Match a single inline markdown link ``[text](url)``. Anchored on both ends
# so a plain-text cell never accidentally loses its text. Whitespace inside
# the brackets is preserved as part of the link text.
_MARKDOWN_LINK_RE = re.compile(r"^\[([^\]]*)\]\(([^)]*)\)$")


def _sniff_upload_kind(content: bytes, suffix: str) -> str:
    """Return the inferred upload kind ("tsv", "csv", "xlsx", "xls", "binary").

    Anything that contains a NUL byte in the first 1 KiB is treated as
    binary and rejected unless the extension is an Office format whose
    magic matches. Empty files are rejected.
    """
    head = content[:1024]
    if not head:
        return "binary"
    sfx = suffix.lower()
    if sfx in (".xlsx", ".xls"):
        if sfx == ".xlsx" and head.startswith(_XLSX_MAGIC):
            return "xlsx"
        if sfx == ".xls" and head.startswith(_XLS_MAGIC):
            return "xls"
        return "binary"
    # CSV / TSV: must be plain ASCII / UTF-8 with no NUL bytes and must
    # contain at least one tab (TSV) or comma (CSV).
    try:
        head_text = head.decode("utf-8")
    except UnicodeDecodeError:
        return "binary"
    if "\x00" in head_text:
        return "binary"
    if sfx == ".csv":
        return "csv" if "," in head_text else "tsv"
    if "\t" in head_text:
        return "tsv"
    if "," in head_text:
        return "csv"
    return "tsv"


def _genome_chromosome_keys() -> list[str]:
    """Return the chromosome keys of the preloaded genome, or [] if not ready.

    Blocks briefly on the genome preload completion so the check runs
    after the FASTA has been parsed. Returns an empty list if the genome
    failed to preload, no genome is configured, or the preload is still
    running after the short wait; in any of those cases the chromosome
    check is skipped downstream so the UI flow is not blocked.
    """
    try:
        genome = genome_cache.wait_until_ready(timeout=5.0)
    except (RuntimeError, TimeoutError):
        return []
    if not genome:
        return []
    return list(genome.keys())


def _render_template(name, **kwargs) -> str:
    """Render a Jinja2 template to HTML string."""
    base_dir = Path(__file__).resolve().parent
    env = Environment(loader=FileSystemLoader(str(base_dir / "templates")))
    template = env.get_template(name)
    return template.render(**kwargs)


def _parse_source_link(value: object) -> dict:
    """Parse a TSV ``source`` cell into ``{text, url}``.

    Accepts the inline markdown form ``[text](url)`` and keeps the URL only
    when it uses an ``http(s)`` scheme (anything else — e.g. ``javascript:``
    — collapses to ``url=None`` so the frontend renders plain text).
    Plain-text cells return ``{text: <value>, url: None}``.
    """
    text = "" if value is None else str(value)
    match = _MARKDOWN_LINK_RE.match(text.strip())
    if not match:
        return {"text": text, "url": None}
    link_text, url = match.group(1), match.group(2).strip()
    if not url:
        return {"text": link_text, "url": None}
    parsed = urlparse(url)
    if parsed.scheme.lower() not in ("http", "https"):
        return {"text": text, "url": None}
    return {"text": link_text, "url": url}


def _df_to_table_rows(df: pd.DataFrame, source: str) -> list[dict]:
    """Convert a profiles DataFrame to the condensed row format used by
    the API. Each row includes an ``origin`` field ("default" / "custom")
    and a ``source`` field with ``{text, url}`` parsed from the TSV's
    optional ``source`` column (empty when the column is absent, e.g. for
    uploaded custom profiles)."""
    has_source_col = "source" in df.columns
    rows = []
    for _, r in df.iterrows():
        if has_source_col:
            source_link = _parse_source_link(r["source"])
        else:
            source_link = {"text": "", "url": None}
        rows.append({
            "name": str(r["name"]),
            "editing": str(r["editing"]),
            "general": f"{int(r['general_edit_start'])}-{int(r['general_edit_end'])}",
            "high": f"{int(r['high_edit_start'])}-{int(r['high_edit_end'])}",
            "pam": str(r["pam"]),
            "cas": str(r.get("cas", "")) if "cas" in df.columns else "",
            "origin": source,
            "source": source_link,
        })
    return rows


def _decode_targets(raw_b64: str) -> pd.DataFrame:
    """Decode a base64-encoded TSV into a DataFrame. Raises HTTPException(400)."""
    try:
        content = base64.b64decode(raw_b64.encode("ascii"))
    except Exception as exc:
        raise HTTPException(400, f"Invalid base64 payload: {exc}") from exc
    try:
        df = pd.read_csv(io.BytesIO(content), sep="\t")
    except Exception as exc:
        raise HTTPException(400, f"Failed to parse targets TSV: {exc}") from exc
    df.columns = df.columns.str.strip()
    return df


def _partition_names_by_source(names: list[str]):
    """Split a list of editor names into the subset that lives in the
    custom-profiles store vs the default profiles file.

    Returns (custom_names, default_names, missing_names).
    """
    store = get_profiles_store()
    custom_names = set(store.names()) if store.is_loaded() else set()
    try:
        default_df = pd.read_csv(
            get_config().base_editor_profiles_file, sep="\t"
        )
        default_names = set(default_df["name"].astype(str).tolist())
    except Exception:
        default_names = set()

    selected = set(names)
    custom_subset = sorted(selected & custom_names)
    default_subset = sorted(selected & default_names)
    missing = sorted(selected - custom_names - default_names)
    return custom_subset, default_subset, missing


def _build_merged_profiles_b64(custom_names: list[str], default_names: list[str]):
    """Build a base64-encoded TSV that contains the union of custom + default
    profiles, restricted to the requested names. Returns the base64 string,
    or None if no merge is needed (only one source in play)."""
    if not custom_names and not default_names:
        return None
    if not custom_names:
        return None  # default-only: caller uses the default file directly
    if not default_names:
        return None  # custom-only: caller uses the custom bytes directly

    store = get_profiles_store()
    custom_df = pd.read_csv(io.BytesIO(store.get_bytes()), sep="\t")
    default_df = pd.read_csv(
        get_config().base_editor_profiles_file, sep="\t"
    )
    custom_rows = custom_df[custom_df["name"].isin(custom_names)]
    default_rows = default_df[default_df["name"].isin(default_names)]
    merged = pd.concat([custom_rows, default_rows], ignore_index=True)
    return base64.b64encode(merged.to_csv(sep="\t", index=False).encode()).decode("ascii")


def _load_profiles_for_check(base_editor: str):
    """Resolve a base_editor dropdown value to (profiles, names, error_dict).

    ``profiles`` is the result of ``load_base_editor_profiles``; ``names`` is
    the list of editor names the user actually selected. On failure the
    third element is a JSON error body and the first two are ``None``.

    Accepted forms for ``base_editor``:
      * ``"custom"`` — all uploaded custom profile names.
      * A group key (``"all"``, ``"all ABEs"``, ``"all CBEs"``, or any single
        editor name) — resolved via ``load_base_editor_groups`` (backward
        compatibility for the single-select dropdown).
      * A comma-separated list of individual editor names — each name is
        looked up in the custom store first, then the default file.
        When both sources are involved, a merged profiles dict is loaded
        from a temporary TSV built on the fly.
    """
    if base_editor == CUSTOM_EDITOR_KEY:
        store = get_profiles_store()
        if not store.is_loaded():
            return None, None, {"ok": False, "error": "No custom profiles uploaded."}
        content = store.get_bytes()
        try:
            profiles = load_base_editor_profiles(io.BytesIO(content))
        except (ValueError, KeyError) as exc:
            return None, None, {"ok": False, "error": f"Invalid custom profiles: {exc}"}
        return profiles, store.names(), None

    cfg = get_config()
    groups = load_base_editor_groups()
    if base_editor in groups:
        try:
            profiles = load_base_editor_profiles(cfg.base_editor_profiles_file)
        except (ValueError, KeyError) as exc:
            return None, None, {"ok": False, "error": f"Failed to load profiles: {exc}"}
        return profiles, groups[base_editor], None

    # Comma-separated list of individual editor names (multi-select UI).
    requested = [n.strip() for n in base_editor.split(",") if n.strip()]
    if not requested:
        return None, None, {"ok": False, "error": "No base editor selected."}

    custom_subset, default_subset, missing = _partition_names_by_source(requested)
    if missing:
        return None, None, {
            "ok": False,
            "error": f"Unknown base editor(s): {', '.join(missing)}",
        }

    merged_b64 = _build_merged_profiles_b64(custom_subset, default_subset)
    if merged_b64 is not None:
        try:
            profiles = load_base_editor_profiles(
                io.BytesIO(base64.b64decode(merged_b64.encode("ascii")))
            )
        except (ValueError, KeyError) as exc:
            return None, None, {"ok": False, "error": f"Failed to load merged profiles: {exc}"}
    elif custom_subset:
        try:
            profiles = load_base_editor_profiles(
                io.BytesIO(get_profiles_store().get_bytes())
            )
        except (ValueError, KeyError) as exc:
            return None, None, {"ok": False, "error": f"Invalid custom profiles: {exc}"}
    else:
        try:
            profiles = load_base_editor_profiles(cfg.base_editor_profiles_file)
        except (ValueError, KeyError) as exc:
            return None, None, {"ok": False, "error": f"Failed to load profiles: {exc}"}

    return profiles, requested, None


def _resolve_base_editor_for_run(base_editor: str):
    """Resolve a base_editor dropdown value for a pipeline run.

    Returns (payload_dict, error_dict). Exactly one is None.

    Accepted forms for ``base_editor``: see ``_load_profiles_for_check``.
    For comma-separated individual editor names spanning both the custom
    store and the default file, a ``merged_profiles_b64`` is included so
    the pipeline can load a single merged TSV.
    """
    if base_editor == CUSTOM_EDITOR_KEY:
        store = get_profiles_store()
        if not store.is_loaded():
            return None, {"ok": False, "error": "No custom profiles uploaded."}
        return {
            "base_editor_names": store.names(),
            "use_custom_profiles": True,
            "custom_profiles_b64": store.get_b64(),
        }, None

    groups = load_base_editor_groups()
    if base_editor in groups:
        return {
            "base_editor_names": groups[base_editor],
            "use_custom_profiles": False,
            "custom_profiles_b64": None,
            "merged_profiles_b64": None,
        }, None

    requested = [n.strip() for n in base_editor.split(",") if n.strip()]
    if not requested:
        return None, {"ok": False, "error": "No base editor selected."}

    custom_subset, default_subset, missing = _partition_names_by_source(requested)
    if missing:
        return None, {
            "ok": False,
            "error": f"Unknown base editor(s): {', '.join(missing)}",
        }

    merged_b64 = _build_merged_profiles_b64(custom_subset, default_subset)
    if merged_b64 is None and custom_subset:
        # custom-only selection: write the full custom TSV
        merged_b64 = get_profiles_store().get_b64()

    return {
        "base_editor_names": requested,
        "use_custom_profiles": False,
        "custom_profiles_b64": None,
        "merged_profiles_b64": merged_b64,
    }, None


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    cfg = get_config()
    app = FastAPI(title="BEguided Web")

    # Load legal-page content once at startup. Missing or empty values
    # are tolerated (the pages render with blanks + a startup warning);
    # malformed YAML fails fast. The result is attached to ``app.state``
    # so route handlers can read it without re-reading the file.
    legal_cfg = load_legal_config(
        Path(cfg.legal_config_path) if cfg.legal_config_path else None
    )
    app.state.legal_config = legal_cfg
    install_legal_config(legal_cfg)

    # Cap request body size to keep memory bounded.
    app.add_middleware(MaxBodySizeMiddleware, max_bytes=cfg.max_upload_mb * 1024 * 1024)
    # Add the fixed security response headers (CSP, X-Content-Type-Options,
    # Referrer-Policy, Permissions-Policy, Cross-Origin-Opener-Policy;
    # HSTS only on HTTPS). Stripped uvicorn ``Server`` header.
    app.add_middleware(SecurityHeadersMiddleware)
    # Resolve the real client IP from the trusted reverse proxy so the
    # per-IP cooldown + access logs work once Caddy is in front.
    app.add_middleware(
        ProxyHeadersMiddleware,
        trusted_hosts=cfg.trusted_proxy_ips or None,
    )

    # Static files
    base_dir = Path(__file__).resolve().parent
    app.mount("/static", StaticFiles(directory=str(base_dir / "static")), name="static")

    @app.on_event("startup")
    async def startup() -> None:
        setup_logging(level=logging.INFO)
        attach_sse_handler_to_beguided_logs()
        # Validate every operator-configured path; warn but don't fail if
        # the genome FASTA / GTF / .npz are absent (genome-free reruns).
        validate_operator_paths(cfg)
        # Kick off the genome preload in the background. If it fails, the
        # preload thread calls os._exit(1) (fail-fast) so this process never
        # serves jobs against a missing/broken genome.
        if cfg.preload_genome and cfg.genome_fasta and Path(cfg.genome_fasta).exists():
            genome_cache.preload_async(cfg.genome_fasta, label=cfg.genome_label)
            logger.info(
                "Genome preload started in background for %s (label=%r)",
                cfg.genome_fasta,
                cfg.genome_label or "<auto>",
            )
        else:
            logger.info(
                "Genome preload disabled (preload_genome=%s, genome_fasta=%s)",
                cfg.preload_genome,
                cfg.genome_fasta,
            )
        # Same fail-fast policy for the GTF annotation: parse it once at
        # startup so step 2 (bystander analysis) does not reload per job.
        if cfg.preload_gtf and cfg.gtf_path and Path(cfg.gtf_path).exists():
            gtf_cache.preload_async(cfg.gtf_path)
            logger.info("GTF preload started in background for %s", cfg.gtf_path)
        else:
            logger.info(
                "GTF preload disabled (preload_gtf=%s, gtf_path=%s)",
                cfg.preload_gtf,
                cfg.gtf_path,
            )
        # Install and start the scheduler.
        scheduler = Scheduler(
            db_path=cfg.db_path,
            queue_max=cfg.queue_max,
            cooldown_seconds=cfg.cooldown_seconds,
            job_retention_seconds=cfg.job_retention_hours * 3600,
            log_ip_salt=cfg.log_ip_salt,
        )
        install_scheduler(scheduler)
        scheduler.start()
        logger.info(
            "BEguided web server started. queue_max=%d cooldown=%ds db=%s retention=%dh",
            cfg.queue_max,
            int(cfg.cooldown_seconds),
            cfg.db_path,
            int(cfg.job_retention_hours),
        )

    @app.get("/", response_class=HTMLResponse)
    def index() -> HTMLResponse:
        cfg = get_config()
        base_editor_names = get_base_editor_display_names()
        html = _render_template(
            "index.html",
            base_editor_names=base_editor_names,
            filter_modes=list(VALID_FILTER_MODES),
            beguidedoff_available=cfg.has_beguidedoff,
            preload_genome=cfg.preload_genome,
            preload_gtf=cfg.preload_gtf,
        )
        return HTMLResponse(content=html)

    def _render_legal(template_name: str) -> HTMLResponse:
        cfg = get_config()
        html = _render_template(
            template_name,
            legal=current_legal_config(),
            retention_hours=int(cfg.job_retention_hours),
        )
        return HTMLResponse(content=html)

    @app.get("/impressum", response_class=HTMLResponse)
    def impressum() -> HTMLResponse:
        """Impressum (TMG/DDG §5). Mandatory for any German-public web presence."""
        return _render_legal("impressum.html")

    @app.get("/datenschutz", response_class=HTMLResponse)
    def datenschutz() -> HTMLResponse:
        """Datenschutzerklärung (GDPR Art. 13)."""
        return _render_legal("datenschutz.html")

    @app.get("/agb", response_class=HTMLResponse)
    def agb() -> HTMLResponse:
        """Research-only terms of use."""
        return _render_legal("agb.html")

    @app.post("/admin/purge-all")
    async def admin_purge_all(request: Request) -> JSONResponse:
        """Manually purge every job from the DB.

        Required to honour GDPR Art. 17 (right to erasure) requests that
        arrive after the retention sweeper has already pruned the user's
        data, or to wipe the table after a security incident. Gated by a
        static ``BEGUIDED_ADMIN_TOKEN``; the endpoint refuses every
        request when the token is not configured.
        """
        cfg = get_config()
        if not cfg.admin_token:
            return JSONResponse(
                {"ok": False, "error": "Admin endpoint not configured."},
                status_code=403,
            )
        provided = request.headers.get("x-admin-token", "")
        # Constant-time comparison so a timing oracle cannot leak the token.
        if not provided or not secrets.compare_digest(provided, cfg.admin_token):
            return JSONResponse(
                {"ok": False, "error": "Invalid admin token."},
                status_code=403,
            )
        scheduler = current_scheduler()
        deleted = scheduler.purge_all_finished()
        logger.warning(
            "Admin purge-all invoked: deleted %d job(s).", deleted
        )
        return JSONResponse({"ok": True, "deleted": deleted})

    @app.get("/config")
    def server_config() -> dict:
        """Return server configuration (safe for frontend consumption)."""
        cfg = get_config()
        return {
            "ok": True,
            "beguidedoff_available": cfg.has_beguidedoff,
            "max_input_targets": cfg.max_input_targets,
            "queue_max": cfg.queue_max,
            "cooldown_seconds": int(cfg.cooldown_seconds),
            "max_upload_mb": cfg.max_upload_mb,
            "preload_genome": cfg.preload_genome,
            "preload_gtf": cfg.preload_gtf,
        }

    @app.get("/genome/status")
    def genome_status() -> dict:
        """Return the preload state of the genome FASTA."""
        cfg = get_config()
        status = genome_cache.get_status()
        return {
            "ok": True,
            "preload_enabled": cfg.preload_genome,
            **status,
        }

    @app.get("/gtf/status")
    def gtf_status() -> dict:
        """Return the preload state of the GTF annotation."""
        cfg = get_config()
        status = gtf_cache.get_status()
        return {
            "ok": True,
            "preload_enabled": cfg.preload_gtf,
            **status,
        }

    @app.get("/base-editor-table")
    def base_editor_table(source: str = "default") -> dict:
        """Return a condensed view of the base editor profiles TSV.

        ``source=default`` reads the operator-configured default file.
        ``source=custom`` reads the in-memory uploaded TSV (rejected if no
        upload has happened). ``source=combined`` returns the union of
        default + custom profiles in one list, with each row tagged via its
        ``origin`` field ("default" / "custom").

        Merges general_edit_start/end into a ``general`` column and
        high_edit_start/end into a ``high`` column, each formatted as
        ``"start-end"``.
        """
        try:
            default_df = pd.read_csv(get_config().base_editor_profiles_file, sep="\t")
        except Exception as exc:
            return {"ok": False, "error": f"Failed to load default profiles: {exc}"}

        if source == "custom":
            store = get_profiles_store()
            if not store.is_loaded():
                return {"ok": False, "error": "No custom profiles uploaded."}
            try:
                custom_df = pd.read_csv(io.BytesIO(store.get_bytes()), sep="\t")
            except Exception as exc:
                return {"ok": False, "error": f"Failed to load custom profiles: {exc}"}
            rows = _df_to_table_rows(custom_df, "custom")
            return {"ok": True, "rows": rows, "source": "custom"}

        if source == "combined":
            rows = _df_to_table_rows(default_df, "default")
            store = get_profiles_store()
            if store.is_loaded():
                try:
                    custom_df = pd.read_csv(io.BytesIO(store.get_bytes()), sep="\t")
                    # Custom profiles are listed first so the multi-select UI
                    # can render them in their own group at the top.
                    rows = _df_to_table_rows(custom_df, "custom") + rows
                except Exception:
                    pass
            return {"ok": True, "rows": rows, "source": "combined"}

        rows = _df_to_table_rows(default_df, "default")
        return {"ok": True, "rows": rows, "source": "default"}

    @app.post("/base-editor-profiles/upload")
    async def upload_profiles(file: UploadFile = File(...)) -> JSONResponse:
        """Validate and store a custom base editor profiles spreadsheet.

        Accepts ``.xlsx`` (preferred, matches the downloaded template)
        and ``.tsv`` for backward compatibility. ``.csv`` is also
        accepted defensively. The store normalises everything to TSV
        bytes so the downstream pipeline keeps reading TSV.
        """
        content = await file.read()
        suffix = Path(file.filename or "").suffix.lower()
        if not content:
            return JSONResponse(
                {"ok": False, "error": "Uploaded file is empty."},
                status_code=400,
            )
        cfg = get_config()
        try:
            df = validate_uploaded_profiles(
                content,
                suffix,
                cfg.max_profile_upload_mb * 1024 * 1024,
                default_profiles_file=cfg.base_editor_profiles_file,
            )
        except ProfilesUploadError as exc:
            return JSONResponse(
                {"ok": False, "error": str(exc)},
                status_code=exc.status_code,
            )

        try:
            load_base_editor_profiles(io.BytesIO(
                df.to_csv(sep="\t", index=False).encode("utf-8")
            ))
        except (ValueError, KeyError) as exc:
            return JSONResponse(
                {"ok": False, "error": f"Invalid profiles: {exc}"},
                status_code=400,
            )

        store = get_profiles_store()
        store.set_dataframe(df)
        names = store.names()
        return JSONResponse({
            "ok": True,
            "names": names,
            "count": len(names),
            "label": CUSTOM_EDITOR_LABEL,
            "key": CUSTOM_EDITOR_KEY,
        })

    @app.get("/base-editor-profiles/download")
    def download_profiles() -> Response:
        """Serve a base editor profiles XLSX template.

        The template contains one placeholder row (``example``) and one
        real example row (ABE8e) so users can see the expected schema.
        Columns: name, editing, general_edit_start, general_edit_end,
        high_edit_start, high_edit_end, pam, cas. The ``cas`` column is
        optional on upload; everything else is required. Users can still
        upload a plain TSV (or CSV) — the server accepts both — but the
        XLSX template is the recommended starting point.
        """
        df = pd.DataFrame(
            [
                {
                    "name": "example",
                    "editing": "ABE",
                    "general_edit_start": 2,
                    "general_edit_end": 13,
                    "high_edit_start": 3,
                    "high_edit_end": 8,
                    "pam": "NGG",
                    "cas": "SpCas9",
                },
                {
                    "name": "ABE8e",
                    "editing": "ABE",
                    "general_edit_start": 2,
                    "general_edit_end": 13,
                    "high_edit_start": 3,
                    "high_edit_end": 8,
                    "pam": "NGG",
                    "cas": "SpCas9",
                },
            ]
        )
        buf = io.BytesIO()
        df.to_excel(buf, index=False, engine="openpyxl")
        return Response(
            buf.getvalue(),
            media_type=(
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            ),
            headers={
                "Content-Disposition":
                    "attachment; filename=base_editor_profiles_template.xlsx"
            },
        )

    @app.get("/base-editor-profiles/status")
    def profiles_status() -> dict:
        """Report whether custom profiles are loaded."""
        store = get_profiles_store()
        return {
            "ok": True,
            "loaded": store.is_loaded(),
            "names": store.names(),
            "label": CUSTOM_EDITOR_LABEL,
            "key": CUSTOM_EDITOR_KEY,
        }

    @app.post("/upload")
    async def upload(file: UploadFile = File(...)) -> dict:
        """Parse uploaded TSV/CSV/XLSX and return preview."""
        content = await file.read()
        suffix = Path(file.filename or "").suffix.lower()

        if not content:
            return {"ok": False, "error": "Uploaded file is empty."}

        # Reject uploads whose extension lies about their contents. A
        # ``.tsv`` with a ZIP/OLE magic byte would otherwise be parsed as
        # TSV, surface garbled columns, and reach the chromosome validator
        # below with junk values.
        kind = _sniff_upload_kind(content, suffix)
        if kind == "binary":
            return {
                "ok": False,
                "error": (
                    f"Uploaded file contents do not match the "
                    f"{suffix or '<no extension>'} extension. "
                    "Allowed: .tsv, .csv, .xlsx, .xls."
                ),
            }

        try:
            if suffix == ".csv" or (kind == "csv" and suffix not in (".xlsx", ".xls")):
                df = pd.read_csv(io.BytesIO(content))
            elif suffix in (".xlsx", ".xls"):
                df = pd.read_excel(io.BytesIO(content))
            else:
                df = pd.read_csv(io.BytesIO(content), sep="\t")
        except Exception as exc:
            return {"ok": False, "error": f"Failed to parse file: {exc}"}

        columns = set(df.columns.str.strip())
        missing = REQUIRED_COLUMNS - columns
        if missing:
            return {
                "ok": False,
                "error": (
                    f"Missing required columns: {', '.join(sorted(missing))}. "
                    f"Found columns: {', '.join(df.columns.tolist())}"
                ),
            }

        cfg = get_config()
        if len(df) > cfg.max_input_targets:
            return {
                "ok": False,
                "error": (
                    f"Too many targets: {len(df)} (max {cfg.max_input_targets})"
                ),
            }

        df.columns = df.columns.str.strip()
        # Validate every chromosome against the loaded genome if available
        # (best-effort: if the genome is still preloading or absent the
        # check is skipped so the UI flow is not blocked).
        if genome_cache.is_ready():
            invalid = validate_chromosomes(
                df["chromosome"].astype(str).tolist(),
                _genome_chromosome_keys(),
            )
            if invalid:
                sample = ", ".join(sorted(set(invalid))[:5])
                more = "" if len(set(invalid)) <= 5 else f" (+{len(set(invalid)) - 5} more)"
                return {
                    "ok": False,
                    "error": (
                        f"Unknown chromosome(s) for the loaded genome: "
                        f"{sample}{more}. Check spelling (e.g. 'chr1' vs '1')."
                    ),
                }

        preview = df.head(100).to_dict(orient="records")
        # Normalize to a 4-column TSV before encoding. /check and /run both
        # call _decode_targets(), which parses raw_b64 as TSV (sep="\t"), so
        # passing through the original CSV/XLSX bytes would break the
        # downstream column lookups (KeyError on "target_base").
        targets_tsv = df[["chromosome", "position", "target_base", "desired_base"]]
        targets_b64 = base64.b64encode(
            targets_tsv.to_csv(sep="\t", index=False).encode()
        ).decode("ascii")
        return {
            "ok": True,
            "row_count": len(df),
            "columns": df.columns.tolist(),
            "preview": preview,
            "raw_b64": targets_b64,
        }

    @app.post("/sequence/parse")
    async def sequence_parse(text: Annotated[str, Form()] = "") -> dict:
        """Parse raw sequence text with ``[X>Y]`` brackets into targets.

        Returns the same shape as ``/upload`` (so the existing ``/check`` and
        ``/run`` flows keep working) plus a ``custom_genome_b64`` payload
        containing the FASTA representation of the parsed sequences. The
        caller passes ``custom_genome_b64`` to ``/run`` to switch the job
        into sequence-input mode.
        """
        try:
            records = parse_sequence_text(text)
        except SequenceParseError as exc:
            prefix = (
                f"Sequence #{exc.record_index + 1}: "
                if exc.record_index >= 0
                else ""
            )
            return {"ok": False, "error": prefix + str(exc)}

        cfg = get_config()
        if len(records) > cfg.max_input_targets:
            return {
                "ok": False,
                "error": (
                    f"Too many sequences: {len(records)} "
                    f"(max {cfg.max_input_targets})"
                ),
            }

        df = records_to_targets_dataframe(records)
        fasta_text = records_to_fasta(records)

        # The standard preview columns are the four core targets columns;
        # append a sequence snippet so the user can spot-check the edit
        # location without scrolling through the full raw text.
        preview_columns = [
            "chromosome",
            "position",
            "target_base",
            "desired_base",
            "sequence_preview",
        ]
        preview = df[preview_columns].head(100).to_dict(orient="records")

        # Targets DataFrame is base64-encoded for parity with /upload. We
        # serialise only the core columns so the downstream pipeline sees
        # exactly the same TSV it would have seen from the spreadsheet path.
        targets_tsv = df[["chromosome", "position", "target_base", "desired_base"]]
        targets_b64 = base64.b64encode(
            targets_tsv.to_csv(sep="\t", index=False).encode()
        ).decode("ascii")

        return {
            "ok": True,
            "row_count": len(df),
            "columns": preview_columns,
            "preview": preview,
            "raw_b64": targets_b64,
            "custom_genome_b64": base64.b64encode(fasta_text.encode()).decode("ascii"),
        }

    @app.post("/check")
    async def check(
        raw_b64: Annotated[str, Form()],
        base_editor: Annotated[str, Form()],
    ) -> dict:
        """Pre-run compatibility check: report target->base editor matching."""
        df = _decode_targets(raw_b64)

        profiles, base_editor_names, error = _load_profiles_for_check(base_editor)
        if error is not None:
            return error

        expected_keys = get_supported_editor_keys(profiles, base_editor_names)

        matched = 0
        mismatched = 0
        detail = []
        for _, row in df.iterrows():
            key = (str(row["target_base"]), str(row["desired_base"]))
            if key in expected_keys:
                matched += 1
                detail.append(
                    f"{row['chromosome']}:{row['position']} "
                    f"({row['target_base']}→{row['desired_base']}) → matched"
                )
            else:
                mismatched += 1
                detail.append(
                    f"{row['chromosome']}:{row['position']} "
                    f"({row['target_base']}→{row['desired_base']}) → NO matching base editor"
                )

        can_run = mismatched == 0
        return {
            "ok": True,
            "can_run": can_run,
            "total": len(df),
            "matched": matched,
            "mismatched": mismatched,
            "editor_types": sorted({f"{k[0]}-to-{k[1]}" for k in expected_keys}),
            "detail": detail,
        }

    @app.post("/run")
    async def run(
        request: Request,
        raw_b64: Annotated[str, Form()],
        base_editor: Annotated[str, Form()],
        mode: Annotated[str, Form()] = "spreadsheet",
        custom_genome_b64: Annotated[str, Form()] = "",
    ) -> JSONResponse:
        """Enqueue a pipeline run. Returns immediately with job_id + revoke_token.

        The web UI always runs the pipeline with all guides (no bystander /
        score filtering) and skips BEguidedOFF. The user filters the result
        table dynamically in the Results tab.

        When ``mode="sequence"`` (and ``custom_genome_b64`` is supplied) the
        job runs in sequence-input mode: the supplied custom FASTA is used
        instead of the operator-configured genome, bystander analysis is
        skipped, and off-target search is disabled.
        """
        if mode not in VALID_INPUT_MODES:
            return JSONResponse(
                {"ok": False, "error": f"Unknown input mode: {mode!r}"},
                status_code=400,
            )

        df = _decode_targets(raw_b64)

        resolved, error = _resolve_base_editor_for_run(base_editor)
        if error is not None:
            return JSONResponse(error, status_code=400)
        base_editor_names = resolved["base_editor_names"]
        use_custom_profiles = resolved["use_custom_profiles"]
        custom_profiles_b64 = resolved["custom_profiles_b64"]
        merged_profiles_b64 = resolved.get("merged_profiles_b64")

        custom_genome_b64 = custom_genome_b64.strip()
        if mode == "sequence":
            if not custom_genome_b64:
                return JSONResponse(
                    {
                        "ok": False,
                        "error": (
                            "Sequence mode requires custom_genome_b64; "
                            "parse the input via /sequence/parse first."
                        ),
                    },
                    status_code=400,
                )
        else:
            # Spreadsheet mode: validate chromosome names against the
            # loaded genome before queueing the job so the user gets an
            # immediate 400 instead of a per-target pipeline error.
            chrom_keys = _genome_chromosome_keys()
            if chrom_keys:
                invalid = validate_chromosomes(
                    df["chromosome"].astype(str).tolist(), chrom_keys
                )
                if invalid:
                    sample = ", ".join(sorted(set(invalid))[:5])
                    more = "" if len(set(invalid)) <= 5 else f" (+{len(set(invalid)) - 5} more)"
                    return JSONResponse(
                        {
                            "ok": False,
                            "error": (
                                f"Unknown chromosome(s) for the loaded genome: "
                                f"{sample}{more}. Check spelling (e.g. 'chr1' vs '1')."
                            ),
                        },
                        status_code=400,
                    )

        scheduler = current_scheduler()
        client_ip = request.client.host if request.client else "unknown"
        try:
            job: Job = scheduler.submit(
                client_ip,
                {
                    "targets_df": df,
                    "base_editor_names": base_editor_names,
                    "filter_mode": "all",
                    "search_offtargets": False,
                    "guide_score_threshold": 0.0,
                    "passthrough_extras": False,
                    "mode": mode,
                    "custom_genome_b64": custom_genome_b64 or None,
                    "use_custom_profiles": use_custom_profiles,
                    "custom_profiles_b64": custom_profiles_b64,
                    "merged_profiles_b64": merged_profiles_b64,
                },
            )
        except Cooldown as exc:
            return JSONResponse(
                {
                    "ok": False,
                    "error": (
                        f"You can submit another job in {exc.retry_after:.0f} seconds. "
                        "This prevents accidental duplicate runs."
                    ),
                    "reason": "cooldown",
                    "retry_after": int(exc.retry_after),
                },
                status_code=429,
                headers={"Retry-After": str(int(exc.retry_after))},
            )
        except QueueFull as exc:
            return JSONResponse(
                {
                    "ok": False,
                    "error": str(exc),
                    "reason": "queue_full",
                },
                status_code=503,
                headers={"Retry-After": "30"},
            )

        return JSONResponse(
            {
                "ok": True,
                "job_id": job.job_id,
                "revoke_token": job.revoke_token,
                "status": "queued",
            }
        )
        # The SSE queue is registered lazily by the /logs/stream endpoint on
        # first connect (which is fine — it's idempotent).

    @app.post("/cancel")
    async def cancel(
        job_id: Annotated[str, Form()],
        revoke_token: Annotated[str, Form()],
    ) -> JSONResponse:
        """Cancel a queued or running job. Token-gated."""
        scheduler = current_scheduler()
        ok = scheduler.cancel(job_id, revoke_token)
        if not ok:
            # Distinguish: unknown job vs bad token. Both are non-2xx.
            if scheduler.get_status(job_id) == "unknown":
                return JSONResponse(
                    {"ok": False, "error": "Unknown job_id."}, status_code=404
                )
            return JSONResponse(
                {"ok": False, "error": "Invalid revoke_token."}, status_code=403
            )
        return JSONResponse({"ok": True, "job_id": job_id, "status": scheduler.get_status(job_id)})

    @app.get("/logs/stream")
    async def log_stream(job_id: str) -> StreamingResponse:
        """SSE endpoint for real-time pipeline logs."""
        sse = get_sse_handler()
        q = sse.register_job(job_id)

        async def event_generator():
            try:
                while True:
                    try:
                        msg = await asyncio.wait_for(q.get(), timeout=30.0)
                        yield f"data: {msg}\n\n"
                    except asyncio.TimeoutError:
                        yield ":keepalive\n\n"
                        if current_scheduler().get_status(job_id) != "running":
                            break
            finally:
                sse.unregister_job(job_id)

        return StreamingResponse(
            event_generator(),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache", "Connection": "keep-alive"},
        )

    @app.get("/status")
    async def status(job_id: str) -> dict:
        """Poll endpoint to check job status."""
        return {"status": current_scheduler().get_status(job_id)}

    @app.get("/results")
    async def results(job_id: str) -> dict:
        """Return finished result records as JSON."""
        result = current_scheduler().get_result(job_id)
        if result is None:
            return {"ok": False, "error": "Unknown job_id."}
        if not result.get("records"):
            return {"ok": False, "message": result.get("message", "No results.")}
        return {
            "ok": True,
            "records": result["records"],
            "message": result.get("message"),
        }

    @app.get("/download/template")
    def download_template() -> PlainTextResponse:
        """Serve a TSV template for target input."""
        template = "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
        return PlainTextResponse(
            template,
            media_type="text/tab-separated-values",
            headers={"Content-Disposition": "attachment; filename=targets_template.tsv"},
        )

    @app.get("/download/results")
    async def download_results(
        job_id: str,
        filter_mode: Annotated[
            str, Query(description="Bystander filter: all | no_bystanders | silent")
        ] = "all",
        guide_score_threshold: Annotated[
            float,
            Query(ge=0.0, le=1.0, description="Minimum on-target guide score"),
        ] = 0.0,
        gc_content_filter: Annotated[
            str,
            Query(description="GC content tier: off | good | acceptable"),
        ] = "off",
        filter_homopolymers: Annotated[
            bool,
            Query(description="Drop guides with >=4-mer homopolymer runs"),
        ] = False,
        trinuc_edit_score_min: Annotated[
            Optional[float],
            Query(
                ge=0.0, le=1.0,
                description="Minimum Marquart et al. 2021 ABE8e trinucleotide context score (ABE only; CBE rows kept)",
            ),
        ] = None,
        condense: Annotated[
            bool,
            Query(description="Collapse rows that differ only in base_editor"),
        ] = False,
        format: Annotated[
            str,
            Query(description="Output format: tsv | xlsx"),
        ] = "tsv",
    ) -> Response:
        """Serve the final result in TSV or XLSX, applying the requested dynamic filters."""
        if filter_mode not in VALID_FILTER_MODES:
            return PlainTextResponse(
                f"Invalid filter_mode: {filter_mode!r}. "
                f"Expected one of: {', '.join(VALID_FILTER_MODES)}.",
                status_code=400,
            )
        if gc_content_filter not in VALID_GC_FILTERS:
            return PlainTextResponse(
                f"Invalid gc_content_filter: {gc_content_filter!r}. "
                f"Expected one of: {', '.join(VALID_GC_FILTERS)}.",
                status_code=400,
            )
        if format not in ("tsv", "xlsx"):
            return PlainTextResponse(
                f"Invalid format: {format!r}. Expected one of: tsv, xlsx.",
                status_code=400,
            )
        result = current_scheduler().get_result(job_id)
        if not result or not result.get("records"):
            return PlainTextResponse("Results not available.", status_code=404)
        df = pd.DataFrame.from_records(result["records"])
        df = filter_guides(df, {
            "mode": filter_mode,
            "guide_score_threshold": guide_score_threshold,
            "gc_content_filter": gc_content_filter,
            "filter_homopolymers": filter_homopolymers,
            "trinuc_edit_score_min": trinuc_edit_score_min,
        })
        if condense:
            df = condense_df(df)
        suffix_parts = []
        if filter_mode != "all" or guide_score_threshold != 0.0:
            suffix_parts.append(f"{filter_mode}_s{guide_score_threshold:g}")
        if gc_content_filter != "off":
            suffix_parts.append(f"gc-{gc_content_filter}")
        if filter_homopolymers:
            suffix_parts.append("no-homopolymers")
        if trinuc_edit_score_min is not None:
            suffix_parts.append(f"trinuc-{trinuc_edit_score_min:g}")
        if condense:
            suffix_parts.append("condensed")
        suffix = "" if not suffix_parts else "_" + "_".join(suffix_parts)
        filename = f"beguided_results{suffix}.{format}"
        if format == "xlsx":
            buf = io.BytesIO()
            df.to_excel(buf, index=False, engine="openpyxl")
            return Response(
                buf.getvalue(),
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                headers={
                    "Content-Disposition": f"attachment; filename={filename}"
                },
            )
        tsv = df.to_csv(sep="\t", index=False)
        return PlainTextResponse(
            tsv,
            media_type="text/tab-separated-values",
            headers={
                "Content-Disposition": f"attachment; filename={filename}"
            },
        )

    return app
