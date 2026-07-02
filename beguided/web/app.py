"""FastAPI application for the BEguided web interface."""

from __future__ import annotations

import asyncio
import base64
import io
import logging
from pathlib import Path
from typing import Annotated

import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, Query, Request, Response, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from jinja2 import Environment, FileSystemLoader

from beguided.filtering.filters import filter_guides
from beguided.logging_config import setup_logging
from beguided.web.config import (
    get_config,
    get_preset_display_names,
    load_preset_groups,
)
from beguided.editing.presets import load_presets, get_supported_editor_keys
from beguided.web import genome_cache, gtf_cache
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
    validate_operator_paths,
)
from beguided.web.sse_logger import attach_sse_handler_to_beguided_logs, get_sse_handler

logger = logging.getLogger("beguided.web.app")

REQUIRED_COLUMNS = {"chromosome", "position", "target_base", "desired_base"}
VALID_FILTER_MODES = ("all", "no_bystanders", "silent")


def _render_template(name, **kwargs) -> str:
    """Render a Jinja2 template to HTML string."""
    base_dir = Path(__file__).resolve().parent
    env = Environment(loader=FileSystemLoader(str(base_dir / "templates")))
    template = env.get_template(name)
    return template.render(**kwargs)


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


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    cfg = get_config()
    app = FastAPI(title="BEguided Web")

    # Cap request body size to keep memory bounded.
    app.add_middleware(MaxBodySizeMiddleware, max_bytes=cfg.max_upload_mb * 1024 * 1024)

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
        )
        install_scheduler(scheduler)
        scheduler.start()
        logger.info(
            "BEguided web server started. queue_max=%d cooldown=%ds db=%s",
            cfg.queue_max,
            int(cfg.cooldown_seconds),
            cfg.db_path,
        )

    @app.get("/", response_class=HTMLResponse)
    def index() -> HTMLResponse:
        cfg = get_config()
        preset_names = get_preset_display_names()
        html = _render_template(
            "index.html",
            preset_names=preset_names,
            filter_modes=list(VALID_FILTER_MODES),
            beguidedoff_available=cfg.has_beguidedoff,
            preload_genome=cfg.preload_genome,
            preload_gtf=cfg.preload_gtf,
        )
        return HTMLResponse(content=html)

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

    @app.get("/preset-table")
    def preset_table() -> dict:
        """Return a condensed view of the presets TSV for the Guide tab.

        Drops guide_length since most are 20 nt. Merges
        general_edit_start/end into a ``general`` column and
        high_edit_start/end into a ``high`` column, each formatted as
        ``"start-end"``.
        """
        try:
            df = pd.read_csv(get_config().presets_file, sep="\t")
        except Exception as exc:
            return {"ok": False, "error": f"Failed to load presets: {exc}"}
        rows = []
        for _, r in df.iterrows():
            rows.append({
                "name": str(r["name"]),
                "editing": str(r["editing"]),
                "general": f"{int(r['general_edit_start'])}-{int(r['general_edit_end'])}",
                "high": f"{int(r['high_edit_start'])}-{int(r['high_edit_end'])}",
                "pam": str(r["pam"]),
                "cas": str(r["cas"]),
            })
        return {"ok": True, "rows": rows}

    @app.post("/upload")
    async def upload(file: UploadFile = File(...)) -> dict:
        """Parse uploaded TSV/CSV/XLSX and return preview."""
        content = await file.read()
        suffix = Path(file.filename or "").suffix.lower()

        try:
            if suffix == ".csv":
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
        preview = df.head(100).to_dict(orient="records")
        return {
            "ok": True,
            "row_count": len(df),
            "columns": df.columns.tolist(),
            "preview": preview,
            "raw_b64": base64.b64encode(content).decode("ascii"),
        }

    @app.post("/check")
    async def check(
        raw_b64: Annotated[str, Form()],
        preset: Annotated[str, Form()],
    ) -> dict:
        """Pre-run compatibility check: report target->preset matching."""
        df = _decode_targets(raw_b64)

        cfg = get_config()
        groups = load_preset_groups()
        preset_names = groups.get(preset, [preset])
        presets = load_presets(cfg.presets_file)
        expected_keys = get_supported_editor_keys(presets, preset_names)

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
                    f"({row['target_base']}→{row['desired_base']}) → NO matching preset"
                )

        can_run = mismatched == 0
        return {
            "ok": True,
            "can_run": can_run,
            "total": len(df),
            "matched": matched,
            "mismatched": mismatched,
            "preset_types": sorted({f"{k[0]}-to-{k[1]}" for k in expected_keys}),
            "detail": detail,
        }

    @app.post("/run")
    async def run(
        request: Request,
        raw_b64: Annotated[str, Form()],
        preset: Annotated[str, Form()],
    ) -> JSONResponse:
        """Enqueue a pipeline run. Returns immediately with job_id + revoke_token.

        The web UI always runs the pipeline with all guides (no bystander /
        score filtering) and skips BEguidedOFF. The user filters the result
        table dynamically in the Results tab.
        """
        df = _decode_targets(raw_b64)

        groups = load_preset_groups()
        if preset not in groups:
            return JSONResponse(
                {"ok": False, "error": f"Unknown preset group: {preset!r}"},
                status_code=400,
            )
        preset_names = groups[preset]

        scheduler = current_scheduler()
        client_ip = request.client.host if request.client else "unknown"
        try:
            job: Job = scheduler.submit(
                client_ip,
                {
                    "targets_df": df,
                    "preset_names": preset_names,
                    "filter_mode": "all",
                    "search_offtargets": False,
                    "guide_score_threshold": 0.0,
                    "passthrough_extras": False,
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
        })
        suffix = "" if filter_mode == "all" and guide_score_threshold == 0.0 else (
            f"_{filter_mode}_s{guide_score_threshold:g}"
        )
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
