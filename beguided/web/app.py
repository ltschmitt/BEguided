"""FastAPI application for the BEguided web interface."""

import asyncio
import base64
import io
import logging
import uuid
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, File, Form, UploadFile
from fastapi.responses import HTMLResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from jinja2 import Environment, FileSystemLoader

from beguided.logging_config import setup_logging
from beguided.web.config import (
    get_config,
    get_preset_display_names,
    load_preset_groups,
)
from beguided.editing.presets import load_presets, get_supported_editor_keys
from beguided.web.sse_logger import attach_sse_handler_to_beguided_logs, get_sse_handler
from beguided.web.tasks import run_pipeline_task

logger = logging.getLogger("beguided.web.app")

# In-memory job storage
_jobs = {}
_jobs_lock = asyncio.Lock()

# Expected input columns
REQUIRED_COLUMNS = {"chromosome", "position", "target_base", "desired_base"}


def _render_template(name, **kwargs):
    """Render a Jinja2 template to HTML string."""
    base_dir = Path(__file__).resolve().parent
    env = Environment(loader=FileSystemLoader(str(base_dir / "templates")))
    template = env.get_template(name)
    return template.render(**kwargs)


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    app = FastAPI(title="BEguided Web")

    # Setup static files
    base_dir = Path(__file__).resolve().parent
    app.mount("/static", StaticFiles(directory=str(base_dir / "static")), name="static")

    @app.on_event("startup")
    async def startup():
        setup_logging(level=logging.INFO)
        attach_sse_handler_to_beguided_logs()
        logger.info("BEguided web server started.")

    @app.get("/", response_class=HTMLResponse)
    def index():
        cfg = get_config()
        preset_names = get_preset_display_names()
        html = _render_template(
            "index.html",
            preset_names=preset_names,
            filter_modes=["all", "no_bystanders", "silent"],
            beguidedoff_available=cfg.has_beguidedoff,
        )
        return HTMLResponse(content=html)

    @app.get("/config")
    def server_config():
        """Return server configuration (safe for frontend consumption)."""
        cfg = get_config()
        return {
            "ok": True,
            "beguidedoff_available": cfg.has_beguidedoff,
            "max_input_targets": cfg.max_input_targets,
        }

    @app.post("/upload")
    async def upload(file: UploadFile = File(...)):
        """Parse uploaded TSV/CSV/XLSX and return preview."""
        content = await file.read()
        suffix = Path(file.filename or "").suffix.lower()

        try:
            if suffix == ".csv":
                df = pd.read_csv(io.BytesIO(content))
            elif suffix in (".xlsx", ".xls"):
                df = pd.read_excel(io.BytesIO(content))
            else:
                # Default to TSV
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

        # Sanitize column names
        df.columns = df.columns.str.strip()

        # Return first 100 rows for preview
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
        raw_b64: str = Form(...),
        preset: str = Form(...),
    ):
        """Pre-run compatibility check: report target->preset matching."""
        content = base64.b64decode(raw_b64.encode("ascii"))
        df = pd.read_csv(io.BytesIO(content), sep="\t")
        df.columns = df.columns.str.strip()

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
        raw_b64: str = Form(...),
        preset: str = Form(...),
        filter_mode: str = Form(...),
        guide_score_threshold: str = Form("0.5"),
        search_offtargets: str = Form("true"),
    ):
        """Start pipeline run in background and return job_id."""
        content = base64.b64decode(raw_b64.encode("ascii"))
        df = pd.read_csv(io.BytesIO(content), sep="\t")
        df.columns = df.columns.str.strip()

        try:
            score_threshold = float(guide_score_threshold)
        except ValueError:
            return {"ok": False, "error": f"Invalid guide_score_threshold: {guide_score_threshold!r}"}

        job_id = str(uuid.uuid4())
        groups = load_preset_groups()
        preset_names = groups.get(preset, [preset])
        run_offtargets = search_offtargets.lower() in ("true", "1", "yes", "on")

        # Fast-fail if off-targets requested but unavailable
        from beguided.web.tasks import _validate_offtarget_config
        try:
            _validate_offtarget_config(run_offtargets)
        except RuntimeError as exc:
            return {"ok": False, "error": str(exc)}

        task = asyncio.create_task(
            run_pipeline_task(
                job_id=job_id,
                targets_df=df,
                preset_names=preset_names,
                filter_mode=filter_mode,
                guide_score_threshold=score_threshold,
                sse_handler=get_sse_handler(),
                search_offtargets=run_offtargets,
            )
        )

        async with _jobs_lock:
            _jobs[job_id] = {
                "task": task,
                "status": "running",
                "result": None,
            }

        # Attach completion callback
        def on_done(t):
            asyncio.create_task(_store_result(job_id, t))

        task.add_done_callback(on_done)

        return {"ok": True, "job_id": job_id}

    async def _store_result(job_id, task):
        try:
            result = task.result()
        except Exception as exc:
            result = {
                "job_id": job_id,
                "status": "error",
                "message": str(exc),
                "records": [],
            }
        async with _jobs_lock:
            if job_id in _jobs:
                _jobs[job_id]["status"] = result.get("status", "error")
                _jobs[job_id]["result"] = result

    @app.get("/logs/stream")
    async def log_stream(job_id: str):
        """SSE endpoint for real-time pipeline logs."""
        q = get_sse_handler().register_job(job_id)

        async def event_generator():
            try:
                while True:
                    try:
                        msg = await asyncio.wait_for(q.get(), timeout=30.0)
                        yield f"data: {msg}\n\n"
                    except asyncio.TimeoutError:
                        yield ":keepalive\n\n"
                        # Check if job is done
                        async with _jobs_lock:
                            job = _jobs.get(job_id)
                        if job and job["status"] != "running":
                            break
            finally:
                get_sse_handler().unregister_job(job_id)

        return StreamingResponse(
            event_generator(),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache", "Connection": "keep-alive"},
        )

    @app.get("/status")
    async def status(job_id: str):
        """Poll endpoint to check job status."""
        async with _jobs_lock:
            job = _jobs.get(job_id)
        if not job:
            return {"status": "unknown"}
        return {"status": job["status"]}

    @app.get("/results")
    async def results(job_id: str):
        """Return finished result records as JSON."""
        async with _jobs_lock:
            job = _jobs.get(job_id)
        if not job or not job["result"]:
            return {"ok": False, "error": "Results not available."}
        return {
            "ok": True,
            "records": job["result"]["records"],
            "message": job["result"]["message"],
        }

    @app.get("/download/template")
    def download_template():
        """Serve a TSV template for target input."""
        template = "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
        return PlainTextResponse(
            template,
            media_type="text/tab-separated-values",
            headers={"Content-Disposition": "attachment; filename=targets_template.tsv"},
        )

    @app.get("/download/results")
    async def download_results(job_id: str):
        """Serve the final result TSV."""
        async with _jobs_lock:
            job = _jobs.get(job_id)
        if not job or not job["result"]:
            return PlainTextResponse("Results not available.", status_code=404)
        records = job["result"]["records"]
        if not records:
            return PlainTextResponse("No records.", status_code=404)
        df = pd.DataFrame.from_records(records)
        tsv = df.to_csv(sep="\t", index=False)
        return PlainTextResponse(
            tsv,
            media_type="text/tab-separated-values",
            headers={
                "Content-Disposition": "attachment; filename=beguided_results.tsv"
            },
        )

    return app
