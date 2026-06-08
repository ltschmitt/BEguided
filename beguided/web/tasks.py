"""Async task wrapper for running the BEguided pipeline in the background."""

import asyncio
import logging
import tempfile
import traceback
from pathlib import Path

import pandas as pd

from beguided.pipeline import run_pipeline
from beguided.web.config import get_config

logger = logging.getLogger("beguided.web.tasks")


def _build_config(targets_df, preset_names, filter_mode, work_dir, search_offtargets=True):
    """Build a pipeline config dict from web UI parameters.

    Args:
        targets_df: DataFrame with target rows.
        preset_names: List of preset name strings.
        filter_mode: One of 'all', 'no_bystanders', 'silent'.
        work_dir: Temporary directory for outputs.
        search_offtargets: Whether to run BEguidedOFF off-target steps.

    Returns:
        dict: Compatible with run_pipeline().
    """
    cfg = get_config()
    targets_path = Path(work_dir) / "targets.tsv"
    targets_df.to_csv(targets_path, sep="\t", index=False)

    return {
        "input": {
            "targets": str(targets_path),
            "genome": cfg.genome_fasta,
            "presets": cfg.presets_file,
            "gff_file": cfg.gtf_path,
        },
        "output_prefix": str(work_dir),
        "use_presets": ",".join(preset_names),
        "keep_bystanders": filter_mode,
        "search_offtargets": search_offtargets,
        "offtargets": {
            "beguidedoff_genome": cfg.beguidedoff_genome,
            "mismatch_tolerance": cfg.mismatch_tolerance,
            "count_threshold_score": cfg.count_threshold_score,
            "remove_guides_score": cfg.remove_guides_score,
            "remove_guides_n_offtargets": cfg.remove_guides_n_offtargets,
        },
    }


async def run_pipeline_task(job_id, targets_df, preset_names, filter_mode, sse_handler, search_offtargets=True):
    """Run the pipeline asynchronously and stream logs via SSE.

    Args:
        job_id: Unique job identifier.
        targets_df: Input targets DataFrame.
        preset_names: Resolved list of preset names.
        filter_mode: Filter mode string.
        sse_handler: AsyncQueueHandler instance.
        search_offtargets: Whether to run BEguidedOFF off-target steps.

    Returns:
        dict: Result metadata including path to final output.
    """
    loop = asyncio.get_running_loop()
    q = sse_handler.register_job(job_id)

    with tempfile.TemporaryDirectory(prefix="beguided_") as work_dir:
        config = _build_config(targets_df, preset_names, filter_mode, work_dir, search_offtargets)
        final_output = Path(work_dir) / "beguided_results.tsv"

        try:
            await loop.run_in_executor(
                None,
                lambda: run_pipeline(config, start_step=1, end_step=5)
            )

            if final_output.exists():
                result_df = pd.read_csv(final_output, sep="\t")
                result_records = [
                    {k: (None if pd.isna(v) else v) for k, v in row.items()}
                    for row in result_df.to_dict(orient="records")
                ]
                status = "complete"
                message = f"Pipeline complete: {len(result_df)} rows."
            else:
                result_records = []
                status = "error"
                message = "Pipeline finished but output file was not found."
        except Exception as exc:
            error_text = traceback.format_exc()
            for line in error_text.splitlines():
                try:
                    q.put_nowait(line)
                except asyncio.QueueFull:
                    pass
            result_records = []
            status = "error"
            message = str(exc)
        finally:
            sse_handler.unregister_job(job_id)

    return {
        "job_id": job_id,
        "status": status,
        "message": message,
        "records": result_records,
    }
