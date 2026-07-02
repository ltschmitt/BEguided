"""Sync job execution helpers.

The web server runs heavy pipeline work on a single dedicated worker thread
(see ``beguided.web.scheduler``). This module owns the *what* of one job:
build the pipeline config from web parameters, run it, and collect the
result records back from disk into a JSON-serialisable dict.
"""

from __future__ import annotations

import logging
import tempfile
import threading
import traceback
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from beguided.pipeline import PipelineCancelled, run_pipeline
from beguided.web.config import get_config

logger = logging.getLogger("beguided.web.tasks")


def _build_config(
    targets_df: pd.DataFrame,
    preset_names: list[str],
    filter_mode: str,
    work_dir: str,
    search_offtargets: bool = True,
    guide_score_threshold: float = 0.5,
    passthrough_extras: bool = False,
) -> dict:
    """Build a pipeline config dict from web UI parameters.

    Args:
        targets_df: DataFrame with target rows.
        preset_names: List of preset name strings.
        filter_mode: One of 'all', 'no_bystanders', 'silent'.
        work_dir: Working directory for outputs.
        search_offtargets: Whether to run BEguidedOFF off-target steps.
        guide_score_threshold: Minimum on-target guide score (0-1).
        passthrough_extras: Whether non-core columns from the input TSV should
            be carried through to the final output. Defaults to ``False`` in
            the web app.

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
        "guide_score_threshold": guide_score_threshold,
        "search_offtargets": search_offtargets,
        "passthrough_extras": passthrough_extras,
        "offtargets": {
            "beguidedoff_genome": cfg.beguidedoff_genome,
            "mismatch_tolerance": cfg.mismatch_tolerance,
            "count_threshold_score": cfg.count_threshold_score,
            "remove_guides_score": cfg.remove_guides_score,
            "remove_guides_n_offtargets": cfg.remove_guides_n_offtargets,
        },
    }


def _validate_offtarget_config(search_offtargets: bool) -> None:
    """Raise a clear error if off-target search is requested but unavailable.

    Args:
        search_offtargets: Whether the UI/user requested off-target search.

    Raises:
        RuntimeError: If BEguidedOFF is not installed or not configured.
    """
    if not search_offtargets:
        return
    cfg = get_config()
    if not cfg.has_beguidedoff:
        raise RuntimeError(
            "Off-target search is unavailable because BEguidedOFF is not installed. "
            "Install it from ssh://git@gitlab.mn.tu-dresden.de/buchholzlab/beguidedoff.git "
            "or uncheck the Off-target option."
        )
    if not cfg.beguidedoff_genome:
        raise RuntimeError(
            "Off-target search requires --beguidedoff-genome or BEGUIDED_BEGUIDEDOFF_GENOME."
        )


def _df_to_records(path: Path) -> list[dict]:
    """Read a TSV from ``path`` and return JSON-safe records."""
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    out = []
    for row in df.to_dict(orient="records"):
        record = {k: (None if pd.isna(v) else v) for k, v in row.items()}
        # Limit guide_score to 3 decimal places for the web UI and
        # downloads; downstream filtering compares against the rounded value.
        score = record.get("guide_score")
        if isinstance(score, (int, float)) and not isinstance(score, bool):
            record["guide_score"] = round(float(score), 3)
        out.append(record)
    return out


def execute_job(
    job_id: str,
    targets_df: pd.DataFrame,
    preset_names: list[str],
    filter_mode: str,
    search_offtargets: bool = True,
    guide_score_threshold: float = 0.5,
    cancel_event: Optional[threading.Event] = None,
    passthrough_extras: bool = False,
) -> dict:
    """Run a single job synchronously. Called from the scheduler worker thread.

    Builds a per-job tempdir, runs ``run_pipeline`` with cooperative cancel,
    and reads the final TSV back into records. On exception, the traceback is
    forwarded to the SSE handler (if any) so the operator can see it in the
    browser.

    Returns:
        dict with keys: job_id, status, message, records.
    """
    from beguided.web.sse_logger import get_sse_handler
    from beguided.web.safety import validate_chromosomes

    sse = get_sse_handler()

    _validate_offtarget_config(search_offtargets)

    # Cheap pre-flight: if the targets DataFrame is empty, fail early without
    # loading the genome.
    if len(targets_df) == 0:
        return {
            "job_id": job_id,
            "status": "error",
            "message": "No targets supplied.",
            "records": [],
        }

    # Validate chromosomes against the configured genome if it's available.
    # Genome-free reruns (steps 4-5) skip this check.
    cfg = get_config()
    genome: Optional[dict] = None
    if cfg.genome_fasta and Path(cfg.genome_fasta).exists():
        try:
            if cfg.preload_genome:
                # Block the worker thread until the background preload
                # completes (or fails). The single-worker scheduler means
                # only the very first job waits; subsequent jobs see the
                # already-cached genome.
                from beguided.web.genome_cache import wait_until_ready

                genome = wait_until_ready(timeout=None)
            else:
                from beguided.genome.io import load_genome

                genome = load_genome(cfg.genome_fasta)
            unknown = validate_chromosomes(
                targets_df["chromosome"].astype(str).tolist(), genome.keys()
            )
            if unknown:
                return {
                    "job_id": job_id,
                    "status": "error",
                    "message": (
                        f"Unknown chromosome(s): {', '.join(sorted(set(unknown)))}. "
                        f"Available: {', '.join(sorted(genome.keys()))}"
                    ),
                    "records": [],
                }
        except Exception as exc:
            return {
                "job_id": job_id,
                "status": "error",
                "message": f"Failed to load genome for validation: {exc}",
                "records": [],
            }

    # Load (or wait for) the GTF interval trees if a GTF is configured.
    # Step 2 reuses the preloaded dict; otherwise it falls back to parsing
    # the file itself (the pipeline will still work when preload_gtf=False).
    interval_trees: Optional[dict] = None
    if cfg.gtf_path and Path(cfg.gtf_path).exists():
        try:
            if cfg.preload_gtf:
                from beguided.web.gtf_cache import wait_until_ready

                interval_trees = wait_until_ready(timeout=None)
            else:
                from beguided.bystander.annotation import load_gtf_to_interval_tree

                interval_trees = load_gtf_to_interval_tree(cfg.gtf_path)
        except Exception as exc:
            return {
                "job_id": job_id,
                "status": "error",
                "message": f"Failed to load GTF annotation: {exc}",
                "records": [],
            }

    with tempfile.TemporaryDirectory(prefix="beguided_") as work_dir:
        config = _build_config(
            targets_df,
            preset_names,
            filter_mode,
            work_dir,
            search_offtargets,
            guide_score_threshold,
            passthrough_extras,
        )
        final_output = Path(work_dir) / "beguided_results.tsv"

        try:
            run_pipeline(
                config,
                start_step=1,
                end_step=5,
                cancel_event=cancel_event,
                genome=genome,
                interval_trees=interval_trees,
            )
            records = _df_to_records(final_output)
            if not records:
                return {
                    "job_id": job_id,
                    "status": "error",
                    "message": "Pipeline finished but output file was not found.",
                    "records": [],
                }
            return {
                "job_id": job_id,
                "status": "complete",
                "message": f"Pipeline complete: {len(records)} rows.",
                "records": records,
            }
        except PipelineCancelled:
            return {
                "job_id": job_id,
                "status": "cancelled",
                "message": "Cancelled by user.",
                "records": [],
            }
        except Exception as exc:
            tb = traceback.format_exc()
            for line in tb.splitlines():
                sse.push(job_id, line)
            logger.exception("Job %s failed: %s", job_id, exc)
            return {
                "job_id": job_id,
                "status": "error",
                "message": str(exc),
                "records": [],
            }
