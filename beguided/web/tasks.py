"""Sync job execution helpers.

The web server runs heavy pipeline work on a single dedicated worker thread
(see ``beguided.web.scheduler``). This module owns the *what* of one job:
build the pipeline config from web parameters, run it, and collect the
result records back from disk into a JSON-serialisable dict.
"""

from __future__ import annotations

import base64
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

# Columns that are dropped from the records returned to the web UI for
# sequence-mode jobs. ``genome_base`` is redundant with the user-supplied
# input (``genome_base == target_base`` by construction). Everything else —
# including ``position`` and ``strand_used`` — is kept so the records stay
# structurally identical to spreadsheet-mode records and can be screened
# for off-targets with the shared row-key logic; the sequence-mode results
# table hides the sequence-relative columns client-side. The on-disk TSV is
# left untouched.
SEQUENCE_MODE_DROP_COLUMNS: tuple[str, ...] = (
    "genome_base",
)


def _build_config(
    targets_df: pd.DataFrame,
    base_editor_names: list[str],
    filter_mode: str,
    work_dir: str,
    search_offtargets: bool = True,
    guide_score_threshold: float = 0.5,
    passthrough_extras: bool = False,
    genome_fasta: Optional[str] = None,
    gtf_path: Optional[str] = None,
    skip_bystander: bool = False,
    use_custom_profiles: bool = False,
    custom_profiles_b64: Optional[str] = None,
    merged_profiles_b64: Optional[str] = None,
    gc_content_filter: Optional[str] = None,
    filter_homopolymers: Optional[bool] = None,
    trinuc_edit_score_min: Optional[float] = None,
    keep_targets_without_guides: bool = False,
    genome_key: str = "human",
) -> dict:
    """Build a pipeline config dict from web UI parameters.

    Args:
        targets_df: DataFrame with target rows.
        base_editor_names: List of base editor name strings.
        filter_mode: One of 'all', 'none', 'synonymous'.
        work_dir: Working directory for outputs.
        search_offtargets: Whether to run BEguidedOFF off-target steps.
        guide_score_threshold: Minimum on-target guide score (0-1).
        passthrough_extras: Whether non-core columns from the input TSV should
            be carried through to the final output. Defaults to ``False`` in
            the web app.
        genome_fasta: Path to the FASTA file to use as the genome. Defaults to
            the operator-configured ``cfg.genome_fasta``. Pass a custom path
            (e.g. a sequence-mode job's temp FASTA) to override it.
        gtf_path: Path to the GTF annotation to use. Defaults to the
            operator-configured ``cfg.gtf_path``; the per-genome selection
            passes the chosen genome's GTF here.
        skip_bystander: When True, skip the bystander analysis step AND skip
            the bystander filter in step 3. Used by sequence-mode jobs
            (no GTF annotation available).
        use_custom_profiles: When True, decode ``custom_profiles_b64`` and
            write it to ``work_dir / custom_profiles.tsv``; that file replaces
            the operator-configured default for this job and every profile it
            contains is considered (via ``use_base_editors="all"``).
        custom_profiles_b64: Base64-encoded TSV bytes of the user-uploaded
            profiles. Required when ``use_custom_profiles`` is True.
        merged_profiles_b64: Base64-encoded TSV bytes of a merged profiles
            file combining custom + default editors. When provided, written
            to ``work_dir / merged_profiles.tsv`` and used as the profiles
            file for this job (with ``base_editor_names`` controlling which
            profiles are considered).
        gc_content_filter: One of ``"off"``, ``"high"`` (40-60%) or
            ``"medium"`` (30-70%); the deprecated names ``"good"`` /
            ``"acceptable"`` are accepted as aliases. When ``None`` the
            operator-configured default (``cfg.gc_content_filter``) is used.
        filter_homopolymers: When True, drop guides with a >=4-mer
            homopolymer run. When ``None`` the operator-configured default
            (``cfg.filter_homopolymers``) is used.
        trinuc_edit_score_min: Minimum Marquart et al. 2021 normalised ABE8e
            trinucleotide context score (0..1) to keep an ABE guide. CBE
            rows are always kept (NaN). When ``None`` the operator-configured
            default (``cfg.trinuc_edit_score_min``) is used.
        keep_targets_without_guides: When True the final results retain one
            placeholder row per input target that produced no guide, so the
            web UI can always show a row per target even when no guide RNA
            could be designed (or every guide was filtered out). Off by
            default; the web server always sets this.

    Returns:
        dict: Compatible with run_pipeline().
    """
    cfg = get_config()
    targets_path = Path(work_dir) / "targets.tsv"
    targets_df.to_csv(targets_path, sep="\t", index=False)

    # BEguidedOFF configuration for the selected genome: per-genome env vars
    # (BEGUIDED_BEGUIDEDOFF_GENOME_<KEY> / ..._INDEX_<KEY>) win over the
    # global setting, which only falls back for the default genome (and for
    # single-genome setups) so no species is screened against another
    # species' index.
    bgo = cfg.beguidedoff_resolved(genome_key)

    if merged_profiles_b64:
        profiles_path = Path(work_dir) / "merged_profiles.tsv"
        profiles_path.write_bytes(base64.b64decode(merged_profiles_b64))
        profiles_file = str(profiles_path)
        use_base_editors = ",".join(base_editor_names)
    elif use_custom_profiles:
        if not custom_profiles_b64:
            raise ValueError(
                "use_custom_profiles=True requires custom_profiles_b64."
            )
        profiles_path = Path(work_dir) / "custom_profiles.tsv"
        profiles_path.write_bytes(base64.b64decode(custom_profiles_b64))
        profiles_file = str(profiles_path)
        use_base_editors = "all"
    else:
        profiles_file = cfg.base_editor_profiles_file
        use_base_editors = ",".join(base_editor_names)

    if gc_content_filter is None:
        gc_content_filter = cfg.gc_content_filter
    if filter_homopolymers is None:
        filter_homopolymers = cfg.filter_homopolymers
    if trinuc_edit_score_min is None:
        trinuc_edit_score_min = cfg.trinuc_edit_score_min

    return {
        "input": {
            "targets": str(targets_path),
            "genome": genome_fasta if genome_fasta is not None else cfg.genome_fasta,
            "base_editor_profiles": profiles_file,
            "gff_file": gtf_path if gtf_path is not None else cfg.gtf_path,
        },
        "output_prefix": str(work_dir),
        "use_base_editors": use_base_editors,
        "keep_bystanders": filter_mode,
        "guide_score_threshold": guide_score_threshold,
        "gc_content_filter": gc_content_filter,
        "filter_homopolymers": bool(filter_homopolymers),
        "trinuc_edit_score_min": trinuc_edit_score_min,
        "keep_targets_without_guides": bool(keep_targets_without_guides),
        "search_offtargets": search_offtargets,
        "passthrough_extras": passthrough_extras,
        "skip_bystander": skip_bystander,
        "offtargets": {
            "mode": bgo["mode"],
            "beguidedoff_genome": bgo["genome"],
            "index_path": bgo["index"],
            # The pipeline key is max_mismatches; ``mismatch_tolerance`` was
            # the old web-side name under which the value was silently
            # ignored by the pipeline.
            "max_mismatches": cfg.mismatch_tolerance,
            "count_threshold_score": cfg.count_threshold_score,
            "remove_guides_score": cfg.remove_guides_score,
            "remove_guides_n_offtargets": cfg.remove_guides_n_offtargets,
        },
    }


def _validate_offtarget_config(search_offtargets: bool, genome_key: str = "human") -> None:
    """Raise a clear error if off-target search is requested but unavailable.

    Args:
        search_offtargets: Whether the UI/user requested off-target search.
        genome_key: Which selectable genome the search will run against; the
            encoded genome / FM-index are resolved per genome.

    Raises:
        RuntimeError: If BEguidedOFF is not installed or not configured.
    """
    if not search_offtargets:
        return
    cfg = get_config()
    if not cfg.has_beguidedoff:
        raise RuntimeError(
            "Off-target search is unavailable because BEguidedOFF is not installed. "
            "Install it from https://github.com/ltschmitt/BEguidedOFF "
            "or uncheck the Off-target option."
        )
    bgo = cfg.beguidedoff_resolved(genome_key)
    if bgo["mode"] == "extensive" and not bgo["genome"]:
        raise RuntimeError(
            "Off-target search requires --beguidedoff-genome or "
            f"BEGUIDED_BEGUIDEDOFF_GENOME_{genome_key.upper()} "
            "(or the global BEGUIDED_BEGUIDEDOFF_GENOME)."
        )
    if bgo["mode"] == "fast":
        index = bgo["index"]
        if not index:
            raise RuntimeError(
                "fast-mode off-target search requires an FM-index: set "
                "--beguidedoff-index / BEGUIDED_BEGUIDEDOFF_INDEX (or "
                "--beguidedoff-genome / BEGUIDED_BEGUIDEDOFF_GENOME to "
                "derive the index path from the encoded .npz)."
            )
        if not Path(index).exists():
            raise RuntimeError(
                "fast-mode off-target search requires an FM-index (.fmi); "
                f"not found: {index}. Build one with `python -m beguidedoff.cli "
                "build-index <genome.fa> <genome.fmi>` and point "
                "--beguidedoff-index / BEGUIDED_BEGUIDEDOFF_INDEX at it, or "
                "switch to --beguidedoff-mode extensive (GPU)."
            )
    for name in (bgo["genome"], bgo["index"]):
        p = Path(name) if name else None
        if p is not None and ("mini_genome" in p.name or "tests" in p.parts):
            logger.warning(
                "Off-target search is configured against the bundled mini "
                "test genome (%s), not a real reference. Set "
                "BEGUIDED_BEGUIDEDOFF_INDEX (fast mode) or "
                "BEGUIDED_BEGUIDEDOFF_GENOME (extensive mode) for "
                "production use.",
                p,
            )
            break


def _df_to_records(path: Path, drop_columns: tuple[str, ...] = ()) -> list[dict]:
    """Read a TSV from ``path`` and return JSON-safe records.

    ``drop_columns`` lists columns that should be removed from each record
    before it is returned. Used by sequence-mode jobs to hide the one column
    that carries no information beyond what the user supplied on the input
    side (``genome_base`` is equal to the target base by construction). The
    on-disk TSV is left untouched.
    """
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    if drop_columns:
        df = df.drop(columns=[c for c in drop_columns if c in df.columns])
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
    base_editor_names: list[str],
    filter_mode: str,
    search_offtargets: bool = True,
    guide_score_threshold: float = 0.5,
    cancel_event: Optional[threading.Event] = None,
    passthrough_extras: bool = False,
    keep_targets_without_guides: bool = False,
    genome_key: str = "human",
    mode: str = "spreadsheet",
    custom_genome_b64: Optional[str] = None,
    use_custom_profiles: bool = False,
    custom_profiles_b64: Optional[str] = None,
    merged_profiles_b64: Optional[str] = None,
) -> dict:
    """Run a single job synchronously. Called from the scheduler worker thread.

    Builds a per-job tempdir, runs ``run_pipeline`` with cooperative cancel,
    and reads the final TSV back into records. On exception, the traceback is
    forwarded to the SSE handler (if any) so the operator can see it in the
    browser.

    Args:
        mode: ``"spreadsheet"`` (default), ``"single"`` (one manually-entered
            genomic target; behaves exactly like spreadsheet mode) or
            ``"sequence"``. In sequence mode the caller must also supply
            ``custom_genome_b64``; the job then uses the supplied custom
            FASTA instead of the configured genome, bypasses GTF loading,
            and forces bystander analysis + off-target search off.
        genome_key: Which configured genome to use in spreadsheet mode
            (``"human"`` / ``"mouse"``). Ignored in sequence mode.
        custom_genome_b64: Base64-encoded FASTA text containing the user's
            sequences. Required when ``mode="sequence"``.
        use_custom_profiles: When True, ``custom_profiles_b64`` is decoded
            into ``work_dir / custom_profiles.tsv`` and used as the profiles
            file for this job (with all profiles considered).
        custom_profiles_b64: Base64-encoded profiles TSV. Required when
            ``use_custom_profiles`` is True.

    Returns:
        dict with keys: job_id, status, message, records.
    """
    from beguided.web.sse_logger import get_sse_handler
    from beguided.web.safety import validate_chromosomes

    sse = get_sse_handler()

    sequence_mode = mode == "sequence"
    if sequence_mode and not custom_genome_b64:
        return {
            "job_id": job_id,
            "status": "error",
            "message": "Sequence mode requires custom_genome_b64.",
            "records": [],
        }

    # Sequence mode has no GTF annotation and no pre-encoded off-target
    # genome, so we (a) skip the bystander step entirely (via the explicit
    # ``skip_bystander`` pipeline flag, NOT by setting filter_mode to
    # "none", which would also drop every guide with bystanders)
    # and (b) disable BEguidedOFF.
    skip_bystander = sequence_mode
    if sequence_mode:
        search_offtargets = False

    _validate_offtarget_config(search_offtargets, genome_key=genome_key)

    # Cheap pre-flight: if the targets DataFrame is empty, fail early without
    # loading the genome.
    if len(targets_df) == 0:
        return {
            "job_id": job_id,
            "status": "error",
            "message": "No targets supplied.",
            "records": [],
        }

    cfg = get_config()
    genome: Optional[dict] = None
    custom_fasta_path: Optional[str] = None

    # Resolve the selected genome. ``genome_option`` falls back to the
    # configured default when the key is unknown; a truly empty genome set
    # means nothing can run in spreadsheet mode.
    genome_option = cfg.genome_option(genome_key)
    if genome_option is None and not sequence_mode:
        return {
            "job_id": job_id,
            "status": "error",
            "message": f"Unknown genome {genome_key!r}.",
            "records": [],
        }

    if sequence_mode:
        # Decode the custom FASTA and write it to a temp file. The pipeline
        # accepts a path via ``cfg["input"]["genome"]`` and the existing
        # chromosome validation keys off ``genome.keys()`` — both work with
        # a freshly-written temp file.
        try:
            fasta_text = base64.b64decode(custom_genome_b64).decode("ascii")
        except Exception as exc:
            return {
                "job_id": job_id,
                "status": "error",
                "message": f"Failed to decode custom genome payload: {exc}",
                "records": [],
            }
        try:
            from beguided.genome.io import load_genome as _load_genome

            genome = _load_genome_from_string(fasta_text, _load_genome)
        except Exception as exc:
            return {
                "job_id": job_id,
                "status": "error",
                "message": f"Failed to parse custom genome FASTA: {exc}",
                "records": [],
            }
        # Chromosome validation: every chromosome in the targets DF must
        # appear in the custom genome.
        unknown = validate_chromosomes(
            targets_df["chromosome"].astype(str).tolist(), genome.keys()
        )
        if unknown:
            return {
                "job_id": job_id,
                "status": "error",
                "message": (
                    f"Unknown chromosome(s) in sequence input: "
                    f"{', '.join(sorted(set(unknown)))}."
                ),
                "records": [],
            }
    elif (
        genome_option is not None
        and genome_option.genome_fasta
        and Path(genome_option.genome_fasta).exists()
    ):
        try:
            if cfg.preload_genome:
                # Block the worker thread until the background preload
                # completes (or fails). The single-worker scheduler means
                # only the very first job waits; subsequent jobs see the
                # already-cached genome.
                from beguided.web.genome_cache import wait_until_ready

                genome = wait_until_ready(genome_option.key, timeout=None)
            else:
                from beguided.genome.io import load_genome

                genome = load_genome(genome_option.genome_fasta)
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
    # Sequence mode skips this entirely — there's no annotation for custom
    # sequences.
    interval_trees: Optional[dict] = None
    if (
        not sequence_mode
        and genome_option is not None
        and genome_option.gtf_path
        and Path(genome_option.gtf_path).exists()
    ):
        try:
            if cfg.preload_gtf:
                from beguided.web.gtf_cache import wait_until_ready

                interval_trees = wait_until_ready(genome_option.key, timeout=None)
            else:
                from beguided.bystander.annotation import load_gtf_to_interval_tree

                interval_trees = load_gtf_to_interval_tree(
                    genome_option.gtf_path, select_tag=genome_option.gtf_tag
                )
        except Exception as exc:
            return {
                "job_id": job_id,
                "status": "error",
                "message": f"Failed to load GTF annotation: {exc}",
                "records": [],
            }

    with tempfile.TemporaryDirectory(prefix="beguided_") as work_dir:
        if sequence_mode:
            # Persist the custom FASTA in the work dir so run_pipeline can
            # load it through its normal path. We retain the in-memory
            # ``genome`` dict and pass it via the ``genome=`` kwarg so
            # guidefinder doesn't reparse the file.
            custom_fasta_path = str(Path(work_dir) / "custom_genome.fa")
            Path(custom_fasta_path).write_text(
                base64.b64decode(custom_genome_b64).decode("ascii")
            )

        config = _build_config(
            targets_df,
            base_editor_names,
            filter_mode,
            work_dir,
            search_offtargets,
            guide_score_threshold,
            passthrough_extras,
            genome_fasta=(
                custom_fasta_path
                if sequence_mode
                else (genome_option.genome_fasta if genome_option is not None else None)
            ),
            gtf_path=genome_option.gtf_path if genome_option is not None else None,
            skip_bystander=skip_bystander,
            use_custom_profiles=use_custom_profiles,
            custom_profiles_b64=custom_profiles_b64,
            merged_profiles_b64=merged_profiles_b64,
            keep_targets_without_guides=keep_targets_without_guides,
            genome_key=genome_key,
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
            records = _df_to_records(
                final_output,
                drop_columns=SEQUENCE_MODE_DROP_COLUMNS if sequence_mode else (),
            )
            # Sequence mode keeps the spreadsheet record shape (chromosome /
            # position / target_base / desired_base) minus ``genome_base``:
            # the off-target screening row keys and the BEguidedOFF step-4
            # input need those fields. The sequence-specific display (the
            # "sequence" first column, the merged desired-edit column, the
            # hidden sequence-relative coordinates) is applied client-side.
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
                # Metadata a follow-up off-target screening job needs to
                # reconstruct this job's inputs (profiles file, genome
                # selection). Rides along in result_json; never exposed by
                # the /results endpoint, which returns only ``records``.
                "_screen_meta": {
                    "mode": mode,
                    "genome_key": genome_key,
                    "genome_fasta": (
                        genome_option.genome_fasta if genome_option is not None else ""
                    ),
                    "gtf_path": (
                        genome_option.gtf_path if genome_option is not None else ""
                    ),
                    "base_editor_names": list(base_editor_names),
                    "use_custom_profiles": bool(use_custom_profiles),
                    "custom_profiles_b64": custom_profiles_b64,
                    "merged_profiles_b64": merged_profiles_b64,
                },
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


# Placeholder-row markers mirroring the pipeline's no-guide rows (and the
# client-side PLACEHOLDER_MARKERS in app.js). Rows like these carry no guide
# sequence and can never be screened.
_SCREEN_PLACEHOLDER_MARKERS = (
    "**No guide possible**",
    "**Base change incompatible**",
    "**Guides removed by filter**",
)


def _screen_row_key(record: dict) -> tuple:
    """Identity key of one result row for screening requests.

    Matches the key the frontend builds per table row: chromosome, position,
    target_base, desired_base and base_editor uniquely identify a guide row
    of a completed job.
    """

    def _norm(value) -> str:
        if value is None:
            return ""
        if isinstance(value, float) and float(value).is_integer():
            return str(int(value))
        return str(value)

    return (
        _norm(record.get("chromosome")),
        _norm(record.get("position")),
        _norm(record.get("target_base")),
        _norm(record.get("desired_base")),
        _norm(record.get("base_editor")),
    )


def execute_screen_job(
    job_id: str,
    parent_job_id: str,
    rows: list,
    parent_records: list,
    screen_meta: dict,
    genome_key_override: Optional[str] = None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """Run BEguidedOFF off-target screening on a filtered subset of a job.

    Called from the scheduler worker thread for jobs whose payload has
    ``kind="screen"``. The parent job's pipeline results are *immutable*;
    this job writes its own artifact set (step 4 + merged final TSV) with
    the removal rule disabled — every screened row is kept and the client
    renders the off-target verdict from the operator-configured thresholds.

    Sequence-mode parents are supported with an explicit acknowledgment
    from the client (validated in the /offtarget endpoint): the on-target
    site cannot be excluded because the input sequence's location in the
    searched reference genome is unknown, so a perfect-match hit may be
    the on-target itself.

    Args:
        job_id: This (child) job's id.
        parent_job_id: The completed pipeline job the rows belong to.
        rows: List of row-key dicts (chromosome, position, target_base,
            desired_base, base_editor) as selected by the user.
        parent_records: The parent job's result records (fetched by the
            worker from the scheduler's result store).
        screen_meta: The parent job's ``_screen_meta`` (profiles, genome key).
        genome_key_override: Reference genome to search instead of the
            parent's ``_screen_meta.genome_key``. Sequence-mode jobs let
            the user pick it next to the search button; when set, the
            genome FASTA/GTF paths are re-resolved for the chosen genome.
        cancel_event: Cooperative cancel event.

    Returns:
        dict with keys: job_id, status, message, records.
    """
    from beguided.web.sse_logger import get_sse_handler

    sse = get_sse_handler()
    cfg = get_config()
    # Work on a copy: the parent result's ``_screen_meta`` is shared state
    # and must not be mutated when an override re-resolves the paths.
    meta = dict(screen_meta or {})
    genome_key = meta.get("genome_key") or "human"
    if genome_key_override:
        genome_key = genome_key_override
        option = cfg.genome_option(genome_key)
        meta["genome_fasta"] = option.genome_fasta if option is not None else ""
        meta["gtf_path"] = option.gtf_path if option is not None else ""

    if not parent_records:
        return {
            "job_id": job_id,
            "status": "error",
            "message": (
                "Parent job has no stored results to screen; re-run the "
                "guide search first."
            ),
            "records": [],
        }
    if not rows:
        return {
            "job_id": job_id,
            "status": "error",
            "message": "No guides selected for off-target screening.",
            "records": [],
        }

    try:
        _validate_offtarget_config(True, genome_key=genome_key)
    except RuntimeError as exc:
        return {"job_id": job_id, "status": "error", "message": str(exc), "records": []}

    # Match the requested keys against the parent records. Placeholder rows
    # (no guide) can never be screened and are ignored even if requested.
    keyset = {_screen_row_key(r) for r in rows}
    placeholders = set(_SCREEN_PLACEHOLDER_MARKERS)
    matched = [
        r for r in parent_records
        if _screen_row_key(r) in keyset
        and str(r.get("guide_rna_sequence") or "") not in placeholders
    ]
    matched_keys = {_screen_row_key(r) for r in matched}
    missing = keyset - matched_keys
    if missing:
        return {
            "job_id": job_id,
            "status": "error",
            "message": (
                f"{len(missing)} of {len(keyset)} requested guides were not "
                "found in the parent job's results; refresh the results and "
                "try again."
            ),
            "records": [],
        }

    resolved = cfg.beguidedoff_resolved(genome_key)

    with tempfile.TemporaryDirectory(prefix="beguided_screen_") as work_dir:
        # The filtered rows are written as the step-3 input the pipeline's
        # off-target steps expect; column order matches the parent's final
        # TSV, which is a superset of what steps 4-5 need.
        step3_path = Path(work_dir) / "step3_filtered.tsv"
        pd.DataFrame(matched).to_csv(step3_path, sep="\t", index=False)

        # Reconstruct the profiles file the parent job used, so BEguidedOFF
        # resolves every row's editor (PAM / editing window) identically.
        if meta.get("merged_profiles_b64"):
            profiles_path = Path(work_dir) / "merged_profiles.tsv"
            profiles_path.write_bytes(base64.b64decode(meta["merged_profiles_b64"]))
            profiles_file = str(profiles_path)
        elif meta.get("use_custom_profiles") and meta.get("custom_profiles_b64"):
            profiles_path = Path(work_dir) / "custom_profiles.tsv"
            profiles_path.write_bytes(base64.b64decode(meta["custom_profiles_b64"]))
            profiles_file = str(profiles_path)
        else:
            profiles_file = cfg.base_editor_profiles_file

        config = {
            "input": {
                "targets": str(step3_path),
                "genome": meta.get("genome_fasta", ""),
                "base_editor_profiles": profiles_file,
                "gff_file": meta.get("gtf_path", ""),
            },
            "output_prefix": work_dir,
            "use_base_editors": "all",
            "keep_bystanders": "all",
            "search_offtargets": True,
            "keep_targets_without_guides": False,
            "offtargets": {
                "mode": resolved["mode"],
                "beguidedoff_genome": resolved["genome"],
                "index_path": resolved["index"],
                "max_mismatches": cfg.mismatch_tolerance,
                "count_threshold_score": cfg.count_threshold_score,
                # Removal disabled: the web UI keeps every screened row and
                # shows the verdict (None is the pipeline's "rule disabled"
                # sentinel; the fast search cutoff then falls back to the
                # count threshold, clamped to the sound maximum).
                "remove_guides_score": None,
                "remove_guides_n_offtargets": None,
            },
        }

        # Reuse the preloaded FM-index handle when one is cached for this
        # path (startup preload or an earlier screening job).
        index_handle = None
        if resolved["mode"] == "fast":
            from beguided.web.offtarget_index_cache import get_handle

            index_handle = get_handle(resolved["index"])

        final_output = Path(work_dir) / "beguided_results.tsv"
        try:
            run_pipeline(
                config,
                start_step=4,
                end_step=5,
                cancel_event=cancel_event,
                offtarget_index_handle=index_handle,
            )
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
            logger.exception("Screening job %s failed: %s", job_id, exc)
            return {
                "job_id": job_id,
                "status": "error",
                "message": str(exc),
                "records": [],
            }

        records = _df_to_records(final_output)
        if not records:
            return {
                "job_id": job_id,
                "status": "error",
                "message": "Screening finished but no output rows were produced.",
                "records": [],
            }
        n_hit = sum(
            1 for r in records if r.get("worst_offtarget_score") not in (None, "")
        )
        return {
            "job_id": job_id,
            "status": "complete",
            "message": (
                f"Off-target screening complete: {len(records)} guides "
                f"screened, {n_hit} with reported off-targets."
            ),
            "records": records,
        }


def _load_genome_from_string(fasta_text: str, load_genome_fn) -> dict:
    """Load a FASTA text blob into a ``{chrom: sequence}`` dict.

    BioPython's ``SeqIO.parse`` consumes a file path; we write to a
    NamedTemporaryFile so we can reuse the existing ``load_genome`` loader.
    The temp file is closed immediately (and thus OS-cleaned) once the
    loader has read it.
    """
    import tempfile as _tf

    tmp = _tf.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False, encoding="ascii"
    )
    try:
        tmp.write(fasta_text)
        tmp.flush()
        tmp.close()
        return load_genome_fn(tmp.name)
    finally:
        try:
            Path(tmp.name).unlink()
        except OSError:
            pass
