"""Pipeline orchestration: direct function calls, no subprocess, no state files."""

import logging
import threading
from pathlib import Path

import pandas as pd

from beguided.genome.coordinates import guide_coordinates
from beguided.genome.io import load_genome
from beguided.editing.context import EditContextTable, load_edit_context_table
from beguided.editing.presets import (
    load_base_editor_profiles,
    get_available_base_editor_names,
    get_base_editor_names_by_type,
    filter_compatible_targets,
)
from beguided.guidefinder.core import run_genome_wide as run_guidefinder, results_to_dataframe
from beguided.bystander.annotation import load_gtf_to_interval_tree
from beguided.bystander.analysis import run_all as run_bystander
from beguided.filtering.filters import load_filter_criteria, filter_guides
from beguided.off_target.merger import aggregate_beguidedoff_results, merge_with_input
from beguided.logging_config import get_logger
from beguided.models import (
    BASE_CHANGE_INCOMPATIBLE,
    CORE_COLUMN_ORDER,
    NO_GUIDE_FOUND,
)

logger = get_logger("pipeline")

# Largest CCTOP score cutoff the fast (FM-index) mode can search soundly.
# The fast engine anchors on the exact "cold" region — every position whose
# single-mismatch cost (1.2 ** position, 5'-proximal = position 1) reaches the
# cutoff. A cutoff of 1.2**9 (~5.16) keeps the 12 bp PAM-proximal anchor
# (positions 9-20), which occurs a manageable number of times genome-wide.
# Beyond it the anchor shrinks (9 bp at 7.4, 6 bp at 15), hit counts blow past
# the engine's per-anchor cap, and results would be silently incomplete.
FAST_MAX_SEARCH_CUTOFF = 1.2 ** 9


class PipelineError(Exception):
    """Custom exception for pipeline errors."""
    pass


class PipelineCancelled(PipelineError):
    """Raised when a cooperative cancel_event trips mid-run."""
    pass


def _check_cancel(cancel_event):
    """Raise PipelineCancelled if the event is set."""
    if cancel_event is not None and cancel_event.is_set():
        raise PipelineCancelled("Run cancelled by user.")


def _apply_final_column_order(df):
    """Reorder DataFrame columns so CORE_COLUMN_ORDER columns come first.

    Columns in ``CORE_COLUMN_ORDER`` that are missing are silently skipped
    (e.g. ``bystander_effect`` when the bystander step was skipped). All other
    columns keep their existing relative order and follow the core columns.
    User-provided ``extra_columns`` from the input TSV therefore always land
    at the very end.
    """
    seen = set()
    ordered = []
    for c in CORE_COLUMN_ORDER:
        if c in df.columns and c not in seen:
            ordered.append(c)
            seen.add(c)
    tail = [c for c in df.columns if c not in seen]
    return df[ordered + tail]


def _resolve_base_editor_names(profiles, base_editor_arg):
    """Resolve base editor argument to list of names.

    Supports:
      - "all"        → every available base editor profile
      - "all ABEs"   → all ABE (A-to-G) profiles
      - "all CBEs"   → all CBE (C-to-T) profiles
      - comma-separated list of individual base editor names
    """
    available = get_available_base_editor_names(profiles)
    arg = base_editor_arg.lower()
    if arg == "all":
        return sorted(available)
    if arg == "all abes":
        return sorted(get_base_editor_names_by_type(profiles)["abe"])
    if arg == "all cbes":
        return sorted(get_base_editor_names_by_type(profiles)["cbe"])
    names = [p.strip() for p in base_editor_arg.split(",")]
    for name in names:
        if name not in available:
            raise PipelineError(
                f"Unknown base editor '{name}'. Available: {', '.join(sorted(available))}"
            )
    return names


def _parse_output_prefix(output_prefix):
    """Parse output_prefix into base directory.
    """
    if not output_prefix:
        return Path(".")
    return Path(output_prefix)


def _read_new_style_config(config):
    """Normalise config into internal dict.

    Only the flat format is supported.
    """
    result = {}

    if "output_prefix" not in config:
        raise PipelineError("Missing required key: 'output_prefix'")

    result["base_dir"] = _parse_output_prefix(config["output_prefix"])

    # --- inputs ---
    input_cfg = config.get("input", {})
    if "targets" not in input_cfg:
        raise PipelineError("Missing required key: input.targets")
    if "genome" not in input_cfg:
        raise PipelineError("Missing required key: input.genome")

    result["targets"] = input_cfg.get("targets")
    result["genome"] = input_cfg.get("genome")
    result["base_editor_profiles_file"] = input_cfg.get(
        "base_editor_profiles", "configs/base_editor_profiles.tsv"
    )
    result["gff_file"] = input_cfg.get("gff_file")

    # --- use_base_editors (which base editor profiles to apply) ---
    result["base_editor_arg"] = config.get("use_base_editors", "all")

    # --- keep_bystanders ---
    result["filter_mode"] = config.get("keep_bystanders", "all")

    # --- guide_score_threshold ---
    guide_score_threshold = config.get("guide_score_threshold", 0.5)
    try:
        guide_score_threshold = float(guide_score_threshold)
    except (TypeError, ValueError) as exc:
        raise PipelineError(
            f"guide_score_threshold must be a number, got {guide_score_threshold!r}"
        ) from exc
    if not (0 <= guide_score_threshold <= 1):
        raise PipelineError(
            f"guide_score_threshold must be between 0 and 1, got {guide_score_threshold}"
        )
    result["guide_score_threshold"] = guide_score_threshold

    # --- gc_content_filter ---
    gc_content_filter = config.get("gc_content_filter", "off")
    if gc_content_filter not in ("off", "good", "acceptable"):
        raise PipelineError(
            f"gc_content_filter must be one of off|good|acceptable, "
            f"got {gc_content_filter!r}"
        )
    result["gc_content_filter"] = gc_content_filter

    # --- filter_homopolymers ---
    result["filter_homopolymers"] = bool(config.get("filter_homopolymers", False))

    # --- trinuc_edit_score_min ---
    # Optional minimum trinucleotide context score (Marquart et al. 2021,
    # normalised ABE8e) below which ABE guides are dropped in step 3.
    # ``None`` disables the filter. CBE rows are always kept (their
    # ``trinuc_edit_score`` is NaN by construction).
    raw_trinuc_min = config.get("trinuc_edit_score_min", None)
    if raw_trinuc_min is None or raw_trinuc_min == "":
        result["trinuc_edit_score_min"] = None
    else:
        try:
            trinuc_min = float(raw_trinuc_min)
        except (TypeError, ValueError) as exc:
            raise PipelineError(
                f"trinuc_edit_score_min must be a number, got {raw_trinuc_min!r}"
            ) from exc
        if not (0.0 <= trinuc_min <= 1.0):
            raise PipelineError(
                f"trinuc_edit_score_min must be between 0.0 and 1.0, "
                f"got {trinuc_min}"
            )
        result["trinuc_edit_score_min"] = trinuc_min

    # --- edit_context_table ---
    # Optional path to a custom trinucleotide context TSV overriding the
    # bundled copy. ``None`` (default) uses the bundled table. An empty
    # string disables the lookup entirely (every guide's
    # ``trinuc_edit_score`` will be NaN).
    raw_table_path = config.get("edit_context_table", None)
    if raw_table_path == "":
        raw_table_path = ""  # explicit opt-out
    result["edit_context_table"] = raw_table_path

    # --- search_offtargets ---
    result["search_offtargets"] = config.get("search_offtargets", True)

    # --- offtargets block ---
    offtargets = config.get("offtargets", {})
    result["beguidedoff_genome"] = offtargets.get("beguidedoff_genome")

    # --- offtargets.mode ---
    # "extensive" (GPU brute force, default) or "fast" (CPU FM-index filter,
    # sound at the removal score cutoff, no GPU required). The web server
    # always sends an explicit value; plain CLI configs default to the
    # historical extensive behaviour.
    mode = str(offtargets.get("mode", "extensive")).strip().lower()
    if mode not in ("extensive", "fast"):
        raise PipelineError(
            f"offtargets.mode must be 'extensive' or 'fast', "
            f"got {offtargets.get('mode')!r}"
        )
    result["beguidedoff_mode"] = mode
    # FM-index path for fast mode. Empty means "derive from the .npz path".
    result["beguidedoff_index"] = offtargets.get("index_path") or ""

    # --- offtargets.max_mismatches ---
    # Mismatch budget of the extensive (GPU) search; ignored by fast mode,
    # which is driven by the score cutoff instead. The web UI historically
    # sent this value under the name ``mismatch_tolerance`` (where it was
    # silently ignored); accept that key as a fallback alias.
    result["max_mismatches"] = offtargets.get(
        "max_mismatches", offtargets.get("mismatch_tolerance", 5)
    )
    result["count_threshold_score"] = offtargets.get("count_threshold_score", 15.0)
    result["remove_guides_score"] = offtargets.get("remove_guides_score", 5.0)
    result["remove_guides_n_offtargets"] = offtargets.get("remove_guides_n_offtargets", 100)

    # --- keep_intermediates ---
    result["keep_intermediates"] = config.get("keep_intermediates", True)

    # --- passthrough_extras ---
    result["passthrough_extras"] = config.get("passthrough_extras", True)

    # --- skip_bystander ---
    # Explicit flag to skip step 2 (bystander analysis). This is the sole
    # gate for whether step 2 runs; ``keep_bystanders`` only governs step 3
    # filtering. Used by the sequence-input web path where no GTF annotation
    # is available — the step is skipped, and the downstream filter must NOT
    # apply a "none" cut on ``bystander_positions`` either (which
    # would drop every guide on A-rich sequences).
    result["skip_bystander"] = bool(config.get("skip_bystander", False))

    # --- keep_targets_without_guides ---
    # When True, the final results table keeps one row per input target even
    # if no guide could be designed (or every guide was filtered out). The
    # placeholder rows carry the input identity columns plus a
    # ``guide_rna_sequence`` marker (``**No guide possible**`` or
    # ``**Base change incompatible**``); all other result fields are empty.
    # Placement rows are added ONLY in the very last step so they never
    # interfere with guidefinder, bystander analysis, filtering or off-target
    # search. Off by default (CLI); the web server always enables it.
    result["keep_targets_without_guides"] = bool(
        config.get("keep_targets_without_guides", False)
    )

    return result


def _fill_targets_without_guides(df, cfg):
    """Append one placeholder row per input target missing from ``df``.

    Called in the very last pipeline step when ``keep_targets_without_guides``
    is enabled. Targets are re-loaded from the original input file and split
    by base-editor compatibility:

    * targets whose mutation type is not supported by the selected editors
      get ``BASE_CHANGE_INCOMPATIBLE`` as their ``guide_rna_sequence``;
    * targets that were compatible but produced no surviving guide (no valid
      PAM, or every guide was filtered out) get ``NO_GUIDE_FOUND``.

    Every other result column is left empty/``None`` so downstream consumers
    (incl. the web UI) can distinguish real guides from placeholder rows.
    Any extra columns the original target carried are copied through when
    they also exist in ``df`` (``passthrough_extras``-style metadata).

    Args:
        df: the final, already-column-ordered results DataFrame.
        cfg: the normalised pipeline config dict (needs ``targets``,
            ``base_editor_profiles_file`` and ``base_editor_arg``).

    Returns:
        The input ``df`` with placeholder rows appended (converted to the
        same column set). Row order is finalised by the caller, which
        re-orders the table to match the input file order.
    """
    if df.empty and "chromosome" not in df.columns:
        return df

    targets_df = pd.read_csv(cfg["targets"], sep="\t")
    # Re-run the same step-1 compatibility split so placeholders can tell
    # "mutation type not supported" apart from "no guide possible".
    profiles = load_base_editor_profiles(cfg["base_editor_profiles_file"])
    base_editor_names = _resolve_base_editor_names(profiles, cfg["base_editor_arg"])
    compatible, incompatible = filter_compatible_targets(
        targets_df, profiles, base_editor_names
    )

    # Keys already present in the final table (real guides or earlier
    # placeholder rows).
    present = set()
    if "chromosome" in df.columns and "position" in df.columns:
        present = set(
            zip(
                df["chromosome"].astype(str),
                df["position"].astype(int),
            )
        )

    placeholder_rows = []
    for _, row in incompatible.iterrows():
        placeholder_rows.append((row, BASE_CHANGE_INCOMPATIBLE))
    for _, row in compatible.iterrows():
        key = (str(row["chromosome"]), int(row["position"]))
        if key not in present:
            placeholder_rows.append((row, NO_GUIDE_FOUND))

    if not placeholder_rows:
        return df

    def _placeholder(input_row, marker):
        record = {}
        for col in df.columns:
            if col in ("chromosome", "position", "target_base", "desired_base"):
                record[col] = input_row.get(col)
            elif col in input_row.index:
                record[col] = input_row[col]
            else:
                record[col] = None
            if col == "guide_rna_sequence":
                record[col] = marker
        return record

    new_rows = pd.DataFrame(
        [_placeholder(row, marker) for row, marker in placeholder_rows],
        columns=df.columns,
    )
    merged_df = pd.concat([df, new_rows], ignore_index=True)
    return merged_df


def _reorder_by_input_order(df, targets_path):
    """Stable-sort *df* rows to follow the operator's input target order.

    The pipeline computes results in lexicographic ``chromosome``/``position``
    order (so ``chr10`` sorts before ``chr2``); this re-orders the final table
    so targets appear in the order they have in the input file. All rows of
    one target stay grouped, keeping their current relative order (editors
    alphabetically, guides by editing position). Targets missing from the
    input file (should not happen) keep their relative order at the end.

    Returns *df* unchanged when the position columns are missing or the
    targets file cannot be read/parsed.
    """
    if df.empty or not {"chromosome", "position"} <= set(df.columns):
        return df
    try:
        targets_df = pd.read_csv(targets_path, sep="\t")
        order: dict = {}
        for _, row in targets_df.iterrows():
            key = (str(row["chromosome"]), int(row["position"]))
            if key not in order:
                order[key] = len(order)
        keys = list(zip(df["chromosome"].astype(str), df["position"].astype(int)))
    except (ValueError, TypeError, KeyError, OSError) as exc:
        logger.warning(
            "Could not derive input target order from %s (%s); "
            "keeping current row order", targets_path, exc,
        )
        return df
    if not order:
        return df
    fallback = len(order)
    ordered = df.assign(_input_order=[order.get(k, fallback) for k in keys])
    ordered = ordered.sort_values("_input_order", kind="stable")
    return ordered.drop(columns="_input_order").reset_index(drop=True)


def _prepare_beguidedoff_input(df, output_path):
    """Prepare a TSV for BEguidedOFF with correct on-target positions.

    BEguidedOFF expects ``position`` to be the 1-based leftmost guide
    position for ``+`` strand guides and the 1-based rightmost guide
    position for ``-`` strand guides (used for on-target exclusion).
    We compute it from the target base position, editing position
    within the guide, and strand.
    """
    valid = df[df["guide_rna_sequence"] != "**No guide possible**"].copy()
    if valid.empty:
        valid.to_csv(output_path, sep="\t", index=False)
        return

    positions = []
    for _, row in valid.iterrows():
        gs, ge = guide_coordinates(
            row["position"],
            row["editing_position"],
            row["strand_used"],
            guide_length=20,
        )
        if row["strand_used"] == "+":
            positions.append(gs + 1)  # leftmost 1-based
        else:
            positions.append(ge)      # rightmost 1-based

    out = pd.DataFrame({
        "chromosome": valid["chromosome"].astype(str),
        "position": positions,
        "guide_rna_sequence": valid["guide_rna_sequence"],
        "base_editor": valid["base_editor"],
        "strand_used": valid["strand_used"],
    })
    out.to_csv(output_path, sep="\t", index=False)


def _run_beguidedoff(filter_file, presets_file, genome_npz, max_mismatches,
                       count_threshold_score, output_file, batch_size=1500,
                       cancel_event=None, mode="extensive", index_path="",
                       search_cutoff=None, index_handle=None):
    """Run BEguidedOFF off-target search on the filtered guide file.

    Guides are loaded explicitly and processed in batches so the pipeline
    stays within BEguidedOFF's internal limits regardless of input size.

    Args:
        mode: ``"extensive"`` (GPU brute force; searches every site with at
            most ``max_mismatches`` mismatches) or ``"fast"`` (CPU FM-index;
            enumerates exactly the sites with CCTOP score <=
            ``search_cutoff`` — same remove/keep verdicts, no GPU needed).
        index_path: FM-index path (fast mode only).
        search_cutoff: CCTOP score cutoff (fast mode only). Must be at least
            ``remove_guides_score`` so the removal rule sees every site it
            could act on; the call site enforces this.
        index_handle: Optional preloaded FM-index object (fast mode only).
            When given it is used directly instead of opening ``index_path``,
            letting long-lived processes (e.g. the web server) share one warm
            handle across jobs.
        cancel_event: Optional ``threading.Event``. Checked between batches;
            raising ``PipelineCancelled`` exits before the next batch starts.
            A single batch cannot be interrupted.

    Returns:
        The BEguidedOFF result DataFrame.
    """
    try:
        from beguidedoff import beguidedoff
        from beguidedoff.io import load_guides_from_tsv
    except ImportError as exc:
        raise PipelineError(
            "BEguidedOFF is not installed. Off-target search is unavailable.\n"
            "Install it from: https://github.com/ltschmitt/BEguidedOFF\n"
            "Or set search_offtargets=false in your config to skip off-target search."
        ) from exc

    _check_cancel(cancel_event)

    if mode == "fast":
        if not index_path:
            raise PipelineError(
                "fast mode requires an FM-index path (offtargets.index_path). "
                "Build one from a reference FASTA with: "
                "python -m beguidedoff.cli build-index <genome.fa> <genome.fmi>"
            )
        finder = beguidedoff(
            mode="fast",
            index_path=index_path,
            search_cutoff=search_cutoff,
            n_top_offtargets=1,
            scoring_method="CCTOP",
            index=index_handle,
        )
        # Fast mode enumerates exactly the sites with score <= search_cutoff.
        # Counting at a *higher* threshold would silently widen the search to
        # a shorter exact anchor whose genome-wide hit count routinely
        # exceeds the engine's per-anchor cap (incomplete counts); clamp the
        # count threshold to the search cutoff instead.
        effective_count_threshold = count_threshold_score
        if (effective_count_threshold is not None
                and float(effective_count_threshold) > float(search_cutoff)):
            logger.warning(
                "count_threshold_score=%s exceeds the fast-mode search cutoff "
                "%s; counting off-targets at %s instead.",
                effective_count_threshold, search_cutoff, search_cutoff,
            )
            effective_count_threshold = float(search_cutoff)
    elif mode == "extensive":
        finder = beguidedoff(
            genome_path=genome_npz,
            max_mismatches=max_mismatches,
            n_top_offtargets=1,
            scoring_method="CCTOP",
        )
        effective_count_threshold = count_threshold_score
    else:
        raise PipelineError(f"Unknown BEguidedOFF mode: {mode!r}")

    # Load guides explicitly so we can batch them.
    guides = load_guides_from_tsv(str(filter_file), presets_path=str(presets_file))

    total_guides = len(guides)
    n_batches = (total_guides + batch_size - 1) // batch_size
    dfs = []
    for i in range(0, total_guides, batch_size):
        _check_cancel(cancel_event)
        batch_num = i // batch_size + 1
        batch = guides[i:i + batch_size]
        logger.info(
            "BEguidedOFF batch %d/%d (%d guides)", batch_num, n_batches, len(batch)
        )
        ot_df = finder.search(
            guides=batch,
            return_all=False,
            count_score_threshold=effective_count_threshold,
            verbose=False,
            logger=logger,
        )
        if ot_df is not None and not ot_df.empty:
            dfs.append(ot_df)

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
    else:
        combined = pd.DataFrame()

    if output_file:
        combined.to_csv(str(output_file), sep="\t", index=False)

    return combined


def _resolve_fast_search_cutoff(cfg):
    """Resolve the FM-index search cutoff for fast mode.

    The removal rule drops guides whose worst off-target score is below
    ``remove_guides_score``; the search must therefore cover every site below
    that threshold, i.e. the cutoff must be at least the removal threshold.
    When the removal rule is disabled, fall back to the count threshold
    (clamped to the largest sound cutoff) so counting stays complete.

    Raises:
        PipelineError: If the resolved cutoff exceeds
            :data:`FAST_MAX_SEARCH_CUTOFF` (a shorter exact anchor would
            routinely overflow the engine's per-anchor hit cap and make
            results silently incomplete).
    """
    remove_thr = cfg.get("remove_guides_score")
    count_thr = cfg.get("count_threshold_score")
    if remove_thr is not None:
        cutoff = float(remove_thr)
    else:
        cutoff = float(count_thr) if count_thr is not None else 5.0
        cutoff = min(cutoff, FAST_MAX_SEARCH_CUTOFF)
    if cutoff < 0:
        raise PipelineError(
            "fast mode requires a non-negative search cutoff; "
            f"the resolved cutoff is {cutoff}. Check remove_guides_score / "
            "count_threshold_score in the offtargets config."
        )
    if cutoff > FAST_MAX_SEARCH_CUTOFF:
        raise PipelineError(
            f"fast mode supports search cutoffs up to "
            f"{FAST_MAX_SEARCH_CUTOFF:.2f} (12 bp exact anchor); "
            f"the resolved cutoff {cutoff:.2f} would need a shorter anchor "
            f"whose genome-wide hit count routinely exceeds the engine's "
            f"per-anchor cap and makes results incomplete. Lower "
            f"remove_guides_score or use mode='extensive' (GPU)."
        )
    return cutoff


def run_pipeline(config, start_step=1, end_step=5, cancel_event=None, genome=None,
                interval_trees=None, offtarget_index_handle=None):
    """Run pipeline steps as direct function calls.

    Output layout::

        <base_dir>/
            step1_guidefinder.tsv
            step2_bystander.tsv
            step3_filtered.tsv
            step4_offtargets.tsv     (if search_offtargets)
            beguided_results.tsv

    Args:
        config: dict loaded from YAML.
        start_step: First step to run (1-5).
        end_step: Last step to run (1-5).
        cancel_event: Optional ``threading.Event``. Checked at the entry of
            each step and between BEguidedOFF batches. Setting it raises
            ``PipelineCancelled`` at the next checkpoint.
        genome: Optional pre-loaded genome dict (``{chrom_id: sequence}``).
            When supplied, Step 1 and Step 2 reuse this dict instead of calling
            ``load_genome()`` on the configured FASTA. The caller is responsible
            for loading it (e.g. ``beguided.web.genome_cache``).
        interval_trees: Optional pre-loaded GTF interval trees
            (``{chrom_id: IntervalTree}``). When supplied, Step 2 reuses this
            dict instead of calling ``load_gtf_to_interval_tree()`` on the
            configured GTF. The caller is responsible for loading it (e.g.
            ``beguided.web.gtf_cache``).
        offtarget_index_handle: Optional pre-loaded FM-index object for
            fast-mode Step 4 (e.g. from ``beguided.web.offtarget_index_cache``).
            When supplied it is used instead of opening
            ``offtargets.index_path``; the path is still required for logging
            and existence checks.

    Returns:
        Path to final output file.
    """
    import shutil

    cfg = _read_new_style_config(config)

    base_dir = cfg["base_dir"]
    base_dir.mkdir(parents=True, exist_ok=True)

    guidefinder_output = base_dir / "step1_guidefinder.tsv"
    bystander_output = base_dir / "step2_bystander.tsv"
    filter_file = base_dir / "step3_filtered.tsv"
    offtarget_output = base_dir / "step4_offtarget.tsv"
    final_output = base_dir / "beguided_results.tsv"
    screen_offtargets = cfg["search_offtargets"]

    # Step 1: Guidefinder
    _check_cancel(cancel_event)
    if start_step <= 1 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 1/5: Guidefinder - Generating guide RNAs")
        logger.info("=" * 60)
        presets_file = cfg["base_editor_profiles_file"]
        logger.info(f"Loading base editor profiles from: {presets_file}")
        profiles = load_base_editor_profiles(presets_file)

        base_editor_names = _resolve_base_editor_names(profiles, cfg["base_editor_arg"])
        logger.info(f"Using base editors: {', '.join(base_editor_names)}")

        # Load the trinucleotide edit-context table (Marquart et al. 2021,
        # ABE8e, normalised). Empty string => explicit opt-out; None =>
        # bundled default. CBE rows leave ``trinuc_edit_score`` at NaN.
        edit_context_table: EditContextTable | None
        table_path = cfg["edit_context_table"]
        if table_path == "":
            edit_context_table = None
            logger.info("Edit-context table disabled (trinuc_edit_score will be NaN)")
        else:
            try:
                edit_context_table = load_edit_context_table(
                    table_path if table_path is not None else None
                )
                logger.info(
                    "Loaded edit-context table (%s). "
                    "Note: scores are ABE8e-derived and approximate for other ABE editors; "
                    "CBE rows will be NaN.",
                    edit_context_table.source_label,
                )
            except FileNotFoundError as exc:
                raise PipelineError(str(exc)) from exc

        genome_file = cfg["genome"]
        if genome is None:
            logger.info(f"Loading genome from: {genome_file}")
            genome = load_genome(genome_file)
            logger.info(f"Loaded {len(genome)} chromosomes")
        else:
            logger.info(f"Using pre-loaded genome ({len(genome)} chromosomes)")

        targets_df = pd.read_csv(cfg["targets"], sep="\t")
        original_target_count = len(targets_df)
        logger.info(f"Loaded {original_target_count} targets")

        targets_df, skipped_df = filter_compatible_targets(targets_df, profiles, base_editor_names)
        if len(skipped_df) > 0:
            logger.warning(
                f"{len(skipped_df)}/{original_target_count} target(s) were skipped "
                f"because their mutation type is not supported by the selected base editors. "
                f"Skipped: "
                f"{', '.join([f'{r.chromosome}:{r.position}' for r in skipped_df.itertuples()][:5])}"
                f"{' ...' if len(skipped_df) > 5 else ''}"
            )
        if len(targets_df) == 0:
            if cfg["keep_targets_without_guides"]:
                # With the flag on we want a row per target in the final
                # output even though every target was dropped for editor
                # incompatibility. Continue with an empty frame: guidefinder
                # returns no results, later steps see an empty table and the
                # final step fills one placeholder row per original target.
                logger.warning(
                    "No compatible targets remain; final output will contain "
                    "one '**Base change incompatible**' row per target."
                )
            else:
                raise PipelineError(
                    "No targets remain after filtering for base editor compatibility."
                )
        logger.info(f"Processing {len(targets_df)} compatible targets")

        results = run_guidefinder(
            genome, targets_df, profiles, base_editor_names,
            progress_callback=lambda c, p: logger.debug(f"Processing {c}:{p}"),
            passthrough_extras=cfg["passthrough_extras"],
            edit_context_table=edit_context_table,
        )
        df = results_to_dataframe(
            results, passthrough_extras=cfg["passthrough_extras"],
        )
        df.to_csv(guidefinder_output, sep="\t", index=False)

        invalid_mask = df["guide_rna_sequence"].fillna("").str.startswith("**")
        if invalid_mask.any():
            target_guide_counts = (
                df.groupby(["chromosome", "position"])
                .agg(total=("guide_rna_sequence", "size"),
                     invalid=("guide_rna_sequence", lambda s: s.str.startswith("**").all()))
                .reset_index()
            )
            no_valid_guides = target_guide_counts[target_guide_counts["invalid"]]
            total_targets = len(target_guide_counts)

            if len(no_valid_guides) > 0:
                pure_mismatch = []
                pure_pam_fail = []
                all_rows = df.groupby(["chromosome", "position"])
                for (chrom, pos), group in all_rows:
                    seqs = group["guide_rna_sequence"].fillna("").tolist()
                    if all(s.startswith("**") for s in seqs):
                        if all("doesn't match base editor profile" in s for s in seqs):
                            pure_mismatch.append(f"{chrom}:{pos}")
                        else:
                            pure_pam_fail.append(f"{chrom}:{pos}")

                if pure_mismatch:
                    logger.warning(
                        f"{len(pure_mismatch)}/{total_targets} target(s) use a mutation "
                        f"type not supported by the selected base editors: "
                        f"{', '.join(pure_mismatch[:5])}"
                        f"{' ...' if len(pure_mismatch) > 5 else ''}"
                    )
                if pure_pam_fail:
                    logger.warning(
                        f"{len(pure_pam_fail)}/{total_targets} target(s) had no valid PAM"
                    )

        total_result_rows = len(df)
        total_invalid_rows = invalid_mask.sum()
        if cfg["keep_targets_without_guides"]:
            if total_result_rows > 0 and total_invalid_rows == total_result_rows:
                logger.warning(
                    "No valid guides for any compatible target; final output "
                    "will contain one '**No guide possible**' row per target."
                )
        elif total_invalid_rows == total_result_rows:
            # Historic behaviour: raise when not a single real guide was
            # produced (kept for the CLI default so an empty result set is
            # surfaced as an error rather than a silent empty table).
            raise PipelineError(
                "No valid guides found for the given base editors."
            )

        logger.info(
            f"Step 1 complete: {len(df)} results written to {guidefinder_output}"
        )

    # Step 2: Bystander (conditional)
    _check_cancel(cancel_event)
    filter_mode = cfg["filter_mode"]
    skip_bystander = cfg["skip_bystander"]
    needs_bystander = not skip_bystander

    if start_step <= 2 <= end_step:
        logger.info("=" * 60)
        if needs_bystander:
            logger.info("Step 2/5: Bystander analysis")
        else:
            logger.info("Step 2/5: Bystander analysis (SKIPPED)")
        logger.info("=" * 60)

        if needs_bystander:
            # Load the genome if not pre-loaded — required so the
            # codon-aligned extract can splice across exon junctions.
            if genome is None:
                genome_file = cfg["genome"]
                if not genome_file:
                    raise PipelineError(
                        "Genome FASTA is required for bystander analysis "
                        "(when starting from step 2 or later, supply it via "
                        "the genome= kwarg)."
                    )
                logger.info(f"Loading genome from: {genome_file}")
                genome = load_genome(genome_file)
                logger.info(f"Loaded {len(genome)} chromosomes")
            if interval_trees is None:
                gtf_file = cfg["gff_file"]
                if not gtf_file:
                    raise PipelineError(
                        "GFF/GTF file is required for bystander analysis."
                    )
                logger.info(f"Loading GTF annotation from: {gtf_file}")
                interval_trees = load_gtf_to_interval_tree(gtf_file)
                _chrom_count = len(interval_trees) - (
                    1 if "__tx_codons__" in interval_trees else 0
                )
                logger.info(f"Loaded annotations for {_chrom_count} chromosomes")
            else:
                _chrom_count = len(interval_trees) - (
                    1 if "__tx_codons__" in interval_trees else 0
                )
                logger.info(
                    f"Using pre-loaded GTF annotation ({_chrom_count} chromosomes)"
                )

            df = pd.read_csv(guidefinder_output, sep="\t")
            df = df[df["guide_rna_sequence"] != "**No guide possible**"]
            logger.info(f"Analyzing {len(df)} guides for bystander effects")

            if len(df) > 0:
                df = run_bystander(
                    df, interval_trees,
                    progress_callback=lambda idx, total: (
                        logger.info(f"  Bystander analysis: {idx + 1}/{total} guides")
                        if idx % 100 == 0 else None
                    ),
                    genome=genome,
                )
                logger.info("Bystander analysis complete")
            df.to_csv(bystander_output, sep="\t", index=False)
        else:
            shutil.copy(str(guidefinder_output), str(bystander_output))
            logger.info(f"Copied guidefinder output to {bystander_output}")

    # Step 3: Filter
    _check_cancel(cancel_event)
    if start_step <= 3 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 3/5: Filtering guides")
        logger.info(f"Filter mode: {filter_mode}")
        logger.info(f"On-target score threshold: {cfg['guide_score_threshold']}")
        if cfg["trinuc_edit_score_min"] is not None:
            logger.info(
                f"Trinucleotide edit-context min: {cfg['trinuc_edit_score_min']} "
                f"(ABE only; CBE rows kept)"
            )
        logger.info("=" * 60)

        df = pd.read_csv(bystander_output, sep="\t")
        original_count = len(df)
        logger.info(f"Input guides: {original_count}")

        criteria = load_filter_criteria(
            mode=filter_mode,
            guide_score_threshold=cfg["guide_score_threshold"],
            gc_content_filter=cfg["gc_content_filter"],
            filter_homopolymers=cfg["filter_homopolymers"],
            trinuc_edit_score_min=cfg["trinuc_edit_score_min"],
        )
        df = filter_guides(df, criteria)
        logger.info(f"After filtering: {len(df)} guides (from {original_count})")

        df = _apply_final_column_order(df)
        df.to_csv(filter_file, sep="\t", index=False)
        logger.info(f"Filtered results written to {filter_file}")

    # Step 4: BEguidedOFF off-target search (single step)
    _check_cancel(cancel_event)
    if screen_offtargets and start_step <= 4 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 4/5: BEguidedOFF off-target search")
        logger.info("=" * 60)

        genome_npz = cfg["beguidedoff_genome"]
        mode = cfg["beguidedoff_mode"]
        if mode == "extensive" and not genome_npz:
            raise PipelineError(
                "beguidedoff_genome is required for extensive (GPU) mode. "
                "Provide the path to the encoded .npz genome file."
            )

        logger.info(f"BEguidedOFF genome: {genome_npz or '(not used in fast mode)'}")
        logger.info(f"BEguidedOFF mode: {mode}")

        index_path = ""
        search_cutoff = None
        if mode == "fast":
            index_path = cfg.get("beguidedoff_index") or ""
            if not index_path:
                if not genome_npz:
                    raise PipelineError(
                        "fast mode requires an FM-index: set "
                        "BEGUIDED_BEGUIDEDOFF_INDEX (or offtargets.index_path), "
                        "or BEGUIDED_BEGUIDEDOFF_GENOME (or "
                        "offtargets.beguidedoff_genome) to derive the index "
                        "path from the encoded .npz."
                    )
                # Default: FM-index next to the .npz, same basename.
                index_path = str(Path(str(genome_npz)).with_suffix(".fmi"))
            if not Path(index_path).exists():
                raise PipelineError(
                    f"FM-index not found: {index_path}. Build it from the "
                    f"reference FASTA with `python -m beguidedoff.cli "
                    f"build-index <genome.fa> {index_path}` or point "
                    f"offtargets.index_path at an existing index."
                )
            logger.info(f"BEguidedOFF index: {index_path}")
            search_cutoff = _resolve_fast_search_cutoff(cfg)
            logger.info(f"Search cutoff (score): {search_cutoff:.2f}")
        else:
            logger.info(f"Max mismatches: {cfg['max_mismatches']}")
        logger.info(f"Score threshold (count): {cfg['count_threshold_score']}")

        # Prepare input file with guide-start positions for on-target exclusion
        bgo_input = base_dir / "step4_offtarget_input.tsv"
        _prepare_beguidedoff_input(
            pd.read_csv(filter_file, sep="\t"),
            bgo_input,
        )
        logger.info(f"Prepared BEguidedOFF input: {bgo_input}")

        ot_df = _run_beguidedoff(
            bgo_input,
            cfg["base_editor_profiles_file"],
            genome_npz,
            cfg["max_mismatches"],
            cfg["count_threshold_score"],
            offtarget_output,
            cancel_event=cancel_event,
            mode=mode,
            index_path=index_path,
            search_cutoff=search_cutoff,
            index_handle=offtarget_index_handle,
        )
        logger.info(
            f"BEguidedOFF complete: {len(ot_df)} off-target rows written to "
            f"{offtarget_output}"
        )

    _check_cancel(cancel_event)
    if start_step <= 5 <= end_step:
        input_df = pd.read_csv(filter_file, sep="\t")

        if screen_offtargets:
            logger.info("=" * 60)
            logger.info("Step 5/5: Merging off-target summary")
            logger.info("=" * 60)

            if not offtarget_output.exists():
                raise PipelineError(
                    f"Off-target output not found: {offtarget_output}. "
                    f"Run step 4 first."
                )

            try:
                ot_df = pd.read_csv(offtarget_output, sep="\t")
            except pd.errors.EmptyDataError:
                ot_df = pd.DataFrame()
            summary_df = aggregate_beguidedoff_results(ot_df)
            logger.info(
                f"Loaded off-target summary: {len(summary_df)} unique guides"
            )

            logger.info(f"Merging with {len(input_df)} filtered guides")
            merged = merge_with_input(input_df, summary_df)
            logger.info(f"Merged: {len(merged)} total rows")

            # Remove verbose intermediate columns
            drop_cols = [
                "genomic_context_upstream",
                "genomic_context_downstream",
                "simulated_outcome",
                "aa_wanted",
                "aa_bystander",
            ]
            merged = merged.drop(columns=drop_cols, errors="ignore")

            worst_threshold = cfg["remove_guides_score"]
            max_offtargets_remove = cfg["remove_guides_n_offtargets"]

            if worst_threshold is not None:
                before = len(merged)
                # Keep guides with no off-targets (NaN) or with score >= threshold
                mask = (
                    merged["worst_offtarget_score"].isna()
                    | (merged["worst_offtarget_score"] >= worst_threshold)
                )
                merged = merged[mask]
                logger.info(
                    f"Filtered by worst_offtarget_score >= {worst_threshold}: "
                    f"{len(merged)}/{before} guides remain"
                )

            if max_offtargets_remove is not None:
                before = len(merged)
                merged = merged[merged["offtarget_count"] <= max_offtargets_remove]
                logger.info(
                    f"Filtered by max off-targets <= {max_offtargets_remove}: "
                    f"{len(merged)}/{before} guides remain"
                )

            merged = _apply_final_column_order(merged)
        else:
            logger.info("=" * 60)
            logger.info("Step 5/5: Finalising results (off-target screening disabled)")
            logger.info("=" * 60)

            logger.info(f"Using {len(input_df)} filtered guides (no off-target data)")

            drop_cols = [
                "genomic_context_upstream",
                "genomic_context_downstream",
                "simulated_outcome",
                "aa_wanted",
                "aa_bystander",
            ]
            merged = input_df.drop(columns=drop_cols, errors="ignore")

            merged = _apply_final_column_order(merged)

        # keep_targets_without_guides: add one placeholder row per input target
        # that has no row in the final table. Done here, at the very end, so it
        # cannot interfere with any earlier processing step.
        if cfg["keep_targets_without_guides"]:
            merged = _fill_targets_without_guides(merged, cfg)
            merged = _apply_final_column_order(merged)

        # Present the rows in the operator's input order (all rows of a target
        # stay grouped) instead of the lexicographic chromosome/position order
        # the pipeline computes in.
        merged = _reorder_by_input_order(merged, cfg["targets"])

        merged.to_csv(final_output, sep="\t", index=False)
        size_mb = final_output.stat().st_size / (1024 * 1024)
        logger.info(
            f"Final output written: {final_output} "
            f"({len(merged)} rows, {size_mb:.1f} MB)"
        )

    if not cfg["keep_intermediates"]:
        logger.info("keep_intermediates=false (cleanup not yet implemented)")

    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETE")
    logger.info(f"Output: {final_output}")
    logger.info("=" * 60)

    return final_output
