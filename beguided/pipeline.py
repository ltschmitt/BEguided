"""Pipeline orchestration: direct function calls, no subprocess, no state files."""

import logging
from pathlib import Path

import pandas as pd

from beguided.genome.coordinates import guide_coordinates
from beguided.genome.io import load_genome
from beguided.editing.presets import load_presets, get_available_preset_names, get_preset_names_by_type, filter_compatible_targets
from beguided.guidefinder.core import run_genome_wide as run_guidefinder, results_to_dataframe
from beguided.bystander.annotation import load_gtf_to_interval_tree
from beguided.bystander.analysis import run_all as run_bystander
from beguided.filtering.filters import load_filter_criteria, filter_guides
from beguided.off_target.merger import aggregate_beguidedoff_results, merge_with_input
from beguided.logging_config import get_logger

logger = get_logger("pipeline")


class PipelineError(Exception):
    """Custom exception for pipeline errors."""
    pass


def _resolve_preset_names(presets, preset_arg):
    """Resolve preset argument to list of names.

    Supports:
      - "all"        → every available preset
      - "all ABEs"   → all A-to-G presets
      - "all CBEs"   → all C-to-T presets
      - comma-separated list of individual preset names
    """
    available = get_available_preset_names(presets)
    arg = preset_arg.lower()
    if arg == "all":
        return sorted(available)
    if arg == "all abes":
        return sorted(get_preset_names_by_type(presets)["abe"])
    if arg == "all cbes":
        return sorted(get_preset_names_by_type(presets)["cbe"])
    names = [p.strip() for p in preset_arg.split(",")]
    for name in names:
        if name not in available:
            raise PipelineError(
                f"Unknown preset '{name}'. Available: {', '.join(sorted(available))}"
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
    result["presets_file"] = input_cfg.get("presets", "be_presets.tsv")
    result["gff_file"] = input_cfg.get("gff_file")

    # --- use_presets ---
    result["preset_arg"] = config.get("use_presets", "all")

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

    # --- search_offtargets ---
    result["search_offtargets"] = config.get("search_offtargets", True)

    # --- offtargets block ---
    offtargets = config.get("offtargets", {})
    result["beguidedoff_genome"] = offtargets.get("beguidedoff_genome")
    result["max_mismatches"] = offtargets.get("max_mismatches", 5)
    result["count_threshold_score"] = offtargets.get("count_threshold_score", 15.0)
    result["remove_guides_score"] = offtargets.get("remove_guides_score", 5.0)
    result["remove_guides_n_offtargets"] = offtargets.get("remove_guides_n_offtargets", 100)

    # --- keep_intermediates ---
    result["keep_intermediates"] = config.get("keep_intermediates", True)

    return result


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
        "be_preset": valid["be_preset"],
        "strand_used": valid["strand_used"],
    })
    out.to_csv(output_path, sep="\t", index=False)


def _run_beguidedoff(filter_file, presets_file, genome_npz, max_mismatches,
                       count_threshold_score, output_file, batch_size=1500):
    """Run BEguidedOFF off-target search on the filtered guide file.

    Guides are loaded explicitly and processed in batches so the pipeline
    stays within BEguidedOFF's internal limits regardless of input size.

    Returns the BEguidedOFF result DataFrame.
    """
    try:
        from beguidedoff import beguidedoff
        from beguidedoff.io import load_guides_from_tsv
    except ImportError as exc:
        raise PipelineError(
            "BEguidedOFF is not installed. Off-target search is unavailable.\n"
            "Install it from: ssh://git@gitlab.mn.tu-dresden.de/buchholzlab/beguidedoff.git\n"
            "Or set search_offtargets=false in your config to skip off-target search."
        ) from exc

    finder = beguidedoff(
        genome_path=genome_npz,
        max_mismatches=max_mismatches,
        n_top_offtargets=1,
        scoring_method="CCTOP",
    )

    # Load guides explicitly so we can batch them.
    guides = load_guides_from_tsv(str(filter_file), presets_path=str(presets_file))

    total_guides = len(guides)
    n_batches = (total_guides + batch_size - 1) // batch_size
    dfs = []
    for i in range(0, total_guides, batch_size):
        batch_num = i // batch_size + 1
        batch = guides[i:i + batch_size]
        logger.info(
            "BEguidedOFF batch %d/%d (%d guides)", batch_num, n_batches, len(batch)
        )
        ot_df = finder.search(
            guides=batch,
            return_all=False,
            count_score_threshold=count_threshold_score,
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


def run_pipeline(config, start_step=1, end_step=5):
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
    if start_step <= 1 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 1/5: Guidefinder - Generating guide RNAs")
        logger.info("=" * 60)
        presets_file = cfg["presets_file"]
        logger.info(f"Loading presets from: {presets_file}")
        presets = load_presets(presets_file)

        preset_names = _resolve_preset_names(presets, cfg["preset_arg"])
        logger.info(f"Using presets: {', '.join(preset_names)}")

        genome_file = cfg["genome"]
        logger.info(f"Loading genome from: {genome_file}")
        genome = load_genome(genome_file)
        logger.info(f"Loaded {len(genome)} chromosomes")

        targets_df = pd.read_csv(cfg["targets"], sep="\t")
        original_target_count = len(targets_df)
        logger.info(f"Loaded {original_target_count} targets")

        targets_df, skipped_df = filter_compatible_targets(targets_df, presets, preset_names)
        if len(skipped_df) > 0:
            logger.warning(
                f"{len(skipped_df)}/{original_target_count} target(s) were skipped "
                f"because their mutation type is not supported by the selected presets. "
                f"Skipped: "
                f"{', '.join([f'{r.chromosome}:{r.position}' for r in skipped_df.itertuples()][:5])}"
                f"{' ...' if len(skipped_df) > 5 else ''}"
            )
        if len(targets_df) == 0:
            raise PipelineError(
                "No targets remain after filtering for preset compatibility."
            )
        logger.info(f"Processing {len(targets_df)} compatible targets")

        results = run_guidefinder(
            genome, targets_df, presets, preset_names,
            progress_callback=lambda c, p: logger.debug(f"Processing {c}:{p}")
        )
        df = results_to_dataframe(results)
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
                        if all("doesn't match preset" in s for s in seqs):
                            pure_mismatch.append(f"{chrom}:{pos}")
                        else:
                            pure_pam_fail.append(f"{chrom}:{pos}")

                if pure_mismatch:
                    logger.warning(
                        f"{len(pure_mismatch)}/{total_targets} target(s) use a mutation "
                        f"type not supported by the selected presets: "
                        f"{', '.join(pure_mismatch[:5])}"
                        f"{' ...' if len(pure_mismatch) > 5 else ''}"
                    )
                if pure_pam_fail:
                    logger.warning(
                        f"{len(pure_pam_fail)}/{total_targets} target(s) had no valid PAM"
                    )

        total_result_rows = len(df)
        total_invalid_rows = invalid_mask.sum()
        if total_invalid_rows == total_result_rows:
            raise PipelineError(
                "No valid guides found for the given presets."
            )

        logger.info(
            f"Step 1 complete: {len(df)} results written to {guidefinder_output}"
        )

    # Step 2: Bystander (conditional)
    filter_mode = cfg["filter_mode"]
    needs_bystander = filter_mode != "no_bystanders"

    if start_step <= 2 <= end_step:
        logger.info("=" * 60)
        if needs_bystander:
            logger.info("Step 2/5: Bystander analysis")
        else:
            logger.info("Step 2/5: Bystander analysis (SKIPPED)")
        logger.info("=" * 60)

        if needs_bystander:
            gtf_file = cfg["gff_file"]
            if not gtf_file:
                raise PipelineError(
                    "GFF/GTF file is required for bystander analysis."
                )
            logger.info(f"Loading GTF annotation from: {gtf_file}")
            interval_trees = load_gtf_to_interval_tree(gtf_file)
            logger.info(f"Loaded annotations for {len(interval_trees)} chromosomes")

            df = pd.read_csv(guidefinder_output, sep="\t")
            df = df[df["guide_rna_sequence"] != "**No guide possible**"]
            logger.info(f"Analyzing {len(df)} guides for bystander effects")

            if len(df) > 0:
                if "genomic_context_upstream" not in df.columns:
                    genome_file = cfg["genome"]
                    logger.info(f"Loading genome from: {genome_file}")
                    genome = load_genome(genome_file)
                    logger.info(f"Loaded {len(genome)} chromosomes")
                else:
                    genome = None
                    logger.info("Using guidefinder context columns")

                df = run_bystander(
                    df, interval_trees,
                    progress_callback=lambda idx, total: (
                        logger.info(f"  Bystander analysis: {idx + 1}/{total} guides")
                        if idx % 100 == 0 else None
                    )
                )
                logger.info("Bystander analysis complete")
            df.to_csv(bystander_output, sep="\t", index=False)
        else:
            shutil.copy(str(guidefinder_output), str(bystander_output))
            logger.info(f"Copied guidefinder output to {bystander_output}")

    # Step 3: Filter
    if start_step <= 3 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 3/5: Filtering guides")
        logger.info(f"Filter mode: {filter_mode}")
        logger.info(f"On-target score threshold: {cfg['guide_score_threshold']}")
        logger.info("=" * 60)

        df = pd.read_csv(bystander_output, sep="\t")
        original_count = len(df)
        logger.info(f"Input guides: {original_count}")

        criteria = load_filter_criteria(
            mode=filter_mode,
            guide_score_threshold=cfg["guide_score_threshold"],
        )
        df = filter_guides(df, criteria)
        logger.info(f"After filtering: {len(df)} guides (from {original_count})")

        df.to_csv(filter_file, sep="\t", index=False)
        logger.info(f"Filtered results written to {filter_file}")

    # Step 4: BEguidedOFF off-target search (single step)
    if screen_offtargets and start_step <= 4 <= end_step:
        logger.info("=" * 60)
        logger.info("Step 4/5: BEguidedOFF off-target search")
        logger.info("=" * 60)

        genome_npz = cfg["beguidedoff_genome"]
        if not genome_npz:
            raise PipelineError(
                "beguidedoff_genome is required when search_offtargets=true. "
                "Provide the path to the encoded .npz genome file."
            )

        logger.info(f"BEguidedOFF genome: {genome_npz}")
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
            cfg["presets_file"],
            genome_npz,
            cfg["max_mismatches"],
            cfg["count_threshold_score"],
            offtarget_output,
        )
        logger.info(
            f"BEguidedOFF complete: {len(ot_df)} off-target rows written to "
            f"{offtarget_output}"
        )

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

            merged.to_csv(final_output, sep="\t", index=False)
            size_mb = final_output.stat().st_size / (1024 * 1024)
            logger.info(
                f"Final output written: {final_output} "
                f"({len(merged)} rows, {size_mb:.1f} MB)"
            )
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
