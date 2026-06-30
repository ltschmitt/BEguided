"""Single CLI entry point for BEguided with subcommands."""

import argparse
import sys
from pathlib import Path

import pandas as pd

from beguided.genome.io import load_genome
from beguided.editing.presets import load_presets, get_available_preset_names, get_preset_names_by_type, filter_compatible_targets
from beguided.guidefinder.core import run_genome_wide as run_guidefinder, results_to_dataframe
from beguided.bystander.annotation import load_gtf_to_interval_tree
from beguided.bystander.analysis import run_all as run_bystander
from beguided.filtering.filters import load_filter_criteria, filter_guides
from beguided.pipeline import run_pipeline, PipelineError


def _add_common_io(parser, input_help, output_help):
    parser.add_argument("--input", "-i", required=True, help=input_help)
    parser.add_argument("--output", "-o", required=True, help=output_help)


def cmd_guidefinder(args):
    """Run guidefinder subcommand."""
    try:
        presets = load_presets(args.presets_file)
    except FileNotFoundError:
        print(f"Error: Presets file not found: {args.presets_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading presets: {e}", file=sys.stderr)
        sys.exit(1)

    available = get_available_preset_names(presets)
    type_map = get_preset_names_by_type(presets)
    preset_arg = args.presets
    if preset_arg.lower() == "all":
        preset_names = sorted(available)
    elif preset_arg.lower() == "all abes":
        preset_names = type_map["abe"]
    elif preset_arg.lower() == "all cbes":
        preset_names = type_map["cbe"]
    else:
        preset_names = [p.strip() for p in preset_arg.split(",")]
        for name in preset_names:
            if name not in available:
                print(
                    f"Error: Unknown preset '{name}'. "
                    f"Available: {', '.join(sorted(available))}",
                    file=sys.stderr,
                )
                sys.exit(1)

    genome = load_genome(args.genome)
    df = pd.read_csv(args.input, sep="\t")
    original_count = len(df)

    df, skipped = filter_compatible_targets(df, presets, preset_names)
    if len(skipped) > 0:
        print(
            f"Warning: {len(skipped)}/{original_count} target(s) were skipped "
            f"because their mutation type is not supported by the selected presets "
            f"(supports both sense and antisense complement directions). "
            f"Skipped targets: "
            f"{', '.join([f'{r.chromosome}:{r.position}' for r in skipped.itertuples()][:5])}"
            f"{' ...' if len(skipped) > 5 else ''}",
            file=sys.stderr,
        )
    if len(df) == 0:
        print(
            "Error: No targets remain after filtering for preset compatibility. "
            "Check that your target mutations match the selected preset editing directions "
            "(e.g., ABE presets handle A→G / T→C; CBE presets handle C→T / G→A).",
            file=sys.stderr,
        )
        sys.exit(1)

    results = run_guidefinder(
        genome, df, presets, preset_names,
        progress_callback=lambda c, p: print(f"{c}:{p}")
    )
    result_df = results_to_dataframe(results)
    result_df.to_csv(args.output, sep="\t", index=False)


def cmd_bystander(args):
    """Run bystander scanner subcommand."""
    df = pd.read_csv(args.input, sep="\t")
    df = df[df["guide_rna_sequence"] != "**No guide possible**"]

    if len(df) == 0:
        print("No valid guides to analyze")
        df.to_csv(args.output, sep="\t", index=False)
        return

    has_context = (
        "genomic_context_upstream" in df.columns
        and "genomic_context_downstream" in df.columns
    )

    if has_context:
        if args.genome:
            print("Note: Context columns present; --genome is ignored.")
        genome = None
    else:
        if not args.genome:
            print(
                "Error: --genome is required when input lacks context columns. "
                "Run guidefinder first to produce context columns.",
                file=sys.stderr,
            )
            sys.exit(1)
        genome = load_genome(args.genome)

    interval_trees = load_gtf_to_interval_tree(args.gtf_annotation)

    input_chroms = set(df["chromosome"].astype(str))
    gtf_chroms = set(interval_trees.keys())
    unmatched = input_chroms - gtf_chroms
    if unmatched:
        stripped = {c.replace("chr", "") for c in unmatched}
        if stripped.issubset(gtf_chroms):
            df["chromosome"] = df["chromosome"].astype(str).str.replace("chr", "")
        elif {f"chr{c}" for c in unmatched if not c.startswith("chr")}.issubset(gtf_chroms):
            df["chromosome"] = "chr" + df["chromosome"].astype(str)
        else:
            print("Error: Chromosome names don't match between input and GTF", file=sys.stderr)
            sys.exit(1)

    df = run_bystander(
        df, interval_trees,
        progress_callback=lambda idx, total: print(f"  Processed {idx}/{total} guides")
    )
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nBystander effect summary:")
    print(df["bystander_effect"].value_counts())


def cmd_filter(args):
    """Run filter subcommand."""
    df = pd.read_csv(args.input, sep="\t")
    criteria = load_filter_criteria(config_file=args.config, mode=args.mode)
    df = filter_guides(df, criteria)
    df.to_csv(args.output, sep="\t", index=False)
    
    # Report filtering summary
    threshold = criteria.get("guide_score_threshold")
    if threshold is not None:
        print(f"On-target score threshold: {threshold}")
    print(f"Filtered results written to {args.output}")


def cmd_web(args):
    """Start the BEguided web server."""
    try:
        import uvicorn  # noqa: F401
        from beguided.web.app import create_app
        from beguided.web.config import WebConfig, set_config
    except ImportError as exc:
        print(
            "Error: The web server requires extra dependencies.\n"
            "Install them with one of the following:\n"
            "  pip install -e \".[web]\"\n"
            "  pixi install --features web",
            file=sys.stderr,
        )
        sys.exit(1)

    kwargs = {}
    if args.presets_file is not None:
        kwargs["presets_file"] = args.presets_file
    if args.genome_fasta is not None:
        kwargs["genome_fasta"] = args.genome_fasta
    if args.gtf_path is not None:
        kwargs["gtf_path"] = args.gtf_path
    if args.beguidedoff_genome is not None:
        kwargs["beguidedoff_genome"] = args.beguidedoff_genome
    cfg = WebConfig(**kwargs)
    set_config(cfg)

    app = create_app()
    uvicorn.run(
        app,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


def cmd_run_pipeline(args):
    """Run complete pipeline subcommand."""
    import yaml
    from beguided.logging_config import setup_logging

    level = "DEBUG" if args.verbose else "INFO"
    setup_logging(level=level)

    if not Path(args.config).exists():
        print(f"Error: Config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    for key in ("input", "output_prefix"):
        if key not in config:
            print(f"Error: Missing required config key: {key}", file=sys.stderr)
            sys.exit(1)

    start_step = args.from_step if args.from_step else 1

    try:
        final = run_pipeline(config, start_step=start_step, end_step=5)
        print(f"\nPipeline completed: {final}")
    except PipelineError as e:
        print(f"Pipeline failed: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        prog="beguided",
        description="Cas9 Base Editor Guide RNA Generator",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # guidefinder
    p_gf = subparsers.add_parser("guidefinder", help="Generate guide RNAs")
    _add_common_io(p_gf, "Input TSV file", "Output TSV file")
    p_gf.add_argument("--genome", "-g", required=True, help="Genome FASTA file")
    p_gf.add_argument("--presets", default="ABE8e", help="Comma-separated presets, or 'all', 'all ABEs', 'all CBEs'")
    p_gf.add_argument("--presets-file", default="be_presets.tsv", help="Presets TSV file")
    p_gf.set_defaults(func=cmd_guidefinder)

    # bystander
    p_bs = subparsers.add_parser("bystander", help="Analyze bystander effects")
    _add_common_io(p_bs, "Input TSV from guidefinder", "Output TSV file")
    p_bs.add_argument("--genome", "-g", help="Genome FASTA file (only if input lacks context columns)")
    p_bs.add_argument("--gtf-annotation", required=True, help="GTF annotation file")
    p_bs.set_defaults(func=cmd_bystander)

    # filter
    p_f = subparsers.add_parser("filter", help="Filter guides")
    _add_common_io(p_f, "Input TSV", "Output TSV")
    p_f.add_argument("--config", "-c", help="Pipeline config YAML file")
    p_f.add_argument("--mode", "-m", choices=["all", "no_bystanders", "silent"], help="Filter mode")
    p_f.set_defaults(func=cmd_filter)

    # web
    p_web = subparsers.add_parser("web", help="Start the BEguided web server")
    p_web.add_argument("--host", default="0.0.0.0", help="Bind host (default: 0.0.0.0)")
    p_web.add_argument("--port", type=int, default=8000, help="Bind port (default: 8000)")
    p_web.add_argument("--reload", action="store_true", help="Enable auto-reload for development")
    p_web.add_argument("--presets-file", default=None, help="Presets TSV file (default: env BEGUIDED_PRESETS_FILE or be_presets.tsv)")
    p_web.add_argument("--genome-fasta", default=None, help="Genome FASTA file (default: env BEGUIDED_GENOME_FASTA or genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa)")
    p_web.add_argument("--gtf-path", default=None, help="GTF annotation file (default: env BEGUIDED_GTF_PATH or genome/Homo_sapiens.GRCh38.112.gtf)")
    p_web.add_argument("--beguidedoff-genome", default=None, help="BEguidedOFF encoded genome .npz (default: env BEGUIDED_BEGUIDEDOFF_GENOME)")
    p_web.set_defaults(func=cmd_web)

    # run-pipeline
    p_pipe = subparsers.add_parser("run-pipeline", help="Run complete pipeline")
    p_pipe.add_argument("--config", "-c", required=True, help="Pipeline config YAML")
    p_pipe.add_argument("--from-step", type=int, choices=range(1, 6), help="Start from step (1-5)")
    p_pipe.add_argument("--verbose", "-v", action="store_true", help="Enable DEBUG logging")
    p_pipe.set_defaults(func=cmd_run_pipeline)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
