"""CLI entry point for the BEguided web server."""

import argparse
import sys


def main():
    """Start the BEguided web server."""
    try:
        import uvicorn  # noqa: F401
        from beguided.web.app import create_app
        from beguided.web.config import WebConfig, set_config
        from beguided.web.safety import validate_operator_paths
    except ImportError as exc:
        print(
            "Error: The web server requires extra dependencies.\n"
            "Install them with one of the following:\n"
            '  pip install -e ".[web]"\n'
            "  pixi install --features web",
            file=sys.stderr,
        )
        sys.exit(1)

    parser = argparse.ArgumentParser(
        prog="beguided-web",
        description="Start the BEguided web server.",
    )
    parser.add_argument("--host", default="0.0.0.0", help="Bind host (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8000, help="Bind port (default: 8000)")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload for development")
    parser.add_argument(
        "--base-editor-profiles-file",
        default=None,
        help="Base editor profiles TSV file (default: env BEGUIDED_BASE_EDITOR_PROFILES_FILE or base_editor_profiles.tsv)",
    )
    parser.add_argument(
        "--genome-fasta",
        default=None,
        help="Genome FASTA file (default: env BEGUIDED_GENOME_FASTA or genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa)",
    )
    parser.add_argument(
        "--gtf-path",
        default=None,
        help="GTF annotation file (default: env BEGUIDED_GTF_PATH or genome/Homo_sapiens.GRCh38.112.gtf)",
    )
    parser.add_argument(
        "--beguidedoff-genome",
        default=None,
        help="BEguidedOFF encoded genome .npz (default: env BEGUIDED_BEGUIDEDOFF_GENOME)",
    )
    parser.add_argument(
        "--db-path",
        default=None,
        help="SQLite path for the job store (default: env BEGUIDED_DB_PATH or /var/lib/beguided/jobs.db)",
    )
    parser.add_argument(
        "--queue-max",
        type=int,
        default=None,
        help="Max pending jobs in queue (default: 8)",
    )
    parser.add_argument(
        "--cooldown-seconds",
        type=float,
        default=None,
        help="Per-IP cooldown between submits (default: 60)",
    )
    parser.add_argument(
        "--max-upload-mb",
        type=int,
        default=None,
        help="Cap on request body size in MB (default: 10)",
    )
    parser.add_argument(
        "--preload-genome",
        dest="preload_genome",
        action="store_true",
        default=None,
        help="Preload the genome FASTA at startup (background thread). Default: enabled.",
    )
    parser.add_argument(
        "--no-preload-genome",
        dest="preload_genome",
        action="store_false",
        default=None,
        help="Disable genome preload. The genome is then parsed per-job instead.",
    )
    parser.add_argument(
        "--genome-label",
        default=None,
        help=(
            "Display label for the genome (e.g. 'GRCh38.112'). "
            "Empty = derive from FASTA filename. "
            "Env: BEGUIDED_GENOME_LABEL."
        ),
    )
    parser.add_argument(
        "--preload-gtf",
        dest="preload_gtf",
        action="store_true",
        default=None,
        help="Preload the GTF annotation at startup (background thread). Default: enabled.",
    )
    parser.add_argument(
        "--no-preload-gtf",
        dest="preload_gtf",
        action="store_false",
        default=None,
        help="Disable GTF preload. The GTF is then parsed per-job instead.",
    )

    args = parser.parse_args()

    kwargs = {}
    if args.base_editor_profiles_file is not None:
        kwargs["base_editor_profiles_file"] = args.base_editor_profiles_file
    if args.genome_fasta is not None:
        kwargs["genome_fasta"] = args.genome_fasta
    if args.gtf_path is not None:
        kwargs["gtf_path"] = args.gtf_path
    if args.beguidedoff_genome is not None:
        kwargs["beguidedoff_genome"] = args.beguidedoff_genome
    if args.db_path is not None:
        kwargs["db_path"] = args.db_path
    if args.queue_max is not None:
        kwargs["queue_max"] = args.queue_max
    if args.cooldown_seconds is not None:
        kwargs["cooldown_seconds"] = args.cooldown_seconds
    if args.max_upload_mb is not None:
        kwargs["max_upload_mb"] = args.max_upload_mb
    if args.preload_genome is not None:
        kwargs["preload_genome"] = args.preload_genome
    if args.genome_label is not None:
        kwargs["genome_label"] = args.genome_label
    if args.preload_gtf is not None:
        kwargs["preload_gtf"] = args.preload_gtf
    cfg = WebConfig(**kwargs)
    set_config(cfg)

    # Fail fast on missing presets file (the genome FASTA / GTF / .npz are
    # best-effort: some deployments only do steps 4-5 where they aren't needed).
    try:
        validate_operator_paths(cfg)
    except (FileNotFoundError, PermissionError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    # Reminder: do NOT run with --workers N > 1. Each worker would claim the
    # GPU independently and the FIFO guarantee would break.
    print(
        f"BEguided web: queue_max={cfg.queue_max} cooldown={int(cfg.cooldown_seconds)}s "
        f"db={cfg.db_path} upload_cap={cfg.max_upload_mb}MB",
        file=sys.stderr,
    )

    app = create_app()
    uvicorn.run(
        app,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
