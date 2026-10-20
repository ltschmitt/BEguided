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
        help="Base editor profiles TSV file (default: env BEGUIDED_BASE_EDITOR_PROFILES_FILE or configs/base_editor_profiles.tsv)",
    )
    parser.add_argument(
        "--genome-fasta",
        default=None,
        help="Human genome FASTA file (default: env BEGUIDED_GENOME_FASTA or genome/Homo_sapiens.GRCh38.116.fa)",
    )
    parser.add_argument(
        "--gtf-path",
        default=None,
        help="Human GTF annotation file (default: env BEGUIDED_GTF_PATH or genome/Homo_sapiens.GRCh38.116.gtf)",
    )
    parser.add_argument(
        "--mouse-genome-fasta",
        default=None,
        help="Mouse genome FASTA file (default: env BEGUIDED_MOUSE_GENOME_FASTA or genome/Mus_musculus.GRCm39.2026_04.fa)",
    )
    parser.add_argument(
        "--mouse-gtf-path",
        default=None,
        help="Mouse GTF annotation file (default: env BEGUIDED_MOUSE_GTF_PATH or genome/Mus_musculus.GRCm39.2026_04.gtf)",
    )
    parser.add_argument(
        "--mouse-genome-label",
        default=None,
        help="Assembly label for the mouse genome, shown in the UI as 'Mus musculus (<label>)' (default: 'GRCm39'). Env: BEGUIDED_MOUSE_GENOME_LABEL.",
    )
    parser.add_argument(
        "--default-genome",
        default=None,
        choices=["human", "mouse"],
        help="Genome selected by default in the UI: human|mouse (default: human). Env: BEGUIDED_DEFAULT_GENOME.",
    )
    parser.add_argument(
        "--beguidedoff-genome",
        default=None,
        help="BEguidedOFF encoded genome .npz (default: env BEGUIDED_BEGUIDEDOFF_GENOME)",
    )
    parser.add_argument(
        "--beguidedoff-mode",
        default=None,
        choices=["fast", "extensive"],
        help="BEguidedOFF search mode: fast (CPU FM-index filter, default) or "
        "extensive (GPU brute force). Env: BEGUIDED_BEGUIDEDOFF_MODE.",
    )
    parser.add_argument(
        "--beguidedoff-index",
        default=None,
        help="BEguidedOFF FM-index .fmi for fast mode (default: env "
        "BEGUIDED_BEGUIDEDOFF_INDEX, or the .npz path with an .fmi suffix)",
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
        "--max-input-targets",
        type=int,
        default=None,
        help=(
            "Max target rows per upload (default: env BEGUIDED_MAX_INPUT_TARGETS "
            "or 500)"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "Number of parallel pipeline worker threads (default: env "
            "BEGUIDED_WORKERS or 1). Threads share the preloaded genome/GTF, so "
            "more workers cost CPU time but not RAM; keep it <= CPU cores."
        ),
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
            "Assembly label for the human genome, shown in the UI as "
            "'Homo sapiens (<label>)' (e.g. 'GRCh38'). "
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
    parser.add_argument(
        "--job-retention-hours",
        type=float,
        default=None,
        help=(
            "Retention window for finished jobs in jobs.db (default: 24). "
            "Sweeper deletes any finished row older than this. "
            "Set to 0 to disable. Env: BEGUIDED_JOB_RETENTION_HOURS."
        ),
    )
    parser.add_argument(
        "--trusted-proxy-ips",
        default=None,
        help=(
            "Comma-separated list of upstream proxy IPs trusted to set "
            "X-Forwarded-For / X-Forwarded-Proto (default: '127.0.0.1,::1'). "
            "Env: BEGUIDED_TRUSTED_PROXY_IPS."
        ),
    )
    parser.add_argument(
        "--log-ip-salt",
        default=None,
        help=(
            "Salt used to hash client IPs in log lines (GDPR hygiene). "
            "Empty = no hashing (development only). "
            "Env: BEGUIDED_LOG_IP_SALT."
        ),
    )
    parser.add_argument(
        "--admin-token",
        default=None,
        help=(
            "Static token required by POST /admin/purge-all. "
            "Empty = endpoint refuses every request. "
            "Env: BEGUIDED_ADMIN_TOKEN."
        ),
    )
    parser.add_argument(
        "--legal-config",
        default=None,
        help=(
            "Path to the YAML file with legal-page content (institution "
            "name, address, DPO, competent DPA, etc.). Default: "
            "/etc/beguided/legal.yaml. Empty = serve empty legal pages "
            "with a startup warning. The file is gitignored; only the "
            "schema reference at configs/web-legal.yaml is committed. "
            "Env: BEGUIDED_LEGAL_CONFIG."
        ),
    )
    parser.add_argument(
        "--usage-stats-path",
        default=None,
        help=(
            "Base path for the anonymous usage-statistics CSV. Written "
            "monthly as <stem>_YYYY-MM.csv. Default: "
            "/var/lib/beguided/usage_stats.csv. Env: BEGUIDED_USAGE_STATS_PATH."
        ),
    )
    parser.add_argument(
        "--no-usage-stats",
        dest="usage_stats_enabled",
        action="store_false",
        default=None,
        help=(
            "Disable anonymous usage-statistics recording. Default: "
            "enabled. Env: BEGUIDED_USAGE_STATS_ENABLED."
        ),
    )

    args = parser.parse_args()

    kwargs = {}
    if args.base_editor_profiles_file is not None:
        kwargs["base_editor_profiles_file"] = args.base_editor_profiles_file
    if args.genome_fasta is not None:
        kwargs["genome_fasta"] = args.genome_fasta
    if args.gtf_path is not None:
        kwargs["gtf_path"] = args.gtf_path
    if args.mouse_genome_fasta is not None:
        kwargs["mouse_genome_fasta"] = args.mouse_genome_fasta
    if args.mouse_gtf_path is not None:
        kwargs["mouse_gtf_path"] = args.mouse_gtf_path
    if args.mouse_genome_label is not None:
        kwargs["mouse_genome_label"] = args.mouse_genome_label
    if args.default_genome is not None:
        kwargs["default_genome"] = args.default_genome
    if args.beguidedoff_genome is not None:
        kwargs["beguidedoff_genome"] = args.beguidedoff_genome
    if args.beguidedoff_mode is not None:
        kwargs["beguidedoff_mode"] = args.beguidedoff_mode
    if args.beguidedoff_index is not None:
        kwargs["beguidedoff_index"] = args.beguidedoff_index
    if args.db_path is not None:
        kwargs["db_path"] = args.db_path
    if args.queue_max is not None:
        kwargs["queue_max"] = args.queue_max
    if args.cooldown_seconds is not None:
        kwargs["cooldown_seconds"] = args.cooldown_seconds
    if args.max_upload_mb is not None:
        kwargs["max_upload_mb"] = args.max_upload_mb
    if args.max_input_targets is not None:
        kwargs["max_input_targets"] = args.max_input_targets
    if args.workers is not None:
        kwargs["workers"] = args.workers
    if args.preload_genome is not None:
        kwargs["preload_genome"] = args.preload_genome
    if args.genome_label is not None:
        kwargs["genome_label"] = args.genome_label
    if args.preload_gtf is not None:
        kwargs["preload_gtf"] = args.preload_gtf
    if args.job_retention_hours is not None:
        kwargs["job_retention_hours"] = float(args.job_retention_hours)
    if args.trusted_proxy_ips is not None:
        kwargs["trusted_proxy_ips"] = [
            s.strip() for s in args.trusted_proxy_ips.split(",") if s.strip()
        ]
    if args.log_ip_salt is not None:
        kwargs["log_ip_salt"] = args.log_ip_salt
    if args.admin_token is not None:
        kwargs["admin_token"] = args.admin_token
    if args.legal_config is not None:
        kwargs["legal_config_path"] = args.legal_config
    if args.usage_stats_path is not None:
        kwargs["usage_stats_path"] = args.usage_stats_path
    if args.usage_stats_enabled is not None:
        kwargs["usage_stats_enabled"] = args.usage_stats_enabled
    cfg = WebConfig(**kwargs)
    set_config(cfg)

    # Fail fast on missing presets file (the genome FASTA / GTF / .npz are
    # best-effort: some deployments only do steps 4-5 where they aren't needed).
    try:
        validate_operator_paths(cfg)
    except (FileNotFoundError, PermissionError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    # Reminder: do NOT pass uvicorn a worker count > 1. That would spawn
    # separate processes, each reloading the genome (wasting GBs of RAM) and
    # each with its own queue/cooldown. Use --workers instead, which adds
    # in-process scheduler threads that share the preloaded genome.
    print(
        f"BEguided web: workers={cfg.workers} queue_max={cfg.queue_max} "
        f"cooldown={int(cfg.cooldown_seconds)}s db={cfg.db_path} "
        f"upload_cap={cfg.max_upload_mb}MB max_targets={cfg.max_input_targets}",
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
