"""CLI entry point for BEguided pipeline execution."""

import argparse
import sys
from pathlib import Path

from beguided.pipeline import run_pipeline, PipelineError


def main():
    parser = argparse.ArgumentParser(
        prog="beguided",
        description="Cas9 Base Editor Guide RNA Generator — run the full pipeline.",
    )
    parser.add_argument("--config", "-c", required=True, help="Pipeline config YAML")
    parser.add_argument("--from-step", type=int, choices=range(1, 6), help="Start from step (1-5)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable DEBUG logging")
    parser.add_argument(
        "--passthrough-extras",
        dest="passthrough_extras",
        action="store_true",
        default=None,
        help="Pass non-core columns from the input TSV through to the output (default: enabled).",
    )
    parser.add_argument(
        "--no-passthrough-extras",
        dest="passthrough_extras",
        action="store_false",
        default=None,
        help="Drop non-core columns from the input TSV.",
    )

    args = parser.parse_args()

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

    if args.passthrough_extras is not None:
        config["passthrough_extras"] = args.passthrough_extras

    try:
        final = run_pipeline(config, start_step=start_step, end_step=5)
        print(f"\nPipeline completed: {final}")
    except PipelineError as e:
        print(f"Pipeline failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
