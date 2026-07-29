# CLI usage

The CLI entry point is `beguided`. The pipeline is configured entirely via a YAML file (see [`configs/cli.yaml`](../../configs/cli.yaml) for the full template).

## Running the full pipeline

```bash
beguided --config configs/cli.yaml
```

This runs every step end-to-end: guide generation → bystander analysis → filtering → (optional) off-target search → final merge.

## Running a subset of steps

Use `--from-step` to start at a later step, useful when you've already produced intermediate files and only want to re-run the tail:

```bash
# Start from step 3 (re-run filtering + off-target + merge)
beguided --config configs/cli.yaml --from-step 3
```

`-v` switches the log level to DEBUG. `--no-passthrough-extras` drops every column from your input TSV that is not one of the four core columns (`chromosome`, `position`, `target_base`, `desired_base`).

## Web server

```bash
# Default: http://localhost:8000
beguided-web

# Custom address
beguided-web --host 0.0.0.0 --port 8080

# Development with auto-reload
beguided-web --reload
```

For full deployment instructions (reverse proxy, TLS, compliance), see [`web-deployment.md`](web-deployment.md).

## Configuration file

The pipeline is configured via a YAML file. Key sections:

```yaml
input:
  targets: "path/to/targets.tsv"
  genome: "path/to/genome.fa"
  base_editor_profiles: "configs/base_editor_profiles.tsv"
  gff_file: "path/to/annotation.gtf"

output_prefix: "./results"
use_base_editors: "all"
keep_bystanders: "all"
guide_score_threshold: 0.5
search_offtargets: true

offtargets:
  beguidedoff_genome: "path/to/genome_otbe.npz"
  mismatch_tolerance: 5
  count_threshold_score: 15.0
  remove_guides_score: 5.0
  remove_guides_n_offtargets: null
```

The full annotated template at [`configs/cli.yaml`](../../configs/cli.yaml) walks through every key. The most common knobs:

| Key | Default | Description |
|---|---|---|
| `guide_score_threshold` | `0.5` | CRISPRater **minimum** score to keep a guide. Range 0–1 (higher = better). `0` disables. |
| `gc_content_filter` | `off` | `off` / `good` (40–60%) / `acceptable` (30–70%). Guides with missing GC are always kept. |
| `filter_homopolymers` | `false` | Drop guides with a run of ≥4 identical bases (Pol III termination signal). |
| `trinuc_edit_score_min` | `null` | Minimum Marquart et al. 2021 ABE8e trinuc context score (0–1) to keep an ABE guide. CBE guides are always kept (NaN by construction). |
| `keep_bystanders` | `all` | `all` / `no_bystanders` / `silent` |
| `use_base_editors` | `all` | `all` / `all ABEs` / `all CBEs` / comma-separated profile names |
| `offtargets.mismatch_tolerance` | `5` | Max mismatches in 20-nt guide |
| `offtargets.count_threshold_score` | `15.0` | CCTOP **maximum** for counting off-targets |
| `offtargets.remove_guides_score` | `5.0` | CCTOP **minimum** to keep a guide |
| `offtargets.remove_guides_n_offtargets` | `null` | Max off-target count to keep a guide |

CCTOP semantics: **lower** score = more concerning (mismatches near PAM). A guide with no off-targets gets `worst_offtarget_score = inf`.

See [`scoring.md`](scoring.md) for the detailed semantics of every score.
