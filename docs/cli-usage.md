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
  # "extensive" = GPU brute force (CLI default); "fast" = CPU FM-index filter
  mode: "extensive"
  beguidedoff_genome: "path/to/genome_otbe.npz"
  # FM-index for fast mode; empty = derived from the .npz (same basename, .fmi)
  index_path: ""
  max_mismatches: 5              # extensive mode only (legacy alias: mismatch_tolerance)
  count_threshold_score: 15.0
  remove_guides_score: 5.0
  remove_guides_n_offtargets: 100
```

The full annotated template at [`configs/cli.yaml`](../../configs/cli.yaml) walks through every key. The most common knobs:

| Key | Default | Description |
|---|---|---|
| `guide_score_threshold` | `0.5` | CRISPRater **minimum** score to keep a guide. Range 0–1 (higher = better). `0` disables. |
| `gc_content_filter` | `off` | `off` / `high` (40–60%) / `medium` (30–70%). The deprecated names `good` (= `high`) and `acceptable` (= `medium`) are still accepted. Guides with missing GC are always kept. |
| `filter_homopolymers` | `false` | Drop guides with a run of ≥4 identical bases (Pol III termination signal). |
| `trinuc_edit_score_min` | `null` | Minimum Marquart et al. 2021 ABE8e trinuc context score (0–1) to keep an ABE guide. CBE guides are always kept (NaN by construction). |
| `keep_targets_without_guides` | `false` | Keep one row per input target in `beguided_results.tsv` even when no guide could be designed or the mutation type is unsupported. Placeholder rows keep the input identity/bases columns; `guide_rna_sequence` is `**No guide possible**` or `**Base change incompatible**`, all other result fields empty. Added only in the final step. (The web UI's filtered TSV/XLSX download can additionally show `**Guides removed by filter**` placeholder rows for targets whose guides were all removed by the active filters.) |
| `keep_bystanders` | `all` | `all` / `none` / `synonymous` |
| `use_base_editors` | `all` | `all` / `all ABEs` / `all CBEs` / comma-separated profile names |
| `offtargets.mode` | `extensive` | `extensive` (GPU brute force) or `fast` (CPU FM-index filter; needs an `.fmi`, see `offtargets.index_path`). The web server always sends an explicit value (`fast`). |
| `offtargets.beguidedoff_genome` | — | Encoded `.npz` genome. Required for `extensive`; in `fast` mode only used to derive the default `.fmi` path. |
| `offtargets.index_path` | `""` | FM-index `.fmi` for `fast` mode. Empty = derived from `beguidedoff_genome` (same basename, `.fmi` suffix). |
| `offtargets.max_mismatches` | `5` | Max mismatches in 20-nt guide (`extensive` mode only; `mismatch_tolerance` is accepted as a legacy alias) |
| `offtargets.count_threshold_score` | `15.0` | CCTOP **maximum** for counting off-targets |
| `offtargets.remove_guides_score` | `5.0` | CCTOP **minimum** to keep a guide |
| `offtargets.remove_guides_n_offtargets` | `100` | Max off-target count to keep a guide. Set to `null` to disable the count limit. |

CCTOP semantics: **lower** score = more concerning (mismatches near PAM). A guide with no off-targets gets `worst_offtarget_score = NaN` (blank in the TSV), which always passes the removal filters.

See [`scoring.md`](scoring.md) for the detailed semantics of every score.
