# Agent Guidelines for BEguided

This document provides guidelines for AI agents working on the BEguided codebase.

## Project Overview

BEguided is a Python package for generating Cas9 base editor guide RNA sequences. It is managed with Pixi (single `pyproject.toml` manifest) and uses BEguidedOFF (standalone sibling repo) for GPU-accelerated off-target search. All functionality lives in the `beguided/` package.

## Architecture

```
beguided/                  # Main package
├── cli.py                   # Single CLI entry with argparse subcommands
├── pipeline.py              # Pipeline orchestrator (function calls, no subprocess)
├── models.py                # Shared dataclasses (BaseEditorConfig, GuideResult, etc.)
├── genome/                  # Shared genome I/O and sequence operations
│   ├── io.py                # FASTA loading
│   ├── sequence.py          # Reverse complement, mutation helpers
│   ├── coordinates.py       # Coordinate math (0-based vs 1-based)
│   └── context.py           # Genomic region reconstruction
├── editing/                 # Core base-editing logic
│   ├── presets.py           # Base editor profile loading, PAM regex, validation
│   ├── search.py            # PAM regex search
│   ├── mutations.py         # Editing simulation, bystander detection
│   ├── context.py           # Trinucleotide edit-context lookup (ABE, Marquart 2021)
│   └── data/
│       └── edit_context_norm_proportion.tsv  # Bundled X-A-Y proportion table
├── guidefinder/             # Step 1: Guide generation
│   └── core.py              # process_target, run_genome_wide, results_to_dataframe
├── bystander/               # Step 2: Bystander analysis
│   ├── annotation.py        # GTF parsing, interval trees
│   ├── protein.py           # Translation, genetic code
│   └── analysis.py          # analyze_guide, run_all
├── off_target/              # Step 4-6: BEguidedOFF integration
│   └── merger.py            # Aggregate BEguidedOFF results + merge with input
├── filtering/               # Step 3
│   └── filters.py           # Guide filtering logic
└── web/                     # FastAPI web interface
    ├── app.py               # FastAPI routes
    ├── config.py            # Web config and base editor groups
    ├── sse_logger.py        # SSE log streaming handler
    ├── tasks.py             # Async pipeline task wrapper
    ├── sequence_input.py    # Sequence-mode parser ([X>Y] bracket syntax)
    ├── templates/
    │   ├── base.html
    │   └── index.html
    └── static/
        ├── style.css
        └── app.js
```

### Removed modules

- `beguided/off_target/casoff_input.py` — deleted (Cas-OFFinder input generation)
- `beguided/off_target/cfd.py` — deleted (CFD scoring)
- Root-level scripts (`beguided.py`, `bystander_scanner.py`, etc.) — removed

## Build and Run Commands

### Environment Setup

```bash
# Pixi is configured inside pyproject.toml (no separate pixi.toml)

# Core environment (CLI + tests)
pixi install

# + Web server dependencies
pixi install --features web

# Install BEguided in editable mode (always required)
pip install -e .

# Run the pipeline
beguided --config config.yaml
```

### CLI Commands

```
beguided                 --config config.yaml [--from-step 3] [-v] \
                         [--passthrough-extras | --no-passthrough-extras]
beguided-web             [--host 0.0.0.0] [--port 8000] [--reload] \
                         [--base-editor-profiles-file ...] [--genome-fasta ...] \
                         [--gtf-path ...] [--beguidedoff-genome ...] \
                         [--preload-genome | --no-preload-genome] \
                         [--preload-gtf | --no-preload-gtf] \
                         [--genome-label ...] \
                         [--job-retention-hours 24] \
                         [--allowed-hosts www.buchholz.org,buchholz.org] \
                         [--trusted-proxy-ips 127.0.0.1,::1] \
                         [--log-ip-salt ...] [--admin-token ...] \
                         [--legal-config /etc/beguided/legal.yaml]
                         # Requires web extra: pixi run --environment web beguided-web
```

### Deployment configuration

The project uses **no `.env` file**. Every deployment value is read from
shell environment variables and `docker-compose.yml` ships inline defaults
for each one, so `docker compose up` works with zero configuration. Override
per host via `docker-compose.override.yml` or by exporting the variable in
your shell before `docker compose up`. See `docker-compose.yml` comments
for the full list.

Legal-page content (institution name, address, DPO, competent DPA, …)
lives in a separate gitignored `legal.yaml` file. The committed template
is `configs/web-legal.yaml`. The container path is set via
`BEGUIDED_LEGAL_CONFIG` (default `/etc/beguided/legal.yaml`); if the file
is missing the legal pages render with empty values + a startup warning.
See `compliance_safety.md` for the full pre-launch checklist.

### Running Tests

```bash
# Pixi task (uses pytest defined in pyproject.toml)
pixi run test

# Or directly
python -m pytest tests/ -v

# Specific module
tests/test_guidefinder/ -v

# Specific test
python -m pytest tests/test_editing/test_presets.py::TestPAMToRegex -v
```

## Packaging: Single `pyproject.toml`

Pixi configuration lives in `[tool.pixi.*]` sections inside `pyproject.toml`. There is **no `pixi.toml`**.

Dependency sources:
- `[project.dependencies]` — standard Python dependencies, shared with pip/setuptools.
- `[project.optional-dependencies]` — optional groups (`test = ["pytest>=9.0.3"]`).
- `[tool.pixi.dependencies]` — Conda-specific overrides/additions.
- `[tool.pixi.feature.<name>]` — optional Pixi features (e.g. `cuda`, `test`).

### Pixi features

| Feature | Adds | Install |
|---|---|---|
| (default) | Core runtime deps + test | `pixi install` |
| `web` | `fastapi`, `uvicorn`, `python-multipart`, `jinja2` | `pixi install --features web` |
| `cuda` | `pytorch-cuda>=12.4` | `pixi install --features cuda` |
| `test` | `pytest>=9.0.3`, `httpx` | Already in default |

### Docker targets

| Target | Base image | CUDA | Web |
|---|---|---|---|
| `base` | `python:3.12-slim-bookworm` | No | No |
| `web` | extends `base` | No | Yes |
| `base-gpu` | `nvidia/cuda:12.4.1-runtime-ubuntu22.04` | Yes | No |
| `cli-cuda` | extends `base-gpu` | Yes | No |
| `web-cuda` | extends `base-gpu` | Yes | Yes |

Pixi resolves the environment by combining `[project.dependencies]` + `[tool.pixi.dependencies]`, locking the result in `pixi.lock`.

### Environments

| Environment | Features | Use case |
|---|---|---|
| `default` | `test` | CLI / pipeline / testing |
| `web` | `web`, `test` | Web server + testing |
| `cuda` | `cuda`, `test` | GPU off-target search |

Run a command in a specific environment:
```bash
pixi run --environment web beguided-web
```

### CUDA status

- `pytorch-cuda` is **not** included in the default environment.
- Add it with `pixi install --features cuda` for GPU-accelerated off-target search.
- PyTorch with CUDA reports `torch.cuda.is_available() == True` and CUDA 12.4.

## Code Style Guidelines

### Python Style
- Follow **PEP 8**
- 4 spaces for indentation
- 100 character line limit
- Two blank lines between top-level functions
- One blank line between methods in a class

### Naming Conventions
- Constants: `UPPER_CASE`
- Functions/variables: `snake_case`
- Preset keys: string names (e.g. `'ABE8e'`)

### Imports
- Group: stdlib, third-party, local
- Use explicit imports (avoid `from module import *`)

### Error Handling
- Validate CLI arguments early
- Use `sys.exit(1)` for fatal errors
- Print errors to stderr
- Raise `ValueError` for invalid business logic

## Testing Approach

### Unit Tests (pytest)

```bash
# Full suite
pixi run test

# Specific module
python -m pytest tests/test_pipeline/ -v

# Single test
python -m pytest tests/test_editing/test_presets.py::TestPAMToRegex -v
```

### End-to-End Testing
1. Use sample files to verify output
2. Test edge cases (chromosome not found, no valid PAMs)
3. Verify both A→G and T→C conversions
4. Test with different base editor profiles

### Test Fixtures
- `tests/fixtures/sample_targets.tsv`: Small test input
- `tests/fixtures/sample_genome.fasta`: Test genome
- `tests/fixtures/mini_genome.npz`: BEguidedOFF-encoded mini genome for off-target tests
- `configs/base_editor_profiles.tsv`: Bundled base editor profile configurations (20-nt guides)

## Important Implementation Details

### Position Conventions
- Input positions are **1-based**
- Internal calculations use **0-based** indexing
- Editing positions in output are **1-based from 5' end**

### BEguidedOFF Integration

- BEguidedOFF is a **separate sibling repository** (`../BEguidedOFF` from BEguided).
- Installed automatically into the `cuda` Pixi env via a pypi-dependencies path source (`beguidedoff = { path = "../BEguidedOFF", editable = true }` under `[tool.pixi.feature.cuda]`); the `default` and `web` envs do **not** include it.
- Outside the `cuda` env: editable `pip install` or add to `PYTHONPATH`.
- The genome must be **pre-encoded** to `.npz` format once via `python -m beguidedoff.cli encode-genome`.
- Pipeline Step 4 calls `BEguidedOFF.search()` via Python API (no subprocess).
- `n_top_offtargets=1` — only the single worst off-target per guide is returned.
- Scoring is **CCTOP** (lower = more concerning).
- On-target exclusion uses the 1-based leftmost guide position for ``+`` strand guides and the 1-based rightmost guide position for ``-`` strand guides, computed from the target coordinates.

### Preload Caches (Web)

The web server mirrors the same two caches used by the pipeline:

- `beguided/web/genome_cache.py` — parses the genome FASTA once at startup; exposes `/genome/status`; `os._exit(1)` on failure (fail-fast).
- `beguided/web/gtf_cache.py` — parses the GTF annotation once at startup; exposes `/gtf/status`; same fail-fast policy.

Both back the `run_pipeline(..., genome=..., interval_trees=...)` kwargs in `beguided.pipeline`. When `preload_genome=False` or `preload_gtf=False`, the pipeline falls back to parsing the file per job.

### Config Keys (offtargets block)

| Key | Default | Description |
|-----|---------|-------------|
| `beguidedoff_genome` | — | Path to `.npz` genome |
| `mismatch_tolerance` | 5 | Max mismatches in guide |
| `count_threshold_score` | 15.0 | CCTOP max for counting |
| `remove_guides_score` | 5.0 | CCTOP min to keep guide |
| `remove_guides_n_offtargets` | null | Max count to keep guide |

### Pipeline Output Layout

```
<base_dir>/
    step1_guidefinder.tsv
    step2_bystander.tsv
    step3_filtered.tsv
    step4_offtargets.tsv      (if search_offtargets)
    beguided_results.tsv
```

### Pipeline Steps

1. **guidefinder** — Generate guide RNAs
2. **bystander** — Analyze bystander effects
3. **filter** — Filter guides
4. **BEguidedOFF** — Off-target search (Python API call)
5. *(integrated)* — BEguidedOFF returns scored results
6. **merge** — Aggregate BEguidedOFF output and merge into final TSV

### Output Column Names (final TSV)

- `worst_offtarget_score` — lowest CCTOP score (most dangerous; `inf` if none)
- `offtarget_count` — number of off-targets ≤ `count_threshold_score`
- `worst_offtarget_chr`, `worst_offtarget_pos`, `worst_offtarget_sequence`
- `trinuc_edit_score` — Marquart et al. 2021 normalised ABE8e proportion for the local `X-A-Y` trinucleotide context (0..1, `NaN` for CBE guides). Bundled table at `beguided/editing/data/edit_context_norm_proportion.tsv`; override with `input.edit_context_table`. Filterable via `trinuc_edit_score_min`.

### CCTOP Score Semantics
- **Lower = more concerning** (mismatches near PAM)
- **Higher = less concerning** (mismatches far from PAM)
- Perfect match = 0.0
- No off-targets = `inf` (always passes `worst_offtarget_score >= X` filter)

### Web Interface

- FastAPI-based, started via `beguided-web`
- Config via CLI flags (override env variables):
  - `--base-editor-profiles-file`
  - `--genome-fasta`
  - `--gtf-path`
  - `--beguidedoff-genome`
- Or via environment variables:
  - `BEGUIDED_BASE_EDITOR_PROFILES_FILE`
  - `BEGUIDED_GENOME_FASTA`
  - `BEGUIDED_GTF_PATH`
  - `BEGUIDED_BEGUIDEDOFF_GENOME`
- `BEGUIDED_PRELOAD_GENOME` / `--preload-genome` / `--no-preload-genome` — preload the genome FASTA at startup (default: enabled).
- `BEGUIDED_PRELOAD_GTF` / `--preload-gtf` / `--no-preload-gtf` — preload the GTF annotation at startup (default: enabled).
- `BEGUIDED_QUEUE_MAX`, `BEGUIDED_COOLDOWN_SECONDS`, `BEGUIDED_MAX_UPLOAD_MB` — scheduler caps (env-only; CLI overrides also exist).
- `BEGUIDED_ALLOWED_HOSTS` / `BEGUIDED_TRUSTED_PROXY_IPS` / `BEGUIDED_JOB_RETENTION_HOURS` / `BEGUIDED_LOG_IP_SALT` / `BEGUIDED_ADMIN_TOKEN` — compliance knobs (see `compliance_safety.md`).
- `BEGUIDED_LEGAL_CONFIG` / `--legal-config` — path to the operator's `web-legal.yaml` with institution / DPO / DPA. Default `/etc/beguided/legal.yaml`. See `configs/web-legal.yaml` (template, committed empty).

### Sequence input mode

The web UI has two input modes, switched at the top of the page:

- **Spreadsheet** (default): upload a TSV/CSV/XLSX with columns
  `chromosome, position, target_base, desired_base`.
- **Sequence**: paste raw sequence text. Each sequence must contain
  exactly one edit bracket, e.g. `ATAGAGT...[A>G]...CCTATCC` to mark an
  A→G edit.

Accepted bracket directions: `[A>G]`, `[T>C]`, `[C>T]`, `[G>A]`
(the four directions actually supported by BEguided base editor profiles including
sense + antisense complement). Anything else (including IUPAC codes
inside the bracket) is rejected.

Format details:

- A single sequence with no FASTA header gets the chromosome name `custom`.
- With a FASTA header (anything after `>` up to the first whitespace), the
  header token is used as the chromosome name.
- The bracket is **replaced** by the base on the left side of `>` (the
  target_base) at the edit position; the text originally after `]` shifts
  one position later in the resulting sequence. For example
  `ATAG[T>C]GAC` becomes `ATAGTGAC` with edit position 5 (the `T` from
  the bracket).
- Edit position is 1-based and points at the base represented by the
  bracket.
- Per-sequence length cap: 1000 nt (`MAX_SEQUENCE_LENGTH` in
  `beguided/web/sequence_input.py`).
- Per-input cap: `cfg.max_input_targets` sequences.

In sequence mode:

- The configured reference genome and GTF annotation are not used.
- Bystander analysis is silently skipped (no GTF for custom sequences).
  This is implemented with the explicit `skip_bystander` pipeline flag so
  that step 3's filter does NOT apply a `mode="no_bystanders"` cut on
  `bystander_positions` either — doing so would drop every guide on
  A-rich sequences.
- BEguidedOFF off-target search is disabled (requires a pre-encoded `.npz`).

Implementation:

- Parser: `beguided/web/sequence_input.py`.
- Endpoint: `POST /sequence/parse` returns the same `{raw_b64, preview, ...}`
  shape as `POST /upload` plus a `custom_genome_b64` payload.
- `POST /run` accepts `mode=spreadsheet|sequence` and (in sequence mode)
  `custom_genome_b64`.
