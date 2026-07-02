# Agent Guidelines for BEguided

This document provides guidelines for AI agents working on the BEguided codebase.

## Project Overview

BEguided is a Python package for generating Cas9 base editor guide RNA sequences. It is managed with Pixi (single `pyproject.toml` manifest) and uses BEguidedOFF (standalone sibling repo) for GPU-accelerated off-target search. All functionality lives in the `beguided/` package.

## Architecture

```
beguided/                  # Main package
├── cli.py                   # Single CLI entry with argparse subcommands
├── pipeline.py              # Pipeline orchestrator (function calls, no subprocess)
├── models.py                # Shared dataclasses (PresetConfig, GuideResult, etc.)
├── genome/                  # Shared genome I/O and sequence operations
│   ├── io.py                # FASTA loading
│   ├── sequence.py          # Reverse complement, mutation helpers
│   ├── coordinates.py       # Coordinate math (0-based vs 1-based)
│   └── context.py           # Genomic region reconstruction
├── editing/                 # Core base-editing logic
│   ├── presets.py           # Preset loading, PAM regex, validation
│   ├── search.py            # PAM regex search
│   └── mutations.py         # Editing simulation, bystander detection
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
    ├── config.py            # Web config and preset groups
    ├── sse_logger.py        # SSE log streaming handler
    ├── tasks.py             # Async pipeline task wrapper
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

# Run guidefinder
python -m beguided.cli guidefinder \
    --input sample_targets.tsv \
    --output results.tsv \
    --genome sample_genome.fasta \
    --presets ABE8e
```

### Single CLI with Subcommands

```
beguided guidefinder     --input ... --genome ... --output ...
beguided bystander       --input ... --genome ... --gtf-annotation ... --output ...
beguided filter          --input ... --output ... --mode no_bystanders
beguided run-pipeline    --config config.yaml [--from-step 3]
beguided web             [--host 0.0.0.0] [--port 8000] [--reload] \
                         [--presets-file ...] [--genome-fasta ...] \
                         [--gtf-path ...] [--beguidedoff-genome ...] \
                         [--preload-genome | --no-preload-genome] \
                         [--preload-gtf | --no-preload-gtf] \
                         [--genome-label ...]
                         # Requires web extra: pixi run --environment web beguided web
```

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
pixi run --environment web beguided web
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
4. Test with different presets

### Test Fixtures
- `tests/fixtures/sample_targets.tsv`: Small test input
- `tests/fixtures/sample_genome.fasta`: Test genome
- `tests/fixtures/mini_genome.npz`: BEguidedOFF-encoded mini genome for off-target tests
- `be_presets.tsv`: Base editor preset configurations

## Important Implementation Details

### Position Conventions
- Input positions are **1-based**
- Internal calculations use **0-based** indexing
- Editing positions in output are **1-based from 5' end**

### BEguidedOFF Integration

- BEguidedOFF is a **separate sibling repository** (`../BEguidedOFF` from BEguided).
- Install: editable pip install or add to `PYTHONPATH`.
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

### CCTOP Score Semantics
- **Lower = more concerning** (mismatches near PAM)
- **Higher = less concerning** (mismatches far from PAM)
- Perfect match = 0.0
- No off-targets = `inf` (always passes `worst_offtarget_score >= X` filter)

### Web Interface

- FastAPI-based, started via `beguided web`
- Config via CLI flags (override env variables):
  - `--presets-file`
  - `--genome-fasta`
  - `--gtf-path`
  - `--beguidedoff-genome`
- Or via environment variables:
  - `BEGUIDED_PRESETS_FILE`
  - `BEGUIDED_GENOME_FASTA`
  - `BEGUIDED_GTF_PATH`
  - `BEGUIDED_BEGUIDEDOFF_GENOME`
- `BEGUIDED_PRELOAD_GENOME` / `--preload-genome` / `--no-preload-genome` — preload the genome FASTA at startup (default: enabled).
- `BEGUIDED_PRELOAD_GTF` / `--preload-gtf` / `--no-preload-gtf` — preload the GTF annotation at startup (default: enabled).
