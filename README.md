# BEguided

Cas9 Base Editor Guide RNA Generator

BEguided is a Python package for designing guide RNAs for Cas9 base editors. It generates candidate guides around target mutations, analyses bystander editing effects on protein coding sequences, filters guides by biological criteria, and scores off-target activity using BEguidedOFF with CCTOP scoring.

## Features

- **Guide RNA generation** – Search around target positions for valid PAM sites, select guide RNAs that place the target base in the high-efficiency editing window, and insert the desired edit.
- **On-target scoring** – CRISPRater score (adapted from CCTOP standalone) is computed for every 20-nt guide to estimate predicted on-target activity.
- **Preset system** – Support multiple base editors with custom PAM patterns, guide lengths and editing windows (loaded from a TSV file).
- **Bystander scanning** – Predict amino-acid consequences of unintended edits within the sequence by integrating GTF annotation, codon checks, and positional analysis.
- **Guide filtering** – Filter by bystander effect (no bystanders, silent only, or accept all).
- **Off-target scoring** – Search for off-target sites using BEguidedOFF (GPU-accelerated or CPU fallback) and score them with the CCTOP penalty model.
- **CCTOP scoring** – Lower CCTOP scores indicate more concerning off-targets (mismatches near the PAM). The score is unbiased across PAM families.
- **Column passthrough** – Any additional columns in your target input file (e.g. `mutation_occurance`, `genomic_mutation_id`) are preserved through all pipeline steps to the final output.
- **Web interface** – A browser-based GUI for uploading targets, selecting presets, streaming pipeline logs, and downloading results. Runs locally or on a server.

## Installation

Requires Python >= 3.10 and [Pixi](https://pixi.sh) (recommended).

```bash
# Clone the repository
cd BEguided

# Install dependencies (including BEguidedOFF from GitLab)
pixi install

# Install BEguided in editable mode
pip install -e .
```

The `beguided` command will then be available:

```bash
beguided --help
```

Or run directly via Pixi:

```bash
pixi run beguided --help
```

### Dependencies

- numpy, biopython, pandas, intervaltree, pyyaml
- fastapi, uvicorn, jinja2 (web interface)
- torch, numba, pysam (BEguidedOFF off-target search)

For GPU-accelerated off-target search, additionally install:

```bash
pixi add pytorch-cuda
```

## Input format

The target input file is a tab-separated file with these required columns:

| Column         | Description                                              |
|----------------|----------------------------------------------------------|
| `chromosome`   | Chromosome name (e.g. `chr1`, `12`, `X`)                 |
| `position`     | 1-based genomic coordinate of the target mutation        |
| `target_base`  | Base currently present at the target (`A`, `T`, `C`, `G`)|
| `desired_base` | Desired base after editing (`A`, `T`, `C`, `G`)          |

You may add any number of extra columns (e.g. `gene`, `mutation_occurance`, `genomic_mutation_id`); they will be carried through to every output file.

Example (`tests/fixtures/sample_targets.tsv`):

```tsv
chromosome	position	target_base	desired_base
chr1	50	T	C
chr1	75	A	G
chr2	30	T	C
chr2	80	A	G
```

## Usage

BEguided provides a single CLI entrypoint with subcommands. You can run individual steps for development/debugging or execute the full pipeline in one go.

### 1. Generate guide RNAs

```bash
beguided guidefinder \
    --input tests/fixtures/sample_targets.tsv \
    --output results_guides.tsv \
    --genome tests/fixtures/sample_genome.fasta \
    --presets all \
    --presets-file be_presets.tsv
```

Preset options:
- `--presets all` — run every available preset
- `--presets all ABEs` — run only A→G / T→C presets (e.g. ABE8e, ABE9, saABE8e)
- `--presets all CBEs` — run only C→T / G→A presets (e.g. CBE6d, TadCBEd)
- `--presets ABE8e,ABE9` — comma-separated list of specific presets

### 2. Analyse bystander effects

```bash
beguided bystander \
    --input results_guides.tsv \
    --output results_bystander.tsv \
    --genome genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf-annotation genome/Homo_sapiens.GRCh38.112.gtf
```

### 3. Filter guides

```bash
beguided filter \
    --input results_bystander.tsv \
    --output results_filtered.tsv \
    --mode silent
```

Filter modes:
- `all` – keep everything (default)
- `no_bystanders` – only guides with no bystander edits (skips bystander analysis when used in the pipeline)
- `silent` – guides that are silent or have no bystanders

### Run the complete pipeline

All steps are executed in one command via a YAML config file:

```bash
beguided run-pipeline --config config.yaml
```

You can also run a subset of steps, e.g. steps 3-6:

```bash
beguided run-pipeline --config config.yaml --from-step 3
```

See `template.yaml` for a commented configuration template.

### Off-target search setup

Before running the pipeline with off-targets, encode the reference genome once:

```bash
python -m beguidedoff.cli encode-genome \
    --fasta genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --output genome_otbe.npz \
    --include-rc \
    --include-n-mask
```

Then reference `genome_otbe.npz` in your pipeline config under `offtargets.beguidedoff_genome`.

### Web interface

BEguided includes a FastAPI-based web interface that can be run locally or served on a remote machine.

**Start the server:**

```bash
# Default: http://localhost:8000
beguided web

# Custom address
beguided web --host 0.0.0.0 --port 8080

# Development with auto-reload
beguided web --reload

# Override file paths (recommended for production)
beguided web \
  --presets-file be_presets.tsv \
  --genome-fasta genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gtf-path genome/Homo_sapiens.GRCh38.112.gtf \
  --beguidedoff-genome /path/to/hg38_otbe.npz
```

The web server looks for four paths.  They can be set via **CLI flags** (shown above) or **environment variables**:

| CLI flag | Environment variable | Default |
|---|---|---|
| `--presets-file` | `BEGUIDED_PRESETS_FILE` | `be_presets.tsv` (repo root) |
| `--genome-fasta` | `BEGUIDED_GENOME_FASTA` | `genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| `--gtf-path` | `BEGUIDED_GTF_PATH` | `genome/Homo_sapiens.GRCh38.112.gtf` |
| `--beguidedoff-genome` | `BEGUIDED_BEGUIDEDOFF_GENOME` | *(none – off-target search is disabled if missing)* |

CLI flags take priority over environment variables.

**Interface overview:**

- **Top bar** – Upload target files (`.tsv`, `.csv`, `.xlsx`, up to 1000 rows), download a template, select a base editor preset, and choose the filter mode.
- **Input tab** – Preview the uploaded target table.
- **Logging tab** – Real-time pipeline logs streamed via Server-Sent Events.
- **Results tab** – Interactive final results table with a download button for `beguided_results.tsv`.

**Configuration (environment variables):**

| Variable | Default | Description |
|----------|---------|-------------|
| `BEGUIDED_PRESETS_FILE` | `be_presets.tsv` | Preset definitions |
| `BEGUIDED_GENOME_FASTA` | `genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | Genome FASTA |
| `BEGUIDED_GTF_PATH` | `genome/Homo_sapiens.GRCh38.112.gtf` | Gene annotation |
| `BEGUIDED_BEGUIDEDOFF_GENOME` | `tests/fixtures/mini_genome.npz` | BEguidedOFF encoded genome |

## Pipeline configuration

The pipeline is configured via a YAML file. Key sections:

```yaml
input:
  targets: "path/to/targets.tsv"
  genome: "path/to/genome.fa"
  presets: "be_presets.tsv"
  gff_file: "path/to/annotation.gtf"

output_prefix: "./results"
use_presets: "all"
keep_bystanders: "all"
guide_score_threshold: 0.5
search_offtargets: true

offtargets:
  beguidedoff_genome: "path/to/genome_otbe.npz"
  mismatch_tolerance: 5
  count_threshold_score: 15.0     # max CCTOP to count
  remove_guides_score: 5.0        # min CCTOP to keep guide
  remove_guides_n_offtargets: null  # max count to keep guide
```

### Config key reference

| Key | Default | Description |
|-----|---------|-------------|
| `guide_score_threshold` | 0.5 | CRISPRater **minimum** score to keep a guide. Guides with a score below this value are removed in step 3. Range: 0–1 (higher = better). Set to 0 to disable. |
| `offtargets.beguidedoff_genome` | — | Path to BEguidedOFF `.npz` genome |
| `offtargets.mismatch_tolerance` | 5 | Max mismatches in 20-nt guide |
| `offtargets.count_threshold_score` | 15.0 | CCTOP **maximum** for counting off-targets |
| `offtargets.remove_guides_score` | 5.0 | CCTOP **minimum** to keep a guide (guides with worse off-targets are removed) |
| `offtargets.remove_guides_n_offtargets` | null | Max off-target count to keep a guide |

**CCTOP semantics:** Lower score = more concerning (mismatches near PAM). A guide with no off-targets gets `worst_offtarget_score = inf`.

## Output layout

When running the full pipeline, files are written directly into the `output_prefix` directory:

```
output_prefix/
├── step1_guidefinder.tsv     (step 1 – guide RNAs)
├── step2_bystander.tsv       (step 2 – bystander analysis)
├── step3_filtered.tsv        (step 3 – filtered results)
├── step4_offtarget.tsv       (step 4 – BEguidedOFF off-target results)
└── beguided_results.tsv    (step 6 – final merged output)
```

The final `beguided_results.tsv` contains all guide results merged with off-target summary columns (`worst_offtarget_score`, `offtarget_count`, `worst_offtarget_chr`, `worst_offtarget_pos`, `worst_offtarget_sequence`). Every row also carries a `guide_score` column — the CRISPRater on-target activity score computed in step 1 (requires a 20-nt guide; non-20-nt presets are skipped with a warning). Any extra columns from your input target file are preserved in the final output.

## Base editor presets

Base editor configurations are stored in a TSV file (default: `be_presets.tsv`):

| Column               | Description                                  |
|----------------------|----------------------------------------------|
| `name`               | Preset name (e.g. `ABE8e`)                   |
| `editing`            | Conversion type (e.g. `A-to-G`)              |
| `general_edit_start` | Start of the general editing window          |
| `general_edit_end`   | End of the general editing window            |
| `high_edit_start`    | Start of the high-efficiency window          |
| `high_edit_end`      | End of the high-efficiency window            |
| `guide_length`       | Length of the guide RNA (e.g. 20)            |
| `pam`                | PAM pattern using IUPAC codes (e.g. `NGG`)   |
| `cas`                | Cas nuclease name (e.g. `spCas9`)            |
| `source`             | Publication or source                        |

Both sense and antisense variants are automatically generated from each preset definition (e.g. `A-to-G` → sense; `T-to-C` → antisense).

## BEguidedOFF integration

BEguidedOFF is an independent Python package (developed alongside BEguided) that searches the genome for off-target sites with GPU acceleration. It is installed automatically via Pixi from its GitLab repository.

- **Repos**: BEguided depends on BEguidedOFF via a Git SSH source in `pyproject.toml`.
- **Genome encoding**: One-time `encode-genome` command produces a `.npz` file used by the pipeline.
- **Top-N results**: BEguidedOFF returns only the single worst off-target per guide (`n_top_offtargets=1`) to save disk space.

## Testing

```bash
# Run the full test suite
pixi run test

# Or directly with pytest
python -m pytest tests/ -v

# Run a specific module
python -m pytest tests/test_guidefinder/ -v

# Run a single test
python -m pytest tests/test_editing/test_presets.py::TestPAMToRegex -v
```

All core logic lives in the `beguided.*` namespace, so it can also be imported and used directly from Python scripts or web frameworks (e.g. FastAPI or Flask) without CLI overhead.
