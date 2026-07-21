# BEguided

Cas9 Base Editor Guide RNA Generator and Evaluator

BEguided is a Python package for designing guide RNAs for Cas9 base editors. It generates candidate guides around target mutations, analyses bystander editing effects on protein coding sequences, filters guides by biological criteria, and scores off-target activity using BEguidedOFF with CCTOP scoring.

## Features

- **Guide RNA generation** – Search around target positions for valid PAM sites, select guide RNAs that place the target base in the high-efficiency editing window, and insert the desired edit.
- **On-target scoring** – CRISPRater score (adapted from CCTOP standalone) is computed for every 20-nt guide to estimate predicted on-target activity.
- **Preset system** – Support multiple base editors with custom PAM patterns, guide lengths and editing windows (loaded from a TSV file).
- **Bystander scanning** – Predict amino-acid consequences of unintended edits within the sequence by integrating GTF annotation, codon checks, and positional analysis.
- **Guide filtering** – Filter by bystander effect (no bystanders, silent only, or accept all), by CRISPRater on-target score, by GC content (40–60% "good" or 30–70% "acceptable"), by homopolymer runs (≥4-mer), and by ABE trinucleotide edit-context score (Marquart et al. 2021, ABE8e-derived).
- **Edit-context scoring** – Every ABE guide gets a `trinuc_edit_score` from a bundled, normalised Marquart et al. 2021 (ABE8e) `X-A-Y` trinucleotide context table (1.0 = `CAC`, 0.02 = `AAA`). CBE guides are `NaN`. The `-` strand is reverse-complemented so the lookup uses the mRNA-sense context.
- **Off-target scoring** – Search for off-target sites using BEguidedOFF (GPU-accelerated or CPU fallback) and score them with the CCTOP penalty model.
- **CCTOP scoring** – Lower CCTOP scores indicate more concerning off-targets (mismatches near the PAM). The score is unbiased across PAM families.
- **Column passthrough** – Any additional columns in your target input file (e.g. `mutation_occurance`, `genomic_mutation_id`) are preserved through all pipeline steps to the final output.
- **Web interface** – A browser-based GUI for uploading targets, selecting base editor profiles, streaming pipeline logs, and downloading results. Runs locally or on a server.

## Installation

Requires Python >= 3.10. BEguided can be installed via **pip** (universal) or **Pixi** (recommended for TU Dresden/internal users who have access to BEguidedOFF).

### pip (any system)

| What you need | Command |
|---|---|
| Core CLI & pipeline only | `pip install -e .` |
| + Web server | `pip install -e ".[web]"` |
| + Tests | `pip install -e ".[test]"` |
| Everything | `pip install -e ".[web,test]"` |

```bash
cd BEguided

# 1) Install core dependencies
pip install -r requirements.txt

# 2) Install BEguided
pip install -e .
```

Core functionality works out of the box. Off-target search requires the separate BEguidedOFF package (install manually if you have access).

### Pixi (TU Dresden / internal)

```bash
cd BEguided

# Core environment (CLI + tests)
pixi install

# + Web server
pixi install --features web
# or use the dedicated web environment:
# pixi run --environment web beguided-web

# + CUDA support for GPU-accelerated off-target search
pixi install --features cuda

# Install BEguided in editable mode
pip install -e .
```

### Platform support

| Platform | Core features | Off-target search (BEguidedOFF) | GPU acceleration |
|---|---|---|---|
| Linux | Yes | Yes (manual install) | Yes |
| macOS | Yes | Yes (manual install) | CPU only |
| Windows | Yes | Yes (manual install) | CPU only |

### Optional: Install BEguidedOFF (off-target search)

When using Pixi, BEguidedOFF is installed automatically — editable, from the sibling repo at `../BEguidedOFF` — into the `cuda` environment:

```bash
pixi install --features cuda
```

Then provide the path to an encoded genome `.npz` via `--beguidedoff-genome` or `BEGUIDED_BEGUIDEDOFF_GENOME`.

If you are not using Pixi (or the sibling repo is not at `../BEguidedOFF`), install it manually:

```bash
pip install "git+ssh://git@gitlab.mn.tu-dresden.de/buchholzlab/beguidedoff.git"
```

### Docker

A multi-stage `Dockerfile` is included. It produces four images:

| Target | Tag | CUDA | Web | Use case |
|---|---|---|---|---|
| `base` | `beguided:cli` | No | No | CLI / pipeline |
| `web` (default) | `beguided:web` | No | Yes | Web interface |
| `cli-cuda` | `beguided:cli-cuda` | Yes | No | CLI + GPU off-target search |
| `web-cuda` | `beguided:web-cuda` | Yes | Yes | Web + GPU off-target search |

**Build:**
```bash
docker build --target base     -t beguided:cli     .
docker build --target web      -t beguided:web     .
docker build --target cli-cuda -t beguided:cli-cuda .
docker build --target web-cuda -t beguided:web-cuda .
```

**Run CPU CLI:**
```bash
docker run -v $(pwd)/data:/app/data beguided:cli \
    --config /app/data/config.yaml
```

**Run CPU web server:**
```bash
docker run -p 8000:8000 \
  -v $(pwd)/base_editor_profiles.tsv:/app/base_editor_profiles.tsv \
  -v $(pwd)/genome:/app/genome \
  beguided:web
```

**Run GPU CLI (requires NVIDIA Container Toolkit):**
```bash
docker run --gpus all \
  -v $(pwd)/data:/app/data \
  beguided:cli-cuda --config /app/data/config.yaml
```

**Run GPU web server:**
```bash
docker run --gpus all -p 8000:8000 \
  -v $(pwd)/base_editor_profiles.tsv:/app/base_editor_profiles.tsv \
  -v $(pwd)/genome:/app/genome \
  beguided:web-cuda
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

BEguided provides two commands: `beguided` (run the pipeline) and `beguided-web` (start the web server).

### Run the complete pipeline

All steps (guide generation, bystander analysis, filtering, off-target search, merge) are executed via a YAML config file:

```bash
beguided --config config.yaml
```

You can also start from a later step, e.g. steps 3-5:

```bash
beguided --config config.yaml --from-step 3
```

Use `-v` for DEBUG logging, and `--no-passthrough-extras` to drop non-core columns from the input TSV.

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

> **Requires the `web` extra.** Install first: `pip install -e ".[web]"` or `pixi install --features web`.

BEguided includes a FastAPI-based web interface that can be run locally or served on a remote machine.

**Start the server:**

```bash
# Default: http://localhost:8000
beguided-web

# Custom address
beguided-web --host 0.0.0.0 --port 8080

# Development with auto-reload
beguided-web --reload

# Override file paths (recommended for production)
beguided-web \
  --base-editor-profiles-file base_editor_profiles.tsv \
  --genome-fasta genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gtf-path genome/Homo_sapiens.GRCh38.112.gtf \
  --beguidedoff-genome /path/to/hg38_otbe.npz
```

The web server looks for four paths.  They can be set via **CLI flags** (shown above) or **environment variables**:

| CLI flag | Environment variable | Default |
|---|---|---|
| `--base-editor-profiles-file` | `BEGUIDED_BASE_EDITOR_PROFILES_FILE` | `base_editor_profiles.tsv` (repo root) |
| `--genome-fasta` | `BEGUIDED_GENOME_FASTA` | `genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| `--gtf-path` | `BEGUIDED_GTF_PATH` | `genome/Homo_sapiens.GRCh38.112.gtf` |
| `--beguidedoff-genome` | `BEGUIDED_BEGUIDEDOFF_GENOME` | *(none – off-target search is disabled if missing)* |

CLI flags take priority over environment variables.

**Interface overview:**

- **Top bar** – Upload target files (`.tsv`, `.csv`, `.xlsx`, up to 1000 rows), download a template, select a base editor profile, and choose the filter mode.
- **Input tab** – Preview the uploaded target table.
- **Logging tab** – Real-time pipeline logs streamed via Server-Sent Events.
- **Results tab** – Interactive final results table with a download button for `beguided_results.tsv`.

**Configuration (environment variables):**

| Variable | Default | Description |
|----------|---------|-------------|
| `BEGUIDED_BASE_EDITOR_PROFILES_FILE` | `base_editor_profiles.tsv` | Base editor profile definitions |
| `BEGUIDED_GENOME_FASTA` | `genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | Genome FASTA |
| `BEGUIDED_GTF_PATH` | `genome/Homo_sapiens.GRCh38.112.gtf` | Gene annotation |
| `BEGUIDED_BEGUIDEDOFF_GENOME` | `tests/fixtures/mini_genome.npz` | BEguidedOFF encoded genome |

## Pipeline configuration

The pipeline is configured via a YAML file. Key sections:

```yaml
input:
  targets: "path/to/targets.tsv"
  genome: "path/to/genome.fa"
  base_editor_profiles: "base_editor_profiles.tsv"
  gff_file: "path/to/annotation.gtf"

output_prefix: "./results"
use_base_editors: "all"
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
| `gc_content_filter` | `off` | Filter guides by GC content of the 20-nt protospacer: `off` (no filter), `good` (keep 40–60%) or `acceptable` (keep 30–70%). Guides with a missing GC value are always kept. |
| `filter_homopolymers` | `false` | When `true`, drop guides that contain a run of ≥4 identical bases (AAAA, TTTT, CCCC, GGGG). Set to `false` to keep everything. |
| `trinuc_edit_score_min` | `null` | Minimum trinucleotide edit-context score (Marquart et al. 2021, ABE8e, normalised 0–1) to keep an ABE guide. The score reflects how likely the local `X-A-Y` context is to be edited (`CAC` = 1.0, `AAA` = 0.02). CBE rows are always kept (their `trinuc_edit_score` is `NaN` because the table is ABE-only). Set to `null` to disable. |
| `input.edit_context_table` | bundled | Optional custom trinucleotide context TSV overriding the bundled copy. Empty string (`""`) disables the `trinuc_edit_score` column entirely. |
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

The final `beguided_results.tsv` contains all guide results merged with off-target summary columns (`worst_offtarget_score`, `offtarget_count`, `worst_offtarget_chr`, `worst_offtarget_pos`, `worst_offtarget_sequence`). Every row also carries a `guide_score` column — the CRISPRater on-target activity score computed in step 1 (requires a 20-nt guide; non-20-nt base editor profiles are skipped with a warning), plus `gc_content` (whole-guide GC percent, 0–100), `homopolymer` (boolean: `true` if the guide contains any run of ≥4 identical bases) and `trinuc_edit_score` (the Marquart et al. 2021 normalised ABE8e proportion for the local `X-A-Y` context, 0–1; `NaN` for CBE guides). The three sequence-property columns are computed for every guide length. Any extra columns from your input target file are preserved in the final output.

> **Note on `trinuc_edit_score`:** The bundled table is derived from ABE8e spCas9 data (Marquart et al. 2021, https://doi.org/10.1038/s41467-021-25375-z) and is applied to **every ABE profile** as an approximation — actual editing efficiencies will differ between ABE8e, ABE9, ABE-Umax, etc. CBE profiles always produce `NaN` because the table is ABE-only. The `-` strand is handled by reverse-complementing the genomic 3-mer before lookup, so the value reflects the mRNA-sense context the paper measured.

## Base editor profiles

Base editor configurations are stored in a TSV file (default: `base_editor_profiles.tsv`):

| Column               | Description                                  |
|----------------------|----------------------------------------------|
| `name`               | Profile name (e.g. `ABE8e`)                  |
| `editing`            | Editor family (`ABE` or `CBE`)               |
| `general_edit_start` | Start of the general editing window          |
| `general_edit_end`   | End of the general editing window            |
| `high_edit_start`    | Start of the high-efficiency window          |
| `high_edit_end`      | End of the high-efficiency window            |
| `pam`                | PAM pattern using IUPAC codes (e.g. `NGG`)   |
| `cas`                | Cas nuclease name (e.g. `SpCas9`)            |
| `source`             | Publication or source                        |

All base editor profiles in the bundled file use a 20-nt guide. The guide RNA length is therefore fixed at 20 bp and is not stored in the TSV. Both sense and antisense variants are automatically generated from each profile definition (e.g. `ABE` → sense `A→G`, antisense `T→C`; `CBE` → sense `C→T`, antisense `G→A`).

> **Note:** The `guide_score` column (CRISPRater on-target activity score) is only computed for **20 bp** guides. Because every base editor profile in `base_editor_profiles.tsv` uses a 20-nt guide, the score is always populated. If you add a custom profile with a different guide length, the score will be skipped with a warning and the `guide_score` column will be empty.

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
