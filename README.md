# BEguided

**Cas9 base editor guide RNA design.** Generates candidate guides for multiple base editors and multiple target mutations, predicts bystander editing, filters by biological criteria, and scores off-target activity.

Two entry points:

- **`beguided`** — CLI pipeline driven by a YAML config. **Primary entry point** — this README focuses on it.
- **`beguided-web`** — FastAPI server with a browser UI. For interactive use; see [`docs/web-deployment.md`](docs/web-deployment.md).

## Install

Pick one of the following.

### CLI — pip

```bash
git clone <repo-url> BEguided
cd BEguided
pip install -e .
beguided --help
```

### CLI — Pixi

```bash
git clone <repo-url> BEguided
cd BEguided
pixi install
```

### CLI + CUDA off-target search

The CUDA path uses the separate BEguidedOFF package for GPU-accelerated off-target search.

```bash
git clone <repo-url> BEguided
cd BEguided

# Option A — Pixi with sibling BEguidedOFF at ../BEguidedOFF
pixi install --features cuda

# Option B — pip + manual BEguidedOFF
pip install -e .
pip install "git+ssh://git@gitlab.mn.tu-dresden.de/buchholzlab/beguidedoff.git"
```

See [`docs/installation.md`](docs/installation.md#beguidedoff-integration) for the full CUDA setup, including encoding the genome once to `.npz`.

### Web server (optional, interactive use only)

```bash
pip install -e ".[web]"        # add [test] for the test suite
beguided-web                   # → http://localhost:8000
# or via Pixi:
pixi install --features web
```

See [`docs/web-deployment.md`](docs/web-deployment.md) for public deployment (Caddy, TLS, GDPR).

### Docker

The repo ships four image targets: `beguided:cli`, `beguided:cli-cuda`, `beguided:web`, `beguided:web-cuda`. See [`docs/installation.md#docker`](docs/installation.md#docker).

## Get the reference genome

BEguided needs a reference FASTA and a matching GTF annotation. A helper
script downloads both from Ensembl into `./genome/`:

```bash
scripts/download_genome.sh                       # default: ./genome, release 116
scripts/download_genome.sh --target-dir /data/refs/GRCh38.116
scripts/download_genome.sh --release 115
scripts/download_genome.sh --force               # re-download even if present
```

After the script finishes you will have:

```
./genome/Homo_sapiens.GRCh38.116.fa
./genome/Homo_sapiens.GRCh38.116.gtf
```

Point your CLI config at these paths (see `configs/cli.yaml`) or set
`BEGUIDED_GENOME_FASTA` / `BEGUIDED_GTF_PATH` for the web server. Re-running
the script is safe: existing files are kept unless `--force` is passed.

## Quick start — CLI

```bash
# 1. Configure (copy the template, edit paths)
cp configs/cli.yaml my_run.yaml

# 2. Encode the genome once for off-target search (optional)
python -m beguidedoff.cli encode-genome \
    --fasta genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --output genome_otbe.npz --include-rc --include-n-mask

# 3. Run the pipeline
beguided --config my_run.yaml
# → results/step1_guidefinder.tsv, ..., beguided_results.tsv
```

See [`docs/cli-usage.md`](docs/cli-usage.md) for the full config reference, [`docs/input-output.md`](docs/input-output.md) for file formats, and [`docs/pipeline-steps.md`](docs/pipeline-steps.md) for what each of the six steps does.

## Quick start — web UI

```bash
# Start the server
beguided-web
# → http://localhost:8000
```

Open the URL, paste sequences (or upload a targets TSV), pick base editor profiles, hit **Run**. Results stream in via Server-Sent Events and can be downloaded as TSV or XLSX with dynamic filters.

For public deployment (Caddy, TLS, GDPR), see [`docs/web-deployment.md`](docs/web-deployment.md).

## Input / output

- **Input**: TSV with columns `chromosome`, `position`, `target_base`, `desired_base` (extras preserved).
- **Output**: `beguided_results.tsv` with every guide, its bystander analysis, on-target score, off-target summary, and your passthrough columns.

Full schema: [`docs/input-output.md`](docs/input-output.md).

## Configuration files

```
configs/
├── cli.yaml                   ← annotated template for the CLI pipeline
├── base_editor_profiles.tsv   ← bundled ABE / CBE profiles (20-nt guides)
└── web-legal.yaml             ← empty template for the legal pages (operator fills in their own file)
```

The operator's real legal YAML is gitignored and pointed at via `BEGUIDED_LEGAL_CONFIG` — see [`docs/web-deployment.md`](docs/web-deployment.md#legal-yaml).

## Documentation

| | |
|---|---|
| [Installation](docs/installation.md) | pip, pixi, docker; platform support; BEguidedOFF |
| [CLI usage](docs/cli-usage.md) | Running the pipeline; config keys; web server |
| [Input / output](docs/input-output.md) | Target file format, output TSV columns |
| [Pipeline steps](docs/pipeline-steps.md) | What each of the six steps does |
| [Scoring](docs/scoring.md) | `guide_score`, `gc_content`, `trinuc_edit_score`, CCTOP |
| [Web deployment](docs/web-deployment.md) | Caddy, Docker Compose, GDPR checklist, legal YAML |
| [AGENTS.md](docs/AGENTS.md) | Notes for AI coding agents working on this repo |

## Testing

```bash
pixi run test                    # full suite
python -m pytest tests/test_pipeline_steps/ -v   # one module
python -m pytest tests/test_pipeline_steps/test_filtering.py::test_x -v   # one test
```

## License

MIT. See [LICENSE](LICENSE).
