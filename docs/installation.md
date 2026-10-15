# Installation

BEguided requires Python ≥ 3.10. Install via **pip** (any system), **Pixi** (recommended when you have access to BEguidedOFF), or **Docker** (no local Python needed).

## pip (any system)

| What you need | Command |
|---|---|
| Core CLI & pipeline only | `pip install -e .` |
| + Web server | `pip install -e ".[web]"` |
| + Tests | `pip install -e ".[test]"` |
| Everything | `pip install -e ".[web,test]"` |

```bash
git clone https://github.com/ltschmitt/BEguided.git BEguided
cd BEguided
pip install -e .
```

Core functionality works out of the box. Off-target search requires the separate `BEguidedOFF` package — see [BEguidedOFF integration](#beguidedoff-integration) below.

## Pixi

```bash
git clone https://github.com/ltschmitt/BEguided.git BEguided
cd BEguided

# Core environment (CLI + tests)
pixi install

# + Web server
pixi install --features web
# or use the dedicated web environment:
# pixi run --environment web beguided-web

# + CUDA support for GPU-accelerated off-target search
pixi install --features cuda
```

The `cuda` feature installs `BEguidedOFF` automatically from the sibling repo at `../BEguidedOFF` (editable). If you don't have that sibling repo, see the manual install instructions below.

## Platform support

| Platform | Core features | Off-target search | GPU acceleration |
|---|---|---|---|
| Linux | Yes | Yes (manual install) | Yes |
| macOS | Yes | Yes (manual install) | CPU only |
| Windows | Yes | Yes (manual install) | CPU only |

## BEguidedOFF integration

Off-target search is optional. When `BEguidedOFF` is missing, the pipeline skips step 4 silently and the off-target columns in `beguided_results.tsv` are `NaN`.

### Via Pixi (recommended)

```bash
pixi install --features cuda
```

This installs the `cuda` pixi environment with `BEguidedOFF` editable from `../BEguidedOFF`. Then provide the path to an encoded genome `.npz` via `--beguidedoff-genome` or `BEGUIDED_BEGUIDEDOFF_GENOME`.

### Manual install

```bash
pip install "git+https://github.com/ltschmitt/BEguidedOFF.git"
```

### Encoding a genome

Before running the pipeline with off-targets, encode the reference genome once:

```bash
python -m beguidedoff.cli encode-genome \
    --fasta genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --output genome_otbe.npz \
    --include-rc \
    --include-n-mask
```

Reference the resulting `.npz` in your pipeline config under `offtargets.beguidedoff_genome`, or pass `--beguidedoff-genome /path/to/genome_otbe.npz` on the CLI.

## Docker

The repo ships a multi-stage `Dockerfile` that produces four images:

| Target | Tag | CUDA | Web | Use case |
|---|---|---|---|---|
| `base` | `beguided:cli` | No | No | CLI / pipeline |
| `web` (default) | `beguided:web` | No | Yes | Web interface |
| `cli-cuda` | `beguided:cli-cuda` | Yes | No | CLI + GPU off-target search |
| `web-cuda` | `beguided:web-cuda` | Yes | Yes | Web + GPU off-target search |

### Build

```bash
docker build --target base     -t beguided:cli     .
docker build --target web      -t beguided:web     .
docker build --target cli-cuda -t beguided:cli-cuda .
docker build --target web-cuda -t beguided:web-cuda .
```

### Run CPU CLI

```bash
docker run -v $(pwd)/data:/app/data beguided:cli \
    --config /app/data/config.yaml
```

### Run CPU web server

```bash
docker run -p 8000:8000 \
  -v $(pwd)/configs/base_editor_profiles.tsv:/app/configs/base_editor_profiles.tsv:ro \
  -v $(pwd)/genome:/app/genome:ro \
  beguided:web
```

### Run GPU web server (requires NVIDIA Container Toolkit)

```bash
docker run --gpus all -p 8000:8000 \
  -v $(pwd)/configs/base_editor_profiles.tsv:/app/configs/base_editor_profiles.tsv:ro \
  -v $(pwd)/genome:/app/genome:ro \
  beguided:web-cuda
```

## Verifying the install

```bash
beguided --help
beguided-web --help
```

If both commands print their usage strings, the install is complete.
