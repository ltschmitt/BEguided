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

Run from a checkout — the recommended setup for development and for the CUDA/web features:

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

The `cuda` feature installs `BEguidedOFF` automatically from its public GitHub repo and builds its native FM-index extension against conda's `sdsl-lite` — no sibling checkout needed. (macOS/Windows: see the manual install instructions below.)

## Pixi global install (no clone)

The repository doubles as a pixi **source package** (the `pixi-build` preview feature). Pixi ≥ 0.81 builds it straight from GitHub and exposes both entry points system-wide:

```bash
pixi global install --git https://github.com/ltschmitt/BEguided.git
beguided --help
```

Scope of the global install:

- Ships the **core CLI** (`beguided` and `beguided-web` entry points) with all runtime dependencies from conda-forge.
- No web-server framework. Add it to the same global environment with:

  ```bash
  pixi global install --git https://github.com/ltschmitt/BEguided.git beguided \
      --with fastapi --with uvicorn --with jinja2 --with python-multipart
  ```

- No off-target search. BEguidedOFF is a separate package with a native extension and build-time requirements (`sdsl-lite`, a C++ toolchain) that don't fit a one-line global install — use the workspace setup below instead.

## BEguided as a dependency in your own Pixi workspace

Both repos are pixi source packages, so a workspace of your own can consume them directly from GitHub — no BEguided checkout. In your `pyproject.toml`:

```toml
[tool.pixi.workspace]
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]
preview = ["pixi-build"]   # source packages require the pixi-build preview flag

[tool.pixi.dependencies]
beguided = { git = "https://github.com/ltschmitt/BEguided.git" }

# Off-target search (fast CPU mode):
beguidedoff = { git = "https://github.com/ltschmitt/BEguidedOFF.git" }
```

or via the CLI: `pixi add --git https://github.com/ltschmitt/BEguided.git beguided` (same for `beguidedoff`). Pixi builds each package on your machine (pure Python for BEguided; C++ FM-index extension against conda's `sdsl-lite` for BEguidedOFF) and resolves all runtime dependencies — `bioconda` is needed for `pysam`.

Notes:

- BEguidedOFF's FM-index binding is built for one CPython ABI; the recipe caps the build at CPython < 3.15 because bioconda does not ship py3.15 `pysam` builds yet. Bump that ceiling in BEguidedOFF's recipe when bioconda catches up.
- Fast mode then needs its FM-index exactly as in a checkout — see [BEguidedOFF integration](#beguidedoff-integration) for `encode-genome` / `build-index` and the `BEGUIDED_BEGUIDEDOFF_INDEX*` variables.

### GPU (extensive mode) in your own workspace

Add a CUDA-enabled PyTorch from conda-forge (the `pytorch` conda channel has stopped publishing Linux builds after the 2.5.1/py3.12 era, so it cannot be combined with current CPython builds):

```toml
[tool.pixi.workspace]
channels = ["conda-forge", "bioconda"]
# pixi's solver only sees the CUDA driver if it is declared here — the
# detected `__cuda` virtual package is not visible to the solve otherwise.
# This makes the workspace GPU-only, which is what you want in a dedicated
# GPU workspace (and impossible to combine with CPU-only use of the same
# workspace).
platforms = [{ platform = "linux-64", __cuda = "12" }]
preview = ["pixi-build"]

[tool.pixi.dependencies]
beguided = { git = "https://github.com/ltschmitt/BEguided.git" }
beguidedoff = { git = "https://github.com/ltschmitt/BEguidedOFF.git" }
pytorch-gpu = ">=2.5.1,<3"   # conda-forge CUDA build + CUDA runtime libraries
numba = ">=0.65.1,<0.66"     # BEguidedOFF's GPU kernels run through numba.cuda
cuda-nvvm = "*"              # NVVM + libdevice that numba.cuda JIT-compiles against
```

Runtime-verified on an RTX 3060 (driver CUDA 13.4): torch 2.14.1 (cuda129 build), `numba.cuda` JIT against conda-forge `cuda-nvvm` 12.9, and the `beguidedoff` FM-index binding — all on CPython 3.14, all packages pulled from Git. `pytorch-gpu` depends on the `__cuda` virtual package, so this workspace only runs on machines with an NVIDIA driver — expected for a GPU workspace, but don't merge these platform lines into a workspace that must also run CPU-only.

BEguided's own `cuda` feature (checkout route above) stays on the `pytorch` channel's pytorch 2.5.1/CUDA 12.4 stack on purpose: those are the only CUDA builds that do not depend on `__cuda`, so the shared workspace keeps working on GPU-less machines (CI, CPU servers). A fresh `pixi lock` of that stack is verified to still solve.


## Platform support

| Platform | Core features | Off-target search | GPU acceleration |
|---|---|---|---|
| Linux | Yes | Yes (auto with Pixi) | Yes |
| macOS | Yes | Yes (manual install) | CPU only |
| Windows | Yes | Yes (manual install) | CPU only |

## BEguidedOFF integration

Off-target search is optional. Set `search_offtargets: false` in the pipeline config to skip steps 4–5's off-target handling entirely — the final `beguided_results.tsv` then contains no off-target columns. When `search_offtargets: true` and `BEguidedOFF` is not installed, the pipeline aborts with an error pointing at the install instructions (it does not skip silently).

BEguidedOFF has two search modes (pipeline config key `offtargets.mode`):

- **`fast`** (default for the web server) — CPU FM-index filter. Enumerates exactly the sites with CCTOP score ≤ the search cutoff (the removal threshold `offtargets.remove_guides_score`), so removal decisions match the GPU mode without a GPU. Needs an `.fmi` index built once per reference genome.
- **`extensive`** — GPU brute force (≤ `offtargets.max_mismatches` mismatches). Requires CUDA; used via the `cuda` pixi environment (not available in docker).

### Via Pixi (recommended)

```bash
# GPU-less environments (web server, CPU-only CLI): fast mode
pixi install -e web

# GPU off-target search: CUDA environment
pixi install --features cuda
```

Both install `BEguidedOFF` from its public GitHub repo. The `web` environment uses the CPU-only base install (no torch); the `cuda` environment additionally provides the GPU (`extensive`) mode.

Then provide an FM-index for fast mode via `--beguidedoff-index` or `BEGUIDED_BEGUIDEDOFF_INDEX` — the FM-index is the only search input fast mode needs. Alternatively, set `--beguidedoff-genome` / `BEGUIDED_BEGUIDEDOFF_GENOME` to the encoded genome `.npz`: fast mode then derives the index path from the `.npz` sibling (same basename, `.fmi` suffix). The `.npz` itself is only read by extensive (GPU) mode, which requires it.

### Per-genome configuration

With several selectable genomes (e.g. human and mouse), each genome can have its own encoded genome and FM-index. The per-genome variables win over the global ones; the global setting is the fallback only for the **default genome** (and for single-genome deployments), so a non-default genome is never silently screened against another species' index:

| Variable | Purpose |
|---|---|
| `BEGUIDED_BEGUIDEDOFF_INDEX_HUMAN` / `..._MOUSE` | FM-index `.fmi` for that genome key (fast mode; primary) |
| `BEGUIDED_BEGUIDEDOFF_GENOME_HUMAN` / `..._MOUSE` | Encoded `.npz` for that genome key — optional in fast mode (index is derived from its sibling), required for extensive mode |
| `BEGUIDED_BEGUIDEDOFF_INDEX` / `_GENOME` / `_MODE` | Global fallback (default genome and single-genome deployments only) |

The `/config` endpoint reports per-genome availability (`offtarget.genomes`) so the UI can enable or disable the off-target search button per selected genome; a configured genome without its own BEguidedOFF index is reported unavailable with the variables to set. For fast mode it also reports the search window (`offtarget.fast_search`) — the exact PAM-proximal seed length and the tolerated PAM-distal mismatches, derived from the effective search cutoff — which the UI shows as a plain-language note in the results toolbar.

### Off-target screening from the web UI

The web UI uses a two-phase flow: a guide search runs the pipeline without BEguidedOFF; the user filters the Guides table; the **Check severe off-targets** button (in the filter bar above the results) screens exactly the visible guides. The screening job re-runs only steps 4–5 on that subset and keeps every row (the removal rule is disabled — the UI shows each guide's verdict from the operator thresholds). The verdict appears directly in the Guides table as a first **Off-target check** column (pass/fail badge; hover for the reason; unscreened guides show "not screened"). With **Condense results** enabled, a target whose editors disagree splits into a failing row and a passing row, so **Hide failing guides** removes only the failing rows and the target stays visible through its passing editors. The **Off-target report** tab holds a compact per-guide specificity report (worst off-target score, count and worst site) plus the full TSV download of the screened subset. The search enumerates only the most severe off-targets (CCTOP score at or below the operator's score threshold): "pass" means none that severe were found, not that the guide has no off-targets.

| Variable | Default | Purpose |
|---|---|---|
| `BEGUIDED_MAX_SCREEN_GUIDES` | `20000` | Hard cap on guides per screening job; at the cap the UI asks the user to filter further. Covers the largest legitimate payload (max targets × profiles). |
| `BEGUIDED_SCREEN_COOLDOWN_SECONDS` | `2` | Per-IP cooldown between screening submits. Tracked separately from the main submit cooldown (60 s), so a screening submit neither inherits nor suppresses it. |
| `BEGUIDED_PRELOAD_OFFTARGET_INDEX` | `1` | Open every configured fast-mode FM-index once at startup and share the handle across jobs (mmap-backed, ~1.7 GB page cache per genome; set `0` to load lazily per job). |

### Manual install

```bash
pip install "git+https://github.com/ltschmitt/BEguidedOFF.git"
```

The `fast` mode additionally needs the compiled FM-index binding (`beguidedoff._fmindex`), which is built automatically when the sdsl-lite headers/library are available at install time.

### Encoding a genome

Before running the pipeline with off-targets, encode the reference genome once:

```bash
python -m beguidedoff.cli encode-genome \
    --fasta genome/Homo_sapiens.GRCh38.116.fa \
    --output genome_otbe.npz \
    --include-rc \
    --include-n-mask
```

For fast mode, build the FM-index from the same reference FASTA:

```bash
python -m beguidedoff.cli build-index \
    --fasta genome/Homo_sapiens.GRCh38.116.fa \
    --output genome_otbe.fmi
```

Reference the resulting `.npz` in your pipeline config under `offtargets.beguidedoff_genome`, or pass `--beguidedoff-genome /path/to/genome_otbe.npz` on the CLI. Fast mode reads the `.fmi` (default location: next to the `.npz`; it cannot load without its `.meta` sidecar and needs the `.nmask` sidecar for the N-mismatch semantics — keep the trio together when copying); `extensive` mode reads the `.npz`.

Note: in fast mode the off-target count (`offtarget_count`) covers sites with score ≤ the search cutoff. A configured `offtargets.count_threshold_score` above that cutoff is clamped to it (with a warning), because deeper thresholds would need a shorter exact anchor whose genome-wide hit count exceeds the engine's per-anchor cap and would make counts incomplete.

## Docker

The repo ships a multi-stage `Dockerfile` that produces two images:

| Target | Tag | CUDA | Web | Use case |
|---|---|---|---|---|
| `base` | `beguided:cli` | No | No | CLI / pipeline |
| `web` (default) | `beguided:web` | No | Yes | Web interface |

The `web` image also installs [BEguidedOFF](https://github.com/ltschmitt/BEguidedOFF) (CPU fast mode) with its FM-index binding, so off-target screening works out of the box once the index files are provided — see [Off-target screening in docker](#off-target-screening-in-docker) below. There are no GPU docker targets: web screening is fast/CPU-only, and GPU off-target search runs from a pixi environment (see the manual install section).

### Build

```bash
docker build --target base -t beguided:cli .
docker build --target web  -t beguided:web .
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

### Off-target screening in docker

The `web` image ships BEguidedOFF (fast/CPU mode) from
`github.com/ltschmitt/BEguidedOFF`, but **not** the genome index files —
they are operator data (the human FM-index alone is ~1.7 GB) and are mounted
at `/app/genome`. Fast mode is configured with the FM-index alone: place
the `.fmi` and its sidecars (`.fmi.meta` chromosome table, required to load;
`.fmi.nmask` N runs, required for the N-mismatch semantics — all three files
must stay together) in the mounted genome directory, e.g. for human:

```
genome/hg38_beguidedoff.fmi        # FM-index (fast mode)
genome/hg38_beguidedoff.fmi.meta   # chromosome-table sidecar
genome/hg38_beguidedoff.fmi.nmask  # N-run sidecar (N-mismatch semantics)
```

Build them once from the reference FASTA with BEguidedOFF's CLI (see
"Encoding a genome" above — the `.npz` is only needed for extensive GPU
mode), or copy an existing pair from a workstation that already built it.
`docker-compose.yml` defaults `BEGUIDED_BEGUIDEDOFF_INDEX` to
`/app/genome/hg38_beguidedoff.fmi`. Until the files exist, `/config`
reports the genome as unavailable and the UI's off-target controls show a
clear error — the rest of the pipeline is unaffected. Mouse stays disabled
until an equivalent `*.fmi` pair and `BEGUIDED_BEGUIDEDOFF_INDEX_MOUSE`
are configured.

## Verifying the install

```bash
beguided --help
beguided-web --help
```

If both commands print their usage strings, the install is complete.
