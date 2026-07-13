# BEguided Docker image
# Multi-stage build with CPU and GPU variants.
#
# Targets:
#   base      – CLI-only, CPU (python:3.12-slim)
#   web       – Web server, CPU (default)
#   base-gpu  – CLI-capable, CUDA 12.4 runtime
#   cli-cuda  – CLI-only with GPU support
#   web-cuda  – Web server with GPU support
#
# Build examples:
#   docker build --target base     -t beguided:cli     .
#   docker build --target web      -t beguided:web     .
#   docker build --target cli-cuda -t beguided:cli-cuda .
#   docker build --target web-cuda -t beguided:web-cuda .
#
# Run GPU targets (requires NVIDIA Container Toolkit):
#   docker run --gpus all beguided:cli-cuda --help
#   docker run --gpus all -p 8000:8000 beguided:web-cuda
#
# Off-target search via BEguidedOFF is optional — install it manually if needed.

# ============================================================================
# 1. CPU BASE — CLI-only, minimal image (no web server deps)
# ============================================================================
FROM python:3.12-slim-bookworm AS base

WORKDIR /app

# Prevent interactive prompts (e.g. tzdata) during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy source and install BEguided (core only)
COPY pyproject.toml README.md ./
COPY beguided/ ./beguided/
COPY base_editor_profiles.tsv ./
RUN pip install --no-cache-dir -e "."

# Create directories for user-provided data
RUN mkdir -p /app/genome /app/data

# Default entrypoint for CLI usage
ENTRYPOINT ["beguided"]
CMD ["--help"]

# ============================================================================
# 2. CPU WEB — FastAPI web server (default target)
# ============================================================================
FROM base AS web

RUN pip install --no-cache-dir -e ".[web]"
RUN mkdir -p /var/lib/beguided

ENV PYTHONUNBUFFERED=1
ENV BEGUIDED_BASE_EDITOR_PROFILES_FILE=/app/base_editor_profiles.tsv

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/config', timeout=5)" || exit 1

ENTRYPOINT []

CMD ["beguided-web", "--host", "0.0.0.0", "--port", "8000"]

# ============================================================================
# 3. GPU BASE — CUDA 12.4 runtime, CLI-capable
# ============================================================================
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04 AS base-gpu

WORKDIR /app

# Prevent interactive prompts (e.g. tzdata) during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.12 and build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update && apt-get install -y --no-install-recommends \
    python3.12 \
    python3.12-venv \
    python3.12-dev \
    gcc \
    g++ \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Create an isolated virtual environment
RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy source and install BEguided (core only, includes torch)
COPY pyproject.toml README.md ./
COPY beguided/ ./beguided/
RUN pip install --no-cache-dir -e "."

# Create directories for user-provided data
RUN mkdir -p /app/genome /app/data

# Default entrypoint for CLI usage
ENTRYPOINT ["beguided"]
CMD ["--help"]

# ============================================================================
# 4. GPU CLI — CUDA-enabled CLI image
# ============================================================================
FROM base-gpu AS cli-cuda

# Inherits base-gpu configuration

# ============================================================================
# 5. GPU WEB — CUDA-enabled web server
# ============================================================================
FROM base-gpu AS web-cuda

RUN pip install --no-cache-dir -e ".[web]"
RUN mkdir -p /var/lib/beguided

ENV PYTHONUNBUFFERED=1
ENV BEGUIDED_BASE_EDITOR_PROFILES_FILE=/app/base_editor_profiles.tsv

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/config', timeout=5)" || exit 1

ENTRYPOINT []

CMD ["beguided-web", "--host", "0.0.0.0", "--port", "8000"]
