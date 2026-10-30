# BEguided Docker image
#
# Targets:
#   base      – CLI-only, CPU (python:3.12-slim)
#   web       – Web server, CPU (default)
#
# The `web` target additionally installs BEguidedOFF
# (https://github.com/ltschmitt/BEguidedOFF) for off-target screening,
# CPU fast mode only (FM-index) — the web server never uses the GPU. The
# FM-index binding is a native extension built against sdsl-lite at image
# build time (libsdsl-dev) and forced on with BEGUIDEDOFF_FAST=1 so a
# missing sdsl-lite fails the build instead of silently shipping without
# it. The .npz/.fmi index files are operator-provided data, not part of
# the image (mounted at /app/genome; see docker-compose.yml and
# docs/installation.md).
#
# Build examples:
#   docker build --target base -t beguided:cli     .
#   docker build --target web  -t beguided:web     .
#
# GPU off-target search (BEguidedOFF extensive mode) is not available in
# docker; run the CLI from a pixi environment with the `cuda` feature.
#
# Package installs use uv (https://docs.astral.sh/uv/) pulled from its
# official image — a single static binary, no wheel cache to clean up.
# `UV_BREAK_SYSTEM_PACKAGES=1` lets uv install into the externally-managed
# Python from `python:3.12-slim`; it's a no-op for the venv-based GPU targets.

# ============================================================================
# 1. CPU BASE — CLI-only, minimal image (no web server deps)
# ============================================================================
FROM python:3.12-slim-bookworm AS base

WORKDIR /app

# Prevent interactive prompts (e.g. tzdata) during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV UV_BREAK_SYSTEM_PACKAGES=1

# Install uv (single static binary) from the official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Install system build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy source and install BEguided (core only)
COPY pyproject.toml README.md ./
COPY beguided/ ./beguided/
COPY configs/base_editor_profiles.tsv ./configs/
RUN uv pip install --system --no-cache -e "."

# Create directories for user-provided data
RUN mkdir -p /app/genome /app/data

# Default entrypoint for CLI usage
ENTRYPOINT ["beguided"]
CMD ["--help"]

# ============================================================================
# 2. CPU WEB — FastAPI web server (default target)
# ============================================================================
FROM base AS web

# sdsl-lite for the BEguidedOFF FM-index binding, and git for uv's
# git+https dependency fetch. Debian's libsdsl-dev puts the library in the
# multiarch directory, which BEguidedOFF's setup.py does not search by
# default, so SDSL_LIB_DIR is set explicitly.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libsdsl-dev \
    && rm -rf /var/lib/apt/lists/*
ENV SDSL_INCLUDE_DIR=/usr/include \
    SDSL_LIB_DIR=/usr/lib/x86_64-linux-gnu \
    BEGUIDEDOFF_FAST=1

RUN uv pip install --system --no-cache -e ".[web]" \
    && uv pip install --system --no-cache git+https://github.com/ltschmitt/BEguidedOFF.git
RUN mkdir -p /var/lib/beguided

ENV PYTHONUNBUFFERED=1
ENV BEGUIDED_BASE_EDITOR_PROFILES_FILE=/app/configs/base_editor_profiles.tsv

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/config', timeout=5)" || exit 1

ENTRYPOINT []

CMD ["beguided-web", "--host", "0.0.0.0", "--port", "8000"]
