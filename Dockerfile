# BEguided Docker image
# Builds a container with core BEguided functionality.
# Off-target search via BEguidedOFF is optional — install it manually if needed.

FROM python:3.12-slim-bookworm AS base

WORKDIR /app

# Install system build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy source and install BEguided
COPY pyproject.toml README.md ./
COPY beguided/ ./beguided/
RUN pip install --no-cache-dir -e "."

# Create directories for user-provided data
RUN mkdir -p /app/genome /app/data

# ---------------------------------------------------------------------------
# Web server stage (default)
# ---------------------------------------------------------------------------
FROM base AS web

ENV PYTHONUNBUFFERED=1
ENV BEGUIDED_PRESETS_FILE=/app/be_presets.tsv

EXPOSE 8000

# Healthcheck
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/config', timeout=5)" || exit 1

# Default: run web server with no external data mounted
# Mount your genome/presets/gtf at runtime:
#   docker run -v $(pwd)/be_presets.tsv:/app/be_presets.tsv \
#              -v $(pwd)/genome:/app/genome \
#              -p 8000:8000 \
#              beguided web --host 0.0.0.0 --port 8000
CMD ["beguided", "web", "--host", "0.0.0.0", "--port", "8000"]
