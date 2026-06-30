#!/bin/bash
# Launch the BEguided web server from a Docker image.
# Mounts host genome/ data into the container.
#
# Usage:
#   ./run_web_docker.sh              # use web-cuda image (GPU)
#   ./run_web_docker.sh web          # use web image (CPU)
#   ./run_web_docker.sh web-cuda --port 8080   # custom port
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=${1:-beguided:web-cuda}
if [ "$IMAGE" != "web-cuda" ] && [ "$IMAGE" != "web" ] && [ "$IMAGE" != "cli-cuda" ] && [ "$IMAGE" != "cli" ]; then
    # First arg was not a known image tag – treat it as an extra docker/CLI arg
    IMAGE="beguided:web-cuda"
else
    shift  # consume the image tag
fi
docker run --gpus all \
    -p 8000:8000 \
    -v "$(pwd)/be_presets_wt.tsv:/app/be_presets_wt.tsv:ro" \
    -v "$(pwd)/be_presets.tsv:/app/be_presets.tsv:ro" \
    -v "$(pwd)/genome:/app/genome:ro" \
    -e BEGUIDED_PRESETS_FILE=/app/be_presets_wt.tsv \
    -e BEGUIDED_GENOME_FASTA=/app/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -e BEGUIDED_GTF_PATH=/app/genome/Homo_sapiens.GRCh38.112.gtf \
    -e BEGUIDED_BEGUIDEDOFF_GENOME=/app/genome/hg38_otbe.npz \
    "$IMAGE" \
    beguided web "$@" \
    --host 0.0.0.0 \
    --port 8000
