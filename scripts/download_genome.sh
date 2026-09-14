#!/bin/bash
# Download the Ensembl human reference genome (FASTA + GTF) into a local
# directory, decompress, and rename the FASTA to a release-tagged name.
#
# Designed for fresh machines with no genome yet. Re-running is safe: existing
# final outputs are skipped unless --force is given.
#
# Usage:
#   scripts/download_genome.sh                       # default: ./genome, release 116
#   scripts/download_genome.sh --target-dir /data/refs/GRCh38.116
#   scripts/download_genome.sh --release 115
#   scripts/download_genome.sh --force               # re-download even if present
#   scripts/download_genome.sh -h | --help
#
# Outputs:
#   <target-dir>/Homo_sapiens.GRCh38.<release>.fa
#   <target-dir>/Homo_sapiens.GRCh38.<release>.gtf

set -euo pipefail

TARGET_DIR="./genome"
RELEASE=116
FORCE=0

usage() {
    sed -n '2,18p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-dir)
            TARGET_DIR="$2"
            shift 2
            ;;
        --release)
            RELEASE="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

FASTA_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
FASTA_GZ_NAME="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_GZ_NAME="Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
FASTA_FINAL="Homo_sapiens.GRCh38.${RELEASE}.fa"
GTF_FINAL="Homo_sapiens.GRCh38.${RELEASE}.gtf"

command -v curl >/dev/null 2>&1 || { echo "Error: 'curl' not found in PATH" >&2; exit 1; }
command -v gunzip >/dev/null 2>&1 || { echo "Error: 'gunzip' not found in PATH" >&2; exit 1; }

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

cleanup() {
    rm -f ./*.tmp 2>/dev/null || true
}
trap cleanup EXIT

download_and_unpack() {
    local url="$1"
    local gz_name="$2"
    local final="$3"
    local stem="${gz_name%.gz}"

    if [[ -f "$final" && $FORCE -eq 0 ]]; then
        echo "Exists, skipping: $final"
        return 0
    fi

    echo "Downloading $url"
    curl -fSL --retry 3 -o "${gz_name}.tmp" "$url"
    mv "${gz_name}.tmp" "$gz_name"

    echo "Decompressing $gz_name -> $final"
    gunzip -- "$gz_name"  # strips the .gz, leaves $stem
    if [[ "$stem" != "$final" ]]; then
        mv -- "$stem" "$final"
    fi
}

download_and_unpack "$FASTA_URL" "$FASTA_GZ_NAME" "$FASTA_FINAL"
download_and_unpack "$GTF_URL" "$GTF_GZ_NAME" "$GTF_FINAL"

echo
echo "Done. Files in $(pwd):"
ls -lh -- "$FASTA_FINAL" "$GTF_FINAL"
