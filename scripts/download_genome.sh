#!/bin/bash
# Download reference genome(s) (FASTA + GTF) into a local directory,
# decompress, and rename the FASTA to a release-tagged name.
#
# Supported species:
#   human  Ensembl GRCh38 (FASTA + GTF)
#   mouse  EBI ensemblorganisms GRCm39 (GCA_000001635.9; unmasked FASTA + GTF)
#
# Re-running is safe: existing final outputs are skipped unless --force is given.
#
# Usage:
#   scripts/download_genome.sh                       # both species, ./genome
#   scripts/download_genome.sh --species human
#   scripts/download_genome.sh --species mouse
#   scripts/download_genome.sh --target-dir /data/refs
#   scripts/download_genome.sh --release 115
#   scripts/download_genome.sh --mouse-release 2026_04
#   scripts/download_genome.sh --force
#   scripts/download_genome.sh -h | --help
#
# Outputs:
#   <target-dir>/Homo_sapiens.GRCh38.<release>.fa
#   <target-dir>/Homo_sapiens.GRCh38.<release>.gtf
#   <target-dir>/Mus_musculus.GRCm39.<mouse-release>.fa
#   <target-dir>/Mus_musculus.GRCm39.<mouse-release>.gtf

set -euo pipefail

TARGET_DIR="./genome"
RELEASE=116
MOUSE_RELEASE="2026_04"
SPECIES="all"
FORCE=0

usage() {
    sed -n '2,26p' "$0" | sed 's/^# \{0,1\}//'
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
        --mouse-release)
            MOUSE_RELEASE="$2"
            shift 2
            ;;
        --species)
            SPECIES="$2"
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

case "$SPECIES" in
    human|mouse|all) ;;
    *)
        echo "Error: --species must be one of: human, mouse, all (got '$SPECIES')" >&2
        exit 1
        ;;
esac

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

    if [[ -f "$final" && $FORCE -eq 0 ]]; then
        echo "Exists, skipping: $final"
        return 0
    fi

    echo "Downloading $url"
    curl -fSL --retry 3 -o "${gz_name}.tmp" "$url"
    mv "${gz_name}.tmp" "$gz_name"

    echo "Decompressing $gz_name -> $final"
    # Use `gunzip -c` rather than in-place `gunzip`: the mouse FASTA is a
    # BGZF file named `unmasked.fa.bgz`, and plain `gunzip` would leave an
    # intermediate `*.fa.b` instead of stripping the suffix.
    gunzip -c -- "$gz_name" > "${final}.tmp"
    mv -- "${final}.tmp" "$final"
    rm -f -- "$gz_name"
}

if [[ "$SPECIES" == "human" || "$SPECIES" == "all" ]]; then
    FASTA_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
    FASTA_GZ_NAME="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GTF_GZ_NAME="Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
    FASTA_FINAL="Homo_sapiens.GRCh38.${RELEASE}.fa"
    GTF_FINAL="Homo_sapiens.GRCh38.${RELEASE}.gtf"

    download_and_unpack "$FASTA_URL" "$FASTA_GZ_NAME" "$FASTA_FINAL"
    download_and_unpack "$GTF_URL" "$GTF_GZ_NAME" "$GTF_FINAL"
fi

if [[ "$SPECIES" == "mouse" || "$SPECIES" == "all" ]]; then
    MOUSE_BASE="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/635/9/ensembl/${MOUSE_RELEASE}"
    FASTA_URL="${MOUSE_BASE}/genome/unmasked.fa.bgz"
    GTF_URL="${MOUSE_BASE}/geneset/genes.gtf.gz"
    FASTA_GZ_NAME="unmasked.fa.bgz"
    GTF_GZ_NAME="genes.gtf.gz"
    FASTA_FINAL="Mus_musculus.GRCm39.${MOUSE_RELEASE}.fa"
    GTF_FINAL="Mus_musculus.GRCm39.${MOUSE_RELEASE}.gtf"

    download_and_unpack "$FASTA_URL" "$FASTA_GZ_NAME" "$FASTA_FINAL"
    download_and_unpack "$GTF_URL" "$GTF_GZ_NAME" "$GTF_FINAL"
fi

echo
echo "Done. Files in $(pwd):"
ls -lh -- ./*.fa ./*.gtf
