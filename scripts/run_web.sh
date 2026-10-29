#!/bin/bash
# Local launcher for the BEguided web server.
#
# Configures the scheduler (queue, cooldown, SQLite job store), the
# operator-supplied paths (presets, genome, GTF, BEguidedOFF npz/index)
# and the off-target screening knobs, then hands off to `beguided-web`.
#
# All `beguided-web` flags are exposed as environment variables below;
# any CLI argument passed to this script is forwarded to the server last,
# so it always wins over the defaults.
#
# Usage:
#   ./run_web.sh
#   ./run_web.sh --port 8080
#   BEGUIDED_DB_PATH=/tmp/beguided.db ./run_web.sh
#
# NOTE: This script is for personal/local use. It is tracked in the repo,
# but its hardcoded paths assume a local checkout; production uses the
# container entrypoint (see docs/web-deployment.md).

set -euo pipefail

# Repo root (the script lives in scripts/): all default paths below are
# relative to the checkout, and the FM-index readability checks must resolve
# them the same way the server does.
cd "$(dirname "$0")/.."

# ---------------------------------------------------------------------------
# Scheduler / safety knobs (override via env vars)
# ---------------------------------------------------------------------------
: "${BEGUIDED_DB_PATH:=/var/lib/beguided/jobs.db}"
: "${BEGUIDED_QUEUE_MAX:=1}"
: "${BEGUIDED_COOLDOWN_SECONDS:=5}"
: "${BEGUIDED_MAX_UPLOAD_MB:=5}"

# Make sure the SQLite directory exists so the scheduler can open the file.
DB_DIR="$(dirname "$BEGUIDED_DB_PATH")"
if [ ! -d "$DB_DIR" ]; then
    mkdir -p "$DB_DIR" 2>/dev/null || {
        echo "Warning: cannot create $DB_DIR (need sudo?). Falling back to /tmp/beguided_jobs.db." >&2
        BEGUIDED_DB_PATH="/tmp/beguided_jobs.db"
    }
fi

# ---------------------------------------------------------------------------
# Data paths (override via env vars or pass --base-editor-profiles-file /
# --genome-fasta etc. on the command line)
# ---------------------------------------------------------------------------
: "${BEGUIDED_BASE_EDITOR_PROFILES_FILE:=configs/base_editor_profiles.tsv}"
: "${BEGUIDED_GENOME_FASTA:=genome/Homo_sapiens.GRCh38.116.fa}"
: "${BEGUIDED_GTF_PATH:=genome/Homo_sapiens.GRCh38.116.gtf}"
: "${BEGUIDED_MOUSE_GENOME_FASTA:=genome/Mus_musculus.GRCm39.2026_04.fa}"
: "${BEGUIDED_MOUSE_GTF_PATH:=genome/Mus_musculus.GRCm39.2026_04.gtf}"
# BEguidedOFF off-target search. genome/hg38_beguidedoff.{npz,fmi} are
# symlinks into the BEguidedOFF checkout's nogit/ artifacts (GRCh38 primary
# assembly; coordinates match this repo's GRCh38.116 reference FASTA).
# The global setting covers the human (default) genome; screening for any
# other genome needs its own BEGUIDED_BEGUIDEDOFF_GENOME_<KEY> /
# BEGUIDED_BEGUIDEDOFF_INDEX_<KEY> (e.g. _MOUSE once a mouse index exists).
: "${BEGUIDED_BEGUIDEDOFF_GENOME:=genome/hg38_beguidedoff.npz}"
: "${BEGUIDED_BEGUIDEDOFF_INDEX:=genome/hg38_beguidedoff.fmi}"
: "${BEGUIDED_BEGUIDEDOFF_MODE:=fast}"

# Off-target screening knobs (the UI's "Search off-targets" button on
# filtered results). Screening uses a separate 2 s cooldown and a hard cap
# on the number of guides per job.
: "${BEGUIDED_MAX_SCREEN_GUIDES:=20000}"
: "${BEGUIDED_SCREEN_COOLDOWN_SECONDS:=2}"
: "${BEGUIDED_PRELOAD_OFFTARGET_INDEX:=1}"

if [ "$BEGUIDED_BEGUIDEDOFF_MODE" = "fast" ]; then
    if [ ! -r "$BEGUIDED_BEGUIDEDOFF_INDEX" ]; then
        echo "Warning: fast-mode FM-index not found: $BEGUIDED_BEGUIDEDOFF_INDEX" >&2
        echo "         off-target screening will show as unavailable (see /config)." >&2
    elif [ ! -r "$BEGUIDED_BEGUIDEDOFF_INDEX.meta" ]; then
        echo "Warning: FM-index sidecar not found: $BEGUIDED_BEGUIDEDOFF_INDEX.meta" >&2
        echo "         an .fmi cannot load without its .meta file; off-target screening will show as unavailable." >&2
    fi
fi

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------
echo "BEguided web launcher" >&2
echo "  db:               $BEGUIDED_DB_PATH" >&2
echo "  queue_max:        $BEGUIDED_QUEUE_MAX" >&2
echo "  cooldown:         ${BEGUIDED_COOLDOWN_SECONDS}s" >&2
echo "  upload cap:       ${BEGUIDED_MAX_UPLOAD_MB}MB" >&2
echo "  base editor file: $BEGUIDED_BASE_EDITOR_PROFILES_FILE" >&2
echo "  genome:           $BEGUIDED_GENOME_FASTA" >&2
echo "  gtf:              $BEGUIDED_GTF_PATH" >&2
echo "  mouse genome:     $BEGUIDED_MOUSE_GENOME_FASTA" >&2
echo "  mouse gtf:        $BEGUIDED_MOUSE_GTF_PATH" >&2
echo "  beguidedoff npz:  $BEGUIDED_BEGUIDEDOFF_GENOME" >&2
echo "  beguidedoff idx:  $BEGUIDED_BEGUIDEDOFF_INDEX" >&2
echo "  beguidedoff mode: $BEGUIDED_BEGUIDEDOFF_MODE" >&2
echo "  screen cap:       $BEGUIDED_MAX_SCREEN_GUIDES guides" >&2
echo "  screen cooldown:  ${BEGUIDED_SCREEN_COOLDOWN_SECONDS}s" >&2
echo "  screen preload:   $BEGUIDED_PRELOAD_OFFTARGET_INDEX" >&2

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
exec pixi run --environment web beguided-web \
    --host 0.0.0.0 \
    --port 8000 \
    --db-path "$BEGUIDED_DB_PATH" \
    --queue-max "$BEGUIDED_QUEUE_MAX" \
    --cooldown-seconds "$BEGUIDED_COOLDOWN_SECONDS" \
    --max-upload-mb "$BEGUIDED_MAX_UPLOAD_MB" \
    --base-editor-profiles-file "$BEGUIDED_BASE_EDITOR_PROFILES_FILE" \
    --genome-fasta "$BEGUIDED_GENOME_FASTA" \
    --gtf-path "$BEGUIDED_GTF_PATH" \
    --mouse-genome-fasta "$BEGUIDED_MOUSE_GENOME_FASTA" \
    --mouse-gtf-path "$BEGUIDED_MOUSE_GTF_PATH" \
    --beguidedoff-genome "$BEGUIDED_BEGUIDEDOFF_GENOME" \
    --beguidedoff-index "$BEGUIDED_BEGUIDEDOFF_INDEX" \
    --beguidedoff-mode "$BEGUIDED_BEGUIDEDOFF_MODE" \
    "$@"
