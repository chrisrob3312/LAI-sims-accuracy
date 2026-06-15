#!/usr/bin/env bash
# ============================================================================
# submit_chunked_phasing.sh
#
# Orchestrate the 3-stage chunked-phasing pipeline with SLURM dependency
# chains so the user runs one command:
#
#   ./submit_chunked_phasing.sh
#
# Steps:
#   0. Build resources/b38_chunks.tsv (8cM phasing windows, 2cM overlap)
#   1. Submit 1b_phase_prep_clm.sh as array 1-22 -> prep_jobid
#   2. Submit 1b_phase_chunks_clm.sh as array 1-<n_chunks>, dep:afterok prep
#   3. Submit 1b_phase_ligate_clm.sh as array 1-22,    dep:afterok chunks
#
# Re-running is safe: each stage's script is idempotent (skips work already
# done). To force a re-run of one chunk, delete its phased BCF file.
# ============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"
RESOURCES="${RESOURCES:-${PROJECT_ROOT}/resources}"
CHUNKS_TSV="${CHUNKS_TSV:-${RESOURCES}/b38_chunks.tsv}"

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$RESOURCES"

# 0. Build chunks file if not already present.
if [[ ! -s "$CHUNKS_TSV" ]]; then
    echo "[$(date +%T)] building $CHUNKS_TSV"
    module load anaconda3/2024.06 2>/dev/null || true
    python3 "${REPO_DIR}/resources/build_b38_chunks.py" "$GMAP_DIR" "$CHUNKS_TSV"
fi
N_CHUNKS=$(( $(wc -l < "$CHUNKS_TSV") - 1 ))
echo "[$(date +%T)] $N_CHUNKS chunks in $CHUNKS_TSV"
echo "  per-chr breakdown:"
awk -F'\t' 'NR>1 {c[$2]++} END {for (k in c) print "    " k ": " c[k] " chunks"}' "$CHUNKS_TSV" | sort

# Pass chunks path to child jobs via env.
export PROJECT_ROOT CHUNKS_TSV

# 1. prep stage
PREP_JID=$(sbatch --parsable "${REPO_DIR}/1b_phase_prep_clm.sh")
echo "[$(date +%T)] submitted prep array (1-22) as $PREP_JID"

# 2. chunks stage (depends on all 22 prep tasks)
CHUNK_JID=$(sbatch --parsable \
    --dependency=afterok:"$PREP_JID" \
    --array=1-"$N_CHUNKS"%50 \
    "${REPO_DIR}/1b_phase_chunks_clm.sh")
echo "[$(date +%T)] submitted chunks array (1-${N_CHUNKS}) as $CHUNK_JID (waits on prep)"
echo "  (%50 throttle = at most 50 concurrent chunk tasks; raise/lower to taste)"

# 3. ligate stage (depends on all chunks)
LIG_JID=$(sbatch --parsable \
    --dependency=afterok:"$CHUNK_JID" \
    "${REPO_DIR}/1b_phase_ligate_clm.sh")
echo "[$(date +%T)] submitted ligate array (1-22) as $LIG_JID (waits on chunks)"

cat <<EOF

Pipeline submitted.
  prep    job $PREP_JID   array 1-22
  chunks  job $CHUNK_JID  array 1-$N_CHUNKS
  ligate  job $LIG_JID    array 1-22

Monitor with:
  squeue -u \$USER
  watch -n 30 'squeue -u \$USER'

Outputs land in:
  $PROJECT_ROOT/01_merged_phased_panel/merged_chr<N>.shapeit5_phased.softunion_maf005.bcf
  $PROJECT_ROOT/01_merged_phased_panel/merged_chr<N>.shapeit5_phased.softunion_maf005.rechr.bcf
EOF
