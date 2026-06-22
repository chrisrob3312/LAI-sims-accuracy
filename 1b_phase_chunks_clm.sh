#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_chunk
#SBATCH --partition=mhgcp
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=06:00:00
#SBATCH --time-min=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=12
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_chunk_%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_chunk_%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Stage 2 of the chunked-phasing pipeline. Array task M reads line M of
# $CHUNKS_TSV (built by resources/build_b38_chunks.py) and runs
# SHAPEIT5_phase_common on that (chr, region) pair against the per-chr
# softunion BCF produced by 1b_phase_prep_clm.sh.
#
# Submit with --array=1-<N> where N = total chunk count (lines in chunks.tsv
# minus the header). The submit_chunked_phasing.sh wrapper figures this out
# automatically.
#
# Output per chunk:
#   ${OUTDIR}/chunks/chr${CHR}_chunk${IDX}.phased.bcf{,.csi}
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"

set -euo pipefail
TASK=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-32}

[[ -s "$CHUNKS_TSV" ]] || { echo "ERROR: missing $CHUNKS_TSV"; exit 1; }

# Read line $TASK+1 (skip header). Columns: chr_num chr chunk_idx phase_region ligate_region
LINE=$(awk -v t="$TASK" 'NR==t+1' "$CHUNKS_TSV")
[[ -n "$LINE" ]] || { echo "ERROR: task $TASK has no chunks.tsv row"; exit 1; }
read -r CHR_NUM CHR CHUNK_IDX PHASE_REGION LIG_REGION <<< "$LINE"

GMAP="${GMAP_DIR}/${CHR}.b38.gmap.gz"
SOFTUNION="${OUTDIR}/merged_${CHR}.softunion.bcf"
CHUNK_DIR="${OUTDIR}/chunks"
PHASED_CHUNK="${CHUNK_DIR}/${CHR}_chunk${CHUNK_IDX}.phased.bcf"
mkdir -p "$CHUNK_DIR"

[[ -s "$SOFTUNION" ]] || { echo "ERROR: missing $SOFTUNION (run 1b_phase_prep first)"; exit 1; }
[[ -s "$GMAP"      ]] || { echo "ERROR: missing $GMAP"; exit 1; }

# Idempotent: skip if already done.
if [[ -s "$PHASED_CHUNK" && -s "${PHASED_CHUNK}.csi" ]]; then
    echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] already done, skipping"
    exit 0
fi

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] phase_common ${PHASE_REGION}"
SHAPEIT5_phase_common \
    --input "$SOFTUNION" \
    --map "$GMAP" \
    --region "$PHASE_REGION" \
    --output "$PHASED_CHUNK" \
    --thread "$THREADS" \
    --filter-maf 0.001

echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] DONE -> $PHASED_CHUNK"
