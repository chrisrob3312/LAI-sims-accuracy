#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_chunk
#SBATCH --partition=atkinson,mhgcp
# DO NOT REMOVE mhgcp-t* FROM THIS EXCLUSION TO GAIN CONCURRENCY.
#
# Those nodes' CPUs lack instructions the SHAPEIT5 / bcftools builds in this
# conda env are compiled with. A task dispatched to one dies instantly with
# SIGILL -- exit 132, ~3 seconds wall clock, ~2 MB RSS, and an EMPTY log,
# because the binary never reaches its first write.
#
# Measured 2026-11-07 on a 594-task chunk array run with the exclusion lifted:
#   270 / 270 tasks that landed on t02,t05,t06,t07,t08,t10,t11  FAILED
#     0 / 57  tasks that landed on a01, m00                     failed
#
# They look attractive because sinfo shows them idle with free CPUs. They are
# idle because nothing that needs these instruction sets can run there. An
# earlier 1d_almarri failure on t02 -- empty logs, 4 MB RSS, never explained
# -- was the same thing.
#
# mhgcp-c02 IS A SEPARATE QUESTION and the evidence above says nothing about
# it. No task landed on c02 in that run, so it is neither condemned nor
# cleared. It has the same specs as c00 and c01 (56 CPUs, 232 GB), which both
# run this software without trouble and are NOT excluded, so the reason for
# excluding c02 is something else -- unknown to this script. Four test tasks
# would settle it and would add 56 CPUs.
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=08:00:00
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
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"

set -euo pipefail
TASK=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-32}

# The environment is set up FIRST, before anything that runs a tool.
#
# It used to sit just above the SHAPEIT5 call, which was fine while nothing
# above it needed a binary. The staleness check added below calls bcftools,
# and with conda still unactivated that is a command not found: exit 127,
# under `set -e`, after the task had already read its chunk row. The array
# then failed as "Mixed [0-127]" -- tasks with no existing chunk skipped the
# check and exited 0, tasks that found a June chunk died -- which reads like a
# data problem and is not one.
module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH after activating ${CONDA_ENV}"; exit 1; }

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
    # "Exists" is not "current". These chunks persist across runs, and the
    # panel they were phased from changes: adding SGDP and Almarri to the
    # merge leaves June's chunks on disk, byte-valid and silently wrong. A
    # rerun would skip every one of them and ligate a panel missing two whole
    # cohorts, with nothing in any log saying so.
    #
    # The sample count is the cheap discriminator: it comes from the index
    # header, and a chunk phased from a different panel has a different one.
    _n_have=$(bcftools query -l "$PHASED_CHUNK" 2>/dev/null | wc -l)
    _n_want=$(bcftools query -l "$SOFTUNION" 2>/dev/null | wc -l)
    if [[ "$_n_have" -eq "$_n_want" && "$_n_have" -gt 0 ]]; then
        echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] already done (${_n_have} samples), skipping"
        exit 0
    fi
    echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] STALE: has ${_n_have} samples,"
    echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] soft-union has ${_n_want} -- rebuilding"
    rm -f "$PHASED_CHUNK" "${PHASED_CHUNK}.csi"
fi

echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] phase_common ${PHASE_REGION}"
SHAPEIT5_phase_common \
    --input "$SOFTUNION" \
    --map "$GMAP" \
    --region "$PHASE_REGION" \
    --output "$PHASED_CHUNK" \
    --thread "$THREADS" \
    --filter-maf 0.001

echo "[$(date +%T)] [${CHR} chunk${CHUNK_IDX}] DONE -> $PHASED_CHUNK"
