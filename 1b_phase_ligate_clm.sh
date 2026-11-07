#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_ligate
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
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=00:15:00
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_ligate_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_ligate_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Stage 3 of the chunked-phasing pipeline. Per-chrom (array 1-22):
#   - Enumerate the per-chunk phased BCFs produced by 1b_phase_chunks_clm.sh
#     in chunk-index order
#   - Run SHAPEIT5_ligate to stitch them into one chr-level phased BCF
#   - Rename chr${CHR} -> ${CHR} for RFMix v1
#
# Outputs:
#   ${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf{,.csi}
#   ${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf{,.csi}
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-16}
CHUNK_DIR="${OUTDIR}/chunks"
PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf"

[[ -s "$CHUNKS_TSV" ]] || { echo "ERROR: missing $CHUNKS_TSV"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# Skip if already done -- but "done" must mean "done FROM THE CURRENT CHUNKS",
# not merely "a file with this name exists".
#
# These outputs persist across runs and the panel behind them changes. The
# June build left 22 .rechr.bcf files on disk; after SGDP and Almarri were
# added to the merge, a rerun would have found them, skipped all 22, and
# handed the next stage a panel missing two entire cohorts -- valid files,
# plausible sizes, nothing in any log to say the phasing never happened.
#
# The sample count separates the two cases and costs one index read.
if [[ -s "$RECHR" && -s "${RECHR}.csi" ]]; then
    _n_have=$(bcftools query -l "$RECHR" 2>/dev/null | wc -l)
    # Compare against a chunk that is about to be ligated, not against the
    # soft-union: the chunks are what this step actually consumes.
    _probe=$(awk -v c="chr${CHR}" -F'\t' 'NR>1 && $2==c {print $3; exit}' "$CHUNKS_TSV")
    _probe_bcf="${CHUNK_DIR}/chr${CHR}_chunk${_probe}.phased.bcf"
    _n_want=0
    [[ -s "$_probe_bcf" ]] && _n_want=$(bcftools query -l "$_probe_bcf" 2>/dev/null | wc -l)
    if [[ "$_n_want" -gt 0 && "$_n_have" -eq "$_n_want" ]]; then
        echo "[$(date +%T)] [chr${CHR}] rechr BCF current (${_n_have} samples), skipping"
        exit 0
    fi
    if [[ "$_n_want" -eq 0 ]]; then
        echo "[$(date +%T)] [chr${CHR}] cannot read a chunk to compare against;"
        echo "[$(date +%T)] [chr${CHR}] re-ligating rather than trusting the existing file"
    else
        echo "[$(date +%T)] [chr${CHR}] STALE: existing rechr has ${_n_have} samples,"
        echo "[$(date +%T)] [chr${CHR}] the chunks have ${_n_want} -- re-ligating"
    fi
    rm -f "$RECHR" "${RECHR}.csi" "$PHASED" "${PHASED}.csi"
fi

# Build the ordered chunk-file list for this chr.
CHUNKS_LIST="${OUTDIR}/tmp/chr${CHR}/chunks_list.txt"
mkdir -p "$(dirname "$CHUNKS_LIST")"
awk -v c="chr${CHR}" -F'\t' 'NR>1 && $2==c {print $3}' "$CHUNKS_TSV" \
    | sort -n \
    | awk -v dir="$CHUNK_DIR" -v c="chr${CHR}" '{print dir"/"c"_chunk"$1".phased.bcf"}' \
    > "$CHUNKS_LIST"
NCHUNK=$(wc -l < "$CHUNKS_LIST")
echo "[$(date +%T)] [chr${CHR}] ligating ${NCHUNK} chunks"

# Sanity: all chunk files must exist.
while read -r f; do
    [[ -s "$f" && -s "${f}.csi" ]] || { echo "ERROR: missing chunk $f"; exit 1; }
done < "$CHUNKS_LIST"

# Ligate.
SHAPEIT5_ligate \
    --input "$CHUNKS_LIST" \
    --output "$PHASED" \
    --thread "$THREADS"

# SHAPEIT5_ligate does NOT emit a .csi -- index the ligated BCF here.
bcftools index -f "$PHASED"
NRECS=$(bcftools view "$PHASED" -H | wc -l)
echo "[$(date +%T)] [chr${CHR}] ligated: ${NRECS} records"

# Strip chr-prefix for RFMix v1.
echo "[$(date +%T)] [chr${CHR}] rename chr${CHR} -> ${CHR}"
RENAME_TXT="${OUTDIR}/tmp/chr${CHR}/rename_chr${CHR}.txt"
echo "chr${CHR} ${CHR}" > "$RENAME_TXT"
bcftools annotate --rename-chrs "$RENAME_TXT" --threads "$THREADS" \
    -Ob -o "$RECHR" "$PHASED"
bcftools index -f "$RECHR"
rm "$RENAME_TXT"
echo "[$(date +%T)] [chr${CHR}] DONE: $RECHR"
