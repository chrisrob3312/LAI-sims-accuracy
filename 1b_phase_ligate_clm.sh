#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_ligate
#SBATCH --partition=mhgcp
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
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

# Skip if already done.
if [[ -s "$RECHR" && -s "${RECHR}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] rechr BCF already exists, skipping"
    exit 0
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
