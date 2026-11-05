#!/usr/bin/env bash
#SBATCH --job-name=strip_almarri
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/strip_almarri_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/strip_almarri_chr%a_%j.err
# ============================================================================
# strip-almarri-to-gt.sh   --   array 1-22
#
# One-off: reduce the ALREADY-BUILT almarri_prepped.chr*.bcf to GT only.
#
# 1d now strips at the merge, but the existing files were produced before that
# and carry Long Ranger's full FORMAT set (PL, AD, GQ, barcodes). That is
# ~100 GB for 137 samples -- chr1 alone is 8.4 GB -- and nothing downstream
# reads any of it. Re-running 1d would redo the 137-sample prep and merge per
# chromosome; this just rewrites the finished output, which is far cheaper.
#
# WHY IT MATTERS BEYOND DISK. 1c runs `bcftools merge` across HGDP+1KG, MXB,
# SGDP and Almarri. bcftools merge emits the UNION of the inputs' FORMAT
# fields, so Almarri's extras would be written as missing for every other
# cohort's samples -- the cost is paid across the whole panel, on every
# chromosome, not just on this cohort.
#
# Safe to re-run: a chromosome already reduced is detected and skipped.
#
# Usage:
#   sbatch --array=1-22%6 scripts/strip-almarri-to-gt.sh
# ============================================================================
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"

CHR=${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-8}

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

SRC="${OUTDIR}/almarri_prepped.chr${CHR}.bcf"
[[ -s "$SRC" && -s "${SRC}.csi" ]] || { echo "ERROR: missing $SRC or its index"; exit 1; }

# Already reduced? The FORMAT header lines say so without reading the records.
# Read the header into a variable, never `bcftools view -h | grep -q`: grep -q
# exits at the first match, bcftools dies of SIGPIPE finishing the header, and
# under `set -o pipefail` the test reports the OPPOSITE of the truth.
HDR=$(bcftools view -h "$SRC") || { echo "ERROR: cannot read header of $SRC"; exit 1; }
N_FMT=$(grep -c '^##FORMAT=' <<<"$HDR" || true)
if [[ "$N_FMT" -le 1 ]]; then
    echo "[$(date +%T)] [chr${CHR}] already GT-only (${N_FMT} FORMAT field) -- skipping"
    exit 0
fi
BEFORE=$(stat -c %s "$SRC")
echo "[$(date +%T)] [chr${CHR}] ${N_FMT} FORMAT fields, $((BEFORE/1024/1024)) MB -- reducing to GT"

TMP_OUT="${SRC}.gtonly.partial"
rm -f "$TMP_OUT" "${TMP_OUT}.csi"
bcftools annotate -x 'INFO,^FORMAT/GT' --threads "$THREADS" -Ob -o "$TMP_OUT" "$SRC"
bcftools index "$TMP_OUT"

# Verify BEFORE replacing the only copy: same samples, same records.
S_OLD=$(bcftools query -l "$SRC" | wc -l)
S_NEW=$(bcftools query -l "$TMP_OUT" | wc -l)
R_OLD=$(bcftools index -n "$SRC")
R_NEW=$(bcftools index -n "$TMP_OUT")
if [[ "$S_OLD" -ne "$S_NEW" || "$R_OLD" -ne "$R_NEW" ]]; then
    echo "ERROR: chr${CHR} mismatch after stripping -- samples ${S_OLD}->${S_NEW}, records ${R_OLD}->${R_NEW}" >&2
    echo "       original left untouched at $SRC" >&2
    rm -f "$TMP_OUT" "${TMP_OUT}.csi"
    exit 2
fi

# Index first, then the BCF, so the final path never holds a pair the resume
# check would accept as finished while it is incomplete.
mv "${TMP_OUT}.csi" "${SRC}.csi"
mv "$TMP_OUT" "$SRC"
AFTER=$(stat -c %s "$SRC")
echo "[$(date +%T)] [chr${CHR}] done: $((BEFORE/1024/1024)) MB -> $((AFTER/1024/1024)) MB "\
     "($(awk -v a="$BEFORE" -v b="$AFTER" 'BEGIN{printf "%.1fx", (b>0?a/b:0)}') smaller), "\
     "${S_NEW} samples, ${R_NEW} records unchanged"
