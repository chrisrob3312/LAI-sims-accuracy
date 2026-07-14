#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1c_prep_sgdp
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --time-min=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1c_prep_sgdp_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1c_prep_sgdp_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- Path C stage 1: SGDP acquisition + normalization.
# ============================================================================
# Downloads SGDP hg38 phased per-chr VCFs from the Reich lab public FTP,
# filters to samples in the anchor lists produced by
# scripts/fetch-sgdp-anchors.py (SAS + OCE + MEN + optionally AMR), splits
# multi-allelics, normalizes, and emits sgdp_prepped.chr${CHR}.bcf ready to
# feed into 1c_phase_prep_clm.sh's bcftools merge alongside HGDP+1KG and MXB.
#
# Why re-phase from scratch: RFMix assumes reference haplotypes are jointly
# phased. Appending an already-phased SGDP release to an already-phased
# HGDP+1KG+MXB panel creates cohort-specific switch-error patterns that RFMix
# would misinterpret as ancestry switches. Re-phasing after merge is the
# only way to keep LAI accuracy honest.
#
# URLs:
#   SGDP metadata (scraped by fetch-sgdp-anchors.py):
#     https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/
#   SGDP phased hg38 per-chr VCFs (Reich lab public share):
#     override SGDP_URL_TPL to match the actual layout. Typical patterns:
#       https://reichdata.hms.harvard.edu/pub/datasets/sgdp/vcf_hg38/chr${CHR}.vcf.gz
#       https://sharehost.hms.harvard.edu/reich_lab/sgdp/phased_hg38/chr${CHR}.vcf.gz
#     The default below is a placeholder -- CONFIRM before first run:
#       curl -kI "${SGDP_URL_TPL//__CHR__/22}"
#     should return HTTP 200.
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"

# Where SGDP prepped BCFs land (analogous to mxb_lifted.chr${CHR}.bcf).
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
SGDP_STAGE="${SGDP_STAGE:-${PROJECT_ROOT}/00_sgdp_raw}"   # cache for downloads

# Reich lab public FTP URL template. __CHR__ is substituted with 1..22.
# CONFIRM this URL before running; the exact layout has drifted a few times.
SGDP_URL_TPL="${SGDP_URL_TPL:-https://reichdata.hms.harvard.edu/pub/datasets/sgdp/vcf_hg38/SGDP_phased_hg38.chr__CHR__.vcf.gz}"

# Combined anchor list (SAS + OCE + MEN + AMR from fetch-sgdp-anchors.py).
# One sample ID per line. This gates who ends up in the merged panel.
SGDP_ANCHORS="${SGDP_ANCHORS:-${REPO_DIR}/reference_ids/sgdp_combined_anchors.txt}"

# hg38 reference used for bcftools norm.
HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg38/hg38.fa}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUTDIR" "$SGDP_STAGE" "${OUTDIR}/tmp/chr${CHR}"
TMPDIR="${OUTDIR}/tmp/chr${CHR}"

[[ -s "$SGDP_ANCHORS" ]] || {
    echo "ERROR: missing $SGDP_ANCHORS"
    echo "       Generate it with:"
    echo "         python3 scripts/fetch-sgdp-anchors.py --outdir reference_ids"
    echo "         cat reference_ids/sgdp_{sas,oceania,men,amr}_anchors.txt \\"
    echo "             | sort -u > reference_ids/sgdp_combined_anchors.txt"
    exit 1
}
[[ -s "$HG38_FA" ]] || { echo "ERROR: missing $HG38_FA"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ---------------------------------------------------------------------------
# 1. Download the per-chr SGDP phased VCF (idempotent -- skip if cached).
# ---------------------------------------------------------------------------
SGDP_URL="${SGDP_URL_TPL//__CHR__/${CHR}}"
RAW_VCF="${SGDP_STAGE}/SGDP_phased_hg38.chr${CHR}.vcf.gz"

if [[ ! -s "$RAW_VCF" ]]; then
    echo "[$(date +%T)] [chr${CHR}] downloading SGDP: $SGDP_URL"
    # -k because cluster's CA bundle doesn't cover the Harvard chain; the
    # payload is public academic data so integrity relies on file size sanity
    # and downstream bcftools header checks rather than TLS.
    curl -kSLf --retry 3 --retry-delay 30 -o "$RAW_VCF" "$SGDP_URL" || {
        echo "ERROR: SGDP download failed. Confirm URL:"
        echo "       $SGDP_URL"
        echo "       curl -kI '$SGDP_URL'   # should return 200"
        rm -f "$RAW_VCF"
        exit 1
    }
    # Sanity: file exists and has plausible size (>1MB for a chr VCF).
    if [[ $(stat -c%s "$RAW_VCF") -lt 1000000 ]]; then
        echo "ERROR: downloaded file suspiciously small ($(stat -c%s "$RAW_VCF") bytes)"
        head -c 200 "$RAW_VCF"
        exit 1
    fi
fi

# Fetch the tabix index if the site publishes one, else build locally.
if [[ ! -s "${RAW_VCF}.tbi" ]]; then
    if ! curl -kSLf --retry 3 -o "${RAW_VCF}.tbi" "${SGDP_URL}.tbi" 2>/dev/null; then
        echo "[$(date +%T)] [chr${CHR}] no remote .tbi, building"
        rm -f "${RAW_VCF}.tbi"
        tabix -p vcf "$RAW_VCF"
    fi
fi

# ---------------------------------------------------------------------------
# 2. Filter to anchor samples, normalize, split multi-allelics.
# ---------------------------------------------------------------------------
PREPPED="${OUTDIR}/sgdp_prepped.chr${CHR}.bcf"

if [[ -s "$PREPPED" && -s "${PREPPED}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] $PREPPED exists, skipping"
    exit 0
fi

# Sample-set filter: keep only those in $SGDP_ANCHORS that actually appear in
# the VCF. --force-samples means bcftools will not fail if some anchor IDs
# aren't in this chr's VCF (SGDP occasionally omits samples with poor
# coverage on a given chr).
KEEP_IDS="${TMPDIR}/sgdp_keep_ids.chr${CHR}.txt"
awk 'NF' "$SGDP_ANCHORS" > "$KEEP_IDS"
N_ANCHORS=$(wc -l < "$KEEP_IDS")
echo "[$(date +%T)] [chr${CHR}] anchor list: ${N_ANCHORS} SGDP samples requested"

echo "[$(date +%T)] [chr${CHR}] filter + norm + split multi-allelics"
bcftools view -S "$KEEP_IDS" --force-samples --threads "$THREADS" -Ou "$RAW_VCF" \
  | bcftools norm --threads "$THREADS" --check-ref w -f "$HG38_FA" -Ou \
  | bcftools norm --threads "$THREADS" -m -any -Ob -o "$PREPPED"
bcftools index "$PREPPED"

# ---------------------------------------------------------------------------
# 3. Report contents.
# ---------------------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "$PREPPED" | wc -l)
N_VARIANTS=$(bcftools index -n "$PREPPED")
echo "[$(date +%T)] [chr${CHR}] DONE. ${N_SAMPLES} samples, ${N_VARIANTS} variants -> $PREPPED"

# Emit a per-chr sample list for provenance / debugging.
bcftools query -l "$PREPPED" > "${PREPPED%.bcf}.samples.txt"
