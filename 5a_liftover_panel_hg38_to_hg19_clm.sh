#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5a_lift19
#SBATCH --partition=mhgcp
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5a_lift19_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5a_lift19_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Lift the joint-phased hg38 panel down to hg19 with Picard LiftoverVcf.
# Per-chrom (array 1-22), lifts BOTH:
#   - the full panel:    merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf
#   - the homog 5-pop:   homog_5pop/merged_chr${CHR}.homog_5pop.bcf  (if exists)
#
# Outputs land under $PROJECT_ROOT/05_panel_hg19/:
#   all_samples/merged_chr${CHR}.hg19.bcf{,.csi}
#   homog_5pop/merged_chr${CHR}.homog_5pop.hg19.bcf{,.csi}
#   *.rejected.vcf.gz  (variants that failed liftover)
#
# Why hg19? The user's bulk RNA-seq was called against hg19; this lift
# lets them genotype-match RNA-seq variants against the reference panel
# without having to lift the RNA-seq data instead.
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
HG19_DIR="${HG19_DIR:-${PROJECT_ROOT}/05_panel_hg19}"
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg38ToHg19.over.chain.gz}"
HG19_FA="${HG19_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg19/hg19.fa}"
PICARD_XMX="${PICARD_XMX:-36g}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMPDIR_LIFT="${HG19_DIR}/tmp/chr${CHR}"
mkdir -p "$HG19_DIR/all_samples" "$HG19_DIR/homog_5pop" "$TMPDIR_LIFT"

[[ -s "$CHAIN"   ]] || { echo "ERROR: missing $CHAIN"; exit 1; }
[[ -s "$HG19_FA" ]] || { echo "ERROR: missing $HG19_FA. Override with HG19_FA=..."; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

# Helper to lift a single BCF.
do_lift () {
    local src="$1" dst="$2" reject="$3"
    if [[ -s "$dst" && -s "${dst}.csi" ]]; then
        echo "[$(date +%T)] [chr${CHR}] $dst exists, skipping"
        return 0
    fi
    [[ -s "$src" ]] || { echo "WARN: missing $src, skipping"; return 0; }
    local lift_vcf="${TMPDIR_LIFT}/$(basename "${dst%.bcf}").lifted.vcf.gz"
    echo "[$(date +%T)] [chr${CHR}] lifting $(basename "$src") -> hg19"
    picard "-Xmx${PICARD_XMX}" LiftoverVcf \
        I="$src" \
        O="$lift_vcf" \
        CHAIN="$CHAIN" \
        R="$HG19_FA" \
        REJECT="$reject" \
        RECOVER_SWAPPED_REF_ALT=true \
        WARN_ON_MISSING_CONTIG=true \
        CREATE_INDEX=false \
        MAX_RECORDS_IN_RAM=500000
    # Convert to BCF and index; sort because Picard's output may need it.
    bcftools sort -m 8G -T "${TMPDIR_LIFT}/sort" -Ob -o "$dst" "$lift_vcf"
    bcftools index "$dst"
    rm -f "$lift_vcf"
    local n=$(bcftools view "$dst" -H | wc -l)
    local nr=$(zcat "$reject" 2>/dev/null | grep -vc '^#' || echo 0)
    echo "[$(date +%T)] [chr${CHR}] lifted ${n} records, rejected ${nr}"
}

# Lift full-panel BCF.
do_lift \
    "${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf" \
    "${HG19_DIR}/all_samples/merged_chr${CHR}.hg19.bcf" \
    "${HG19_DIR}/all_samples/merged_chr${CHR}.hg19.rejected.vcf.gz"

# Lift homog-5pop BCF if present.
do_lift \
    "${HOMOG_DIR}/homog_5pop/merged_chr${CHR}.homog_5pop.bcf" \
    "${HG19_DIR}/homog_5pop/merged_chr${CHR}.homog_5pop.hg19.bcf" \
    "${HG19_DIR}/homog_5pop/merged_chr${CHR}.homog_5pop.hg19.rejected.vcf.gz"

echo "[$(date +%T)] [chr${CHR}] DONE."
