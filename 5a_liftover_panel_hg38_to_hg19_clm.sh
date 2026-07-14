#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5a_lift19
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=00:15:00
#SBATCH --time=06:00:00
#SBATCH --mem=64G
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
PICARD_XMX="${PICARD_XMX:-48g}"
# St Jude / Ensembl-annotated RNA-seq uses '1..22'. UCSC hg19.fa yields 'chr1..chr22'
# after Picard. Set RENAME_TO_ENSEMBL=1 to rename output contigs 'chr1' -> '1'.
RENAME_TO_ENSEMBL="${RENAME_TO_ENSEMBL:-1}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMPDIR_LIFT="${HG19_DIR}/tmp/chr${CHR}"
mkdir -p "$HG19_DIR/all_samples" "$HG19_DIR/homog_5pop" "$TMPDIR_LIFT"

[[ -s "$CHAIN"   ]] || { echo "ERROR: missing $CHAIN"; exit 1; }
[[ -s "$HG19_FA" ]] || { echo "ERROR: missing $HG19_FA. Override with HG19_FA=..."; exit 1; }
# Picard requires a .dict sidecar. Without it, LiftoverVcf fails only AFTER the
# slow BCF -> VCF.gz conversion; check up front.
HG19_DICT="${HG19_FA%.fa}.dict"
[[ -s "$HG19_DICT" ]] || { echo "ERROR: missing Picard sequence dictionary $HG19_DICT"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# Helper to lift a single BCF.
do_lift () {
    local src="$1" dst="$2" reject="$3"
    if [[ -s "$dst" && -s "${dst}.csi" ]]; then
        echo "[$(date +%T)] [chr${CHR}] $dst exists, skipping"
        return 0
    fi
    # Fail loud on missing input -- prior silent 'return 0' hid a real 4a
    # regression that dropped all non-AMR homog samples.
    [[ -s "$src" ]] || { echo "ERROR: missing $src"; exit 2; }
    # Picard's htsjdk chokes on bcftools-emitted BCFs with
    #   "Input stream does not contain a BCF encoded file; BCF magic header info not found"
    # even when the file is a valid BCF2. Feed it a bgzipped VCF instead.
    local in_vcf="${TMPDIR_LIFT}/$(basename "${src%.bcf}").vcf.gz"
    local lift_vcf="${TMPDIR_LIFT}/$(basename "${dst%.bcf}").lifted.vcf.gz"
    echo "[$(date +%T)] [chr${CHR}] BCF -> VCF.gz for Picard: $(basename "$src")"
    bcftools view "$src" -Oz -o "$in_vcf"
    bcftools index -t "$in_vcf"
    echo "[$(date +%T)] [chr${CHR}] lifting $(basename "$src") -> hg19"
    # DISABLE_SORT=true: skip Picard's in-memory sort. On 3448-sample panels the
    # internal sort OOMs at 36g heap (java.lang.OutOfMemoryError). bcftools
    # sort below spills to disk and handles the reorder efficiently.
    # MAX_RECORDS_IN_RAM cut from 500000 -> 100000 to further bound heap use.
    picard "-Xmx${PICARD_XMX}" LiftoverVcf \
        I="$in_vcf" \
        O="$lift_vcf" \
        CHAIN="$CHAIN" \
        R="$HG19_FA" \
        REJECT="$reject" \
        RECOVER_SWAPPED_REF_ALT=true \
        WARN_ON_MISSING_CONTIG=true \
        CREATE_INDEX=false \
        DISABLE_SORT=true \
        MAX_RECORDS_IN_RAM=100000
    rm -f "$in_vcf" "${in_vcf}.tbi"
    # Optionally rename contigs to Ensembl style ('chr1' -> '1') for RNA-seq
    # collaborators using Ensembl-annotated GRCh37 (St Jude STAR pipeline).
    if [[ "${RENAME_TO_ENSEMBL}" == "1" ]]; then
        local rename_tsv="${TMPDIR_LIFT}/rename_chrs.tsv"
        if [[ ! -s "$rename_tsv" ]]; then
            : > "$rename_tsv"
            for c in {1..22} X Y MT; do echo -e "chr${c}\t${c}" >> "$rename_tsv"; done
            # UCSC hg19 uses 'chrM'; Ensembl uses 'MT'.
            echo -e "chrM\tMT" >> "$rename_tsv"
        fi
        local lift_vcf_renamed="${TMPDIR_LIFT}/$(basename "${dst%.bcf}").lifted.renamed.vcf.gz"
        bcftools annotate --rename-chrs "$rename_tsv" -Oz -o "$lift_vcf_renamed" "$lift_vcf"
        rm -f "$lift_vcf"
        lift_vcf="$lift_vcf_renamed"
    fi
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
