#!/usr/bin/env bash
# ============================================================================
# 1a_prep-mxb-liftover_clm.sh
#
# One-shot prep step before 1b_phasing-jointcall_clm.sh:
#   1. Optionally rename MXB contigs 1..22 -> chr1..chr22 (UCSC chain expects chr-prefixed)
#   2. Picard LiftoverVcf hg19 -> hg38 with RECOVER_SWAPPED_REF_ALT=true
#      (handles strand flips AND ref/alt swaps that CrossMap silently mis-encodes)
#   3. bcftools norm -m- -f hg38.fa  (split multiallelics, left-align)
#   4. bcftools +fixref --check-ref ws -- -m flip -d  (catch any remaining mismatches)
#   5. bcftools sort, index, split per-chrom
#
# Output: $OUTDIR/mxb_lifted.chr{1..22}.bcf{,.csi}
# Reject log: $OUTDIR/mxb_lifted.rejected.vcf.gz
# ============================================================================

# ----------------------------------------------------------------------------
# USER CONFIG  --  edit for your cluster / paths
# ----------------------------------------------------------------------------
#SBATCH --job-name=prep_mxb_liftover
#SBATCH --output=logs/prep_mxb_liftover_%j.out
#SBATCH --error=logs/prep_mxb_liftover_%j.err
#SBATCH --partition=medium          # ADJUST: cluster partition
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8

# Conda env (from envs/shapeit5.yml)
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Inputs
MXB_HG19="${MXB_HG19:-/storage/atkinson/shared_resources/reference/mexico_biobank/original_download/tosharemxb50wgs/mexican_50_autosomes.vcf.gz}"
REF_FA="${REF_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg19ToHg38.over.chain.gz}"

# Project root -- all outputs anchored under this so nothing collides with
# other lab members' work.
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"

# Output / scratch
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"
TMPDIR="${TMPDIR:-${OUTDIR}/tmp/mxb_prep}"

# Picard heap size (Xmx). Bump if liftover OOMs on chr1.
PICARD_XMX="${PICARD_XMX:-40g}"

# ----------------------------------------------------------------------------
set -euo pipefail

THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# 1. Detect contig naming; rename to chr-prefixed if needed
first_contig=$(bcftools view -h "$MXB_HG19" | awk -F'[<=,]' '/^##contig=<ID=/ {print $3; exit}')
echo "[$(date +%T)] MXB first contig: $first_contig"
if [[ "$first_contig" != chr* ]]; then
    echo "[$(date +%T)] Renaming numeric contigs to chr-prefixed for liftover"
    : > "${TMPDIR}/rename_to_chr.txt"
    for c in {1..22} X Y MT; do echo "$c chr$c" >> "${TMPDIR}/rename_to_chr.txt"; done
    bcftools annotate --rename-chrs "${TMPDIR}/rename_to_chr.txt" \
        --threads "$THREADS" -Oz \
        -o "${TMPDIR}/mxb.hg19.chr.vcf.gz" "$MXB_HG19"
    bcftools index --threads "$THREADS" -t "${TMPDIR}/mxb.hg19.chr.vcf.gz"
    MXB_HG19_INPUT="${TMPDIR}/mxb.hg19.chr.vcf.gz"
else
    MXB_HG19_INPUT="$MXB_HG19"
fi

# 2. Picard LiftoverVcf hg19 -> hg38
echo "[$(date +%T)] Running Picard LiftoverVcf"
picard "-Xmx${PICARD_XMX}" LiftoverVcf \
    I="$MXB_HG19_INPUT" \
    O="${TMPDIR}/mxb.hg38.lifted.vcf.gz" \
    CHAIN="$CHAIN" \
    R="$REF_FA" \
    REJECT="${OUTDIR}/mxb_lifted.rejected.vcf.gz" \
    RECOVER_SWAPPED_REF_ALT=true \
    WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=false \
    MAX_RECORDS_IN_RAM=500000

# 3. norm: split multiallelics, left-align, check ref
echo "[$(date +%T)] bcftools norm"
bcftools norm -f "$REF_FA" -m- --threads "$THREADS" \
    -Ob -o "${TMPDIR}/mxb.hg38.norm.bcf" \
    "${TMPDIR}/mxb.hg38.lifted.vcf.gz"
bcftools index --threads "$THREADS" "${TMPDIR}/mxb.hg38.norm.bcf"

# 4. fixref: catch any remaining REF/ALT mismatches (flip palindromic, drop unfixable)
echo "[$(date +%T)] bcftools +fixref"
bcftools +fixref "${TMPDIR}/mxb.hg38.norm.bcf" \
    --threads "$THREADS" -Ob \
    -o "${TMPDIR}/mxb.hg38.fixref.bcf" \
    -- -f "$REF_FA" -m flip -d
bcftools index --threads "$THREADS" "${TMPDIR}/mxb.hg38.fixref.bcf"

# 5. sort + per-chrom split
echo "[$(date +%T)] sort + per-chrom split"
bcftools sort -m 16G -T "${TMPDIR}/sort" \
    -Ob -o "${TMPDIR}/mxb.hg38.sorted.bcf" \
    "${TMPDIR}/mxb.hg38.fixref.bcf"
bcftools index --threads "$THREADS" "${TMPDIR}/mxb.hg38.sorted.bcf"

for CHR in {1..22}; do
    bcftools view -r "chr${CHR}" --threads "$THREADS" \
        -Ob -o "${OUTDIR}/mxb_lifted.chr${CHR}.bcf" \
        "${TMPDIR}/mxb.hg38.sorted.bcf"
    bcftools index --threads "$THREADS" "${OUTDIR}/mxb_lifted.chr${CHR}.bcf"
done

echo "[$(date +%T)] # Complete. Outputs in ${OUTDIR}/mxb_lifted.chr*.bcf"
