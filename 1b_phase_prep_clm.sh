#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_prep
#SBATCH --partition=mhgcp
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_prep_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_prep_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Stage 1 of the chunked-phasing pipeline. Per-chrom (array 1-22):
#   1. Resolve per-chrom HGDP+1KG filter1 BCF, drop kinship outliers
#   2. bcftools merge with the lifted MXB chunk (output of 1a)
#   3. plink2 site QC (biallelic SNPs, ACGT, dedup)
#   4. Per-superpop pre-filter -> soft-union site list
#   5. Subset full BCF to soft-union sites -> merged_chr${CHR}.softunion.bcf
# Stops before SHAPEIT5; the softunion BCF is the input for stage 2 (chunks).
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/hgdp1kgp_chr__CHR__.shapeit5_phased.filter1_SNP_maf005.vcf.gz}"
OUTLIERS="${OUTLIERS:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"
LAI_MAF="${LAI_MAF:-0.005}"
MIN_SUBPOP_N="${MIN_SUBPOP_N:-10}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-16}
TMPDIR="${OUTDIR}/tmp/chr${CHR}"
MXB_LIFTED="${OUTDIR}/mxb_lifted.chr${CHR}.bcf"
mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

[[ -s "$MXB_LIFTED"    ]] || { echo "ERROR: missing $MXB_LIFTED"; exit 1; }
[[ -s "$SAMPLE_GROUPS" ]] || { echo "ERROR: missing $SAMPLE_GROUPS"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

MERGED="${TMPDIR}/merged_chr${CHR}.bcf"
SOFTUNION="${OUTDIR}/merged_chr${CHR}.softunion.bcf"
QCED_PREFIX="${TMPDIR}/merged_chr${CHR}.qced"
QCED="${QCED_PREFIX}.bcf"

# Skip if final softunion already exists (idempotent re-submit).
if [[ -s "$SOFTUNION" && -s "${SOFTUNION}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] softunion BCF already exists, skipping prep"
    exit 0
fi

# 1. Per-chr HGDP+1KG: drop outliers, chr-prefix
HGDP1KG_CHR_SRC="${HGDP1KG_PATTERN//__CHR__/$CHR}"
[[ -s "$HGDP1KG_CHR_SRC" ]] || { echo "ERROR: missing $HGDP1KG_CHR_SRC"; exit 1; }
echo "[$(date +%T)] [chr${CHR}] drop kinship outliers"
OUTLIERS_IDS="${TMPDIR}/related_outlier_ids.txt"
awk '{print $2}' "$OUTLIERS" > "$OUTLIERS_IDS"
HGDP1KG_PREP="${TMPDIR}/hgdp1kg_chr${CHR}.prep.bcf"
bcftools view -S "^${OUTLIERS_IDS}" --force-samples --threads "$THREADS" \
    -Ob -o "$HGDP1KG_PREP" "$HGDP1KG_CHR_SRC"
bcftools index "$HGDP1KG_PREP"
NHGDP=$(bcftools view "$HGDP1KG_PREP" -H | wc -l)
NSAMP=$(bcftools query -l "$HGDP1KG_PREP" | wc -l)
echo "[$(date +%T)] [chr${CHR}] HGDP+1KG: ${NSAMP} samples, ${NHGDP} records"
[[ "$NHGDP" -gt 0 ]] || { echo "ERROR: empty HGDP+1KG chr${CHR}"; exit 1; }

if bcftools view -h "$HGDP1KG_PREP" | grep -q "^##contig=<ID=chr${CHR}[,>]"; then
    HGDP1KG_CHR="$HGDP1KG_PREP"
elif bcftools view -h "$HGDP1KG_PREP" | grep -q "^##contig=<ID=${CHR}[,>]"; then
    RENAME_TXT="${TMPDIR}/rename_to_chr.txt"
    : > "$RENAME_TXT"
    for c in {1..22} X Y MT; do echo "$c chr$c" >> "$RENAME_TXT"; done
    HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
    bcftools annotate --rename-chrs "$RENAME_TXT" --threads "$THREADS" \
        -Ob -o "$HGDP1KG_CHR" "$HGDP1KG_PREP"
    bcftools index "$HGDP1KG_CHR"
else
    echo "ERROR: HGDP+1KG has neither chr${CHR} nor ${CHR}"; exit 1
fi

# 2. Merge with MXB
echo "[$(date +%T)] [chr${CHR}] bcftools merge"
bcftools merge --threads "$THREADS" -Ob -o "$MERGED" "$HGDP1KG_CHR" "$MXB_LIFTED"
bcftools index "$MERGED"

# 3. plink2 site QC
echo "[$(date +%T)] [chr${CHR}] plink2 site QC"
plink2 --bcf "$MERGED" --set-missing-var-ids '@:#[b38]' --rm-dup exclude-all \
       --max-alleles 2 --snps-only just-acgt --output-chr chrM \
       --export bcf --threads "$THREADS" --out "$QCED_PREFIX"
bcftools index "$QCED"
bcftools view -h "$QCED" | grep -q "^##contig=<ID=chr${CHR}[,>]" || \
    { echo "ERROR: QCED lost chr prefix"; exit 1; }

# 4. Per-superpop pre-filter -> soft-union site list
echo "[$(date +%T)] [chr${CHR}] per-pop soft-union"
VCF_SAMPLES="${TMPDIR}/vcf_samples.txt"
bcftools query -l "$QCED" > "$VCF_SAMPLES"
SITELIST="${TMPDIR}/softunion_sites.tsv"
: > "$SITELIST"
MAF_EXPR="INFO/AN>0 && INFO/AC>=${LAI_MAF}*INFO/AN && (INFO/AN-INFO/AC)>=${LAI_MAF}*INFO/AN"
for GROUP in $(awk '{print $2}' "$SAMPLE_GROUPS" | sort -u); do
    GK="${TMPDIR}/${GROUP}_keep.txt"
    awk -v g="$GROUP" '$2==g {print $1}' "$SAMPLE_GROUPS" \
        | grep -xFf "$VCF_SAMPLES" > "$GK" || true
    N=$(wc -l < "$GK")
    if (( N < MIN_SUBPOP_N )); then
        echo "[$(date +%T)] [chr${CHR}]   skip ${GROUP}: ${N}<${MIN_SUBPOP_N}"
        continue
    fi
    bcftools view "$QCED" -S "$GK" --force-samples --threads "$THREADS" -Ou \
      | bcftools +fill-tags --threads "$THREADS" -Ou -- -t 'AC,AN' \
      | bcftools view -i "$MAF_EXPR" --threads "$THREADS" -Ou \
      | bcftools query -f '%CHROM\t%POS\n' >> "$SITELIST"
    echo "[$(date +%T)] [chr${CHR}]   ${GROUP}: ${N} samples"
done
SORTED="${TMPDIR}/softunion_sites.sorted.tsv"
sort -k1,1 -k2,2n -u "$SITELIST" > "$SORTED"
NUNION=$(wc -l < "$SORTED")
echo "[$(date +%T)] [chr${CHR}] soft-union sites: ${NUNION}"
[[ "$NUNION" -gt 0 ]] || { echo "ERROR: empty soft-union"; exit 1; }
bgzip -f "$SORTED"
tabix -s1 -b2 -e2 -f "${SORTED}.gz"

# 5. Subset full BCF to soft-union sites -> SOFTUNION
echo "[$(date +%T)] [chr${CHR}] subset full BCF to soft-union sites"
bcftools view "$QCED" -T "${SORTED}.gz" --threads "$THREADS" -Ob -o "$SOFTUNION"
bcftools index "$SOFTUNION"
NSOFT=$(bcftools view "$SOFTUNION" -H | wc -l)
echo "[$(date +%T)] [chr${CHR}] DONE prep. softunion: ${NSOFT} records, ${SOFTUNION}"
