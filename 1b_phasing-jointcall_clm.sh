#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_phase
#SBATCH --partition=mhgcp,atkinson
#SBATCH --time=72:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_phase_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_phase_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Per-chrom (SLURM array 1-22):
#   1. Subset HGDP+1KG postoutlier to chr$CHR, drop kinship outliers.
#   2. bcftools merge with the lifted MXB chunk (output of 1a_prep-mxb-liftover_clm.sh).
#      NOTE: merge (column-wise sample join), NOT concat. Different samples,
#      same/overlapping sites -> 4147 sample columns per site. Missing in one
#      panel becomes ./. and is imputed by SHAPEIT5 during phasing.
#   3. bcftools +fill-tags -S sample_groups.tsv -t 'AF'  (per-superpop AF).
#   4. Soft-union per-superpop MAF >= LAI_MAF in any of {AFR, AMR, EUR, EAS,
#      SAS, CSA, OCE, MEN}. Replaces the global `--maf 0.005` from filter1.
#   5. plink2 site QC: biallelic SNPs, ACGT only, geno<=0.1, dedup IDs.
#   6. SHAPEIT5_phase_common joint phase against the SHAPEIT4-format hg38 gmap.
#   7. Rename chr$CHR -> $CHR for downstream RFMix.
#
# RARE VARIANT PHASING is OFF by default. Set RUN_PHASE_RARE=1 to enable
# phase_rare (adds ~2x wall time per chr).
#
# Run AFTER 1a_prep-mxb-liftover_clm.sh AND make_sample_groups.sh.
#
# Outputs per chrom in $OUTDIR/:
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf{,.csi}
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf{,.csi}
#       (chr$CHR -> $CHR for RFMix v1)
#   merged_chr${CHR}.shapeit5_full_phased.bcf{,.csi}        [only if RUN_PHASE_RARE=1]
# ============================================================================

# ----------------------------------------------------------------------------
# Variable configuration  --  paths and tunables (override via env vars at submit time)
# ----------------------------------------------------------------------------
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Project root -- all outputs anchored here so nothing collides with other
# lab members' work.
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"

# Inputs
HGDP1KG="${HGDP1KG:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/archived/TGP_HGDP_hg38/filtered/hgdp_tgp_filtered_postoutlier.vcf.gz}"
OUTLIERS="${OUTLIERS:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt}"
REF_FA="${REF_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"

# Sample -> superpop TSV produced by make_sample_groups.sh
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"

# Outputs / scratch
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"

# Filter thresholds
LAI_MAF="${LAI_MAF:-0.005}"   # soft-union per-superpop MAF for LAI panel

# Optional: rare-variant phasing (off by default; not needed for LAI -- TOPMed
# / All-of-Us covers imputation)
RUN_PHASE_RARE="${RUN_PHASE_RARE:-0}"

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-32}

GMAP="${GMAP_DIR}/chr${CHR}.b38.gmap.gz"
TMPDIR="${OUTDIR}/tmp/chr${CHR}"
MXB_LIFTED="${OUTDIR}/mxb_lifted.chr${CHR}.bcf"
mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

[[ -s "$MXB_LIFTED"     ]] || { echo "ERROR: missing $MXB_LIFTED. Run 1a_prep-mxb-liftover_clm.sh first."; exit 1; }
[[ -s "$SAMPLE_GROUPS"  ]] || { echo "ERROR: missing $SAMPLE_GROUPS. Run make_sample_groups.sh first."; exit 1; }
[[ -s "$GMAP"           ]] || { echo "ERROR: missing $GMAP."; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
MERGED="${TMPDIR}/merged_chr${CHR}.bcf"
TAGGED="${TMPDIR}/merged_chr${CHR}.tagged.bcf"
SOFTUNION="${TMPDIR}/merged_chr${CHR}.softunion.bcf"
QCED_PREFIX="${TMPDIR}/merged_chr${CHR}.qced"
QCED="${QCED_PREFIX}.bcf"
PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf"
FULL_PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_full_phased.bcf"

# 1. Subset HGDP+1KG to chr$CHR; drop kinship outliers
#    related_outliers.txt is 2 cols (super-pop, sample_id) -- bcftools -S
#    expects one ID per line, so extract col 2 to a temp file first.
echo "[$(date +%T)] [chr${CHR}] subset HGDP+1KG, drop outliers"
OUTLIERS_IDS="${TMPDIR}/related_outlier_ids.txt"
awk '{print $2}' "$OUTLIERS" > "$OUTLIERS_IDS"
bcftools view -r "chr${CHR}" -S "^${OUTLIERS_IDS}" --force-samples \
    --threads "$THREADS" -Ob -o "$HGDP1KG_CHR" "$HGDP1KG"
bcftools index --threads "$THREADS" "$HGDP1KG_CHR"

# 2. Merge HGDP+1KG (chr$CHR) with lifted MXB (chr$CHR) -- column-wise sample join
echo "[$(date +%T)] [chr${CHR}] bcftools merge"
bcftools merge --threads "$THREADS" \
    -Ob -o "$MERGED" \
    "$HGDP1KG_CHR" "$MXB_LIFTED"
bcftools index --threads "$THREADS" "$MERGED"

# 3. plink2 site QC -- biallelic SNPs, ACGT only, dedup, missing-var-ids.
#    NO --geno here: global call rate is the wrong filter for LAI because RFMix
#    only uses AFR/EUR/AMR samples (not all 8 superpops). A global --geno 0.1
#    can drop sites with adequate calls in the 3 LAI pops just because of
#    missingness in unused pops or because MXB samples (50/3454) are missing.
#    Per-superpop missingness filter is applied in step 4 instead.
echo "[$(date +%T)] [chr${CHR}] plink2 site QC (biallelic SNPs, ACGT, dedup)"
plink2 --bcf "$MERGED" \
       --set-missing-var-ids '@:#[b38]' \
       --rm-dup exclude-all \
       --max-alleles 2 --snps-only just-acgt \
       --export bcf \
       --threads "$THREADS" \
       --out "$QCED_PREFIX"
bcftools index --threads "$THREADS" "$QCED"

# 4. Per-superpop AF + F_MISSING tags on the QC'd merged BCF.
echo "[$(date +%T)] [chr${CHR}] +fill-tags per-superpop AF + F_MISSING"
bcftools +fill-tags "$QCED" --threads "$THREADS" \
    -Ob -o "$TAGGED" \
    -- -S "$SAMPLE_GROUPS" -t 'AF,F_MISSING'
bcftools index --threads "$THREADS" "$TAGGED"

# 5. Combined site filter:
#    a) Drop sites with >10% missingness in any of the 3 LAI pops (AFR/EUR/AMR)
#       -- these can't reliably contribute to LAI inference. Other superpops
#       (CSA/EAS/OCE/MEN) aren't used by RFMix for the panels we run, so we
#       don't filter on their missingness.
#    b) Soft-union: keep sites with MAF >= LAI_MAF in AT LEAST ONE superpop.
#       Equivalently, exclude sites where AF<LAI_MAF || AF>(1-LAI_MAF) in ALL.
#    NOTE: gnomAD lumps 1KG-SAS into CSA, so the 7 superpops actually present
#    are AFR/AMR/EUR/EAS/CSA/OCE/MEN.
#    SHAPEIT5_phase_common will impute the small remaining missingness during
#    phasing.
echo "[$(date +%T)] [chr${CHR}] per-LAI-pop F_MISSING + soft-union MAF filter"
HI=$(awk -v m="$LAI_MAF" 'BEGIN{printf "%.6f", 1-m}')
bcftools view "$TAGGED" \
    -e "(INFO/F_MISSING_AFR > 0.1 || INFO/F_MISSING_EUR > 0.1 || INFO/F_MISSING_AMR > 0.1) || \
        ((INFO/AF_AFR<${LAI_MAF} || INFO/AF_AFR>${HI}) && \
         (INFO/AF_AMR<${LAI_MAF} || INFO/AF_AMR>${HI}) && \
         (INFO/AF_EUR<${LAI_MAF} || INFO/AF_EUR>${HI}) && \
         (INFO/AF_EAS<${LAI_MAF} || INFO/AF_EAS>${HI}) && \
         (INFO/AF_CSA<${LAI_MAF} || INFO/AF_CSA>${HI}) && \
         (INFO/AF_OCE<${LAI_MAF} || INFO/AF_OCE>${HI}) && \
         (INFO/AF_MEN<${LAI_MAF} || INFO/AF_MEN>${HI}))" \
    --threads "$THREADS" -Ob -o "$SOFTUNION"
bcftools index --threads "$THREADS" "$SOFTUNION"

# 6. SHAPEIT5 phase_common (joint re-phase: HGDP+1KG + MXB together)
echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_common"
SHAPEIT5_phase_common \
    --input "$SOFTUNION" \
    --map "$GMAP" \
    --region "chr${CHR}" \
    --output "$PHASED" \
    --thread "$THREADS" \
    --filter-maf 0.001
bcftools index --threads "$THREADS" "$PHASED"

# 6b. (optional) phase_rare on the un-MAF-filtered QC'd input, conditional on
#     the soft-union scaffold above. Off by default since imputation here uses
#     TOPMed/AoU. The QCED file (already biallelic-SNP/ACGT/dedup/geno<=0.1
#     filtered) is the right input -- it contains the rare variants too.
if [[ "$RUN_PHASE_RARE" == "1" ]]; then
    echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_rare"
    SHAPEIT5_phase_rare \
        --input "$QCED" \
        --scaffold "$PHASED" \
        --map "$GMAP" \
        --input-region "chr${CHR}" \
        --scaffold-region "chr${CHR}" \
        --output "$FULL_PHASED" \
        --thread "$THREADS"
    bcftools index --threads "$THREADS" "$FULL_PHASED"
fi

# 7. Strip "chr" prefix on the LAI panel (RFMix v1 expects bare numeric contigs)
echo "[$(date +%T)] [chr${CHR}] rename chr${CHR} -> ${CHR}"
echo "chr${CHR} ${CHR}" > "${TMPDIR}/rename_chr${CHR}.txt"
bcftools annotate --rename-chrs "${TMPDIR}/rename_chr${CHR}.txt" \
    --threads "$THREADS" -Ob -o "$RECHR" "$PHASED"
bcftools index --threads "$THREADS" "$RECHR"
rm "${TMPDIR}/rename_chr${CHR}.txt"

echo "[$(date +%T)] [chr${CHR}] # Complete."
echo "[$(date +%T)] [chr${CHR}] # LAI panel (chr*)   : $PHASED"
echo "[$(date +%T)] [chr${CHR}] # LAI panel (rechr)  : $RECHR"
[[ "$RUN_PHASE_RARE" == "1" ]] && \
    echo "[$(date +%T)] [chr${CHR}] # Full phased panel  : $FULL_PHASED"
