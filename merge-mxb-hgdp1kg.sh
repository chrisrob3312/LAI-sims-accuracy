#!/usr/bin/env bash
# ============================================================================
# merge-mxb-hgdp1kg.sh
#
# Per-chrom (SLURM array 1-22):
#   1. Subset HGDP+1KG postoutlier to chr$CHR, drop kinship outliers.
#   2. bcftools merge with the lifted MXB chunk (output of prep-mxb-liftover.sh).
#      NOTE: merge (column-wise sample join), NOT concat. Different samples,
#      same/overlapping sites -> 4147 sample columns per site. Missing in one
#      panel becomes ./. and is imputed by SHAPEIT5 during phasing.
#   3. plink2 site QC (biallelic SNPs, ACGT only, geno<=0.1, dedup IDs).
#      No MAF filter at this stage -- preserves rare variants for the
#      general-purpose phased panel.
#   4. SHAPEIT5_phase_common  -> common-variant scaffold (--filter-maf 0.001).
#   5. SHAPEIT5_phase_rare    -> full phased panel (common + rare),
#      conditional on the scaffold. Useful for imputation, rare-variant work.
#   6. bcftools +fill-tags -S sample_groups.tsv -t 'AF' on the full panel
#      (per-superpop AF for downstream filtering).
#   7. Soft-union per-superpop MAF >= 0.005 filter applied POST-phase to derive
#      the LAI-ready sub-panel from the full phased output.
#   8. Rename chr$CHR -> $CHR for downstream RFMix.
#
# Run AFTER prep-mxb-liftover.sh AND make_sample_groups.sh.
#
# Outputs per chrom in $OUTDIR/:
#   merged_chr${CHR}.shapeit5_common_scaffold.bcf{,.csi}              (intermediate)
#   merged_chr${CHR}.shapeit5_full_phased.bcf{,.csi}                  (common + rare; general-purpose)
#   merged_chr${CHR}.shapeit5_full_phased.softunion_maf005.bcf{,.csi} (LAI-ready)
#   merged_chr${CHR}.shapeit5_full_phased.softunion_maf005.rechr.bcf{,.csi}
#       (LAI-ready, chr$CHR -> $CHR for RFMix v1)
# ============================================================================

# ----------------------------------------------------------------------------
# USER CONFIG  --  edit for your cluster / paths
# ----------------------------------------------------------------------------
#SBATCH --job-name=merge_mxb_hgdp1kg
#SBATCH --output=logs/merge_mxb_hgdp1kg_chr%a_%j.out
#SBATCH --error=logs/merge_mxb_hgdp1kg_chr%a_%j.err
#SBATCH --partition=long            # ADJUST: cluster partition
#SBATCH --time=120:00:00            # phase_common + phase_rare on chr1 ~ 48-72h at 32 cores
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22

# Conda env (from envs/shapeit5.yml)
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Inputs
HGDP1KG="${HGDP1KG:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/archived/TGP_HGDP_hg38/filtered/hgdp_tgp_filtered_postoutlier.vcf.gz}"
OUTLIERS="${OUTLIERS:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt}"
REF_FA="${REF_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"

# Sample -> superpop TSV produced by make_sample_groups.sh
SAMPLE_GROUPS="${SAMPLE_GROUPS:-./sample_groups.tsv}"

# Output / scratch
OUTDIR="${OUTDIR:-merged_mxb_hgdp1kg}"
LOGDIR="${LOGDIR:-logs}"

# SHAPEIT5 thresholds
COMMON_MAF="${COMMON_MAF:-0.001}"   # variants below this go to phase_rare; above -> scaffold
LAI_MAF="${LAI_MAF:-0.005}"         # soft-union threshold for LAI-ready sub-panel

# Skip phase_rare to save compute (e.g. LAI-only run)?  Default: run it.
RUN_PHASE_RARE="${RUN_PHASE_RARE:-1}"

# ----------------------------------------------------------------------------
set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-32}

GMAP="${GMAP_DIR}/chr${CHR}.b38.gmap.gz"
TMPDIR="${OUTDIR}/tmp/chr${CHR}"
MXB_LIFTED="${OUTDIR}/mxb_lifted.chr${CHR}.bcf"
mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

[[ -s "$MXB_LIFTED"     ]] || { echo "ERROR: missing $MXB_LIFTED. Run prep-mxb-liftover.sh first."; exit 1; }
[[ -s "$SAMPLE_GROUPS"  ]] || { echo "ERROR: missing $SAMPLE_GROUPS. Run make_sample_groups.sh first."; exit 1; }
[[ -s "$GMAP"           ]] || { echo "ERROR: missing $GMAP."; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
MERGED="${TMPDIR}/merged_chr${CHR}.bcf"
QCED_PREFIX="${TMPDIR}/merged_chr${CHR}.qced"
QCED="${QCED_PREFIX}.bcf"
SCAFFOLD="${OUTDIR}/merged_chr${CHR}.shapeit5_common_scaffold.bcf"
FULL_PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_full_phased.bcf"
FULL_TAGGED="${TMPDIR}/merged_chr${CHR}.full_phased.tagged.bcf"
LAI_READY="${OUTDIR}/merged_chr${CHR}.shapeit5_full_phased.softunion_maf005.bcf"
LAI_RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_full_phased.softunion_maf005.rechr.bcf"

# 1. Subset HGDP+1KG to chr$CHR; drop kinship outliers
echo "[$(date +%T)] [chr${CHR}] subset HGDP+1KG, drop outliers"
bcftools view -r "chr${CHR}" -S "^${OUTLIERS}" --force-samples \
    --threads "$THREADS" -Ob -o "$HGDP1KG_CHR" "$HGDP1KG"
bcftools index --threads "$THREADS" "$HGDP1KG_CHR"

# 2. Merge HGDP+1KG (chr$CHR) with lifted MXB (chr$CHR) -- column-wise sample join
echo "[$(date +%T)] [chr${CHR}] bcftools merge"
bcftools merge --threads "$THREADS" \
    -Ob -o "$MERGED" \
    "$HGDP1KG_CHR" "$MXB_LIFTED"
bcftools index --threads "$THREADS" "$MERGED"

# 3. plink2 site QC (no MAF filter -- rare variants preserved for general-purpose panel)
echo "[$(date +%T)] [chr${CHR}] plink2 site QC (no MAF filter)"
plink2 --bcf "$MERGED" \
       --set-missing-var-ids '@:#[b38]' \
       --rm-dup exclude-all \
       --geno 0.1 \
       --max-alleles 2 --snps-only just-acgt \
       --export bcf \
       --threads "$THREADS" \
       --out "$QCED_PREFIX"
bcftools index --threads "$THREADS" "$QCED"

# 4. SHAPEIT5 phase_common -- common-variant scaffold
echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_common (--filter-maf ${COMMON_MAF})"
SHAPEIT5_phase_common \
    --input "$QCED" \
    --map "$GMAP" \
    --region "chr${CHR}" \
    --output "$SCAFFOLD" \
    --thread "$THREADS" \
    --filter-maf "$COMMON_MAF"
bcftools index --threads "$THREADS" "$SCAFFOLD"

# 5. SHAPEIT5 phase_rare -- full panel (common + rare), conditional on scaffold
if [[ "$RUN_PHASE_RARE" == "1" ]]; then
    echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_rare"
    SHAPEIT5_phase_rare \
        --input "$QCED" \
        --scaffold "$SCAFFOLD" \
        --map "$GMAP" \
        --input-region "chr${CHR}" \
        --scaffold-region "chr${CHR}" \
        --output "$FULL_PHASED" \
        --thread "$THREADS"
    bcftools index --threads "$THREADS" "$FULL_PHASED"
else
    echo "[$(date +%T)] [chr${CHR}] RUN_PHASE_RARE=0 -- skipping phase_rare; full panel = scaffold"
    cp "$SCAFFOLD" "$FULL_PHASED"
    cp "${SCAFFOLD}.csi" "${FULL_PHASED}.csi"
fi

# 6. Per-superpop AF tags on the full phased panel
echo "[$(date +%T)] [chr${CHR}] +fill-tags per-superpop AF (post-phase)"
bcftools +fill-tags "$FULL_PHASED" --threads "$THREADS" \
    -Ob -o "$FULL_TAGGED" \
    -- -S "$SAMPLE_GROUPS" -t 'AF'
bcftools index --threads "$THREADS" "$FULL_TAGGED"

# 7. Soft-union post-filter -> LAI-ready sub-panel.
#    Keep site if MAF >= LAI_MAF in AT LEAST ONE superpop. Equivalently, exclude
#    sites where AF<LAI_MAF || AF>(1-LAI_MAF) in ALL superpops.
echo "[$(date +%T)] [chr${CHR}] soft-union per-superpop MAF >= ${LAI_MAF} (LAI sub-panel)"
HI=$(awk -v m="$LAI_MAF" 'BEGIN{printf "%.6f", 1-m}')
bcftools view "$FULL_TAGGED" \
    -e "(INFO/AF_AFR<${LAI_MAF} || INFO/AF_AFR>${HI}) && (INFO/AF_AMR<${LAI_MAF} || INFO/AF_AMR>${HI}) && (INFO/AF_EUR<${LAI_MAF} || INFO/AF_EUR>${HI}) && (INFO/AF_EAS<${LAI_MAF} || INFO/AF_EAS>${HI}) && (INFO/AF_SAS<${LAI_MAF} || INFO/AF_SAS>${HI}) && (INFO/AF_CSA<${LAI_MAF} || INFO/AF_CSA>${HI}) && (INFO/AF_OCE<${LAI_MAF} || INFO/AF_OCE>${HI}) && (INFO/AF_MEN<${LAI_MAF} || INFO/AF_MEN>${HI})" \
    --threads "$THREADS" -Ob -o "$LAI_READY"
bcftools index --threads "$THREADS" "$LAI_READY"

# 8. Strip "chr" prefix on the LAI panel (RFMix v1 expects bare numeric contigs)
echo "[$(date +%T)] [chr${CHR}] rename chr${CHR} -> ${CHR} (LAI sub-panel)"
echo "chr${CHR} ${CHR}" > "${TMPDIR}/rename_chr${CHR}.txt"
bcftools annotate --rename-chrs "${TMPDIR}/rename_chr${CHR}.txt" \
    --threads "$THREADS" -Ob -o "$LAI_RECHR" "$LAI_READY"
bcftools index --threads "$THREADS" "$LAI_RECHR"
rm "${TMPDIR}/rename_chr${CHR}.txt"

echo "[$(date +%T)] [chr${CHR}] # Complete."
echo "[$(date +%T)] [chr${CHR}] # Full phased panel  : $FULL_PHASED"
echo "[$(date +%T)] [chr${CHR}] # LAI-ready sub-panel: $LAI_READY"
echo "[$(date +%T)] [chr${CHR}] # LAI-ready (rechr)  : $LAI_RECHR"
