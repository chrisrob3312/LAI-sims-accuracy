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
#   3. bcftools +fill-tags -S sample_groups.tsv -t 'AF'  (per-superpop AF).
#   4. Soft-union per-superpop MAF >= 0.005 in any of {AFR, AMR, EUR, EAS, SAS,
#      CSA, OCE, MEN}. Replaces the global `--maf 0.005` from filter1.
#   5. plink2-style site QC: biallelic SNPs, ACGT only, geno<=0.1, dedup IDs.
#   6. SHAPEIT5_phase_common joint phase against the SHAPEIT4-format hg38 gmap.
#   7. Rename chr$CHR -> $CHR for downstream RFMix.
#
# Run AFTER prep-mxb-liftover.sh AND make_sample_groups.sh.
#
# Outputs per chrom in $OUTDIR/:
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf{,.csi}
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf{,.csi}
# ============================================================================

# ----------------------------------------------------------------------------
# USER CONFIG  --  edit for your cluster / paths
# ----------------------------------------------------------------------------
#SBATCH --job-name=merge_mxb_hgdp1kg
#SBATCH --output=logs/merge_mxb_hgdp1kg_chr%a_%j.out
#SBATCH --error=logs/merge_mxb_hgdp1kg_chr%a_%j.err
#SBATCH --partition=long            # ADJUST: cluster partition
#SBATCH --time=96:00:00
#SBATCH --mem=96G
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
TAGGED="${TMPDIR}/merged_chr${CHR}.tagged.bcf"
SOFTUNION="${TMPDIR}/merged_chr${CHR}.softunion.bcf"
QCED_PREFIX="${TMPDIR}/merged_chr${CHR}.qced"
QCED="${QCED_PREFIX}.bcf"
PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf"

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

# 3. Per-superpop AF tags
echo "[$(date +%T)] [chr${CHR}] +fill-tags per-superpop AF"
bcftools +fill-tags "$MERGED" --threads "$THREADS" \
    -Ob -o "$TAGGED" \
    -- -S "$SAMPLE_GROUPS" -t 'AF'
bcftools index --threads "$THREADS" "$TAGGED"

# 4. Soft-union: keep site if MAF >= 0.005 in AT LEAST ONE superpop.
#    Equivalently, exclude sites where AF<0.005 OR AF>0.995 in ALL superpops.
#    +fill-tags emits AF_<grp>; we test both tails since "MAF" = min(AF, 1-AF).
echo "[$(date +%T)] [chr${CHR}] soft-union per-superpop MAF filter"
bcftools view "$TAGGED" \
    --min-alleles 2 --max-alleles 2 --types snps \
    -e '(INFO/AF_AFR<0.005 || INFO/AF_AFR>0.995) && (INFO/AF_AMR<0.005 || INFO/AF_AMR>0.995) && (INFO/AF_EUR<0.005 || INFO/AF_EUR>0.995) && (INFO/AF_EAS<0.005 || INFO/AF_EAS>0.995) && (INFO/AF_SAS<0.005 || INFO/AF_SAS>0.995) && (INFO/AF_CSA<0.005 || INFO/AF_CSA>0.995) && (INFO/AF_OCE<0.005 || INFO/AF_OCE>0.995) && (INFO/AF_MEN<0.005 || INFO/AF_MEN>0.995)' \
    --threads "$THREADS" -Ob -o "$SOFTUNION"
bcftools index --threads "$THREADS" "$SOFTUNION"

# 5. plink2 site QC: dedup, ACGT-only, geno<=0.1 (matches legacy filter1 chain)
echo "[$(date +%T)] [chr${CHR}] plink2 site QC"
plink2 --bcf "$SOFTUNION" \
       --set-missing-var-ids '@:#[b38]' \
       --rm-dup exclude-all \
       --geno 0.1 \
       --max-alleles 2 --snps-only just-acgt \
       --export bcf \
       --threads "$THREADS" \
       --out "$QCED_PREFIX"
bcftools index --threads "$THREADS" "$QCED"

# 6. SHAPEIT5 phase_common (joint re-phase: HGDP+1KG + MXB together)
#    No --reference: phasing from scratch so MXB and HGDP+1KG end up on identical
#    haplotype scaffolds -- avoids the "wonky" rephase-then-merge pattern from
#    previous attempts.
echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_common"
SHAPEIT5_phase_common \
    --input "$QCED" \
    --map "$GMAP" \
    --region "chr${CHR}" \
    --output "$PHASED" \
    --thread "$THREADS" \
    --filter-maf 0.001
bcftools index --threads "$THREADS" "$PHASED"

# 7. Strip "chr" prefix (RFMix v1 expects bare numeric contigs)
echo "[$(date +%T)] [chr${CHR}] rename chr${CHR} -> ${CHR}"
echo "chr${CHR} ${CHR}" > "${TMPDIR}/rename_chr${CHR}.txt"
bcftools annotate --rename-chrs "${TMPDIR}/rename_chr${CHR}.txt" \
    --threads "$THREADS" -Ob -o "$RECHR" "$PHASED"
bcftools index --threads "$THREADS" "$RECHR"
rm "${TMPDIR}/rename_chr${CHR}.txt"

echo "[$(date +%T)] [chr${CHR}] # Complete. Phased panel: $PHASED"
echo "[$(date +%T)] [chr${CHR}] # Renamed for RFMix:    $RECHR"
