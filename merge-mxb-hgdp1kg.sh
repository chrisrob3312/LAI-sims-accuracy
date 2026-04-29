#!/usr/bin/env bash
#SBATCH --partition=long
#SBATCH --time=96:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --job-name=merge_mxb_hgdp1kg
#SBATCH --output=logs/merge_mxb_hgdp1kg_chr%a_%j.out
#
# Per-chrom (SLURM array 1-22):
#   1. Subset HGDP+1KG postoutlier to chr$CHR, drop kinship outliers (matches
#      SHAPEIT5 panel sample set).
#   2. bcftools merge with the lifted MXB chunk (output of prep-mxb-liftover.sh).
#   3. bcftools +fill-tags -S sample_groups.tsv -t 'AF'  (per-superpop AF).
#   4. Soft-union per-superpop MAF >= 0.005 in any of {AFR, AMR, EUR, EAS, SAS,
#      CSA, OCE, MEN}. Replaces the global `--maf 0.005` from filter1.
#   5. plink2-style site QC: biallelic SNPs, ACGT only, geno<=0.1, dedup IDs.
#   6. SHAPEIT5_phase_common joint phase against the SHAPEIT4-format hg38 gmap.
#   7. Rename chr$CHR -> $CHR for downstream RFMix.
#
# Run AFTER prep-mxb-liftover.sh and AFTER make_sample_groups.sh has produced
# sample_groups.tsv (path via env var SAMPLE_GROUPS, default ./sample_groups.tsv).
#
# Outputs per chrom in merged_mxb_hgdp1kg/:
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf{,.csi}
#   merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf{,.csi}

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?must be run as SLURM array job}
THREADS=${SLURM_CPUS_PER_TASK:-32}

# --- Inputs ---
HGDP1KG="/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/archived/TGP_HGDP_hg38/filtered/hgdp_tgp_filtered_postoutlier.vcf.gz"
OUTLIERS="/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt"
REF_FA="/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta"
GMAP="/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38/chr${CHR}.b38.gmap.gz"

SAMPLE_GROUPS="${SAMPLE_GROUPS:-./sample_groups.tsv}"

# --- Output ---
OUTDIR="merged_mxb_hgdp1kg"
TMPDIR="${OUTDIR}/tmp/chr${CHR}"
MXB_LIFTED="${OUTDIR}/mxb_lifted.chr${CHR}.bcf"
mkdir -p "$OUTDIR" "$TMPDIR" logs

[[ -s "$MXB_LIFTED"     ]] || { echo "ERROR: missing $MXB_LIFTED. Run prep-mxb-liftover.sh first."; exit 1; }
[[ -s "$SAMPLE_GROUPS"  ]] || { echo "ERROR: missing $SAMPLE_GROUPS. Run make_sample_groups.sh first."; exit 1; }

# --- Conda env ---
module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate shapeit5

HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
MERGED="${TMPDIR}/merged_chr${CHR}.bcf"
TAGGED="${TMPDIR}/merged_chr${CHR}.tagged.bcf"
QCED="${TMPDIR}/merged_chr${CHR}.qced.bcf"
PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf"

# 1. Subset HGDP+1KG to chr$CHR; drop kinship outliers
echo "[$(date +%T)] [chr${CHR}] subset HGDP+1KG, drop outliers"
bcftools view -r "chr${CHR}" -S "^${OUTLIERS}" --force-samples \
    --threads "$THREADS" -Ob -o "$HGDP1KG_CHR" "$HGDP1KG"
bcftools index --threads "$THREADS" "$HGDP1KG_CHR"

# 2. Merge HGDP+1KG (chr$CHR) with lifted MXB (chr$CHR)
#    --missing-to-ref=NO: missing in one panel stays ./. so SHAPEIT5 imputes it
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

# 4 + 5. Soft-union MAF filter + biallelic-SNP + ACGT + dedup + missingness
#   MAF_<grp> = min(AF_<grp>, 1 - AF_<grp>); +fill-tags emits AF_<grp>, so we
#   compute MAF inline. Keep site if >=0.005 in ANY superpop.
echo "[$(date +%T)] [chr${CHR}] soft-union MAF filter + site QC"
bcftools view "$TAGGED" \
    --min-alleles 2 --max-alleles 2 --types snps \
    -e 'STRLEN(REF)!=1 || STRLEN(ALT)!=1 || REF!~"[ACGT]" || ALT!~"[ACGT]"' \
    -e '(INFO/AF_AFR<0.005 || INFO/AF_AFR>0.995) && (INFO/AF_AMR<0.005 || INFO/AF_AMR>0.995) && (INFO/AF_EUR<0.005 || INFO/AF_EUR>0.995) && (INFO/AF_EAS<0.005 || INFO/AF_EAS>0.995) && (INFO/AF_SAS<0.005 || INFO/AF_SAS>0.995) && (INFO/AF_CSA<0.005 || INFO/AF_CSA>0.995) && (INFO/AF_OCE<0.005 || INFO/AF_OCE>0.995) && (INFO/AF_MEN<0.005 || INFO/AF_MEN>0.995)' \
    --threads "$THREADS" -Ob -o "${TMPDIR}/merged_chr${CHR}.softunion.bcf"
bcftools index --threads "$THREADS" "${TMPDIR}/merged_chr${CHR}.softunion.bcf"

# Convert to plink2 for missingness/dup filter (matches phasing-jointcall.sh QC chain)
plink2 --bcf "${TMPDIR}/merged_chr${CHR}.softunion.bcf" \
       --set-missing-var-ids '@:#[b38]' \
       --rm-dup exclude-all \
       --geno 0.1 \
       --max-alleles 2 --snps-only just-acgt \
       --export bcf \
       --threads "$THREADS" \
       --out "${TMPDIR}/merged_chr${CHR}.qced"
mv "${TMPDIR}/merged_chr${CHR}.qced.bcf" "$QCED"
bcftools index --threads "$THREADS" "$QCED"

# 6. SHAPEIT5 phase_common (joint re-phase: HGDP+1KG + MXB together)
#    No --reference: phasing from scratch so MXB and HGDP+1KG end up on identical
#    haplotype scaffolds and the previous "wonky" rephase-then-merge issue is avoided.
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
