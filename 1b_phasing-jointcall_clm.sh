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
#   1. Take per-chrom HGDP+1KG filter1 file (HGDP1KG_PATTERN, default points
#      at the colleague's `phased_haplotypes_v2_filter1/` per-chrom outputs),
#      drop kinship outliers from related_outliers.txt, and normalize contig
#      naming to chr-prefixed so it lines up with the lifted MXB BCF and with
#      SHAPEIT5's `--region chr${CHR}`. Prior phase in the filter1 input is
#      irrelevant -- SHAPEIT5 re-phases jointly with MXB in step 6.
#   2. bcftools merge with the lifted MXB chunk (output of 1a_prep-mxb-liftover_clm.sh).
#      NOTE: merge (column-wise sample join), NOT concat. Different samples,
#      same/overlapping sites -> 4147 sample columns per site. Missing in one
#      panel becomes ./. and is imputed by SHAPEIT5 during phasing.
#   3. plink2 site QC: biallelic SNPs, ACGT only, dedup IDs (no --geno here;
#      missingness is enforced per-superpop in step 4).
#   4. Per-superpop pre-filter -- subset the QC'd BCF to each superpop in
#      sample_groups.tsv (>= MIN_SUBPOP_N samples), recompute AC,AN on the
#      subset, keep sites where AN>0 AND MAF >= LAI_MAF; emit each pop's
#      passing sites to a TSV. Soft-union = sort -u of all per-pop TSVs.
#      F_MISSING is NOT filtered: the HGDP+1KG postoutlier joint call is
#      sparse-encoded (homref dropped at MXB-private sites), which makes
#      per-pop missingness uninformative. AN>0 already skips sites a pop
#      has no information on, and SHAPEIT5 imputes residual missing GTs.
#   5. Subset the full-sample QC'd BCF to the soft-union site list -- this is
#      the input to SHAPEIT5 phase_common.
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
# Per-chrom filter1 file pattern -- use the literal token __CHR__ as the
# placeholder. Default points at the colleague's per-chrom phased +
# MAF005-filtered + indexed outputs (chr-prefixed, outliers still present).
# SHAPEIT5 will re-phase them jointly with MXB anyway, so the prior phase
# information in these files isn't relied upon.
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/hgdp1kgp_chr__CHR__.shapeit5_phased.filter1_SNP_maf005.vcf.gz}"
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
MIN_SUBPOP_N="${MIN_SUBPOP_N:-10}"  # skip superpops smaller than this
                              # (MAF estimates too noisy)

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

MERGED="${TMPDIR}/merged_chr${CHR}.bcf"
SOFTUNION="${TMPDIR}/merged_chr${CHR}.softunion.bcf"
QCED_PREFIX="${TMPDIR}/merged_chr${CHR}.qced"
QCED="${QCED_PREFIX}.bcf"
PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
RECHR="${OUTDIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.rechr.bcf"
FULL_PHASED="${OUTDIR}/merged_chr${CHR}.shapeit5_full_phased.bcf"

# 1. Resolve per-chrom HGDP+1KG file, drop kinship outliers, normalize contig
#    naming to chr-prefixed (to match the lifted MXB BCF from step 1a).
#
#    HGDP1KG_PATTERN expands the literal token __CHR__ to the chromosome
#    number to form the per-chrom path. Default points at the colleague's
#    filter1 phased BCFs (already MAF005-filtered, indexed, chr-prefixed
#    contigs). SHAPEIT5 re-phases the merge from genotypes in step 6, so
#    any prior phasing in the input is discarded -- this gives us the
#    joint phasing of HGDP+1KG + MXB that we need.
#
#    These filter1 files are pre-outlier (n=4091); we still drop the kinship
#    outliers from related_outliers.txt. The file is 2 cols (super-pop,
#    sample_id) -- bcftools -S expects one ID per line, so extract col 2 to
#    a temp file first.
HGDP1KG_CHR_SRC="${HGDP1KG_PATTERN//__CHR__/$CHR}"
[[ -s "$HGDP1KG_CHR_SRC" ]] || {
    echo "ERROR: missing per-chrom HGDP+1KG file $HGDP1KG_CHR_SRC"
    echo "  HGDP1KG_PATTERN=$HGDP1KG_PATTERN"
    exit 1
}

echo "[$(date +%T)] [chr${CHR}] drop kinship outliers from $HGDP1KG_CHR_SRC"
OUTLIERS_IDS="${TMPDIR}/related_outlier_ids.txt"
awk '{print $2}' "$OUTLIERS" > "$OUTLIERS_IDS"
HGDP1KG_PREP="${TMPDIR}/hgdp1kg_chr${CHR}.prep.bcf"
bcftools view -S "^${OUTLIERS_IDS}" --force-samples \
    --threads "$THREADS" -Ob -o "$HGDP1KG_PREP" "$HGDP1KG_CHR_SRC"
bcftools index --threads "$THREADS" "$HGDP1KG_PREP"
NHGDP1KG=$(bcftools view "$HGDP1KG_PREP" -H | wc -l)
NSAMP_HGDP1KG=$(bcftools query -l "$HGDP1KG_PREP" | wc -l)
echo "[$(date +%T)] [chr${CHR}] HGDP+1KG: ${NSAMP_HGDP1KG} samples, ${NHGDP1KG} records"
[[ "$NHGDP1KG" -gt 0 ]] || { echo "ERROR: HGDP+1KG chr${CHR} is empty"; exit 1; }

# Detect contig naming -- the non-rechr filter1 files use chr-prefixed
# contigs, but if HGDP1KG_PATTERN is pointed at a bare-numeric source
# (e.g. *.rechr.vcf.gz) we rename to chr-prefixed so the merge with the
# chr-prefixed MXB BCF and SHAPEIT5's --region chr${CHR} both line up.
if bcftools view -h "$HGDP1KG_PREP" | grep -q "^##contig=<ID=chr${CHR}[,>]"; then
    HGDP1KG_CHR="$HGDP1KG_PREP"
elif bcftools view -h "$HGDP1KG_PREP" | grep -q "^##contig=<ID=${CHR}[,>]"; then
    echo "[$(date +%T)] [chr${CHR}] rename HGDP+1KG contigs to chr-prefixed"
    RENAME_TXT="${TMPDIR}/rename_to_chr.txt"
    : > "$RENAME_TXT"
    for c in {1..22} X Y MT; do echo "$c chr$c" >> "$RENAME_TXT"; done
    HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
    bcftools annotate --rename-chrs "$RENAME_TXT" \
        --threads "$THREADS" -Ob -o "$HGDP1KG_CHR" "$HGDP1KG_PREP"
    bcftools index --threads "$THREADS" "$HGDP1KG_CHR"
else
    echo "ERROR: neither 'chr${CHR}' nor '${CHR}' contig in $HGDP1KG_CHR_SRC header"
    exit 1
fi

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

# 4. Per-superpop pre-filter, then soft-union of passing sites.
#    For each superpop in $SAMPLE_GROUPS with N >= MIN_SUBPOP_N:
#      a. subset the QC'd BCF to its samples
#      b. recompute AC, AN on the subset (single-pop tags, no -S quirks)
#      c. keep sites where AN>0 AND MAF >= LAI_MAF, expressed without
#         division as `AC>=MAF*AN AND (AN-AC)>=MAF*AN`
#      d. emit chr/pos/ref/alt of passing sites
#    Soft-union = sort -u across per-pop TSVs.
#
#    NOTE: we deliberately do NOT filter on F_MISSING here. The HGDP+1KG
#    postoutlier joint call is sparse-encoded -- samples that are homref
#    at MXB-private sites are emitted as `.|.`, not `0|0`. So F_MISSING in
#    a non-MXB superpop at a MXB-only site is ~100% even though the WGS
#    coverage is fine. The AN>0 check in the MAF expression already drops
#    sites a pop has no information on; SHAPEIT5_phase_common imputes the
#    residual missing genotypes during phasing.
#
#    Tiny superpops (N < MIN_SUBPOP_N, e.g. OCE/MEN in HGDP+1KG) are
#    skipped -- their MAF estimates are too noisy to constrain the panel.
echo "[$(date +%T)] [chr${CHR}] per-superpop pre-filter (MAF>=${LAI_MAF}, soft-union)"
VCF_SAMPLES="${TMPDIR}/vcf_samples.txt"
bcftools query -l "$QCED" > "$VCF_SAMPLES"

SITELIST="${TMPDIR}/softunion_sites.tsv"
: > "$SITELIST"

# MAF >= m, expressed as: AC>=m*AN AND (AN-AC)>=m*AN, with AN>0 sentinel.
MAF_EXPR="INFO/AN>0 && INFO/AC>=${LAI_MAF}*INFO/AN && (INFO/AN-INFO/AC)>=${LAI_MAF}*INFO/AN"

for GROUP in $(awk '{print $2}' "$SAMPLE_GROUPS" | sort -u); do
    GROUP_KEEP="${TMPDIR}/${GROUP}_keep.txt"
    awk -v g="$GROUP" '$2==g {print $1}' "$SAMPLE_GROUPS" \
        | grep -xFf "$VCF_SAMPLES" > "$GROUP_KEEP" || true
    NSAMP=$(wc -l < "$GROUP_KEEP")
    if (( NSAMP < MIN_SUBPOP_N )); then
        echo "[$(date +%T)] [chr${CHR}]   skip ${GROUP}: ${NSAMP} samples < MIN_SUBPOP_N=${MIN_SUBPOP_N}"
        continue
    fi
    GROUP_SITES="${TMPDIR}/${GROUP}_sites.tsv"
    bcftools view "$QCED" -S "$GROUP_KEEP" --force-samples --threads "$THREADS" -Ou \
      | bcftools +fill-tags --threads "$THREADS" -Ou -- -t 'AC,AN' \
      | bcftools view -i "$MAF_EXPR" --threads "$THREADS" -Ou \
      | bcftools query -f '%CHROM\t%POS\n' > "$GROUP_SITES"
    NSITES=$(wc -l < "$GROUP_SITES")
    echo "[$(date +%T)] [chr${CHR}]   ${GROUP}: ${NSAMP} samples, ${NSITES} passing sites"
    cat "$GROUP_SITES" >> "$SITELIST"
done

SITELIST_SORTED="${TMPDIR}/softunion_sites.sorted.tsv"
sort -k1,1 -k2,2n -u "$SITELIST" > "$SITELIST_SORTED"
NUNION=$(wc -l < "$SITELIST_SORTED")
echo "[$(date +%T)] [chr${CHR}] soft-union sites: ${NUNION}"
[[ "$NUNION" -gt 0 ]] || { echo "ERROR: empty soft-union site list"; exit 1; }

# bgzip + tabix the chr/pos targets file so bcftools view -T can stream it.
# 2-col (chr, pos) format avoids the chr/pos/ref/alt → chr/from/to mis-parse
# that bcftools -T does when columns 3/4 are non-numeric (REF/ALT letters).
# Plink2 --rm-dup exclude-all already made positions unique, so chr/pos is
# sufficient to identify each site.
bgzip -f "$SITELIST_SORTED"
tabix -s1 -b2 -e2 -f "${SITELIST_SORTED}.gz"

# 5. Subset the full-sample QC'd BCF to the soft-union sites for SHAPEIT5.
echo "[$(date +%T)] [chr${CHR}] subset full BCF to soft-union sites"
bcftools view "$QCED" -T "${SITELIST_SORTED}.gz" \
    --threads "$THREADS" -Ob -o "$SOFTUNION"
bcftools index --threads "$THREADS" "$SOFTUNION"
NSOFT=$(bcftools view "$SOFTUNION" -H | wc -l)
echo "[$(date +%T)] [chr${CHR}] soft-union BCF: ${NSOFT} records"
[[ "$NSOFT" -gt 0 ]] || { echo "ERROR: subset produced empty BCF"; exit 1; }

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
