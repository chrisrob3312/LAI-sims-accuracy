#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1c_phase
#SBATCH --partition=atkinson,mhgcp
#SBATCH --time=72:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1c_phase_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1c_phase_chr%a_%j.err
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
# SGDP prepped BCF from 1b_sgdp_liftover_clm.sh. OPTIONAL: if present and
# non-empty it is merged into the panel alongside HGDP+1KG and MXB. If
# absent the merge proceeds with just HGDP+1KG + MXB (backwards compatible).
SGDP_PREPPED="${SGDP_PREPPED:-${OUTDIR}/sgdp_prepped.chr${CHR}.bcf}"
# Almarri prepped BCF from 1d_almarri_prep_clm.sh. Same contract as SGDP:
# merged if present and non-empty, skipped if absent.
ALMARRI_PREPPED="${ALMARRI_PREPPED:-${OUTDIR}/almarri_prepped.chr${CHR}.bcf}"

# ----------------------------------------------------------------------------
# EXCLUDE_IDS -- samples dropped from the panel before anything is merged.
#
# This exists because the unfiltered HGDP+1KG source (HGDP1KG_PATTERN pointed
# at the u235158 archive) restores 60 samples that the default relatedness-
# filtered source had removed, and ONE of them must not come back:
#
#   CHMI_CHMI3_WGS2 -- a complete hydatidiform mole cell line. A CHM genome is
#   androgenetic and effectively haploid: every site is homozygous. It is used
#   as a reference-quality assembly control, not as a population sample. In an
#   ADMIXTURE panel it is a sample with no heterozygosity at all, and as an
#   admix-simu donor it would contribute haplotypes that do not come from a
#   population. It is excluded by name, not by a heterozygosity filter, so the
#   reason is recorded rather than inferred.
#
# The default list lives in the repo so the exclusion travels with the code
# and is visible in the diff, instead of being an env var someone forgets.
# Exclusion is applied to EVERY cohort, not just HGDP+1KG, so an ID appearing
# in two sources cannot slip in through the other one.
# ----------------------------------------------------------------------------
# REPO_DIR must come from SLURM_SUBMIT_DIR, NOT from BASH_SOURCE.
#
# sbatch copies the batch script into Slurm's spool directory and runs it from
# there, so under a job `dirname "${BASH_SOURCE[0]}"` is that spool directory
# and not the checkout. This script had no REPO_DIR line at all, so the
# fallback below resolved EXCLUDE_IDS to <spool>/reference_ids/panel_exclude_
# ids.txt, which does not exist. The exclusion block then took its "no
# exclusion list" branch and CHMI_CHMI3_WGS2 -- a hydatidiform-mole cell line
# -- was merged into the panel. The run succeeded and said so only in a line
# nobody was grepping for.
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
EXCLUDE_IDS="${EXCLUDE_IDS:-${REPO_DIR}/reference_ids/panel_exclude_ids.txt}"
# A missing exclusion list is now fatal rather than a note. Every name on that
# list is there because it must not be in the panel; silently proceeding
# without it produces exactly the panel the list exists to prevent.
[[ -s "$EXCLUDE_IDS" ]] || {
    echo "ERROR: no exclusion list at ${EXCLUDE_IDS}" >&2
    echo "       REPO_DIR=${REPO_DIR} (from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset})" >&2
    echo "       Pass EXCLUDE_IDS=<path> explicitly, or submit from the repo." >&2
    exit 1; }

mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

[[ -s "$MXB_LIFTED"     ]] || { echo "ERROR: missing $MXB_LIFTED. Run 1a_prep-mxb-liftover_clm.sh first."; exit 1; }
[[ -s "$SAMPLE_GROUPS"  ]] || { echo "ERROR: missing $SAMPLE_GROUPS. Run make_sample_groups.sh first."; exit 1; }
[[ -s "$GMAP"           ]] || { echo "ERROR: missing $GMAP."; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

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

# Fold the named exclusions in with the kinship outliers so there is ONE drop
# list and one bcftools pass. Reported separately, because "excluded by name"
# and "excluded as a kinship outlier" are different facts about the panel and
# the methods have to state them differently.
DROP_IDS="${TMPDIR}/drop_ids.chr${CHR}.txt"
cp "$OUTLIERS_IDS" "$DROP_IDS"
N_EXCL=0
if [[ -s "$EXCLUDE_IDS" ]]; then
    awk 'NF && $1 !~ /^#/ {print $1}' "$EXCLUDE_IDS" >> "$DROP_IDS"
    N_EXCL=$(awk 'NF && $1 !~ /^#/' "$EXCLUDE_IDS" | wc -l)
    echo "[$(date +%T)] [chr${CHR}] named exclusions from $(basename "$EXCLUDE_IDS"): ${N_EXCL}"
    # Say which of them were actually present, so an ID that has been renamed
    # upstream shows up as "not found" instead of passing as "excluded".
    while read -r _id; do
        if bcftools query -l "$HGDP1KG_CHR_SRC" | grep -Fxq "$_id"; then
            echo "[$(date +%T)] [chr${CHR}]   excluding ${_id} (present in source)"
        else
            echo "[$(date +%T)] [chr${CHR}]   ${_id} NOT in this source (nothing to exclude)"
        fi
    done < <(awk 'NF && $1 !~ /^#/ {print $1}' "$EXCLUDE_IDS")
else
    echo "[$(date +%T)] [chr${CHR}] no exclusion list at $EXCLUDE_IDS"
fi
sort -u -o "$DROP_IDS" "$DROP_IDS"

HGDP1KG_PREP="${TMPDIR}/hgdp1kg_chr${CHR}.prep.bcf"
bcftools view -S "^${DROP_IDS}" --force-samples \
    --threads "$THREADS" -Ob -o "$HGDP1KG_PREP" "$HGDP1KG_CHR_SRC"
bcftools index "$HGDP1KG_PREP"
NHGDP1KG=$(bcftools view "$HGDP1KG_PREP" -H | wc -l)
NSAMP_HGDP1KG=$(bcftools query -l "$HGDP1KG_PREP" | wc -l)
echo "[$(date +%T)] [chr${CHR}] HGDP+1KG: ${NSAMP_HGDP1KG} samples, ${NHGDP1KG} records"
[[ "$NHGDP1KG" -gt 0 ]] || { echo "ERROR: HGDP+1KG chr${CHR} is empty"; exit 1; }

# Detect contig naming -- the non-rechr filter1 files use chr-prefixed
# contigs, but if HGDP1KG_PATTERN is pointed at a bare-numeric source
# (e.g. *.rechr.vcf.gz) we rename to chr-prefixed so the merge with the
# chr-prefixed MXB BCF and SHAPEIT5's --region chr${CHR} both line up.
# Read the header ONCE into a variable, never `bcftools view -h | grep -q`.
# Under `set -o pipefail` that idiom inverts: grep -q exits at the first match,
# bcftools dies of SIGPIPE writing the rest of the header, pipefail returns
# 141, and the `if` is FALSE even though the contig IS present. It cost 1d a
# whole array (exit 141) before it was spotted. Here it would have silently
# renamed contigs that needed no renaming, or aborted claiming the chromosome
# was absent from a file that contains it.
_HDR=$(bcftools view -h "$HGDP1KG_PREP") || {
    echo "ERROR: cannot read header of $HGDP1KG_PREP" >&2; exit 1; }
if grep -q "^##contig=<ID=chr${CHR}[,>]" <<<"$_HDR"; then
    HGDP1KG_CHR="$HGDP1KG_PREP"
elif grep -q "^##contig=<ID=${CHR}[,>]" <<<"$_HDR"; then
    echo "[$(date +%T)] [chr${CHR}] rename HGDP+1KG contigs to chr-prefixed"
    RENAME_TXT="${TMPDIR}/rename_to_chr.txt"
    : > "$RENAME_TXT"
    for c in {1..22} X Y MT; do echo "$c chr$c" >> "$RENAME_TXT"; done
    HGDP1KG_CHR="${TMPDIR}/hgdp1kg_chr${CHR}.bcf"
    bcftools annotate --rename-chrs "$RENAME_TXT" \
        --threads "$THREADS" -Ob -o "$HGDP1KG_CHR" "$HGDP1KG_PREP"
    bcftools index "$HGDP1KG_CHR"
else
    echo "ERROR: neither 'chr${CHR}' nor '${CHR}' contig in $HGDP1KG_CHR_SRC header"
    exit 1
fi

# 2. Merge HGDP+1KG (chr$CHR) + lifted MXB (chr$CHR) [+ SGDP if present]
#    -- column-wise sample join. SGDP contigs are renamed to chr-prefixed
#    if the prepped BCF came from a bare-numeric source, matching the
#    HGDP+1KG/MXB convention so the merge and SHAPEIT5 --region line up.
MERGE_INPUTS=("$HGDP1KG_CHR" "$MXB_LIFTED")

# Running list of sample IDs already in the merge. Each additional cohort is
# deduplicated against ALL of it, not just against HGDP+1KG: once more than one
# optional cohort can join, "already present" has to mean present in any of
# them, or Almarri could re-introduce an ID that SGDP had just contributed.
SEEN_IDS="${TMPDIR}/seen_ids.chr${CHR}.txt"
{ bcftools query -l "$HGDP1KG_CHR"; bcftools query -l "$MXB_LIFTED"; } \
    | sort -u > "$SEEN_IDS"
echo "[$(date +%T)] [chr${CHR}] base panel: $(wc -l < "$SEEN_IDS") samples (HGDP+1KG + MXB)"

# add_cohort <label> <prepped bcf>
#   Normalises contig naming, drops named exclusions and any ID already in the
#   merge, appends to MERGE_INPUTS and extends SEEN_IDS. A cohort whose file is
#   absent is skipped with a note; one whose file is present but unusable is a
#   hard error, because silently phasing a panel that is missing a whole cohort
#   is the failure that is expensive to discover later.
add_cohort () {
    local label="$1" src="$2" lc
    lc=$(echo "$label" | tr 'A-Z' 'a-z')
    if [[ ! -s "$src" ]]; then
        echo "[$(date +%T)] [chr${CHR}] no ${label} prepped BCF at ${src} -- skipping"
        return 0
    fi
    local n; n=$(bcftools query -l "$src" | wc -l)
    echo "[$(date +%T)] [chr${CHR}] ${label} present: ${n} samples -> merge"

    local cur="$src" hdr
    # Header read once; see the SIGPIPE note on the HGDP+1KG block above.
    hdr=$(bcftools view -h "$cur") || {
        echo "ERROR: cannot read header of $src" >&2; return 1; }
    if grep -q "^##contig=<ID=chr${CHR}[,>]" <<<"$hdr"; then
        :
    elif grep -q "^##contig=<ID=${CHR}[,>]" <<<"$hdr"; then
        echo "[$(date +%T)] [chr${CHR}] rename ${label} contigs to chr-prefixed"
        local rn="${TMPDIR}/${lc}_rename_to_chr.txt"
        : > "$rn"
        for c in {1..22} X Y MT; do echo "$c chr$c" >> "$rn"; done
        local out="${TMPDIR}/${lc}_chr${CHR}.rechr.bcf"
        bcftools annotate --rename-chrs "$rn" --threads "$THREADS" -Ob -o "$out" "$cur"
        bcftools index "$out"
        cur="$out"
    else
        echo "ERROR: neither 'chr${CHR}' nor '${CHR}' contig in ${src} header" >&2
        return 1
    fi

    # bcftools merge aborts on duplicate sample names, so collisions have to go
    # before the merge, not after. Named exclusions are applied here too: an ID
    # on the exclusion list must not enter through a cohort other than the one
    # it was noticed in.
    local drop="${TMPDIR}/${lc}_drop_ids.chr${CHR}.txt"
    comm -12 <(bcftools query -l "$cur" | sort -u) "$SEEN_IDS" > "$drop" || true
    if [[ -s "$EXCLUDE_IDS" ]]; then
        awk 'NF && $1 !~ /^#/ {print $1}' "$EXCLUDE_IDS" >> "$drop"
    fi
    sort -u -o "$drop" "$drop"
    # Only count IDs this cohort actually carries, so the log does not report
    # dropping samples that were never here.
    local ndrop
    ndrop=$(comm -12 <(bcftools query -l "$cur" | sort -u) "$drop" | wc -l)
    if [[ "$ndrop" -gt 0 ]]; then
        echo "[$(date +%T)] [chr${CHR}] dropping ${ndrop} ${label} samples (duplicate or excluded)"
        local ded="${TMPDIR}/${lc}_chr${CHR}.dedup.bcf"
        bcftools view -S "^${drop}" --force-samples --threads "$THREADS" \
            -Ob -o "$ded" "$cur"
        bcftools index "$ded"
        cur="$ded"
    fi

    local kept; kept=$(bcftools query -l "$cur" | wc -l)
    if [[ "$kept" -eq 0 ]]; then
        echo "ERROR: ${label} contributed 0 samples to chr${CHR} after dedup." >&2
        echo "       Its file exists, so this is a real problem (ID scheme" >&2
        echo "       mismatch, or every sample already in the panel), not an" >&2
        echo "       absent cohort. Refusing to phase a panel silently missing it." >&2
        return 1
    fi
    echo "[$(date +%T)] [chr${CHR}] ${label} contributes ${kept} samples"
    # Not `sort -m -o "$SEEN_IDS" "$SEEN_IDS" ...`: a merge sort streams its
    # inputs, so writing back over one of them is not safe the way a plain
    # `sort -o file file` is. Build beside it, then move.
    bcftools query -l "$cur" | sort -u > "${TMPDIR}/${lc}_kept.txt"
    sort -u "$SEEN_IDS" "${TMPDIR}/${lc}_kept.txt" > "${SEEN_IDS}.new"
    mv "${SEEN_IDS}.new" "$SEEN_IDS"
    MERGE_INPUTS+=("$cur")
    return 0
}

# REQUIRE_COHORTS -- fail if a cohort that should be in this run is absent.
# Checked BEFORE any cohort is added, so the run dies in seconds rather than
# after the merge. The default is permissive for backwards compatibility, but
# a publication run sets it: the expensive mistake is not a crash, it is
# joint-phasing a panel quietly missing SGDP or Almarri and only noticing at
# ADMIXTURE, by which point everything has to be re-phased.
REQUIRE_COHORTS="${REQUIRE_COHORTS:-}"
for _req in ${REQUIRE_COHORTS//,/ }; do
    case "$_req" in
        SGDP)    [[ -s "$SGDP_PREPPED"    ]] || { echo "ERROR: REQUIRE_COHORTS lists SGDP but $SGDP_PREPPED is missing" >&2; exit 1; } ;;
        ALMARRI) [[ -s "$ALMARRI_PREPPED" ]] || { echo "ERROR: REQUIRE_COHORTS lists ALMARRI but $ALMARRI_PREPPED is missing" >&2; exit 1; } ;;
        *)       echo "ERROR: REQUIRE_COHORTS has unknown cohort '$_req' (valid: SGDP ALMARRI)" >&2; exit 1 ;;
    esac
done
[[ -n "$REQUIRE_COHORTS" ]] && echo "[$(date +%T)] [chr${CHR}] required cohorts present: ${REQUIRE_COHORTS}"

add_cohort SGDP    "$SGDP_PREPPED"    || exit 1
add_cohort ALMARRI "$ALMARRI_PREPPED" || exit 1

echo "[$(date +%T)] [chr${CHR}] panel after all cohorts: $(wc -l < "$SEEN_IDS") samples"

echo "[$(date +%T)] [chr${CHR}] bcftools merge (${#MERGE_INPUTS[@]} inputs)"
bcftools merge --threads "$THREADS" \
    -Ob -o "$MERGED" \
    "${MERGE_INPUTS[@]}"
bcftools index "$MERGED"

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
       --output-chr chrM \
       --export bcf \
       --threads "$THREADS" \
       --out "$QCED_PREFIX"
bcftools index "$QCED"

# Sanity-check: plink2's default for VCF/BCF export is to strip the "chr"
# prefix from contig names, which then makes the soft-union BCF unreadable
# by SHAPEIT5 when invoked with `--region chr${CHR}` (#sites=0, no error).
# `--output-chr chrM` above keeps them chr-prefixed; assert it stuck.
# Header read once. As a PIPELINE this assertion inverted under pipefail: when
# the prefix WAS present, grep -q matched and exited, bcftools died of SIGPIPE
# finishing the header, the pipeline returned 141, `!` made that true, and the
# task aborted reporting the prefix lost. A check that fires on success is
# worse than no check.
_QHDR=$(bcftools view -h "$QCED") || {
    echo "ERROR: cannot read header of $QCED" >&2; exit 1; }
if ! grep -q "^##contig=<ID=chr${CHR}[,>]" <<<"$_QHDR"; then
    echo "ERROR: QCED BCF lost the chr-prefix on contigs (plink2 --output-chr failed?)"
    grep "^##contig" <<<"$_QHDR" | awk 'NR<=3'
    exit 1
fi
unset _QHDR

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
bcftools index "$SOFTUNION"
NSOFT=$(bcftools view "$SOFTUNION" -H | wc -l)
echo "[$(date +%T)] [chr${CHR}] soft-union BCF: ${NSOFT} records"
[[ "$NSOFT" -gt 0 ]] || { echo "ERROR: subset produced empty BCF"; exit 1; }

# 5b. ADMIXTURE site list -- per cohort x super-pop missingness, then overall.
#     Measured HERE, on the pre-phasing genotypes, because SHAPEIT5 fills every
#     missing genotype and the phased panel has none left to measure. In the
#     November panel every SGDP sample had no call at 39.6% of chr22's
#     soft-union sites and every MXB sample at 61.9%; at the SGDP-missing sites
#     the 1KG/HGDP copy of the SAME person was non-ref 27.8% of the time, so
#     those genotypes were neither called nor safely hom-ref. ADMIXTURE then
#     treated the filled-in values as observed, and three of the eight K=8
#     clusters tracked cohorts instead of ancestry.
#
#     Cohort comes from which FILE a sample was merged from (not its ID);
#     super-pop from sample_groups.tsv / the SGDP manifest (Almarri = MEN,
#     MXB = AMR). Only the ADMIXTURE marker set uses this list -- the soft-union
#     itself, and so the phased LAI panel, is unchanged.
ADMIX_GROUP_MAXMISS="${ADMIX_GROUP_MAXMISS:-0.10}"
ADMIX_TOTAL_MAXMISS="${ADMIX_TOTAL_MAXMISS:-0.07}"
ADMIX_MISS_MIN_N="${ADMIX_MISS_MIN_N:-5}"
SGDP_MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
ADMIX_SITES="${OUTDIR}/merged_chr${CHR}.admix_sites.bed"
ADMIX_SUMMARY="${OUTDIR}/merged_chr${CHR}.admix_sites.summary.tsv"
MISS_GROUPS="${TMPDIR}/miss_groups.chr${CHR}.tsv"
bcftools query -l "$SOFTUNION" > "${TMPDIR}/softunion_samples.chr${CHR}.txt"
bcftools query -l "$MXB_LIFTED" > "${TMPDIR}/cohort_mxb.chr${CHR}.txt"
COHORT_ARGS=(--cohort "MXB=${TMPDIR}/cohort_mxb.chr${CHR}.txt")
[[ -s "${TMPDIR}/sgdp_kept.txt"    ]] && COHORT_ARGS+=(--cohort "SGDP=${TMPDIR}/sgdp_kept.txt")
[[ -s "${TMPDIR}/almarri_kept.txt" ]] && COHORT_ARGS+=(--cohort "ALMARRI=${TMPDIR}/almarri_kept.txt")
echo "[$(date +%T)] [chr${CHR}] ADMIXTURE site list: groups <= ${ADMIX_GROUP_MAXMISS} missing" \
     "(n >= ${ADMIX_MISS_MIN_N}), overall <= ${ADMIX_TOTAL_MAXMISS}"
python3 -I "${REPO_DIR}/scripts/build-missingness-groups.py" \
    --samples "${TMPDIR}/softunion_samples.chr${CHR}.txt" "${COHORT_ARGS[@]}" \
    --default-cohort HGDP1KG --groups "$SAMPLE_GROUPS" --sgdp-manifest "$SGDP_MANIFEST" \
    > "$MISS_GROUPS"
ADMIX_FMT=$(python3 -I "${REPO_DIR}/scripts/admix-site-filter.py" --groups "$MISS_GROUPS" --print-format)
bcftools +fill-tags "$SOFTUNION" -Ou -- -S "$MISS_GROUPS" -t AN \
  | bcftools query -f "$ADMIX_FMT" \
  | python3 -I "${REPO_DIR}/scripts/admix-site-filter.py" --groups "$MISS_GROUPS" \
        --max-group-miss "$ADMIX_GROUP_MAXMISS" --max-total-miss "$ADMIX_TOTAL_MAXMISS" \
        --min-group-n "$ADMIX_MISS_MIN_N" \
        --out "${ADMIX_SITES}.partial" --summary "${ADMIX_SUMMARY}.partial"
[[ -s "${ADMIX_SITES}.partial" ]] || { echo "ERROR: no site passed the ADMIXTURE missingness rule on chr${CHR}"; exit 1; }
mv "${ADMIX_SUMMARY}.partial" "$ADMIX_SUMMARY"
mv "${ADMIX_SITES}.partial" "$ADMIX_SITES"
echo "[$(date +%T)] [chr${CHR}] ADMIXTURE site list: $(wc -l < "$ADMIX_SITES") positions -> ${ADMIX_SITES}"
cat "$ADMIX_SUMMARY"

# ---------------------------------------------------------------------------
# PHASE_MODE -- whole-chromosome (default) or hand off to the chunked path.
#
# phase_common on a whole chromosome of this panel takes about 30 hours, so 22
# chromosomes at a throttle of 15 is roughly two and a half days, and stage 2
# gates every later stage.
#
# The repo already has a tested chunked path (1b_phase_chunks + 1b_phase_ligate,
# 8 cM windows with 2 cM overlap) which splits the same work into a few hundred
# independent tasks. Its ligate step writes EXACTLY the path this script
# writes, so it is a drop-in for step 6 and nothing downstream can tell the
# difference.
#
# The one thing missing was the handover: the chunk script reads
# ${OUTDIR}/merged_chr${CHR}.softunion.bcf, and this script had the identical
# file one directory down in its own tmp. So PHASE_MODE=chunked publishes the
# soft-union where the chunk array looks for it and stops.
#
# Everything ABOVE this point -- the four-cohort merge, the exclusions, the
# cohort checks, the per-super-pop MAF soft union -- runs identically either
# way. Only who does the phasing changes.
PHASE_MODE="${PHASE_MODE:-whole}"
if [[ "$PHASE_MODE" == "chunked" ]]; then
    PUB="${OUTDIR}/merged_chr${CHR}.softunion.bcf"
    echo "[$(date +%T)] [chr${CHR}] PHASE_MODE=chunked -- publishing the soft-union"
    cp "$SOFTUNION" "${PUB}.partial"
    bcftools index --threads "$THREADS" "${PUB}.partial"
    # Index first, then the BCF: the chunk array checks for the BCF, so it must
    # not appear before its index exists.
    mv "${PUB}.partial.csi" "${PUB}.csi"
    mv "${PUB}.partial" "$PUB"
    echo "[$(date +%T)] [chr${CHR}] -> ${PUB} ($(bcftools index -n "$PUB") sites, "\
         "$(bcftools query -l "$PUB" | wc -l) samples)"
    echo "[$(date +%T)] [chr${CHR}] phasing is now 1b_phase_chunks + 1b_phase_ligate."
    echo "[$(date +%T)] [chr${CHR}] ligate writes ${PHASED}, which is what this"
    echo "[$(date +%T)] [chr${CHR}] script would have written, so stages 3+ are unchanged."
    exit 0
elif [[ "$PHASE_MODE" != "whole" ]]; then
    echo "ERROR: PHASE_MODE must be 'whole' or 'chunked', got '${PHASE_MODE}'" >&2
    exit 1
fi

# 6. SHAPEIT5 phase_common (joint re-phase: HGDP+1KG + MXB together)
#    SHAPEIT5 writes its own .csi index after phasing, so we do NOT call
#    `bcftools index` on $PHASED -- a second index would error
#    `[E::main_vcfindex] the index file exists` and the script would die
#    after the ~30 hr phase finished. Trust SHAPEIT5's index.
echo "[$(date +%T)] [chr${CHR}] SHAPEIT5_phase_common"
SHAPEIT5_phase_common \
    --input "$SOFTUNION" \
    --map "$GMAP" \
    --region "chr${CHR}" \
    --output "$PHASED" \
    --thread "$THREADS" \
    --filter-maf 0.001

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
    # SHAPEIT5 writes its own index here too -- skip a redundant bcftools index.
fi

# 7. Strip "chr" prefix on the LAI panel (RFMix v1 expects bare numeric contigs)
echo "[$(date +%T)] [chr${CHR}] rename chr${CHR} -> ${CHR}"
echo "chr${CHR} ${CHR}" > "${TMPDIR}/rename_chr${CHR}.txt"
bcftools annotate --rename-chrs "${TMPDIR}/rename_chr${CHR}.txt" \
    --threads "$THREADS" -Ob -o "$RECHR" "$PHASED"
bcftools index -f "$RECHR"
rm "${TMPDIR}/rename_chr${CHR}.txt"

echo "[$(date +%T)] [chr${CHR}] # Complete."
echo "[$(date +%T)] [chr${CHR}] # LAI panel (chr*)   : $PHASED"
echo "[$(date +%T)] [chr${CHR}] # LAI panel (rechr)  : $RECHR"
[[ "$RUN_PHASE_RARE" == "1" ]] && \
    echo "[$(date +%T)] [chr${CHR}] # Full phased panel  : $FULL_PHASED"
