#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1c_phase
#SBATCH --partition=atkinson,mhgcp
# Same exclusion as every other stage (reason in 1b_sgdp_liftover_clm.sh): the
# mhgcp-t* CPUs lack instructions this conda env's bcftools is built with, and
# a task that lands there dies at once with SIGILL (exit 132, empty log).
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
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
#   1. Take the per-chrom HGDP+1KG file (HGDP1KG_PATTERN; default = Koenig et
#      al.'s UNFILTERED phased_haplotypes_v2 release: 4,091 QC'd genomes, SNVs
#      and indels, no MAF cut -- so the MAF floor is applied per super-pop in
#      step 4 instead of inherited from filter1's global -q 0.005), drop the
#      named exclusions (and, only if OUTLIERS names a list, those IDs --
#      by default nobody is dropped for relatedness; see OUTLIERS below),
#      chr-prefix contigs, split every multi-allelic record and check REF
#      against GRCh38 (norm -m -any -c w, +fixref). Prior phase is irrelevant:
#      the merged panel is re-phased jointly.
#   2. Add MXB (1a), SGDP and Almarri (dedup by ID, exclusions), split any
#      multi-allelic record still in them, log site provenance, fill records
#      absent from a complete callset as 0/0 (FILL_REF_COHORTS), then
#      bcftools merge -m none (column-wise sample join, NOT concat).
#   3. plink2 site QC (biallelic ACGT SNPs, dedup); drop positions that still
#      hold more than one record, and positions where a 0/0 fill is
#      contradicted by another cohort of the same super-pop
#      (scripts/fill-contradictions.py).
#   4. Per-superpop MAF soft-union, one pass: keep a site if, in ANY
#      super-pop group with >= MIN_SUBPOP_N samples, AN>0 and MAF >= LAI_MAF.
#      SOFTUNION_GROUPS=superpop (default) pools every cohort's samples by
#      super-pop; =sample_groups is the July/November rule (sample_groups.tsv
#      only, i.e. HGDP+1KG and MXB -- SGDP and Almarri counted nowhere).
#   5. Subset the full-sample QC'd BCF to the soft-union site list -- this is
#      the input to SHAPEIT5 phase_common. 5b writes the ADMIXTURE site list.
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
# Per-chrom HGDP+1KG file pattern -- the literal token __CHR__ is the
# placeholder. Default: the UNFILTERED Koenig et al. 2024 release on disk,
# original_data_downloads/phased_haplotypes_v2 (chr-prefixed BCFs, 4,091
# QC'd genomes, SNVs + indels, site QC = HWE / F_MISSING <= 0.1 / ExcHet, no
# MAF cut). filter1 (processed_data/phased_haplotypes_v2_filter1) is the same
# file cut to biallelic SNPs at a GLOBAL MAF >= 0.005, which on chr21 keeps
# 201,120 of 923,540 SNPs and removes before the soft-union can see them the
# variants that are common in one super-pop and rare overall. The prior phase
# is not used: the merged panel is re-phased jointly.
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chr__CHR__.filtered.SNV_INDEL.phased.shapeit5.bcf}"
# RELATEDNESS. Relatives are KEPT -- the project's policy, stated in
# 1e_king_dedup_clm.sh: only duplicates / MZ twins are removed (1e, KING >
# 0.354); the GWAS models kinship with a GRM (GENESIS); every donor /
# reference / test split keeps relatives at KING >= 0.177 on ONE side (2b, 4a,
# split-homog-donors-reference.py). Pruning relatives would preferentially
# remove less-studied populations, where background relatedness is higher.
#   OUTLIERS=none    (default) drop nobody for relatedness.
#   OUTLIERS=<file>  drop column 2 of a 2-column (super-pop, sample) list --
#                    e.g. the lab's related_outliers.txt (693 IDs, relatives
#                    and outliers), which the July and November panels dropped:
#   /storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt
OUTLIERS="${OUTLIERS:-none}"
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
if [[ "$OUTLIERS" != none && ! -s "$OUTLIERS" ]]; then
    echo "ERROR: OUTLIERS=${OUTLIERS} does not exist or is empty." >&2
    echo "       OUTLIERS=none keeps every sample (relatives included)." >&2
    exit 1
fi

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

# ----------------------------------------------------------------------------
# IDEMPOTENT RE-RUN (PHASE_MODE=chunked). The published soft-union is reused
# when (a) it is newer than every input file and (b) its stamp records the same
# settings AND the same code (this script plus the helpers it runs, hashed).
# Re-submitting the rebuild after a success -- or after a LATER stage failed --
# then costs seconds here instead of a full merge, which would also leave the
# soft-union newer than every phased chunk and so re-phase all of them.
# FORCE_1C=1 rebuilds regardless. The stamp is written AFTER the soft-union is
# in place and deleted when a rebuild starts, so a stamp never vouches for a
# file it did not see finished.
# ----------------------------------------------------------------------------
PUB="${OUTDIR}/merged_chr${CHR}.softunion.bcf"
_code_md5=$(cat "${BASH_SOURCE[0]}" "${REPO_DIR}/scripts/fill-contradictions.py" \
               "${REPO_DIR}/scripts/softunion-sites.py" "${REPO_DIR}/scripts/build-missingness-groups.py" \
               "${REPO_DIR}/scripts/admix-site-filter.py" "${REPO_DIR}/scripts/site-provenance.py" \
            | md5sum | cut -c1-12)
_excl_md5=$(md5sum < "$EXCLUDE_IDS" | cut -c1-12)
# The relatedness drop is part of the stamp: a soft-union built without
# relatives must never be reused as one built with them, or the reverse.
if [[ "$OUTLIERS" == none ]]; then _outl=none; _outl_in=""
else _outl=$(md5sum < "$OUTLIERS" | cut -c1-12); _outl_in="$OUTLIERS"; fi
PUB_RECIPE="1c code=${_code_md5} excl=${_excl_md5} outliers=${_outl} src=${HGDP1KG_PATTERN}"
for _v in FILL_REF_COHORTS FILL_GUARD FILL_GUARD_MIN_N FILL_GUARD_MIN_AC FILL_GUARD_MAX_P \
          SOFTUNION_GROUPS FIXREF_MODE HGDP1KG_NORMALIZE ADMIX_GROUP_MAXMISS ADMIX_TOTAL_MAXMISS \
          ADMIX_MISS_MIN_N REQUIRE_COHORTS; do
    PUB_RECIPE+=" ${_v}=${!_v-<default>}"
done
PUB_RECIPE+=" LAI_MAF=${LAI_MAF} MIN_SUBPOP_N=${MIN_SUBPOP_N}"
if [[ "${PHASE_MODE:-whole}" == "chunked" && "${FORCE_1C:-0}" != "1" \
      && -s "$PUB" && -s "${PUB}.recipe" && -s "${OUTDIR}/merged_chr${CHR}.admix_sites.bed" ]] \
   && [[ "$(cat "${PUB}.recipe")" == "$PUB_RECIPE" ]] && bcftools view -h "$PUB" > /dev/null; then
    _stale=""
    for _in in "${HGDP1KG_PATTERN//__CHR__/$CHR}" "$MXB_LIFTED" "$SGDP_PREPPED" "$ALMARRI_PREPPED" \
               "$EXCLUDE_IDS" "$SAMPLE_GROUPS" "$_outl_in" "${REPO_DIR}/reference_ids/sgdp_manifest.tsv"; do
        if [[ -e "$_in" && "$_in" -nt "$PUB" ]]; then _stale+=" $(basename "$_in")"; fi
    done
    if [[ -z "$_stale" ]]; then
        echo "[$(date +%T)] [chr${CHR}] ${PUB} is newer than every input and was made by this code" \
             "with these settings -- nothing to do (FORCE_1C=1 to rebuild)"
        exit 0
    fi
    echo "[$(date +%T)] [chr${CHR}] rebuilding: newer input(s):${_stale}"
elif [[ -s "$PUB" ]]; then
    echo "[$(date +%T)] [chr${CHR}] rebuilding ${PUB}: settings/code differ from its stamp, or no stamp" \
         "(stamp: '$(cat "${PUB}.recipe" 2>/dev/null || echo none)')"
fi
rm -f "${PUB}.recipe"

# 1. Resolve per-chrom HGDP+1KG file, drop the named exclusions (and the
#    OUTLIERS list, if one is given), normalize contig
#    naming to chr-prefixed (to match the lifted MXB BCF from step 1a).
#
#    HGDP1KG_PATTERN expands the literal token __CHR__ to the chromosome
#    number to form the per-chrom path. Default points at the unfiltered
#    Koenig v2 BCFs (indexed, chr-prefixed contigs). The merged panel is
#    re-phased from genotypes, so any prior phasing in the input is
#    discarded -- this gives the joint phasing of all cohorts that we need.
#
#    v2 (like filter1, which is cut from it) holds 4,091 samples including
#    relatives, and keeps them (OUTLIERS=none, the default). The July and
#    November panels dropped related_outliers.txt here (693 IDs), which took
#    11 of the 31 amr_rfmix samples -- 4 of 6 Karitiana, 4 of 7 Pima -- out of
#    every RFMix run without an error. An OUTLIERS list is 2 cols (super-pop,
#    sample_id); bcftools -S expects one ID per line, so col 2 goes to a temp
#    file first.
HGDP1KG_CHR_SRC="${HGDP1KG_PATTERN//__CHR__/$CHR}"
[[ -s "$HGDP1KG_CHR_SRC" ]] || {
    echo "ERROR: missing per-chrom HGDP+1KG file $HGDP1KG_CHR_SRC"
    echo "  HGDP1KG_PATTERN=$HGDP1KG_PATTERN"
    exit 1
}

OUTLIERS_IDS="${TMPDIR}/related_outlier_ids.chr${CHR}.txt"
if [[ "$OUTLIERS" == none ]]; then
    : > "$OUTLIERS_IDS"
    echo "[$(date +%T)] [chr${CHR}] relatedness drop: none (OUTLIERS=none -- relatives kept;" \
         "1e removes duplicates only)"
else
    awk 'NF {print $2}' "$OUTLIERS" > "$OUTLIERS_IDS"
    echo "[$(date +%T)] [chr${CHR}] dropping the $(wc -l < "$OUTLIERS_IDS") IDs listed in ${OUTLIERS}"
fi

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
# -G: count records without printing ~3,400 genotypes per line to do it.
NHGDP1KG=$(bcftools view -G -H "$HGDP1KG_PREP" | wc -l)
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

# 1b. Normalise HGDP+1KG the way MXB and SGDP are normalised: split EVERY
#     multi-allelic record into biallelic ones (-m -any, as 1a/1b do with -m-;
#     -m -snps left SNP+indel records whole, so the SNP allele inside one
#     could not join another cohort's biallelic SNP record) and check REF
#     against GRCh38 (+fixref: flip / swap consistently with the genotypes,
#     drop what cannot be fixed). plink2 below removes the indel halves.
#     FIXREF_MODE (shared with 1a, 1b_sgdp_raw, 1d): `-m flip -d` treats A/T
#     and C/G SNPs as strand-ambiguous and discards them even when REF matches
#     (fixref.c, "skip ambiguous pairs"); these are forward-strand WGS calls,
#     so flip-all keeps them. Every cohort must use the same mode, or one
#     cohort is absent at every palindromic site -- checked below by stamp.
HGDP1KG_NORMALIZE="${HGDP1KG_NORMALIZE:-1}"
FIXREF_MODE="${FIXREF_MODE:-flip-all}"
case "$FIXREF_MODE" in flip|flip-all) ;; *) echo "ERROR: FIXREF_MODE must be flip or flip-all" >&2; exit 1 ;; esac
if [[ "$HGDP1KG_NORMALIZE" == "1" ]]; then
    HG_NORM="${TMPDIR}/hgdp1kg_chr${CHR}.norm.bcf"
    echo "[$(date +%T)] [chr${CHR}] normalise HGDP+1KG: norm -m -any -c w + fixref -m ${FIXREF_MODE} -d"
    bcftools norm -m -any -f "$REF_FA" -c w --threads "$THREADS" -Ou "$HGDP1KG_CHR" \
      | bcftools +fixref -Ob -o "$HG_NORM" -- -f "$REF_FA" -m "$FIXREF_MODE" -d
    bcftools index -f "$HG_NORM"
    echo "[$(date +%T)] [chr${CHR}] HGDP+1KG after normalisation: $(bcftools view -G -H "$HG_NORM" | wc -l) records (was ${NHGDP1KG})"
    HGDP1KG_CHR="$HG_NORM"
fi

# 2. Merge HGDP+1KG (chr$CHR) + lifted MXB (chr$CHR) [+ SGDP if present]
#    -- column-wise sample join. SGDP contigs are renamed to chr-prefixed
#    if the prepped BCF came from a bare-numeric source, matching the
#    HGDP+1KG/MXB convention so the merge and SHAPEIT5 --region line up.
MERGE_INPUTS=("$HGDP1KG_CHR" "$MXB_LIFTED")
MERGE_LABELS=(HGDP1KG MXB)

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
    MERGE_LABELS+=("$label")
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

# Every prepped cohort must have been through the same fixref mode as the
# HGDP+1KG step above, read from the stamp its prep wrote. A file built before
# the stamps existed was built with `-m flip -d` (no palindromic SNPs). A
# mismatch is fatal: that cohort would be absent at every A/T and C/G site --
# a whole-cohort gap of exactly the kind that made the November batch clusters.
stamp_fixref () {   # stamp_fixref <stamp file> -> mode, or "flip (no stamp)"
    if [[ -s "$1" ]]; then
        awk '{for (i = 1; i <= NF; i++) if ($i ~ /^fixref=/) {sub(/^fixref=/, "", $i); print $i; f = 1; exit}}
             END {if (!f) print "flip (stamp without fixref=)"}' "$1"
    else
        echo "flip (no stamp)"
    fi
}
_fx_bad=0
for _lab in "${MERGE_LABELS[@]}"; do
    case "$_lab" in
        HGDP1KG) continue ;;
        MXB)     _stamp="${OUTDIR}/mxb_lifted.recipe" ;;
        SGDP)    _stamp="${SGDP_PREPPED}.recipe" ;;
        ALMARRI) _stamp="${ALMARRI_PREPPED}.recipe" ;;
        *)       continue ;;
    esac
    _m=$(stamp_fixref "$_stamp")
    echo "[$(date +%T)] [chr${CHR}] ${_lab}: fixref ${_m} (${_stamp})"
    if [[ "$_m" != "$FIXREF_MODE" ]]; then _fx_bad=1; fi
done
if [[ "$_fx_bad" == "1" && "${FIXREF_MIXED_OK:-0}" != "1" ]]; then
    echo "ERROR: a cohort was prepped with a different fixref mode than FIXREF_MODE=${FIXREF_MODE}." >&2
    echo "       Re-run its prep with the same FIXREF_MODE (MXB: 1a, i.e. RERUN_MXB=1 in the" >&2
    echo "       rebuild driver; SGDP: 1b_sgdp_raw_prep; Almarri: 1d), or FIXREF_MIXED_OK=1." >&2
    exit 1
fi

# 2b. Split every multi-allelic record still in an input (norm -m -any), so
#     each ALT allele is its own biallelic record in every cohort. bcftools
#     merge -m none joins records only when REF and ALT match exactly, so a
#     SNP allele inside a cohort's multi-allelic or SNP+indel record stays a
#     separate record from another cohort's biallelic record at that
#     position. With the 0/0 fill below, the cohort that DOES carry the SNP
#     would then be filled 0/0 on the other cohort's record. After the split,
#     a position with two SNP ALTs is two records and 3b removes it whole, as
#     before. HGDP+1KG was split in 1b; 1a/1b split with -m-, while 1d split
#     only SNPs, so this mostly confirms there is nothing left to split.
for _i in "${!MERGE_INPUTS[@]}"; do
    _lab="${MERGE_LABELS[$_i]}"
    if [[ "$_lab" == "HGDP1KG" && "$HGDP1KG_NORMALIZE" == "1" ]]; then continue; fi
    _nm=$(bcftools view -G -H -m3 "${MERGE_INPUTS[$_i]}" | wc -l)
    if [[ "$_nm" -gt 0 ]]; then
        _out="${TMPDIR}/${_lab}.split.chr${CHR}.bcf"
        bcftools norm -m -any --threads "$THREADS" -Ob -o "$_out" "${MERGE_INPUTS[$_i]}"
        bcftools index -f "$_out"
        echo "[$(date +%T)] [chr${CHR}] ${_lab}: split ${_nm} multi-allelic records (norm -m -any)"
        MERGE_INPUTS[$_i]="$_out"
    else
        echo "[$(date +%T)] [chr${CHR}] ${_lab}: no multi-allelic records"
    fi
done

# SITE PROVENANCE, logged at three stages (raw merge, after QC, soft-union)
# BEFORE anything is decided from it: per cohort, how many of its file's
# records reach the panel, how many are unique to that study, and -- for every
# panel record it has no call at -- whether its file lacks the position or has
# it with alleles the merge could not join (swap / flip / other ALT / other
# REF). TSVs in ${OUTDIR}/provenance/, summary lines in this log. The cohort
# site lists are taken from the ORIGINAL files, before any 0/0 fill below.
PROV_DIR="${OUTDIR}/provenance"
mkdir -p "$PROV_DIR"
PROV_ARGS=()
for _i in "${!MERGE_INPUTS[@]}"; do
    _lab="${MERGE_LABELS[$_i]}"
    _st="${TMPDIR}/sites.${_lab}.chr${CHR}.tsv"
    bcftools query -f '%POS\t%REF\t%ALT\n' "${MERGE_INPUTS[$_i]}" > "$_st"
    PROV_ARGS+=(--cohort "${_lab}=${_st}")
done
provenance () {   # provenance <stage> <bcf>
    local st="${TMPDIR}/sites.panel_$1.chr${CHR}.tsv"
    bcftools query -f '%POS\t%REF\t%ALT\n' "$2" > "$st"
    python3 -I "${REPO_DIR}/scripts/site-provenance.py" --stage "$1" --chrom "chr${CHR}" \
        "${PROV_ARGS[@]}" --panel "$st" \
        --out-cohorts "${PROV_DIR}/chr${CHR}.$1.cohorts.tsv" \
        --out-patterns "${PROV_DIR}/chr${CHR}.$1.patterns.tsv"
}

# ABSENT = HOM-REF, for cohorts whose file is a COMPLETE variant callset.
# A record missing from such a file means no sample in that cohort carries the
# ALT, so its samples are 0/0 there -- the assumption Almarri's per-sample
# merge already makes with missing2ref. Without this, `bcftools merge` leaves
# the whole cohort ./. at every such record and SHAPEIT5 has to guess it.
#
# Only for complete callsets: HGDP+1KG v2 (Koenig: every variant passing
# their site QC, no MAF cut), MXB (GATK joint call, PASS, singletons kept),
# Almarri (per-sample + missing2ref). NOT the current SGDP release
# (GLIMPSE-imputed against 1KG phase 3: a missing site is one 1KG p3 lacks,
# and the same person is non-ref there 28% of the time) and NOT a
# MAF-filtered HGDP+1KG source such as filter1 (missing can mean "rare").
#
# The case where "absent" is NOT "nobody carries it" -- the source's own site
# QC dropped a real variant -- is caught in 3a: if another cohort of the same
# super-pop carries the allele at a frequency at which this cohort's samples
# could not all lack it, the position is dropped instead of filled.
#
# How: a union of every input's records, carrying ONE placeholder sample whose
# GT is ./., is merged with the cohort's file using -0, which fills only the
# records the file lacks; ./. calls inside the file stay missing. Then the
# placeholder is removed. Tested on bcftools 1.19 (cluster: 1.20): with a
# sites-only union (view -G) instead, the records a cohort lacked came out
# with no genotype at all ('.'), not 0/0 -- -0 fills only records at which
# some input has FORMAT/GT. Off unless FILL_REF_COHORTS names the cohorts.
FILL_REF_COHORTS="${FILL_REF_COHORTS:-}"
FILLED_COHORTS=""
if [[ -n "$FILL_REF_COHORTS" ]]; then
    for _c in $FILL_REF_COHORTS; do
        case "$_c" in
            HGDP1KG|MXB|SGDP|ALMARRI) ;;
            *) echo "ERROR: FILL_REF_COHORTS has unknown cohort '${_c}' (valid: HGDP1KG MXB SGDP ALMARRI)" >&2; exit 1 ;;
        esac
        if [[ " ${MERGE_LABELS[*]} " != *" ${_c} "* ]]; then
            echo "[$(date +%T)] [chr${CHR}] FILL_REF_COHORTS names ${_c}, which is not in this run -- nothing to fill"
        fi
    done
    # SGDP is a complete callset only when it comes from the raw per-sample
    # calls (1b_sgdp_raw_prep, stamp "sgdp_raw..."); the imputed phased_data2021
    # release lacks every site 1KG phase 3 lacks, so filling it would invent 0/0s.
    if [[ " ${FILL_REF_COHORTS} " == *" SGDP "* && " ${MERGE_LABELS[*]} " == *" SGDP "* ]]; then
        _sg=$(cat "${SGDP_PREPPED}.recipe" 2>/dev/null || echo "none")
        [[ "$_sg" == sgdp_raw* ]] || { echo "ERROR: FILL_REF_COHORTS includes SGDP but ${SGDP_PREPPED}" \
            "was not built by 1b_sgdp_raw_prep_clm.sh (recipe stamp: '${_sg}')" >&2; exit 1; }
    fi
    FILL_PH="__fill_placeholder__"
    _usites=()
    for _i in "${!MERGE_INPUTS[@]}"; do
        _s="${TMPDIR}/sites_only.${MERGE_LABELS[$_i]}.chr${CHR}.bcf"
        bcftools view -G --threads "$THREADS" -Ob -o "$_s" "${MERGE_INPUTS[$_i]}"
        bcftools index -f "$_s"
        _usites+=("$_s")
    done
    UNION_RAW="${TMPDIR}/union_sites_raw.chr${CHR}.bcf"
    UNION_SITES="${TMPDIR}/union_sites.chr${CHR}.bcf"
    bcftools merge -m none --threads "$THREADS" -Ob -o "$UNION_RAW" "${_usites[@]}"
    _uhdr=$(bcftools view -h "$UNION_RAW") || { echo "ERROR: cannot read header of $UNION_RAW" >&2; exit 1; }
    {
        grep '^##' <<<"$_uhdr"
        grep -q '^##FORMAT=<ID=GT,' <<<"$_uhdr" \
            || echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$FILL_PH"
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t.\t.\t.\tGT\t./.\n' "$UNION_RAW"
    } | bcftools view -Ob -o "$UNION_SITES"
    bcftools index -f "$UNION_SITES"
    unset _uhdr
    echo "[$(date +%T)] [chr${CHR}] union of ${#MERGE_INPUTS[@]} inputs: $(bcftools view -G -H "$UNION_SITES" | wc -l) records (for the 0/0 fill)"
    for _i in "${!MERGE_INPUTS[@]}"; do
        _lab="${MERGE_LABELS[$_i]}"
        [[ " ${FILL_REF_COHORTS} " == *" ${_lab} "* ]] || continue
        _out="${TMPDIR}/${_lab}.absent_as_ref.chr${CHR}.bcf"
        _before=$(bcftools view -G -H "${MERGE_INPUTS[$_i]}" | wc -l)
        bcftools merge -m none -0 --threads "$THREADS" -Ou "$UNION_SITES" "${MERGE_INPUTS[$_i]}" \
          | bcftools view -I -s "^${FILL_PH}" --threads "$THREADS" -Ob -o "$_out"
        bcftools index -f "$_out"
        _after=$(bcftools view -G -H "$_out" | wc -l)
        echo "[$(date +%T)] [chr${CHR}] ${_lab}: ${_before} records in its file; ${_after} after" \
             "filling the $(( _after - _before )) records it lacks as 0/0"
        MERGE_INPUTS[$_i]="$_out"
        FILLED_COHORTS+="${_lab} "
    done
fi

# -m none: never build multi-allelic records. Two cohorts that saw different
# ALTs at one position stay two records, which step 3b then removes together.
echo "[$(date +%T)] [chr${CHR}] bcftools merge -m none (${#MERGE_INPUTS[@]} inputs)"
bcftools merge -m none --threads "$THREADS" \
    -Ob -o "$MERGED" \
    "${MERGE_INPUTS[@]}"
bcftools index "$MERGED"
provenance merged "$MERGED"

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

# 3a. Cohort x super-pop groups (scripts/build-missingness-groups.py): cohort
#     = the FILE a sample was merged from, super-pop from sample_groups.tsv /
#     the SGDP manifest (Almarri = MEN, MXB = AMR). Built once, here -- the
#     sample set does not change after QC -- and used by the fill check below,
#     the soft-union groups (step 4) and the ADMIXTURE site list (5b).
SGDP_MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
MISS_GROUPS="${TMPDIR}/miss_groups.chr${CHR}.tsv"
PANEL_SAMPLES="${TMPDIR}/panel_samples.chr${CHR}.txt"
bcftools query -l "$QCED" > "$PANEL_SAMPLES"
bcftools query -l "$MXB_LIFTED" > "${TMPDIR}/cohort_mxb.chr${CHR}.txt"
COHORT_ARGS=(--cohort "MXB=${TMPDIR}/cohort_mxb.chr${CHR}.txt")
# Membership from THIS run's cohorts, not from whether a *_kept.txt exists:
# those files outlive the run that wrote them.
if [[ " ${MERGE_LABELS[*]} " == *" SGDP "* ]];    then COHORT_ARGS+=(--cohort "SGDP=${TMPDIR}/sgdp_kept.txt"); fi
if [[ " ${MERGE_LABELS[*]} " == *" ALMARRI "* ]]; then COHORT_ARGS+=(--cohort "ALMARRI=${TMPDIR}/almarri_kept.txt"); fi
python3 -I "${REPO_DIR}/scripts/build-missingness-groups.py" \
    --samples "$PANEL_SAMPLES" "${COHORT_ARGS[@]}" \
    --default-cohort HGDP1KG --groups "$SAMPLE_GROUPS" --sgdp-manifest "$SGDP_MANIFEST" \
    > "$MISS_GROUPS"
# bcftools +fill-tags -S rejects sample names shorter than 3 characters
# ("Could not parse the file": an off-by-one in its parser, 1.19 and 1.20).
# Every real ID is far longer; fail here, clearly, if that ever changes.
_short=$(awk 'length($1) < 3' "$MISS_GROUPS" | wc -l)
[[ "$_short" -eq 0 ]] || { echo "ERROR: ${_short} sample name(s) shorter than 3 characters;" \
    "bcftools +fill-tags -S cannot read them:" >&2; awk 'length($1) < 3' "$MISS_GROUPS" >&2; exit 1; }

#     FILL CHECK. For each filled cohort, every record its file lacked (so it
#     now reads 0/0 there) is tested against the OTHER cohorts' groups of the
#     same super-pop: when one carries >= FILL_GUARD_MIN_AC copies at a
#     frequency at which this cohort's haplotypes would all lack the allele
#     with probability < FILL_GUARD_MAX_P (both groups >= FILL_GUARD_MIN_N
#     samples), "absent" cannot mean "nobody carries it" -- the source's own
#     QC dropped a real site -- and the position is dropped with 3b's instead
#     of entering the panel as a cohort-wide 0/0. Counts per cohort and
#     super-pop -> ${PROV_DIR}/chr${CHR}.fill_check.tsv.
FILL_GUARD="${FILL_GUARD:-1}"
FILL_GUARD_MIN_N="${FILL_GUARD_MIN_N:-20}"
FILL_GUARD_MIN_AC="${FILL_GUARD_MIN_AC:-5}"
FILL_GUARD_MAX_P="${FILL_GUARD_MAX_P:-1e-3}"
FILL_DROP="${TMPDIR}/fill_contradictions.chr${CHR}.tsv"
: > "$FILL_DROP"
if [[ -n "$FILLED_COHORTS" && "$FILL_GUARD" == "1" ]]; then
    GUARD_ARGS=()
    for _c in $FILLED_COHORTS; do
        GUARD_ARGS+=(--cohort-sites "${_c}=${TMPDIR}/sites.${_c}.chr${CHR}.tsv")
    done
    GUARD_FMT=$(python3 -I "${REPO_DIR}/scripts/fill-contradictions.py" --groups "$MISS_GROUPS" --print-format)
    echo "[$(date +%T)] [chr${CHR}] fill check: ${FILLED_COHORTS}(min_n ${FILL_GUARD_MIN_N}," \
         "min_ac ${FILL_GUARD_MIN_AC}, max_p ${FILL_GUARD_MAX_P})"
    bcftools +fill-tags "$QCED" -Ou -- -S "$MISS_GROUPS" -t AC,AN \
      | bcftools query -f "$GUARD_FMT" \
      | python3 -I "${REPO_DIR}/scripts/fill-contradictions.py" --groups "$MISS_GROUPS" \
            --filled "$FILLED_COHORTS" "${GUARD_ARGS[@]}" \
            --min-n "$FILL_GUARD_MIN_N" --min-ac "$FILL_GUARD_MIN_AC" --max-p "$FILL_GUARD_MAX_P" \
            --chrom "chr${CHR}" --out "$FILL_DROP" --summary "${PROV_DIR}/chr${CHR}.fill_check.tsv"
    cat "${PROV_DIR}/chr${CHR}.fill_check.tsv"
elif [[ -n "$FILLED_COHORTS" ]]; then
    echo "[$(date +%T)] [chr${CHR}] FILL_GUARD=${FILL_GUARD}: the 0/0 fills of ${FILLED_COHORTS}are NOT checked"
fi
N_FILL_DROP=$(wc -l < "$FILL_DROP")

# 3b. Drop every position that still holds more than one SNP record (a split
#     multi-allelic site: two cohorts saw different ALTs at it, or a cohort's
#     own multi-allelic record was split), together with the fill check's
#     positions. The panel is biallelic by design, and at a split position
#     each cohort is called on one record and missing on the other.
MULTI_POS="${TMPDIR}/multirecord_positions.chr${CHR}.tsv"
bcftools query -f '%CHROM\t%POS\n' "$QCED" | uniq -d > "$MULTI_POS"
N_MULTI=$(wc -l < "$MULTI_POS")
N_MULTI_REC=0
if [[ "$N_MULTI" -gt 0 ]]; then
    N_MULTI_REC=$(bcftools query -f '%CHROM\t%POS\n' "$QCED" \
                  | awk 'NR==FNR {d[$1"\t"$2]=1; next} ($1"\t"$2) in d' "$MULTI_POS" - | wc -l)
fi
DROP_POS="${TMPDIR}/drop_positions.chr${CHR}.tsv"
sort -k1,1 -k2,2n -u "$MULTI_POS" "$FILL_DROP" > "$DROP_POS"
N_DROP=$(wc -l < "$DROP_POS")
if [[ "$N_DROP" -gt 0 ]]; then
    bgzip -f "$DROP_POS"
    tabix -s1 -b2 -e2 -f "${DROP_POS}.gz"
    QCED_BI="${QCED_PREFIX}.biallelic.bcf"
    bcftools view -T "^${DROP_POS}.gz" --threads "$THREADS" -Ob -o "$QCED_BI" "$QCED"
    bcftools index -f "$QCED_BI"
    QCED="$QCED_BI"
fi
echo "[$(date +%T)] [chr${CHR}] dropped ${N_DROP} positions: ${N_MULTI} holding more than one record" \
     "(${N_MULTI_REC} records: split multi-allelic sites); ${N_FILL_DROP} where a 0/0 fill was contradicted"
provenance qc "$QCED"

# 4. Per-superpop MAF soft-union. A site is kept if, in ANY group with
#    >= MIN_SUBPOP_N samples, AN>0 and MAF >= LAI_MAF, written without
#    division as AC >= MAF*AN and (AN-AC) >= MAF*AN -- the rule this step has
#    always applied. What changed is HOW: one `+fill-tags -S` pass computes
#    every group's AC/AN (scripts/softunion-sites.py applies the rule) where
#    the loop it replaces read the whole panel once per group; with the
#    unfiltered HGDP+1KG source the panel has ~4.5x the records it had with
#    filter1. Same result: on a test panel both ways kept the same sites.
#
#    SOFTUNION_GROUPS -- WHOSE samples define a super-pop's MAF:
#      superpop       (default) every cohort's samples, pooled by super-pop
#                     (HGDP1KG_AMR + MXB_AMR + SGDP_AMR -> AMR; Almarri ->
#                     MEN), from the 3a groups; samples with no super-pop are
#                     in no group.
#      sample_groups  the July/November rule: sample_groups.tsv only, i.e.
#                     HGDP+1KG and MXB. SGDP and Almarri samples are in no
#                     group, so a variant only they carry can never pass.
#    With FILL_REF_COHORTS, AN counts the filled 0/0 genotypes too, so a
#    cohort-unique variant's MAF is its frequency in the whole super-pop, not
#    in the one cohort that has data there.
#
#    NOTE: no F_MISSING filter here. AN>0 skips sites a group has no calls at,
#    and SHAPEIT5 imputes residual missing genotypes. The ADMIXTURE marker set
#    gets its own missingness rule in 5b. Groups below MIN_SUBPOP_N are
#    skipped -- their MAF estimates are too noisy to constrain the panel.
SOFTUNION_GROUPS="${SOFTUNION_GROUPS:-superpop}"
SU_GROUPS="${TMPDIR}/softunion_groups.chr${CHR}.tsv"
case "$SOFTUNION_GROUPS" in
    superpop)
        awk -F'\t' '{sp = $2; sub(/^[^_]*_/, "", sp); if (sp != "NA") print $1 "\t" sp}' \
            "$MISS_GROUPS" > "$SU_GROUPS" ;;
    sample_groups)
        awk 'NR==FNR {s[$1] = 1; next} ($1 in s) {print $1 "\t" $2}' \
            "$PANEL_SAMPLES" "$SAMPLE_GROUPS" > "$SU_GROUPS" ;;
    *)  echo "ERROR: SOFTUNION_GROUPS must be superpop or sample_groups, got '${SOFTUNION_GROUPS}'" >&2
        exit 1 ;;
esac
[[ -s "$SU_GROUPS" ]] || { echo "ERROR: no sample has a soft-union group (${SOFTUNION_GROUPS})" >&2; exit 1; }
echo "[$(date +%T)] [chr${CHR}] per-superpop soft-union: MAF >= ${LAI_MAF}, groups = ${SOFTUNION_GROUPS}" \
     "($(wc -l < "$SU_GROUPS") of $(wc -l < "$PANEL_SAMPLES") samples in a group)"
SU_FMT=$(python3 -I "${REPO_DIR}/scripts/softunion-sites.py" --groups "$SU_GROUPS" \
             --min-n "$MIN_SUBPOP_N" --print-format)
SITELIST="${TMPDIR}/softunion_sites.tsv"
bcftools +fill-tags "$QCED" -Ou -- -S "$SU_GROUPS" -t AC,AN \
  | bcftools query -f "$SU_FMT" \
  | python3 -I "${REPO_DIR}/scripts/softunion-sites.py" --groups "$SU_GROUPS" \
        --min-n "$MIN_SUBPOP_N" --maf "$LAI_MAF" --chrom "chr${CHR}" \
        --out "$SITELIST" --counts "${PROV_DIR}/chr${CHR}.softunion_groups.tsv"

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
NSOFT=$(bcftools view -G -H "$SOFTUNION" | wc -l)
echo "[$(date +%T)] [chr${CHR}] soft-union BCF: ${NSOFT} records"
[[ "$NSOFT" -gt 0 ]] || { echo "ERROR: subset produced empty BCF"; exit 1; }

provenance softunion "$SOFTUNION"

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
#     Groups = the 3a cohort x super-pop groups. Only the ADMIXTURE marker set
#     uses this list -- the soft-union itself, and so the phased LAI panel, is
#     unchanged.
# NOTE ON THE THRESHOLD. In the November panel a cohort was either fully
# called or fully absent at a site -- every SGDP sample lacked exactly 39.55%
# of chr22's soft-union, every MXB sample 61.9%, every Almarri sample 30.4%,
# every HGDP/1KG sample 12.68% (sd 0.000) -- because absence came from a
# cohort's SOURCE FILE lacking the site, so any per-group threshold between 0
# and 1 selected the same sites: those every cohort has. With FILL_REF_COHORTS
# a filled cohort has no such gaps left (absent = 0/0, or the site is dropped
# by the 3a fill check), and this rule acts on the cohorts that are NOT filled
# (the GLIMPSE-imputed SGDP release) and on genuinely uncalled genotypes.
ADMIX_GROUP_MAXMISS="${ADMIX_GROUP_MAXMISS:-0.15}"
ADMIX_TOTAL_MAXMISS="${ADMIX_TOTAL_MAXMISS:-0.07}"
ADMIX_MISS_MIN_N="${ADMIX_MISS_MIN_N:-5}"
ADMIX_SITES="${OUTDIR}/merged_chr${CHR}.admix_sites.bed"
ADMIX_SUMMARY="${OUTDIR}/merged_chr${CHR}.admix_sites.summary.tsv"
echo "[$(date +%T)] [chr${CHR}] ADMIXTURE site list: groups <= ${ADMIX_GROUP_MAXMISS} missing" \
     "(n >= ${ADMIX_MISS_MIN_N}), overall <= ${ADMIX_TOTAL_MAXMISS}"
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
    # Stamp AFTER the file is in place (see IDEMPOTENT RE-RUN above).
    echo "$PUB_RECIPE" > "${PUB}.recipe.partial"
    mv "${PUB}.recipe.partial" "${PUB}.recipe"
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
