#!/usr/bin/env bash
#SBATCH --job-name=1d_almarri
#SBATCH --partition=atkinson,mhgcp
# DO NOT REMOVE mhgcp-t* FROM THIS EXCLUSION TO GAIN CONCURRENCY.
#
# Those nodes' CPUs lack instructions the SHAPEIT5 / bcftools builds in this
# conda env are compiled with. A task dispatched to one dies instantly with
# SIGILL -- exit 132, ~3 seconds wall clock, ~2 MB RSS, and an EMPTY log,
# because the binary never reaches its first write.
#
# Measured 2026-11-07 on a 594-task chunk array run with the exclusion lifted:
#   270 / 270 tasks that landed on t02,t05,t06,t07,t08,t10,t11  FAILED
#     0 / 57  tasks that landed on a01, m00                     failed
#
# They look attractive because sinfo shows them idle with free CPUs. They are
# idle because nothing that needs these instruction sets can run there. An
# earlier 1d_almarri failure on t02 -- empty logs, 4 MB RSS, never explained
# -- was the same thing.
#
# mhgcp-c02 WAS on this list and has been removed deliberately. It has the
# same specs as c00 and c01 (56 CPUs, 232 GB), both of which run this software
# without trouble and were never excluded, and no evidence ties c02 to the
# instruction-set problem above -- no task ever landed on it. The exclusion
# looks like it dates from c02 simply being busy. Re-add it only if a task
# actually fails there, and write down what happened.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# Absolute log paths, matching 1b. A relative `logs/%A_%a.out` resolves against
# the submit directory, and if that directory has no logs/ Slurm cannot create
# the file and kills the task before the script runs -- a failure with no log
# to explain it, because the log is what failed.
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1d_almarri_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1d_almarri_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu
# ============================================================================
# 1d_almarri_prep_clm.sh   --   array 1-22, one task per chromosome
#
# Turn the 137 Almarri et al. 2021 per-sample phased VCFs into ONE per-chr BCF
# the 1c merge can take, alongside HGDP+1KG, MXB and SGDP.
#
# THE TRAP THIS SCRIPT EXISTS TO HANDLE. These are PER-SAMPLE calls from Long
# Ranger. A single-sample VCF records only the sites where that sample carries
# a non-reference allele; every hom-ref site is simply absent. `bcftools merge`
# turns "absent" into ./. , so a naive merge of 137 of them yields a panel that
# is almost entirely missing. Almarri et al. say so in their own methods:
#
#   "Because variants were called for each sample separately, homozygous ref
#    alleles were assigned as missing when joined with other samples, causing
#    an artificially high missingness rate. To resolve this, we used
#    `bcftools plugin missing2ref` to switch missing alleles to ref."
#
# So absent must become 0/0. This script does that with `bcftools merge -0`,
# and refuses to emit a BCF whose missingness says it did not happen (step 4).
#
# THE PAPER'S QC, which v1-v3 mirrored and v4 does NOT:
#   bcftools view -f PASS ; GQ > 20 ; then missing2ref
#   bcftools norm -c ws -m -snps -m -indels --fasta-ref <GRCh38>
# Dropping a low-quality call and then filling every gap as 0/0 turns that
# call into a confident hom-ref. Measured on chr22 (12 Nov 2026, 6 samples):
# 24-28% of each sample's ALT calls -- ~97% of them heterozygous -- were
# non-PASS (21-24%) or GQ <= 20 (2-4%), so v3 cut Almarri's heterozygosity
# and ALT frequencies by about a quarter: a cohort signal, not ancestry.
#
# v4: a low-quality call becomes MISSING (./.), never 0/0 -- as SGDP's raw prep
# does with FL < 1 -- and SHAPEIT5 fills it from LD like any other gap:
#   * FILTER not PASS (after removing ALMARRI_PASS_FILTERS, for a filter about
#     the PHASE only -- none of Long Ranger's is; see the setting), or
#   * GQ <= MIN_GQ, or GQ missing.
# Then `merge -0` fills only what is ABSENT from a sample's file (no call at
# all) as 0/0; a missing call stays missing. Site mask, as for SGDP: where
# more Almarri genotypes are missing than ALT, all of Almarri is set missing
# there (ALMARRI_SITE_MASK=1). Per-chromosome mean hets and missing per sample
# are logged from `bcftools stats -s`, so the effect is visible in every run.
#
# The files are already GRCh38, so there is no liftover here -- unlike SGDP
# (1b), which is hg19 and must be lifted first.
#
# WHY THESE ARE RE-PHASED LATER. Long Ranger phased each sample INDEPENDENTLY.
# Appending an independently phased cohort to an independently phased panel
# creates cohort-specific switch errors, which RFMix reads as ancestry
# switches. 1c hands the whole merged panel to SHAPEIT5 for joint re-phasing,
# so the incoming phase is deliberately not relied upon.
#
# Usage:
#   sbatch --array=1-22%8 1d_almarri_prep_clm.sh
#   ALMARRI_DIR=~/almarri_vcfs sbatch --array=1-22%8 1d_almarri_prep_clm.sh
#
# Index the source VCFs ONCE before the array (they ship without indexes and
# this script region-queries every one of them):
#   bash 1d_almarri_prep_clm.sh --index
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
ALMARRI_DIR="${ALMARRI_DIR:-${HOME}/almarri_vcfs}"
VCF_DIR="${VCF_DIR:-${ALMARRI_DIR}/vcfs}"
REF_FA="${REF_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
MIN_GQ="${MIN_GQ:-20}"
# Long Ranger FILTER values to treat as passing (space-separated, beyond PASS):
# only a filter about the PHASE of a call, not its genotype, belongs here --
# SHAPEIT5 re-phases the whole panel. None does, so the default is empty.
# Long Ranger 2.2.2's filters, and one sample's chr22 non-PASS ALT calls
# (16,629; ~96% het; 12 Nov 2026):
#   10X_ALLELE_FRACTION_FILTER          AO < 2 or ALT fraction < 0.15     ~8.2k
#   10X_RESCUED_MOLECULE_HIGH_DIVERSITY reads placed by barcode rescue in
#                                       high-divergence sequence          ~5.6k
#   10X_PHASING_INCONSISTENT            alleles that do not segregate with
#                                       the linked-read haplotypes        ~3.6k
#   10X_QUAL_FILTER                     QUAL <= 15 (or hom-ALT QUAL < 50) ~2.1k
#   (counts include calls carrying two or three of them)
# PHASING_INCONSISTENT is not "phase uncertain": an ALT on both haplotypes'
# molecules is the signature of a false het, typically paralogous sequence in
# segmental duplications. Its GQ is high (3,583 of 3,593 > 20) because GQ comes
# from read counts, which a paralog inflates. All four are genotype filters.
ALMARRI_PASS_FILTERS="${ALMARRI_PASS_FILTERS:-}"
# Site mask: where more Almarri genotypes are missing than ALT, set all of
# Almarri missing at that site (1 = on, as SGDP_SITE_MASK).
ALMARRI_SITE_MASK="${ALMARRI_SITE_MASK:-1}"
# Strand/REF mode of the fixref step, shared with 1a, 1b_sgdp_raw and 1c (see
# 1a: `-m flip -d` discards every A/T and C/G SNP; flip-all keeps them).
FIXREF_MODE="${FIXREF_MODE:-flip-all}"
# Above this, absent sites were plainly not filled (an unfilled merge of
# single-sample files is ~99% missing); see step 4. v4 leaves low-quality calls
# missing on purpose, so the bar is higher than v3's 0.10.
MAX_MISSING="${MAX_MISSING:-0.30}"

set -euo pipefail

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

[[ -d "$VCF_DIR" ]] || { echo "ERROR: missing $VCF_DIR (run scripts/fetch-almarri-vcfs.sh)"; exit 1; }

# --------------------------------------------------------------------------
# PREFLIGHT: the bcftools plugins (setGT, fixref), checked BEFORE the
# expensive loop. Step 2 prepares 137 per-sample BCFs, which is the slow part;
# a plugin found missing there throws that work away on every task. bcftools
# ships plugins as separate .so files and finds them via BCFTOOLS_PLUGINS,
# which a conda env does not always set.
# --------------------------------------------------------------------------
# `bcftools plugin -l | grep -q` has the same SIGPIPE inversion as the contig
# check below: grep -q exits at the match, bcftools dies, pipefail returns 141
# and the test reports the plugin MISSING when it is present. Listed once into
# a variable instead.
PLUGINS=$(bcftools plugin -l 2>/dev/null || true)
for _pl in setGT fixref; do
    if ! grep -qw -- "$_pl" <<<"$PLUGINS"; then
        echo "ERROR: the bcftools ${_pl} plugin is not available." >&2
        echo "       bcftools: $(command -v bcftools)" >&2
        echo "       BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}" >&2
        echo "       Try:  export BCFTOOLS_PLUGINS=\$(dirname \$(command -v bcftools))/../libexec/bcftools" >&2
        exit 4
    fi
done
echo "[almarri] setGT and fixref plugins available"

# --------------------------------------------------------------------------
# --index : build the .tbi/.csi the array tasks need. Run ONCE, not per chr.
# The Sanger share publishes the .vcf.gz without indexes, and without an index
# bcftools cannot do the -r chr<N> region query that makes the per-chr array
# tractable. 137 files, so this is a real step, not a formality.
# --------------------------------------------------------------------------
if [[ "${1:-}" == "--index" ]]; then
    n=0; built=0
    for f in "$VCF_DIR"/*.vcf.gz; do
        n=$((n+1))
        if [[ -s "${f}.tbi" || -s "${f}.csi" ]]; then continue; fi
        echo "[almarri] indexing $(basename "$f")"
        bcftools index -t "$f" || { echo "ERROR: could not index $f" >&2; exit 1; }
        built=$((built+1))
    done
    echo "[almarri] ${n} VCFs, ${built} newly indexed"
    exit 0
fi

CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...), or with --index}
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP="${OUTDIR}/tmp/almarri_chr${CHR}"
OUT="${OUTDIR}/almarri_prepped.chr${CHR}.bcf"
mkdir -p "$OUTDIR" "$TMP"

# Reuse a finished chromosome only if it READS BACK and was built by the
# current per-sample recipe (stamped in ${OUT}.recipe). "Non-empty with an
# index" let a file through that bcftools could no longer read -- chr22 was
# found at 0 bytes beside its Nov-5 index -- and would let a v1-recipe file
# (norm -c ws, stderr discarded) stand in for a v2 rebuild.
# v3 = v2 + FIXREF_MODE (v2 used `-m flip -d`, which dropped every A/T and
# C/G SNP). The stamp carries `fixref=<mode>`; 1c refuses to merge cohorts
# whose stamps disagree on it. v4 = low-quality calls missing, not 0/0 (see
# the header); every setting that changes a genotype is in the stamp.
_PASSF=$(tr ' ' '\n' <<<"$ALMARRI_PASS_FILTERS" | awk 'NF' | sort -u | paste -sd, -)
OUT_RECIPE="v4 fixref=${FIXREF_MODE} lowq=missing minGQ=${MIN_GQ} pass=PASS${_PASSF:+,${_PASSF}} sitemask=${ALMARRI_SITE_MASK}"
if [[ -s "$OUT" && -s "${OUT}.csi" ]] \
   && bcftools view -h "$OUT" > /dev/null \
   && [[ "$(cat "${OUT}.recipe" 2>/dev/null)" == "$OUT_RECIPE" ]]; then
    echo "[$(date +%T)] [chr${CHR}] $OUT reads back and is recipe ${OUT_RECIPE} -- nothing to do"
    exit 0
fi
if [[ -e "$OUT" ]]; then
    echo "[$(date +%T)] [chr${CHR}] rebuilding $OUT: unreadable, unindexed, or not recipe ${OUT_RECIPE}" \
         "(recipe stamp: '$(cat "${OUT}.recipe" 2>/dev/null || echo none)')"
fi
for f in "$REF_FA" "${REF_FA}.fai"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# --------------------------------------------------------------------------
# 1. Which contig name does this cohort use?
# Long Ranger output is normally chr-prefixed, but the rest of the pipeline
# (MXB, HGDP+1KG, SGDP after 1b) is chr-prefixed by construction, so this is
# checked rather than assumed: a silent mismatch yields an EMPTY region query,
# which looks like "this cohort has no variants here" instead of an error.
# --------------------------------------------------------------------------
# SIGPIPE, and why none of this is a pipeline any more.
#
# `bcftools view -h "$f" | grep -q PATTERN` under `set -o pipefail` does not
# mean what it reads like. grep -q exits at the FIRST match, bcftools is still
# writing a header with thousands of contig lines, bcftools dies of SIGPIPE
# (141), and pipefail makes the PIPELINE return 141 -- so the `if` is FALSE
# even though the pattern was found. Both branches below missed, execution
# fell into the error path, and the error path's own
# `| grep '^##contig' | head -5` died the same way, aborting the task with
# exit 141 before it could print why.
#
# The header is read ONCE into a variable and matched with a here-string. No
# pipe, so no SIGPIPE, and one read instead of three.
shopt -s nullglob
_vcfs=("$VCF_DIR"/*.vcf.gz)
shopt -u nullglob
[[ ${#_vcfs[@]} -gt 0 ]] || { echo "ERROR: no .vcf.gz in $VCF_DIR"; exit 1; }
FIRST_VCF="${_vcfs[0]}"
HDR=$(bcftools view -h "$FIRST_VCF") || {
    echo "ERROR: cannot read the header of $FIRST_VCF" >&2; exit 1; }
if grep -q "^##contig=<ID=chr${CHR}[,>]" <<<"$HDR"; then
    REGION="chr${CHR}"
elif grep -q "^##contig=<ID=${CHR}[,>]" <<<"$HDR"; then
    REGION="${CHR}"
else
    echo "ERROR: neither 'chr${CHR}' nor '${CHR}' in $(basename "$FIRST_VCF") header" >&2
    echo "       contigs present:" >&2
    grep '^##contig' <<<"$HDR" | awk 'NR<=5' >&2
    exit 1
fi
unset HDR
echo "[$(date +%T)] [chr${CHR}] source region name: ${REGION}"

# --------------------------------------------------------------------------
# 2. Per sample: low-quality calls set MISSING (FILTER not PASS once the
#    ALMARRI_PASS_FILTERS are removed; GQ <= MIN_GQ or missing), split
#    multiallelics, left-align, SNVs only, then the SAME REF check MXB (1a) and
#    SGDP (1b_sgdp_raw) get:
#    `bcftools +fixref -- -f <GRCh38> -m ${FIXREF_MODE} -d` (swap REF/ALT
#    consistently with the genotypes; flip-all also keeps A/T and C/G SNPs,
#    which `-m flip -d` discarded).
#
# Recipe v1 used `norm -c ws`, which SETS a mismatched REF from the FASTA and
# leaves ALT and the genotypes as they were -- for a swapped site that writes a
# record whose genotypes no longer mean what they say. It also sent all three
# commands' stderr to /dev/null, so how many REFs it rewrote is unknowable from
# the logs. v2 warns (`-c w`), lets fixref resolve or drop, and keeps stderr in
# the task's .err.
#
# The recipe tag is part of the cached per-sample filename, so a BCF built by
# an older recipe is never reused by a newer one.
# --------------------------------------------------------------------------
# File-name-safe form of the per-sample part of OUT_RECIPE.
SUB_RECIPE="v4-${FIXREF_MODE}-gq${MIN_GQ}-$(printf '%s' "$_PASSF" | md5sum | cut -c1-6)"
# Remove the filters that are about phase only, so such a call reads as PASS.
# Only those this file's header defines: `annotate -x FILTER/X` on a file
# without X is a hard error (tested, bcftools 1.19), which would kill the task
# on any sample whose header lacks one of them.
pass_filters () {   # pass_filters <source vcf>  (stdin -> stdout)
    local src="$1" defined use
    if [[ -z "$_PASSF" ]]; then cat; return; fi
    defined=$(bcftools view -h "$src" | awk -F'[=,]' '/^##FILTER=<ID=/ {print $3}')
    use=$(tr ',' '\n' <<<"$_PASSF" | grep -Fx -f <(printf '%s\n' "$defined") | sed 's/^/FILTER\//' | paste -sd, -)
    if [[ -n "$use" ]]; then bcftools annotate -x "$use" -Ou; else cat; fi
}
# Two setGT passes, not one: in a single -t q expression that mixes the
# site-level FILTER with the per-sample GQ, bcftools 1.19 silently applies only
# the GQ part (tested: three non-PASS calls kept their genotypes).
mask_low_quality () {
    bcftools +setGT -Ou -- -t q -n . -i 'FILTER!="PASS" & FILTER!="."' \
      | bcftools +setGT -Ou -- -t q -n . -i "FMT/GQ<=${MIN_GQ} | FMT/GQ=\".\""
}
LIST="${TMP}/merge_list.txt"
: > "$LIST"
NS=0; NEMPTY=0
for f in "$VCF_DIR"/*.vcf.gz; do
    [[ -s "${f}.tbi" || -s "${f}.csi" ]] || {
        echo "ERROR: $(basename "$f") has no index. Run: bash $0 --index" >&2; exit 2; }
    base=$(basename "$f" .vcf.gz)
    sub="${TMP}/${base}.${SUB_RECIPE}.bcf"
    # Reuse a cached per-sample BCF only if it READS BACK. Testing for the
    # index alone accepts a truncated BCF beside a complete .csi, which is
    # what a killed task leaves -- and these are cached precisely so a rerun
    # after a kill is cheap, so that is the common case, not a rare one.
    if [[ -s "$sub" && -s "${sub}.csi" ]] && bcftools view -h "$sub" >/dev/null 2>&1; then
        :
    else
        rm -f "$sub" "${sub}.csi"
        # Build beside the target and move, so this loop never leaves a
        # half-written pair for the next run to find.
        echo "[$(date +%T)] [chr${CHR}] ${base}: low-quality -> missing (not PASS${_PASSF:+/${_PASSF}}, GQ<=${MIN_GQ}) -> norm -> fixref" >&2
        bcftools view -r "$REGION" -Ou "$f" \
          | pass_filters "$f" \
          | mask_low_quality \
          | bcftools norm -c w -m -snps --fasta-ref "$REF_FA" -Ou \
          | bcftools view --types snps --max-alleles 2 -Ob -o "${sub}.norm.bcf"
        bcftools +fixref "${sub}.norm.bcf" -Ob -o "${sub}.part" \
            -- -f "$REF_FA" -m "$FIXREF_MODE" -d
        rm -f "${sub}.norm.bcf"
        bcftools index "${sub}.part"
        mv "${sub}.part.csi" "${sub}.csi"
        mv "${sub}.part" "$sub"
    fi
    # Record count from the INDEX, not by reading the file. `bcftools view -H
    # | wc -l` re-reads every per-sample BCF on every run, 137 of them per
    # chromosome, which on a resumed task is most of the work for no new
    # information.
    NREC=$(bcftools index -n "$sub" 2>/dev/null || echo 0)
    if [[ "$NREC" -eq 0 ]]; then
        # A sample with no surviving variant on this chromosome still has to
        # join the merge, or it silently vanishes from the panel on this
        # chromosome only -- giving a cohort whose membership differs per chr.
        NEMPTY=$((NEMPTY+1))
    fi
    echo "$sub" >> "$LIST"
    NS=$((NS+1))
done
echo "[$(date +%T)] [chr${CHR}] prepared ${NS} per-sample BCFs (${NEMPTY} with no variant here)"
[[ "$NS" -gt 0 ]] || { echo "ERROR: no per-sample BCFs built"; exit 1; }

# --------------------------------------------------------------------------
# 3. Merge with -0: a record ABSENT from a sample's file is 0/0 (that sample
#    has no call there); a call set missing in step 2 stays missing -- which
#    `+missing2ref` (v1-v3) would have turned into 0/0 as well. Then the site
#    mask. -m none keeps records biallelic; they were split in step 2 and the
#    panel is biallelic SNVs throughout.
# --------------------------------------------------------------------------
echo "[$(date +%T)] [chr${CHR}] bcftools merge -0 (${NS} samples) + strip to GT + site mask"
MERGED="${TMP}/almarri.merged.bcf"
M0="${TMP}/almarri.merged.unmasked.bcf"
# STRIP EVERYTHING BUT GT.
#
# Long Ranger writes PL, AD, GQ and barcode fields per genotype. Kept, they
# made the merged chr1 8.4 GB and the whole cohort ~100 GB for 137 samples --
# about 34 bytes per genotype where GT alone is a few.
#
# Nothing downstream reads them. SHAPEIT5 phases from GT, the 1c merge only
# joins samples, and bcftools merge across cohorts would otherwise build the
# UNION of every cohort's FORMAT fields, writing Almarri's extras as missing
# for HGDP+1KG, MXB and SGDP -- paying for them across the whole panel, not
# just this cohort.
#
# INFO is dropped for the same reason: these are per-sample call annotations
# whose meaning does not survive a merge anyway.
bcftools merge -m none -0 -l "$LIST" --threads "$THREADS" -Ou \
  | bcftools annotate -x 'INFO,^FORMAT/GT' --threads "$THREADS" -Ob -o "$M0"
bcftools index -f "$M0"

# --------------------------------------------------------------------------
# 4. Prove the absent genotypes were filled -- on the merge BEFORE the site
# mask. Unfilled, a merge of single-sample files is ~99% missing (a record is
# called only in its carriers), the BCF still writes, and the panel would be
# phased looking superficially fine. Measured before the mask because the mask
# blanks exactly the records an unfilled merge leaves mostly missing, so after
# it the failure could pass; and because masked, wholly-missing records -- in
# the segmental duplications the first 20k records of chr22 are full of --
# would fail a correct merge. Here only the low-quality calls are missing.
# --------------------------------------------------------------------------
NSAMP=$(bcftools query -l "$M0" | wc -l)
NREC=$(bcftools view -H "$M0" | wc -l)
[[ "$NREC" -gt 0 ]] || { echo "ERROR: merged chr${CHR} has no records"; exit 1; }
# The producer is wrapped so its SIGPIPE death (head -20000 closes the pipe
# long before a whole chromosome is written) cannot fail the assignment under
# pipefail. Sampling the first 20k records is deliberate: the question is
# whether the absent genotypes were filled at all, and that is visible at once.
FRAC_MISSING=$( { bcftools query -f '[%GT\t]\n' "$M0" || true; } | head -20000 \
  | awk '{for(i=1;i<=NF;i++){t++; if($i=="./."||$i==".|."||$i==".") m++}} END{printf "%.4f", (t>0?m/t:1)}')
echo "[$(date +%T)] [chr${CHR}] merged: ${NSAMP} samples, ${NREC} records, missing GT fraction before the site mask ${FRAC_MISSING}"
if awk -v f="$FRAC_MISSING" -v m="$MAX_MISSING" 'BEGIN{exit !(f>m)}'; then
    echo "ERROR: ${FRAC_MISSING} of genotypes are missing, above MAX_MISSING=${MAX_MISSING}." >&2
    echo "       That is the signature of a merge whose absent genotypes were not" >&2
    echo "       filled: single-sample VCFs carry no hom-ref sites, so a genotype" >&2
    echo "       absent from a sample is ./. unless merge -0 fills it." >&2
    exit 3
fi

# --------------------------------------------------------------------------
# 5. Site mask, as SGDP's: where more Almarri genotypes are missing than ALT,
# the cohort's evidence there is mostly low-quality calls, so all of Almarri is
# set missing rather than half-trusted.
# --------------------------------------------------------------------------
N_SITES_MASK=$(bcftools query -i 'N_MISSING>0 && N_MISSING>COUNT(GT="alt")' -f '%POS\n' "$M0" | wc -l)
if [[ "$ALMARRI_SITE_MASK" == "1" ]]; then
    bcftools +setGT "$M0" -Ob -o "$MERGED" -- -t q -n . -i 'N_MISSING>0 && N_MISSING>COUNT(GT="alt")'
    rm -f "$M0" "${M0}.csi"
    echo "[$(date +%T)] [chr${CHR}] site mask: all of Almarri set missing at ${N_SITES_MASK} of ${NREC} sites (more missing than ALT)"
else
    mv "$M0" "$MERGED"; rm -f "${M0}.csi"
    echo "[$(date +%T)] [chr${CHR}] ALMARRI_SITE_MASK=0: no site mask (${N_SITES_MASK} of ${NREC} sites would qualify)"
fi
bcftools index -f "$MERGED"
# Per-sample heterozygous and missing genotypes (bcftools stats PSC: $6 nHets,
# $14 nMissing), so a cohort-wide shift in either shows in every run's log.
read -r MEAN_HET MEAN_MISS < <(bcftools stats -s - "$MERGED" \
    | awk '/^PSC/ {n++; h += $6; m += $14} END {printf "%.0f %.0f\n", (n ? h / n : 0), (n ? m / n : 0)}')
echo "[$(date +%T)] [chr${CHR}] per sample (mean): ${MEAN_HET} heterozygous, ${MEAN_MISS} missing genotypes"

# Move the index first, then the BCF: $OUT only ever appears once both exist,
# so a kill here cannot leave a BCF that the resume check accepts as finished.
rm -f "${OUT}.recipe"
mv "${MERGED}.csi" "${OUT}.csi"
mv "$MERGED" "$OUT"
echo "$OUT_RECIPE" > "${OUT}.recipe"

# Drop the 137 per-sample intermediates now that the merge is on disk. Each
# chromosome builds its own set and 8 run at once, so leaving them is tens of
# GB of scratch per run for files nothing reads again. KEEP_TMP=1 retains them
# for debugging a specific chromosome.
if [[ -z "${KEEP_TMP:-}" ]]; then
    _sz=$(du -sh "$TMP" 2>/dev/null | cut -f1)
    rm -rf "$TMP"
    echo "[$(date +%T)] [chr${CHR}] cleaned ${_sz:-?} of per-sample intermediates"
else
    echo "[$(date +%T)] [chr${CHR}] KEEP_TMP set -- leaving $TMP"
fi

echo "[$(date +%T)] [chr${CHR}] DONE -> $OUT"
echo "[$(date +%T)] [chr${CHR}] 1c will merge this alongside HGDP+1KG, MXB and SGDP,"
echo "[$(date +%T)] [chr${CHR}] and SHAPEIT5 will re-phase the whole panel jointly."
