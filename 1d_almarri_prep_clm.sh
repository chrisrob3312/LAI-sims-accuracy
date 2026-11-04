#!/usr/bin/env bash
#SBATCH --job-name=1d_almarri
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
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
# So the merge MUST be followed by missing2ref. This script does that, and
# refuses to emit a BCF whose missingness says it did not happen (see step 4).
#
# THE PAPER'S QC, which this mirrors:
#   bcftools view -f PASS ; GQ > 20
#   bcftools norm -c ws -m -snps -m -indels --fasta-ref <GRCh38>
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
# Above this, the merge plainly did not get missing2ref applied; see step 4.
MAX_MISSING="${MAX_MISSING:-0.10}"

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
# PREFLIGHT: the missing2ref plugin, checked BEFORE the expensive loop.
#
# Step 2 prepares 137 per-sample BCFs, which is the slow part. If the plugin
# turns out to be unavailable at step 3, all of that work is thrown away --
# per chromosome, on every task. bcftools ships plugins as separate .so files
# and finds them via BCFTOOLS_PLUGINS, which a conda env does not always set,
# so "is it there" is a real question rather than a formality.
# --------------------------------------------------------------------------
# `bcftools plugin -l | grep -q` has the same SIGPIPE inversion as the contig
# check below: grep -q exits at the match, bcftools dies, pipefail returns 141
# and the test reports the plugin MISSING when it is present. Listed once into
# a variable instead.
PLUGINS=$(bcftools plugin -l 2>/dev/null || true)
if ! grep -q 'missing2ref' <<<"$PLUGINS"; then
    echo "ERROR: the bcftools missing2ref plugin is not available." >&2
    echo "       It is REQUIRED: these are single-sample VCFs, so every site" >&2
    echo "       absent from a sample becomes ./. on merge and the panel is" >&2
    echo "       almost entirely missing without it." >&2
    echo "       bcftools: $(command -v bcftools)" >&2
    echo "       BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}" >&2
    echo "       Try:  export BCFTOOLS_PLUGINS=\$(dirname \$(command -v bcftools))/../libexec/bcftools" >&2
    echo "       then: bcftools plugin -l | grep missing2ref" >&2
    exit 4
fi
echo "[almarri] missing2ref plugin available"

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

if [[ -s "$OUT" && -s "${OUT}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] $OUT exists -- nothing to do"
    exit 0
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
# 2. Per sample: PASS, GQ > MIN_GQ, split multiallelics, left-align, SNVs only.
# `-c ws` makes bcftools norm WARN on a REF mismatch and SET it from the FASTA
# rather than dropping the record silently.
# --------------------------------------------------------------------------
LIST="${TMP}/merge_list.txt"
: > "$LIST"
NS=0; NEMPTY=0
for f in "$VCF_DIR"/*.vcf.gz; do
    [[ -s "${f}.tbi" || -s "${f}.csi" ]] || {
        echo "ERROR: $(basename "$f") has no index. Run: bash $0 --index" >&2; exit 2; }
    base=$(basename "$f" .vcf.gz)
    sub="${TMP}/${base}.bcf"
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
        bcftools view -r "$REGION" -f PASS -i "FMT/GQ>${MIN_GQ}" -Ou "$f" 2>/dev/null \
          | bcftools norm -c ws -m -snps --fasta-ref "$REF_FA" -Ou 2>/dev/null \
          | bcftools view --types snps --max-alleles 2 -Ob -o "${sub}.part" 2>/dev/null
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
# 3. Merge, then missing2ref -- the step that cannot be skipped.
# -m none keeps records biallelic; they were split in step 2 and the panel is
# biallelic SNVs throughout.
# --------------------------------------------------------------------------
echo "[$(date +%T)] [chr${CHR}] bcftools merge (${NS} samples) + missing2ref"
MERGED="${TMP}/almarri.merged.bcf"
bcftools merge -m none -l "$LIST" --threads "$THREADS" -Ou \
  | bcftools +missing2ref -Ob -o "$MERGED" --
bcftools index "$MERGED"

# --------------------------------------------------------------------------
# 4. Prove missing2ref actually ran.
# If the plugin were missing, misnamed, or silently a no-op, the merge would
# still succeed and emit a BCF that is ~99% missing genotypes. That panel
# would then be phased and carried through every downstream stage looking
# superficially fine. Measure it instead of trusting it.
# --------------------------------------------------------------------------
NSAMP=$(bcftools query -l "$MERGED" | wc -l)
NREC=$(bcftools view -H "$MERGED" | wc -l)
[[ "$NREC" -gt 0 ]] || { echo "ERROR: merged chr${CHR} has no records"; exit 1; }
# The producer is wrapped so its SIGPIPE death (head -20000 closes the pipe
# long before a whole chromosome is written) cannot fail the assignment under
# pipefail. Sampling the first 20k records is deliberate: the question is
# whether missing2ref ran at all, and that is visible immediately.
FRAC_MISSING=$( { bcftools query -f '[%GT\t]\n' "$MERGED" || true; } | head -20000 \
  | awk '{for(i=1;i<=NF;i++){t++; if($i=="./."||$i==".|."||$i==".") m++}} END{printf "%.4f", (t>0?m/t:1)}')
echo "[$(date +%T)] [chr${CHR}] merged: ${NSAMP} samples, ${NREC} records, missing GT fraction ${FRAC_MISSING}"
if awk -v f="$FRAC_MISSING" -v m="$MAX_MISSING" 'BEGIN{exit !(f>m)}'; then
    echo "ERROR: ${FRAC_MISSING} of genotypes are missing, above MAX_MISSING=${MAX_MISSING}." >&2
    echo "       That is the signature of a merge whose missing2ref did not" >&2
    echo "       apply: single-sample VCFs carry no hom-ref sites, so every" >&2
    echo "       site absent from a sample becomes ./. on merge." >&2
    echo "       Check that the plugin is available:" >&2
    echo "         bcftools plugin -l | grep missing2ref" >&2
    echo "       and that BCFTOOLS_PLUGINS points at the plugin directory." >&2
    exit 3
fi

# Move the index first, then the BCF: $OUT only ever appears once both exist,
# so a kill here cannot leave a BCF that the resume check accepts as finished.
mv "${MERGED}.csi" "${OUT}.csi"
mv "$MERGED" "$OUT"

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
