#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_sgdp_raw
#SBATCH --partition=atkinson,mhgcp
# Same exclusion as every other stage (reason in 1b_sgdp_liftover_clm.sh):
# SIGILL on the mhgcp-t* CPUs, exit 132 with an empty log.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_sgdp_raw_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_sgdp_raw_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- SGDP from the RAW per-sample WGS calls (vcf_variants)
# ============================================================================
# Replaces 1b_sgdp_liftover_clm.sh, which used phased_data2021: genotypes
# GLIMPSE-imputed against 1KG phase 3, so a site 1KG p3 lacked was missing
# for all of SGDP (39.6% of chr22's November soft-union) and SHAPEIT5 filled
# it in. Same output contract -- ${OUTDIR}/sgdp_prepped.chrN.bcf, hg38,
# chr-prefixed, GT only -- so 1c is unchanged.
#
# Input: scripts/fetch-sgdp-vcf-variants.sh's verified download. One GATK
# UnifiedGenotyper VCF per sample, hs37d5, contigs 1..22, SNPs only, VARIANT
# SITES ONLY (no 0/0 records), FILTER always '.', and a per-call filter level:
#   FL = "filter level in range 0-9 or no value (non-integer: N,?) with zero
#         being least reliable; to threshold at FL=n, use all levels n-9"
# (first sample, SS6004478: 4,109,729 calls; FL 0 9.9%, N 23.9%, 1-9 66.1%).
#
# Per chromosome:
#   1. Map manifest samples to their files (scripts/sgdp-raw-sample-map.py:
#      manifest sample -> SGDP_ID -> Illumina_ID in the public metadata).
#   2. bcftools merge -0 -m none -r N over all of them: a sample with no
#      record where another sample has one is 0/0 -- in a variant-only file
#      no record means no variant call, the same assumption Almarri's 1d
#      merge makes with missing2ref.
#   3. A variant call with FL below SGDP_MIN_FL becomes ./. (missing), NOT
#      0/0: the call is unreliable, which is not evidence of hom-ref. Per
#      sample (`&`, not `&&`, in the expression: `&&` matches across samples).
#      Default 1, the SGDP authors' own: "At 'filter level 1' which we
#      recommend for most analyses ..." and, for their ADMIXTURE analysis,
#      "We used genotypes of at least filter level 1" (Mallick et al. 2016,
#      Nature 538:201, doi:10.1038/nature18964). N and ? ("no value") fail.
#   4. SITE MASK. A variant-only file has no FL for hom-ref sites, so where a
#      site's calls are mostly unreliable (more of SGDP's carriers fail FL
#      than pass), the samples recorded 0/0 there are no more reliable than
#      the calls just removed, and keeping them would leave SGDP looking
#      monomorphic: carriers missing, everyone else 0/0. Such a site is set
#      ./. for ALL of SGDP (SGDP_SITE_MASK=1). It is still in the panel for
#      the other cohorts; 1c's ADMIXTURE site rule sees SGDP as missing there.
#   5. Keep GT only; rename samples to the manifest's names.
#   6. Lift hg19 -> hg38 with the 1a recipe (Picard LiftoverVcf,
#      RECOVER_SWAPPED_REF_ALT), norm -m- + fixref -m ${FIXREF_MODE} -d
#      against hg38, keep chrN, sort -> sgdp_prepped.chrN.bcf, with a .recipe
#      stamp written first.
# The phased_data2021 output found at that path (no .recipe stamp) is moved
# to sgdp_prepped.chrN.phased_data2021.bcf, not deleted.
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
SGDP_VCF_DIR="${SGDP_VCF_DIR:-${SGDP_RAW_DIR}/vcfs}"
SGDP_METADATA_URL="${SGDP_METADATA_URL:-https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt}"
SGDP_METADATA="${SGDP_METADATA:-${SGDP_RAW_DIR}/$(basename "$SGDP_METADATA_URL")}"
SGDP_MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
SGDP_MIN_FL="${SGDP_MIN_FL:-1}"
SGDP_SITE_MASK="${SGDP_SITE_MASK:-1}"
# Shared with 1a, 1c, 1d: `-m flip -d` discards every A/T and C/G SNP even
# when REF matches; flip-all keeps them (these are forward-strand WGS calls).
FIXREF_MODE="${FIXREF_MODE:-flip-all}"
HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg19ToHg38.over.chain.gz}"
# Picard heap derived from the allocation (see 1b_sgdp_liftover_clm.sh: a
# fixed heap above --mem is killed by the kernel, not reported).
if [[ -z "${PICARD_XMX:-}" && -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    _hdr=$(( SLURM_MEM_PER_NODE / 4 )); [[ "$_hdr" -lt 4096 ]] && _hdr=4096
    _heap=$(( (SLURM_MEM_PER_NODE - _hdr) / 1024 )); [[ "$_heap" -lt 4 ]] && _heap=4
    PICARD_XMX="${_heap}g"
fi
PICARD_XMX="${PICARD_XMX:-8g}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMPDIR="${OUTDIR}/tmp/chr${CHR}"
mkdir -p "$OUTDIR" "$TMPDIR"
log () { echo "[$(date +%T)] [chr${CHR}] $*"; }

# The stamp names everything that changes the output, so changing a knob
# rebuilds instead of reusing.
RECIPE="sgdp_raw_v1 minFL=${SGDP_MIN_FL} sitemask=${SGDP_SITE_MASK} fixref=${FIXREF_MODE}"
[[ "$SGDP_MIN_FL" =~ ^[0-9]$ ]] || { echo "ERROR: SGDP_MIN_FL must be 0-9, got '${SGDP_MIN_FL}'" >&2; exit 1; }
[[ -f "${SGDP_VCF_DIR}/DONE" ]] || { echo "ERROR: no ${SGDP_VCF_DIR}/DONE -- run scripts/fetch-sgdp-vcf-variants.sh first" >&2; exit 1; }
[[ -s "$SGDP_MANIFEST" ]] || { echo "ERROR: missing $SGDP_MANIFEST" >&2; exit 1; }
[[ -s "$HG38_FA" ]] || { echo "ERROR: missing $HG38_FA" >&2; exit 1; }
[[ -s "${HG38_FA%.*}.dict" ]] || { echo "ERROR: missing Picard dictionary ${HG38_FA%.*}.dict" >&2; exit 1; }
[[ -s "$CHAIN" ]] || { echo "ERROR: missing hg19->hg38 chain $CHAIN" >&2; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

PREPPED="${OUTDIR}/sgdp_prepped.chr${CHR}.bcf"
STAMP="${PREPPED}.recipe"
log "recipe: ${RECIPE}"

# 0. Reuse only an output built by THIS recipe that still reads back.
if [[ -s "$PREPPED" && -s "$STAMP" && "$(cat "$STAMP")" == "$RECIPE" ]] \
   && bcftools view -h "$PREPPED" > /dev/null && [[ "$(bcftools index -n "$PREPPED")" -gt 0 ]]; then
    log "${PREPPED} was built by this recipe and reads cleanly -- skipping"
    exit 0
fi
if [[ -e "$PREPPED" && ! -e "$STAMP" ]]; then
    OLD="${OUTDIR}/sgdp_prepped.chr${CHR}.phased_data2021.bcf"
    log "${PREPPED} has no recipe stamp: the phased_data2021 (imputed) output -> moving it to ${OLD}"
    mv "$PREPPED" "$OLD"
    if [[ -e "${PREPPED}.csi" ]]; then mv "${PREPPED}.csi" "${OLD}.csi"; fi
fi

# 1. Which file is which panel sample.
if [[ ! -s "$SGDP_METADATA" ]]; then
    log "downloading the SGDP metadata -> ${SGDP_METADATA}"
    curl -fsS --retry 5 --retry-delay 30 -o "${SGDP_METADATA}.partial.$$" "$SGDP_METADATA_URL"
    mv "${SGDP_METADATA}.partial.$$" "$SGDP_METADATA"
fi
MAP="${TMPDIR}/sgdp_raw_map.chr${CHR}.tsv"
python3 -I "${REPO_DIR}/scripts/sgdp-raw-sample-map.py" --manifest "$SGDP_MANIFEST" \
    --metadata "$SGDP_METADATA" --vcf-dir "$SGDP_VCF_DIR" --out "$MAP"
N_SAMPLES=$(wc -l < "$MAP")
cut -f3 "$MAP" > "${TMPDIR}/sgdp_raw_vcfs.chr${CHR}.txt"
awk -F'\t' '{print $1 " " $2}' "$MAP" > "${TMPDIR}/sgdp_raw_rename.chr${CHR}.txt"
log "${N_SAMPLES} per-sample VCFs to merge"

# 2-3. Merge with absent = 0/0, then unreliable calls -> missing.
M0="${TMPDIR}/sgdp_raw.chr${CHR}.hg19.merged.bcf"
M1="${TMPDIR}/sgdp_raw.chr${CHR}.hg19.fl.bcf"
M2="${TMPDIR}/sgdp_raw.chr${CHR}.hg19.bcf"
bcftools merge -0 -m none -r "$CHR" -l "${TMPDIR}/sgdp_raw_vcfs.chr${CHR}.txt" \
    --threads "$THREADS" -Ob -o "$M0"
bcftools index -f "$M0"
LOW_FL='GT="alt"'
for l in $(seq "$SGDP_MIN_FL" 9); do LOW_FL+=" & FMT/FL!=\"${l}\""; done
log "unreliable call = ${LOW_FL}"
bcftools +setGT "$M0" -Ob -o "$M1" -- -t q -n . -i "$LOW_FL"
bcftools index -f "$M1"
N_REC=$(bcftools index -n "$M0")
N_SITES_LOW=$(bcftools query -i 'N_MISSING>0' -f '%POS\n' "$M1" | wc -l)
N_SITES_MASK=$(bcftools query -i 'N_MISSING>0 && N_MISSING>COUNT(GT="alt")' -f '%POS\n' "$M1" | wc -l)
log "${N_REC} records; ${N_SITES_LOW} have at least one call below FL ${SGDP_MIN_FL};" \
    "${N_SITES_MASK} have more failing carriers than passing ones"

# 4-5. Site mask, GT only, panel names.
if [[ "$SGDP_SITE_MASK" == "1" ]]; then
    bcftools +setGT "$M1" -Ou -- -t q -n . -i 'N_MISSING>0 && N_MISSING>COUNT(GT="alt")' \
      | bcftools annotate -x INFO,^FORMAT/GT -Ou \
      | bcftools reheader -s "${TMPDIR}/sgdp_raw_rename.chr${CHR}.txt" -o "$M2"
    log "site mask: all of SGDP set missing at those ${N_SITES_MASK} sites"
else
    bcftools annotate -x INFO,^FORMAT/GT -Ou "$M1" \
      | bcftools reheader -s "${TMPDIR}/sgdp_raw_rename.chr${CHR}.txt" -o "$M2"
    log "SGDP_SITE_MASK=${SGDP_SITE_MASK}: no site mask"
fi
bcftools index -f "$M2"
log "renamed: $(bcftools query -l "$M2" | head -2 | tr '\n' ' ')..."

# 6. Liftover with the 1a recipe, then normalise against hg38.
RENAME="${TMPDIR}/sgdp_raw_to_chr.txt"
: > "$RENAME"; for c in {1..22} X Y MT; do echo "$c chr$c" >> "$RENAME"; done
LIFT_IN="${TMPDIR}/sgdp_raw.chr${CHR}.hg19.chrnamed.vcf.gz"
bcftools annotate --rename-chrs "$RENAME" --threads "$THREADS" -Oz -o "$LIFT_IN" "$M2"
bcftools index -t -f "$LIFT_IN"
LIFTED="${TMPDIR}/sgdp_raw.chr${CHR}.hg38.lifted.vcf.gz"
REJECT="${OUTDIR}/sgdp_prepped.chr${CHR}.rejected.vcf.gz"
log "Picard LiftoverVcf hg19 -> hg38 (-Xmx${PICARD_XMX})"
picard "-Xmx${PICARD_XMX}" LiftoverVcf \
    I="$LIFT_IN" O="$LIFTED" CHAIN="$CHAIN" R="$HG38_FA" REJECT="$REJECT" \
    RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true CREATE_INDEX=false \
    MAX_RECORDS_IN_RAM=200000
NORM="${TMPDIR}/sgdp_raw.chr${CHR}.hg38.norm.bcf"
bcftools norm -f "$HG38_FA" -m- --threads "$THREADS" -Ob -o "$NORM" "$LIFTED"
bcftools index -f "$NORM"
FIXED="${TMPDIR}/sgdp_raw.chr${CHR}.hg38.fixref.bcf"
bcftools +fixref "$NORM" --threads "$THREADS" -Ob -o "$FIXED" -- -f "$HG38_FA" -m "$FIXREF_MODE" -d
PREP_TMP="${TMPDIR}/sgdp_prepped.chr${CHR}.partial.bcf"
rm -f "$PREP_TMP" "${PREP_TMP}.csi"
# -t (streaming), not -r: the fixref output is not indexed (see 1b_sgdp_liftover).
bcftools view -t "chr${CHR}" -Ou "$FIXED" \
  | bcftools sort -m 8G -T "${TMPDIR}/sort_raw" -Ob -o "$PREP_TMP"
bcftools index "$PREP_TMP"
bcftools view -h "$PREP_TMP" > /dev/null || { echo "ERROR: chr${CHR} output unreadable right after writing" >&2; exit 5; }
# Stamp FIRST, then the index, then the BCF: a BCF at the final path always
# has the stamp of the recipe that made it, so a stamp-less file there can
# only be the old imputed output.
echo "$RECIPE" > "${STAMP}.partial"; mv "${STAMP}.partial" "$STAMP"
mv "${PREP_TMP}.csi" "${PREPPED}.csi"
mv "$PREP_TMP" "$PREPPED"

N_OUT=$(bcftools index -n "$PREPPED")
N_REJ=$(bcftools view -H "$REJECT" | wc -l)
REJ_PCT=$(awk -v r="$N_REJ" -v k="$N_OUT" 'BEGIN { t = r + k; printf "%.2f", (t > 0 ? 100 * r / t : 0) }')
MISS=$(bcftools +fill-tags "$PREPPED" -Ou -- -t F_MISSING | bcftools query -f '%INFO/F_MISSING\n' \
       | awk '{s += $1; n++} END {printf "%.4f", (n ? s / n : 0)}')
log "DONE. $(bcftools query -l "$PREPPED" | wc -l) samples, ${N_OUT} records on chr${CHR}," \
    "${N_REJ} liftover rejects (${REJ_PCT}%), mean per-record missingness ${MISS} -> ${PREPPED}"
rm -f "$M0" "$M0.csi" "$M1" "$M1.csi" "$M2" "$M2.csi" "$LIFT_IN" "${LIFT_IN}.tbi" "$LIFTED" \
      "$NORM" "$NORM.csi" "$FIXED"
