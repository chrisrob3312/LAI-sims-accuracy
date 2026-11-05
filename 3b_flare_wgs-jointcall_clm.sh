#!/usr/bin/env bash
#SBATCH --job-name=3b_flare
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=02:00:00
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-198
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/3b_flare_task%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/3b_flare_task%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu
# ============================================================================
# 3b_flare_wgs-jointcall_clm.sh
#
# FLARE local ancestry on the SAME grid as 3b_wgs-rfmix-jointcall_clm.sh:
# the same (track x panel x chromosome) decomposition, the same simulated
# individuals, the same reference panels, the same genetic map.
#
# WHY THE GRID MUST MATCH EXACTLY. The question is "RFMix versus FLARE", so
# the only thing allowed to differ between the two sweeps is the method. If
# FLARE ran on a different panel set, or on re-simulated individuals, or with
# its own subset of chromosomes, a difference in accuracy would not be
# attributable to the method. The task-ID decomposition below is therefore
# copied from 3b, not reinvented, and both arrays are submitted with the same
# ADMIX_POP, GEN, TRACKS_TO_RUN and PANELS_TO_RUN.
#
# WHAT FLARE NEEDS, AND HOW IT DIFFERS FROM RFMix v1
#   RFMix v1 wants .alleles/.classes/.snp_locations built by shapeit2rfmix.py.
#   FLARE wants VCFs:
#     ref=        phased reference VCF (the panel samples)
#     ref-panel=  two columns, sample <TAB> ancestry  (NAT / EUR / AFR)
#     gt=         phased target VCF (the simulated admixed individuals)
#     map=        PLINK-format genetic map
#   so this script converts rather than reusing 3b's intermediates.
#
# The simulated individuals come out of 2b as SHAPEIT .haps/.sample, which
# bcftools converts to VCF. The reference comes from the merged phased panel
# BCF, subset to the panel's keep-file.
#
# OUTPUT: <track>.<panel>.gen${GEN}_chr${CHR}.anc.vcf.gz
#   FORMAT/AN1 and AN2 per sample per site, plus an ##ANCESTRY= header naming
#   the ancestry codes. That is the layout the accuracy reader expects, and
#   the same layout the recorded truth tracts are written in, so truth and
#   calls are read by one reader.
#
# Usage (mirroring the RFMix submission):
#   ADMIX_POP=Brasa GEN=12 sbatch --array=1-198%15 3b_flare_wgs-jointcall_clm.sh
# ============================================================================
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"

PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.rechr.bcf}"
KEEP_DIR="${KEEP_DIR:-${PROJECT_ROOT}/panel_keep_files}"

# Same reference-ID lists 3b uses, so the panels are identical.
NAT_HGDP_RFMIX="${NAT_HGDP_RFMIX:-${REFS}/amr_rfmix.txt}"
NAT_HGDPMXB_RFMIX="${NAT_HGDPMXB_RFMIX:-${REFS}/amr_hgdpmxb_rfmix.txt}"
NAT_HGDPMXB_FULL_RFMIX="${NAT_HGDPMXB_FULL_RFMIX:-${KEEP_DIR}/panel5_HGDPMXB_FULL_IBS_YRI.keep}"
IBS_RFMIX="${IBS_RFMIX:-${REFS}/eur_rfmix.txt}"
YRI_RFMIX="${YRI_RFMIX:-${REFS}/afr_rfmix.txt}"
PEL_RFMIX="${PEL_RFMIX:-${REFS}/pel_rfmix.txt}"
PEL_EAS_RFMIX="${PEL_EAS_RFMIX:-${REFS}/pel_eas_rfmix.txt}"

# FLARE jar. No conda package, so this is a path to flare.jar.
# FLARE is distributed as a single jar by the Browning lab; there is no conda
# package. Searched in a few plausible places so a working install is picked
# up without an env var, but NOT guessed silently: if none is found the error
# says where to get it.
_flare_default=""
for _p in "${PROJECT_ROOT}/software/flare.jar" \
          "/storage/atkinson/shared_resources/software/flare/flare.jar" \
          "/storage/atkinson/shared_resources/software/flare.jar" \
          "${HOME}/software/flare.jar" \
          "${HOME}/flare.jar"; do
    [[ -s "$_p" ]] && { _flare_default="$_p"; break; }
done
FLARE_JAR="${FLARE_JAR:-$_flare_default}"
FLARE_XMX="${FLARE_XMX:-}"      # derived from the allocation below if unset
FLARE_GEN="${FLARE_GEN:-}"      # defaults to GEN (the simulated generations)
FLARE_MIN_MAC="${FLARE_MIN_MAC:-50}"
FLARE_SEED="${FLARE_SEED:-20261105}"

# PLINK-format genetic map. FLARE wants 4 columns: chrom, ID, cM, bp.
GMAP_DIR="${GMAP_DIR:-${REPO_DIR}/resources/gmap_hg38}"
GMAP_TPL="${GMAP_TPL:-${GMAP_DIR}/chr%s.hg38.gmap.txt}"

# CHIP density. Unset = WGS. Set both to run FLARE at array density, so the
# RFMix-vs-FLARE comparison exists at chip density too and not only on WGS.
#
# The crop is applied to the REFERENCE AND THE TARGET ALIKE. Cropping only the
# target would leave FLARE a WGS reference panel to match array data against,
# which is a different (and easier) problem than the one RFMix is given in 5c,
# so the comparison would not be like for like.
CHIP_NAME="${CHIP_NAME:-}"
CHIP_SNPS_TPL="${CHIP_SNPS_TPL:-}"

ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
SIM_DIR="${SIM_DIR:-${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}}"
FLARE_OUTDIR="${FLARE_OUTDIR:-${PROJECT_ROOT}/03_flare${CHIP_NAME:+_${CHIP_NAME}}}"

TRACKS_TO_RUN="${TRACKS_TO_RUN:-NAT NATMXB}"
DEFAULT_PANELS="NAT_HGDP"
[[ -s "$PEL_RFMIX"     ]] && DEFAULT_PANELS="${DEFAULT_PANELS} NAT_PEL"
[[ -s "$PEL_EAS_RFMIX" ]] && DEFAULT_PANELS="${DEFAULT_PANELS} NAT_PEL_EAS"
DEFAULT_PANELS="${DEFAULT_PANELS} NAT_HGDPMXB NAT_HGDPMXB_FULL"
PANELS_TO_RUN="${PANELS_TO_RUN:-$DEFAULT_PANELS}"

THREADS=${SLURM_CPUS_PER_TASK:-4}

# Heap derived from the allocation, never set independently: a JVM given a
# heap larger than its cgroup grows until the kernel kills it, which is how
# 1b lost a whole array. Quarter of the allocation reserved for non-heap use.
if [[ -z "$FLARE_XMX" && -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    _hdr=$(( SLURM_MEM_PER_NODE / 4 )); [[ "$_hdr" -lt 4096 ]] && _hdr=4096
    _heap=$(( (SLURM_MEM_PER_NODE - _hdr) / 1024 )); [[ "$_heap" -lt 4 ]] && _heap=4
    FLARE_XMX="${_heap}g"
fi
FLARE_XMX="${FLARE_XMX:-8g}"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ---------------------------------------------------------------------------
# Task decomposition -- COPIED from 3b_wgs so the two sweeps line up task for
# task. Changing it here would silently compare different cells.
# ---------------------------------------------------------------------------
combos=()
for track in $TRACKS_TO_RUN; do
    for panel in $PANELS_TO_RUN; do
        if [[ "$track" == "NATMXB" && "$panel" == "NAT_HGDPMXB_FULL" ]]; then continue; fi
        combos+=("${track}|${panel}")
    done
done
N_COMBOS=${#combos[@]}
[[ "$N_COMBOS" -gt 0 ]] || { echo "ERROR: empty combo grid"; exit 1; }

TID=$((SLURM_ARRAY_TASK_ID - 1))
CHR=$((TID / N_COMBOS + 1))
CIDX=$((TID % N_COMBOS))
if [[ "$CHR" -gt 22 ]]; then
    echo "[task ${SLURM_ARRAY_TASK_ID}] chr=${CHR} out of range (N_COMBOS=${N_COMBOS}) -- no-op"
    exit 0
fi
COMBO="${combos[$CIDX]}"
TASK_TRACK="${COMBO%|*}"
TASK_PANEL="${COMBO#*|}"
echo "[task ${SLURM_ARRAY_TASK_ID}] chr=${CHR} track=${TASK_TRACK} panel=${TASK_PANEL}"
echo "[task ${SLURM_ARRAY_TASK_ID}] heap -Xmx${FLARE_XMX}, ${THREADS} threads"

PANEL=$(printf "$PANEL_TPL" "$CHR")
GMAP=$(printf "$GMAP_TPL" "$CHR")
WORKDIR="${FLARE_OUTDIR}/${ADMIX_POP}/gen${GEN}"
TMP="${WORKDIR}/tmp/chr${CHR}.${TASK_TRACK}.${TASK_PANEL}"
mkdir -p "$WORKDIR" "$TMP"

OUT_PREFIX="${WORKDIR}/${TASK_TRACK}.${TASK_PANEL}.gen${GEN}_chr${CHR}"
if [[ -s "${OUT_PREFIX}.anc.vcf.gz" ]]; then
    echo "[$(date +%T)] ${OUT_PREFIX}.anc.vcf.gz exists -- skipping"
    exit 0
fi

if [[ -z "$FLARE_JAR" || ! -s "$FLARE_JAR" ]]; then
    cat >&2 <<EOF
ERROR: flare.jar not found.

  Looked in:
    ${PROJECT_ROOT}/software/flare.jar
    /storage/atkinson/shared_resources/software/flare/flare.jar
    /storage/atkinson/shared_resources/software/flare.jar
    ${HOME}/software/flare.jar
    ${HOME}/flare.jar

  FLARE is one jar from the Browning lab, no conda package:
    mkdir -p ${PROJECT_ROOT}/software
    wget -O ${PROJECT_ROOT}/software/flare.jar https://faculty.washington.edu/browning/flare.jar
    java -jar ${PROJECT_ROOT}/software/flare.jar  # prints usage if it works

  Or point at an existing copy:
    FLARE_JAR=/path/to/flare.jar sbatch ... 3b_flare_wgs-jointcall_clm.sh

  If it is already somewhere else on this cluster:
    find /storage/atkinson -name 'flare*.jar' 2>/dev/null
EOF
    exit 1
fi
for f in "$PANEL" "$GMAP"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done
# FLARE's jar is compiled for Java 11 (class file version 55). The cluster's
# default java is 8 (version 52), which fails with
#     UnsupportedClassVersionError: admix/AdmixMain has been compiled by a
#     more recent version of the Java Runtime
# at class-load time -- after the job has started and the reference VCF has
# been built. So the version is checked up front, and a newer JDK inside the
# conda env is preferred over the system one when present.
JAVA_BIN="${JAVA_BIN:-}"
if [[ -z "$JAVA_BIN" ]]; then
    if [[ -x "${CONDA_PREFIX}/bin/java" ]]; then JAVA_BIN="${CONDA_PREFIX}/bin/java"
    else JAVA_BIN="$(command -v java || true)"; fi
fi
[[ -n "$JAVA_BIN" ]] || { echo "ERROR: no java on PATH (FLARE is a jar)"; exit 1; }
# "1.8.0_x" -> 8 ; "11.0.2" / "17.0.1" -> 11 / 17
JAVA_VER=$("$JAVA_BIN" -version 2>&1 | awk -F'"' '/version/ {print $2; exit}' \
           | awk -F. '{ if ($1=="1") print $2; else print $1 }')
echo "[$(date +%T)] java: ${JAVA_BIN} (major ${JAVA_VER:-?})"
if [[ -z "$JAVA_VER" || "$JAVA_VER" -lt 11 ]]; then
    cat >&2 <<EOF
ERROR: FLARE needs Java 11 or newer; ${JAVA_BIN} is major ${JAVA_VER:-unknown}.

  Install one into the conda env (simplest, no module needed):
    conda install -n ${CONDA_ENV} -c conda-forge openjdk=17
  or point at an existing JDK:
    JAVA_BIN=/path/to/jdk-17/bin/java sbatch ... 3b_flare_wgs-jointcall_clm.sh
  or look for a module:
    module avail java 2>&1 | head
EOF
    exit 1
fi

# ---------------------------------------------------------------------------
# 1. The reference VCF and its ancestry map.
#
# ref-panel is sample <TAB> ancestry. The three ancestries are the three
# components the simulation actually used, so FLARE is asked the same question
# RFMix is: NAT, EUR, AFR.
# ---------------------------------------------------------------------------
case "$TASK_PANEL" in
    NAT_HGDP)           NAT_LIST="$NAT_HGDP_RFMIX" ;;
    NAT_PEL)            NAT_LIST="$PEL_RFMIX" ;;
    NAT_PEL_EAS)        NAT_LIST="$PEL_EAS_RFMIX" ;;
    NAT_HGDPMXB)        NAT_LIST="$NAT_HGDPMXB_RFMIX" ;;
    NAT_HGDPMXB_FULL)   NAT_LIST="$NAT_HGDPMXB_FULL_RFMIX" ;;
    *) echo "ERROR: unknown panel ${TASK_PANEL}"; exit 1 ;;
esac
[[ -s "$NAT_LIST" ]] || { echo "ERROR: missing NAT list $NAT_LIST for ${TASK_PANEL}"; exit 1; }

REF_PANEL="${TMP}/ref_panel.txt"
: > "$REF_PANEL"
awk 'NF {print $1"\tNAT"}' "$NAT_LIST"  >> "$REF_PANEL"
awk 'NF {print $1"\tEUR"}' "$IBS_RFMIX" >> "$REF_PANEL"
awk 'NF {print $1"\tAFR"}' "$YRI_RFMIX" >> "$REF_PANEL"

# A sample in two ancestry lists would be trained as two ancestries at once.
DUPS=$(cut -f1 "$REF_PANEL" | sort | uniq -d | wc -l)
[[ "$DUPS" -eq 0 ]] || {
    echo "ERROR: ${DUPS} sample(s) appear under more than one ancestry in ${TASK_PANEL}" >&2
    cut -f1 "$REF_PANEL" | sort | uniq -d | awk 'NR<=10' >&2
    exit 2; }

REF_IDS="${TMP}/ref_ids.txt"
cut -f1 "$REF_PANEL" > "$REF_IDS"
REF_VCF="${TMP}/ref.chr${CHR}.vcf.gz"
echo "[$(date +%T)] [chr${CHR}] reference: $(wc -l < "$REF_IDS") samples from $(basename "$PANEL")"
bcftools view -S "$REF_IDS" --force-samples --threads "$THREADS" -Oz -o "$REF_VCF" "$PANEL"
bcftools index -t "$REF_VCF"

N_REF_IN_VCF=$(bcftools query -l "$REF_VCF" | wc -l)
if [[ "$N_REF_IN_VCF" -lt "$(wc -l < "$REF_IDS")" ]]; then
    echo "[$(date +%T)] [chr${CHR}] NOTE: $(( $(wc -l < "$REF_IDS") - N_REF_IN_VCF )) reference IDs absent from the panel"
    # Keep ref-panel and the VCF in step; FLARE errors on a named sample it
    # cannot find, and --force-samples above already dropped them silently.
    bcftools query -l "$REF_VCF" | sort > "${TMP}/_present.txt"
    awk 'NR==FNR{k[$1];next} ($1 in k)' "${TMP}/_present.txt" "$REF_PANEL" > "${REF_PANEL}.tmp"
    mv "${REF_PANEL}.tmp" "$REF_PANEL"
fi

# ---------------------------------------------------------------------------
# 2. The target VCF: the simulated admixed individuals.
#    2b writes SHAPEIT .haps/.sample; FLARE wants VCF.
# ---------------------------------------------------------------------------
HAPS="${SIM_DIR}/${TASK_TRACK}.${ADMIX_POP}.chr${CHR}.haps"
SAMP="${SIM_DIR}/${TASK_TRACK}.${ADMIX_POP}.chr${CHR}.sample"
[[ -s "$HAPS" ]] || { echo "ERROR: missing $HAPS (run 2b first)"; exit 1; }
[[ -s "$SAMP" ]] || { echo "ERROR: missing $SAMP"; exit 1; }

GT_VCF="${TMP}/target.chr${CHR}.vcf.gz"
echo "[$(date +%T)] [chr${CHR}] converting simulated haplotypes -> VCF"
bcftools convert --hapsample2vcf "${HAPS},${SAMP}" -Oz -o "$GT_VCF" 2>/dev/null || {
    # bcftools wants the pair as prefix when the names are conventional; fall
    # back to the explicit two-file form name difference.
    bcftools convert --hapsample2vcf "${HAPS}" --samples-file "${SAMP}" \
        -Oz -o "$GT_VCF"; }
bcftools index -t "$GT_VCF"

# --- optional chip crop, applied to BOTH sides --------------------------
if [[ -n "$CHIP_SNPS_TPL" ]]; then
    CHIP_SNPS=$(printf "$CHIP_SNPS_TPL" "$CHR")
    [[ -s "$CHIP_SNPS" ]] || { echo "ERROR: missing chip site list $CHIP_SNPS"; exit 1; }
    # The lists are chr:pos with a bare chromosome; the VCFs are chr-prefixed.
    CHIP_REG="${TMP}/chip_regions.chr${CHR}.txt"
    awk -F: 'NF==2 {print "chr"$1"\t"$2}' "$CHIP_SNPS" | sort -k1,1 -k2,2n > "$CHIP_REG"
    N_SITES=$(wc -l < "$CHIP_REG")
    echo "[$(date +%T)] [chr${CHR}] cropping to ${CHIP_NAME}: ${N_SITES} sites"
    for _v in REF GT; do
        _src="${TMP}/$( [[ $_v == REF ]] && echo ref || echo target ).chr${CHR}.vcf.gz"
        _dst="${_src%.vcf.gz}.chip.vcf.gz"
        # -T, not -R: a streaming targets filter, so this does not depend on
        # the index being usable for random access over ~600k scattered sites.
        bcftools view -T "$CHIP_REG" --threads "$THREADS" -Oz -o "$_dst" "$_src"
        bcftools index -t "$_dst"
        mv "$_dst" "$_src"; mv "${_dst}.tbi" "${_src}.tbi"
    done
    N_REF_SITES=$(bcftools index -n "$REF_VCF" 2>/dev/null || echo 0)
    echo "[$(date +%T)] [chr${CHR}] after crop: ${N_REF_SITES} sites in the reference"
    [[ "$N_REF_SITES" -gt 0 ]] || {
        echo "ERROR: chip crop left no sites. The site list and the panel" >&2
        echo "       probably disagree on contig naming or build." >&2; exit 5; }
fi

N_GT=$(bcftools query -l "$GT_VCF" | wc -l)
echo "[$(date +%T)] [chr${CHR}] target: ${N_GT} simulated individuals"
[[ "$N_GT" -gt 0 ]] || { echo "ERROR: target VCF has no samples"; exit 1; }

# A simulated individual must never also be a reference haplotype: that is the
# circularity the whole benchmark is built to avoid, and it is cheap to assert
# here rather than infer from an implausibly good result.
OVERLAP=$(comm -12 <(bcftools query -l "$GT_VCF" | sort -u) \
                   <(bcftools query -l "$REF_VCF" | sort -u) | wc -l)
[[ "$OVERLAP" -eq 0 ]] || {
    echo "ERROR: ${OVERLAP} simulated individuals are also in the reference panel" >&2
    exit 3; }

# ---------------------------------------------------------------------------
# 3. FLARE.
# ---------------------------------------------------------------------------
GEN_ARG="${FLARE_GEN:-$GEN}"
echo "[$(date +%T)] [chr${CHR}] FLARE: gen=${GEN_ARG} min-mac=${FLARE_MIN_MAC}"
"$JAVA_BIN" "-Xmx${FLARE_XMX}" -jar "$FLARE_JAR" \
    ref="$REF_VCF" \
    ref-panel="$REF_PANEL" \
    gt="$GT_VCF" \
    map="$GMAP" \
    out="$OUT_PREFIX" \
    nthreads="$THREADS" \
    gen="$GEN_ARG" \
    min-mac="$FLARE_MIN_MAC" \
    seed="$FLARE_SEED"

[[ -s "${OUT_PREFIX}.anc.vcf.gz" ]] || {
    echo "ERROR: FLARE produced no ${OUT_PREFIX}.anc.vcf.gz" >&2; exit 4; }
bcftools index -t "${OUT_PREFIX}.anc.vcf.gz" 2>/dev/null || true

# The ancestry codes in the header are what the accuracy reader maps calls
# through. If they are missing or in an unexpected order, every call is
# relabelled silently, so they are recorded here per task.
ANC_HDR=$(bcftools view -h "${OUT_PREFIX}.anc.vcf.gz" | grep '^##ANCESTRY' || true)
echo "[$(date +%T)] [chr${CHR}] ${ANC_HDR:-<no ##ANCESTRY header -- investigate>}"
N_OUT=$(bcftools query -l "${OUT_PREFIX}.anc.vcf.gz" | wc -l)
echo "[$(date +%T)] [chr${CHR}] DONE: ${N_OUT} individuals -> ${OUT_PREFIX}.anc.vcf.gz"

rm -rf "$TMP"
