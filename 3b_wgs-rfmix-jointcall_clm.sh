#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=3b_rfmix
#SBATCH --partition=mhgcp,atkinson
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=12
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/3b_rfmix_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/3b_rfmix_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Modernized replacement for wgs-simulation-rfmix-jointcall.sh. Builds RFMix
# v1 inputs and runs RFMix across the 4 reference panels x 2 sim-tracks grid
# from the merged + jointly phased panel produced by 1b_phasing-jointcall_clm.sh.
#
# Reference panels (inner dim of comparison grid):
#   1. NAT_HGDP        HGDP-NAT-rfmix + IBS + YRI                    (legacy)
#   2. NAT_PEL         1KG PEL + IBS + YRI                            (requires pel_rfmix.txt)
#   3. NAT_PEL_EAS     1KG PEL+EAS + IBS + YRI                        (requires pel_eas_rfmix.txt)
#   4. NAT_HGDPMXB     HGDP-NAT-rfmix + MXB-rfmix + IBS + YRI         (NEW)
#
# Sim tracks (outer dim, produced by 2b_simulation_clm.sh):
#   NAT       HGDP-NAT donors only       (legacy)
#   NATMXB    HGDP-NAT + MXB donors      (NEW)
#
# SLURM array 1-22 (per-chrom). One (ADMIX_POP, GEN) per submission; loops
# over panels and sim tracks within each chr task. RFMix v1.5.4 is
# single-threaded per run, so the inner loop is sequential.
#
# Output structure (under $RFMIX_OUTDIR/$ADMIX_POP/gen$GEN/):
#   <track>.<panel>.gen${GEN}_chr${i}.alleles / .classes / .snp_locations / .map
#   <track>.<panel>.gen${GEN}_chr${i}.rfmix.2.Viterbi.txt    (RFMix output)
#   <track>.<panel>.gen${GEN}_chr${i}.Lat3                   (accuracy-prep ready)
# ============================================================================

# ----------------------------------------------------------------------------
# Variable configuration  --  paths and tunables (override via env vars at submit time)
# ----------------------------------------------------------------------------
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Project root -- all outputs anchored under this so nothing collides with
# other lab members' work.
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"

# Phased panel from 1b_phasing-jointcall_clm.sh (chr-stripped, contigs 1..22)
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.rechr.bcf}"

# Sample-ID lists for the 4 panels (panel keep-files from build-panel-keep-files.sh)
KEEP_DIR="${KEEP_DIR:-panel_keep_files}"
REFS="${REFS:-reference_ids}"
NAT_HGDP_RFMIX="${NAT_HGDP_RFMIX:-${REFS}/amr_rfmix.txt}"
NAT_HGDPMXB_RFMIX="${NAT_HGDPMXB_RFMIX:-${REFS}/amr_hgdpmxb_rfmix.txt}"
IBS_RFMIX="${IBS_RFMIX:-${REFS}/eur_rfmix.txt}"
YRI_RFMIX="${YRI_RFMIX:-${REFS}/afr_rfmix.txt}"
PEL_RFMIX="${PEL_RFMIX:-${REFS}/pel_rfmix.txt}"               # supply if running panel 2
PEL_EAS_RFMIX="${PEL_EAS_RFMIX:-${REFS}/pel_eas_rfmix.txt}"   # supply if running panel 3

# Tools
RFMIX_DIR="${RFMIX_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/RFMix_v1.5.4}"
ANCESTRY_PIPELINE_DIR="${ANCESTRY_PIPELINE_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/ancestry_pipeline-master}"

# RFMix-format genetic map (3 cols: pos chr cM)
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/reference/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"
GMAP_TPL="${GMAP_TPL:-${GMAP_DIR}/chr%s.b38.rfmix.gmap.txt}"

# Simulation context (must match what 2b_simulation_clm.sh produced)
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
SIM_DIR="${SIM_DIR:-${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}}"
NOTREF_FILE="${NOTREF_FILE:-${SIM_DIR}/${ADMIX_POP}.notref}"  # IDs of simulated admixed indivs

# Which panels to run (space-separated; comment out 2/3 if you don't have PEL lists yet)
# If pel_rfmix.txt / pel_eas_rfmix.txt exist (built by build-pel-panels.sh),
# auto-include panels 2 and 3 in addition to 1 and 4.
DEFAULT_PANELS="NAT_HGDP NAT_HGDPMXB"
[[ -s "${REFS}/pel_rfmix.txt"     ]] && DEFAULT_PANELS="NAT_HGDP NAT_PEL ${DEFAULT_PANELS#NAT_HGDP }"
[[ -s "${REFS}/pel_eas_rfmix.txt" ]] && DEFAULT_PANELS="${DEFAULT_PANELS/NAT_PEL /NAT_PEL NAT_PEL_EAS }"
PANELS_TO_RUN="${PANELS_TO_RUN:-$DEFAULT_PANELS}"

# Sim tracks to run (space-separated)
TRACKS_TO_RUN="${TRACKS_TO_RUN:-NAT NATMXB}"

# RFMix params (legacy values)
RFMIX_E="${RFMIX_E:-2}"
RFMIX_W="${RFMIX_W:-0.2}"
RFMIX_N="${RFMIX_N:-5}"

# Outputs
RFMIX_OUTDIR="${RFMIX_OUTDIR:-${PROJECT_ROOT}/03_rfmix}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-12}

# shellcheck disable=SC2059
PANEL=$(printf "$PANEL_TPL" "$CHR")
# shellcheck disable=SC2059
GMAP=$(printf "$GMAP_TPL" "$CHR")

[[ -s "$PANEL"        ]] || { echo "ERROR: missing $PANEL"; exit 1; }
[[ -s "$GMAP"         ]] || { echo "ERROR: missing $GMAP"; exit 1; }
[[ -d "$RFMIX_DIR"    ]] || { echo "ERROR: missing $RFMIX_DIR"; exit 1; }
[[ -s "$NOTREF_FILE"  ]] || { echo "ERROR: missing $NOTREF_FILE (IDs of simulated admixed indivs)"; exit 1; }

WORKDIR="${RFMIX_OUTDIR}/${ADMIX_POP}/gen${GEN}"
mkdir -p "$WORKDIR" "$LOGDIR"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# ----------------------------------------------------------------------------
# Helper: per-pop .haps/.sample for one panel-component (NAT, IBS, YRI, etc.)
# ----------------------------------------------------------------------------
extract_haps () {
    local label="$1" keep_file="$2"
    local out="${WORKDIR}/${label}_chr${CHR}"
    if [[ -s "${out}.haps" ]]; then return 0; fi   # cached
    plink2 --bcf "$PANEL" \
           --keep <(awk '{print 0, $1}' "$keep_file") \
           --export haps \
           --threads "$THREADS" \
           --out "$out" \
           > "${out}.plink.log" 2>&1

    # Rebuild .sample header (plink2 reorders rows; legacy convention used 7-col format)
    head -n 2 "${out}.sample" > "${out}.sample.hdr"
    sed '1,2d' "${out}.sample" \
        | awk '{print $1, $2, $3, 0, 0, 0, -9}' \
        > "${out}.sample.body"
    cat "${out}.sample.hdr" "${out}.sample.body" > "${out}.sample"
    rm "${out}.sample.hdr" "${out}.sample.body"
}

echo "[$(date +%T)] [chr${CHR}] extracting per-pop reference haps"
extract_haps NAT_HGDP        "$NAT_HGDP_RFMIX"
extract_haps NAT_HGDPMXB     "$NAT_HGDPMXB_RFMIX"
extract_haps IBS_1KG         "$IBS_RFMIX"
extract_haps YRI_1KG         "$YRI_RFMIX"
[[ -s "$PEL_RFMIX"     ]] && extract_haps NAT_PEL         "$PEL_RFMIX"     || true
[[ -s "$PEL_EAS_RFMIX" ]] && extract_haps NAT_PEL_EAS     "$PEL_EAS_RFMIX" || true

# ----------------------------------------------------------------------------
# .ref keep-files for shapeit2rfmix (sample order: NAT then EUR then AFR)
# ----------------------------------------------------------------------------
build_ref () {
    local panel="$1" nat_label="$2"
    local ref="${WORKDIR}/REF_${panel}.ref"
    sed '1,2d' "${WORKDIR}/${nat_label}_chr${CHR}.sample" | awk '{print $2}' > "$ref"
    sed '1,2d' "${WORKDIR}/IBS_1KG_chr${CHR}.sample"      | awk '{print $2}' >> "$ref"
    sed '1,2d' "${WORKDIR}/YRI_1KG_chr${CHR}.sample"      | awk '{print $2}' >> "$ref"
}

build_ref NAT_HGDP    NAT_HGDP
build_ref NAT_HGDPMXB NAT_HGDPMXB
[[ -s "${WORKDIR}/NAT_PEL_chr${CHR}.haps"     ]] && build_ref NAT_PEL     NAT_PEL     || true
[[ -s "${WORKDIR}/NAT_PEL_EAS_chr${CHR}.haps" ]] && build_ref NAT_PEL_EAS NAT_PEL_EAS || true

# ----------------------------------------------------------------------------
# For each (track, panel): shapeit2rfmix.py -> RFMix v1 -> Viterbi recoding
# ----------------------------------------------------------------------------
run_panel_track () {
    local track="$1" panel="$2"
    local nat_label
    case "$panel" in
        NAT_HGDP)     nat_label="NAT_HGDP" ;;
        NAT_HGDPMXB)  nat_label="NAT_HGDPMXB" ;;
        NAT_PEL)      nat_label="NAT_PEL" ;;
        NAT_PEL_EAS)  nat_label="NAT_PEL_EAS" ;;
        *) echo "ERROR: unknown panel $panel"; exit 1 ;;
    esac
    local nat_haps="${WORKDIR}/${nat_label}_chr${CHR}.haps"
    [[ -s "$nat_haps" ]] || { echo "[chr${CHR}] [$track/$panel] missing $nat_haps -- skipping"; return 0; }

    local admixed_haps="${SIM_DIR}/${track}.${ADMIX_POP}.chr${CHR}.haps"
    local admixed_sample="${SIM_DIR}/${track}.${ADMIX_POP}.chr${CHR}.sample"
    [[ -s "$admixed_haps" ]] || { echo "ERROR: missing $admixed_haps -- run 2b_simulation_clm.sh first"; exit 1; }

    local ref_keep="${WORKDIR}/REF_${panel}.ref"
    local out_prefix="${WORKDIR}/${track}.${panel}.gen${GEN}_chr${CHR}"

    echo "[$(date +%T)] [chr${CHR}] [$track/$panel] shapeit2rfmix"
    python "${ANCESTRY_PIPELINE_DIR}/shapeit2rfmix.py" \
        --shapeit_hap_ref     "${nat_haps},${WORKDIR}/IBS_1KG_chr${CHR}.haps,${WORKDIR}/YRI_1KG_chr${CHR}.haps" \
        --shapeit_hap_admixed "$admixed_haps" \
        --shapeit_sample_ref     "${WORKDIR}/${nat_label}_chr${CHR}.sample,${WORKDIR}/IBS_1KG_chr${CHR}.sample,${WORKDIR}/YRI_1KG_chr${CHR}.sample" \
        --shapeit_sample_admixed "$admixed_sample" \
        --chr "$CHR" \
        --genetic_map "$GMAP" \
        --ref_keep "$ref_keep" \
        --admixed_keep "$NOTREF_FILE" \
        --out "$out_prefix"

    echo "[$(date +%T)] [chr${CHR}] [$track/$panel] RFMix v1"
    ( cd "$RFMIX_DIR" && python RunRFMix.py \
        -e "$RFMIX_E" -w "$RFMIX_W" -n "$RFMIX_N" -G "$GEN" \
        --num-threads "$THREADS" \
        --use-reference-panels-in-EM \
        --forward-backward \
        TrioPhased \
        "${out_prefix}_chr${CHR}.alleles" \
        "${out_prefix}.classes" \
        "${out_prefix}_chr${CHR}.snp_locations" \
        -o "${out_prefix}.rfmix" )

    # Viterbi recoding: 1->0 (NAT), 2->1 (EUR), 3->2 (AFR)  -- matches accuracy.R input
    echo "[$(date +%T)] [chr${CHR}] [$track/$panel] recode Viterbi -> Lat3"
    sed -e 's/1/0/g' -e 's/2/1/g' -e 's/3/2/g' \
        "${out_prefix}.rfmix.2.Viterbi.txt" > "${out_prefix}.recoded"
    awk '{print $1, $3}' "${out_prefix}_chr${CHR}.map" | sed 's/:.*//g' > "${out_prefix}.mappin"
    paste "${out_prefix}.mappin" "${out_prefix}.recoded" | sed 's/\t/ /g' > "${out_prefix}.Lat3"
    rm "${out_prefix}.recoded" "${out_prefix}.mappin"
}

for track in $TRACKS_TO_RUN; do
    for panel in $PANELS_TO_RUN; do
        run_panel_track "$track" "$panel"
    done
done

echo "[$(date +%T)] [chr${CHR}] # Complete. RFMix outputs in $WORKDIR"
