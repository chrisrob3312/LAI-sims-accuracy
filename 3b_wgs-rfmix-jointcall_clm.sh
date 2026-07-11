#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=3b_rfmix
#SBATCH --partition=mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=04:00:00
#SBATCH --time=72:00:00
#SBATCH --mem=64G
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
# Repo root -- BASH_SOURCE under sbatch points at /var/spool, so we fall back
# through SLURM_SUBMIT_DIR before resolving from BASH_SOURCE.
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"

# Phased panel from 1b_phasing-jointcall_clm.sh (chr-stripped, contigs 1..22)
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.rechr.bcf}"

# Sample-ID lists for the 4 panels (panel keep-files from build-panel-keep-files.sh)
KEEP_DIR="${KEEP_DIR:-panel_keep_files}"
REFS="${REFS:-reference_ids}"
NAT_HGDP_RFMIX="${NAT_HGDP_RFMIX:-${REFS}/amr_rfmix.txt}"
NAT_HGDPMXB_RFMIX="${NAT_HGDPMXB_RFMIX:-${REFS}/amr_hgdpmxb_rfmix.txt}"
NAT_HGDPMXB_FULL_RFMIX="${NAT_HGDPMXB_FULL_RFMIX:-${PROJECT_ROOT}/panel_keep_files/panel5_HGDPMXB_FULL_IBS_YRI.keep}"
# panel 5 keep-file is HGDP-NAT + ALL 50 MXB + IBS + YRI; the AMR portion alone
# is amr_rfmix + mxb_rfmix + mxb_simu (built lazily below if needed).
IBS_RFMIX="${IBS_RFMIX:-${REFS}/eur_rfmix.txt}"
YRI_RFMIX="${YRI_RFMIX:-${REFS}/afr_rfmix.txt}"
PEL_RFMIX="${PEL_RFMIX:-${REFS}/pel_rfmix.txt}"               # supply if running panel 2
PEL_EAS_RFMIX="${PEL_EAS_RFMIX:-${REFS}/pel_eas_rfmix.txt}"   # supply if running panel 3

# Tools
# RFMIX_DIR holds the compiled RFMix v1 binaries (PopPhased/, TrioPhased/).
# ANCESTRY_PIPELINE_PY3_DIR holds our Py3-ported copies of Jessica's Py2
# shapeit2rfmix.py and RunRFMix.py (system python is Py3, originals fail).
# RunRFMix.py must still be invoked from inside RFMIX_DIR because it calls
# ./PopPhased/RFMix_PopPhased via a relative path.
RFMIX_DIR="${RFMIX_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/RFMix_v1.5.4}"
ANCESTRY_PIPELINE_PY3_DIR="${ANCESTRY_PIPELINE_PY3_DIR:-${REPO_DIR}/scripts/ancestry_pipeline_py3}"

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
# Rebuild DEFAULT_PANELS unambiguously so PEL_EAS is not silently dropped when
# only pel_eas_rfmix.txt exists (the previous string-substitution was a no-op
# in that case).
DEFAULT_PANELS="NAT_HGDP"
[[ -s "${REFS}/pel_rfmix.txt"     ]] && DEFAULT_PANELS="${DEFAULT_PANELS} NAT_PEL"
[[ -s "${REFS}/pel_eas_rfmix.txt" ]] && DEFAULT_PANELS="${DEFAULT_PANELS} NAT_PEL_EAS"
DEFAULT_PANELS="${DEFAULT_PANELS} NAT_HGDPMXB NAT_HGDPMXB_FULL"
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
set +u
conda activate "$CONDA_ENV"
set -u

# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ----------------------------------------------------------------------------
# Helper: per-pop .haps/.sample for one panel-component (NAT, IBS, YRI, etc.)
# ----------------------------------------------------------------------------
extract_haps () {
    local label="$1" keep_file="$2"
    local out="${WORKDIR}/${label}_chr${CHR}"
    # Extract .haps only if not already cached (expensive plink2 pass).
    if [[ ! -s "${out}.haps" ]]; then
        plink2 --bcf "$PANEL" \
               --keep <(awk '{print 0, $1}' "$keep_file") \
               --export haps \
               --threads "$THREADS" \
               --out "$out" \
               > "${out}.plink.log" 2>&1
    fi

    # Rebuild .sample in SHAPEIT 7-column format UNCONDITIONALLY. Older cached
    # runs may have left a 3-col plink2-native header that shapeit2rfmix.py
    # rejects with 'Shapeit sample file appears to be incorrect'. Skip only if
    # the file already looks 7-col.
    local first
    first=$(head -1 "${out}.sample" 2>/dev/null)
    if [[ "$first" != "ID_1 ID_2 missing father mother sex plink_pheno" ]]; then
        {
            echo "ID_1 ID_2 missing father mother sex plink_pheno"
            echo "0 0 0 D D D B"
            sed '1,2d' "${out}.sample" | awk '{print $1, $2, $3, 0, 0, 0, -9}'
        } > "${out}.sample.rewrite"
        mv "${out}.sample.rewrite" "${out}.sample"
    fi
}

echo "[$(date +%T)] [chr${CHR}] extracting per-pop reference haps"
extract_haps NAT_HGDP        "$NAT_HGDP_RFMIX"
extract_haps NAT_HGDPMXB     "$NAT_HGDPMXB_RFMIX"
extract_haps IBS_1KG         "$IBS_RFMIX"
extract_haps YRI_1KG         "$YRI_RFMIX"
[[ -s "$PEL_RFMIX"     ]] && extract_haps NAT_PEL         "$PEL_RFMIX"     || true
[[ -s "$PEL_EAS_RFMIX" ]] && extract_haps NAT_PEL_EAS     "$PEL_EAS_RFMIX" || true

# Panel 5 (HGDP-NAT + ALL 50 MXB) -- AMR portion = amr_rfmix + mxb_rfmix + mxb_simulation.
# Build a temporary keep-file inline since the union isn't stored separately in REFS/.
{
    cat "${REFS}/amr_rfmix.txt" "${REFS}/mxb_rfmix.txt" "${REFS}/mxb_simulation.txt"
} > "${WORKDIR}/_NAT_HGDPMXB_FULL_keep.chr${CHR}.txt"
extract_haps NAT_HGDPMXB_FULL "${WORKDIR}/_NAT_HGDPMXB_FULL_keep.chr${CHR}.txt"

# ----------------------------------------------------------------------------
# .ref keep-files for shapeit2rfmix (sample order: NAT then EUR then AFR)
# ----------------------------------------------------------------------------
build_ref () {
    local panel="$1" nat_label="$2"
    # Per-chr suffix to avoid concurrent write-vs-read races between array tasks
    # sharing the same $WORKDIR (per-pop, per-generation, not per-chr).
    local ref="${WORKDIR}/REF_${panel}.chr${CHR}.ref"
    sed '1,2d' "${WORKDIR}/${nat_label}_chr${CHR}.sample" | awk '{print $2}' > "$ref"
    sed '1,2d' "${WORKDIR}/IBS_1KG_chr${CHR}.sample"      | awk '{print $2}' >> "$ref"
    sed '1,2d' "${WORKDIR}/YRI_1KG_chr${CHR}.sample"      | awk '{print $2}' >> "$ref"
}

build_ref NAT_HGDP         NAT_HGDP
build_ref NAT_HGDPMXB      NAT_HGDPMXB
build_ref NAT_HGDPMXB_FULL NAT_HGDPMXB_FULL
[[ -s "${WORKDIR}/NAT_PEL_chr${CHR}.haps"     ]] && build_ref NAT_PEL     NAT_PEL     || true
[[ -s "${WORKDIR}/NAT_PEL_EAS_chr${CHR}.haps" ]] && build_ref NAT_PEL_EAS NAT_PEL_EAS || true

# ----------------------------------------------------------------------------
# For each (track, panel): shapeit2rfmix.py -> RFMix v1 -> Viterbi recoding
# ----------------------------------------------------------------------------
run_panel_track () {
    local track="$1" panel="$2"

    # Panel 5 (NAT_HGDPMXB_FULL) puts ALL 50 MXB in the reference. Pairing it
    # with the NATMXB sim track would put the 25 MXB-simu donors in BOTH the
    # simulated haplotypes AND the reference, biasing TPR upward. Skip that
    # combination.
    if [[ "$panel" == "NAT_HGDPMXB_FULL" && "$track" == "NATMXB" ]]; then
        echo "[chr${CHR}] [$track/$panel] skipping: donor/reference overlap -- panel 5 only valid against NAT track"
        return 0
    fi

    local nat_label
    case "$panel" in
        NAT_HGDP)          nat_label="NAT_HGDP" ;;
        NAT_HGDPMXB)       nat_label="NAT_HGDPMXB" ;;
        NAT_HGDPMXB_FULL)  nat_label="NAT_HGDPMXB_FULL" ;;
        NAT_PEL)           nat_label="NAT_PEL" ;;
        NAT_PEL_EAS)       nat_label="NAT_PEL_EAS" ;;
        *) echo "ERROR: unknown panel $panel"; exit 1 ;;
    esac
    local nat_haps="${WORKDIR}/${nat_label}_chr${CHR}.haps"
    [[ -s "$nat_haps" ]] || { echo "[chr${CHR}] [$track/$panel] missing $nat_haps -- skipping"; return 0; }

    local admixed_haps="${SIM_DIR}/${track}.${ADMIX_POP}.chr${CHR}.haps"
    local admixed_sample_raw="${SIM_DIR}/${track}.${ADMIX_POP}.chr${CHR}.sample"
    [[ -s "$admixed_haps" ]] || { echo "ERROR: missing $admixed_haps -- run 2b_simulation_clm.sh first"; exit 1; }

    # 2b emits a 3-col .sample header (ID_1 ID_2 missing) but shapeit2rfmix.py
    # requires the SHAPEIT 7-col header. Normalize into a work copy per (track,chr).
    local admixed_sample="${WORKDIR}/${track}.${ADMIX_POP}.chr${CHR}.7col.sample"
    if [[ ! -s "$admixed_sample" ]]; then
        {
            echo "ID_1 ID_2 missing father mother sex plink_pheno"
            echo "0 0 0 D D D B"
            sed '1,2d' "$admixed_sample_raw" | awk '{print $1, $2, $3, 0, 0, 0, -9}'
        } > "$admixed_sample"
    fi

    local ref_keep="${WORKDIR}/REF_${panel}.chr${CHR}.ref"
    local out_prefix="${WORKDIR}/${track}.${panel}.gen${GEN}_chr${CHR}"

    echo "[$(date +%T)] [chr${CHR}] [$track/$panel] shapeit2rfmix"
    python -u "${ANCESTRY_PIPELINE_PY3_DIR}/shapeit2rfmix.py" \
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
    ( cd "$RFMIX_DIR" && python -u "${ANCESTRY_PIPELINE_PY3_DIR}/RunRFMix.py" \
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
        "${out_prefix}.rfmix.${RFMIX_E}.Viterbi.txt" > "${out_prefix}.recoded"
    awk '{print $1, $3}' "${out_prefix}_chr${CHR}.map" | sed 's/:.*//g' > "${out_prefix}.mappin"
    paste "${out_prefix}.mappin" "${out_prefix}.recoded" | sed 's/\t/ /g' > "${out_prefix}.Lat3"
    rm "${out_prefix}.recoded" "${out_prefix}.mappin"
}

# -----------------------------------------------------------------------------
# Parallelize the (track x panel) combo grid in fixed-size batches. Wall
# clock per task = ceil(N_combos / BATCH) x slowest-combo-in-batch. On 12
# CPUs, BATCH=3 with COMBO_THREADS=4 saturates all cores.
#
# Fixed-size batches (not a rolling wait -n) because bash on this cluster is
# older than 4.3 and doesn't support 'wait -n'. Batched wait is bash 3-safe.
# Skip the (NATMXB, NAT_HGDPMXB_FULL) combo -- it double-dips donors.
# -----------------------------------------------------------------------------
BATCH_SIZE="${COMBOS_MAX_PARALLEL:-3}"
COMBO_THREADS="${COMBO_THREADS:-4}"
THREADS="$COMBO_THREADS"    # picked up by run_panel_track via --num-threads

combos=()
for track in $TRACKS_TO_RUN; do
    for panel in $PANELS_TO_RUN; do
        if [[ "$track" == "NATMXB" && "$panel" == "NAT_HGDPMXB_FULL" ]]; then
            echo "[$(date +%T)] [chr${CHR}] skip invalid combo $track/$panel (donors overlap)"
            continue
        fi
        combos+=("${track}|${panel}")
    done
done

echo "[$(date +%T)] [chr${CHR}] launching ${#combos[@]} combos in batches of ${BATCH_SIZE}, ${COMBO_THREADS} threads each"

batch_num=0
for ((i=0; i<${#combos[@]}; i+=BATCH_SIZE)); do
    batch_num=$((batch_num + 1))
    echo "[$(date +%T)] [chr${CHR}] --- batch ${batch_num}: ${combos[@]:i:BATCH_SIZE} ---"
    for combo in "${combos[@]:i:BATCH_SIZE}"; do
        track="${combo%|*}"; panel="${combo#*|}"
        combo_log="${WORKDIR}/_combo_${track}_${panel}_chr${CHR}.log"
        (
            set +e
            echo "[$(date +%T)] [chr${CHR}] [$track/$panel] START (log: $(basename "$combo_log"))"
            run_panel_track "$track" "$panel" >> "$combo_log" 2>&1
            rc=$?
            echo "[$(date +%T)] [chr${CHR}] [$track/$panel] END rc=${rc}"
            exit "$rc"
        ) &
    done
    # Wait for all subshells in this batch; don't let set -e kill us on a
    # failed combo (the failed combo's log preserves the traceback).
    wait || true
done

echo "[$(date +%T)] [chr${CHR}] # Complete. RFMix outputs in $WORKDIR"
