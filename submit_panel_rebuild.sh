#!/usr/bin/env bash
# =============================================================================
# submit_panel_rebuild.sh -- rebuild the reference panel from the cohort preps
# through 4a, then (CHAIN_PUBLICATION=1) the publication stages, as ONE afterok
# chain that can be submitted before its first input exists:
#
#   [SGDP download job, if still running]
#        |
#   [1a MXB liftover]  1b_sgdp_raw SGDP prep  1d Almarri   (in parallel)
#        \                 |          /
#         1c  merge + QC + soft-union + ADMIXTURE site lists   (PHASE_MODE=chunked)
#                          |
#         1b_phase_chunks  (main array)  +  heavy chunks (own array, longer limit)
#                          |
#         1b_phase_ligate  -> 1e KING (called sites + metadata duplicates)
#                          |
#         4a K=8  (site-filtered markers, unsupervised, data-first discovery,
#                  supervised, homogeneity)
#                          |
#         publication stages 4-6: 2b simulation -> RFMix/FLARE arms -> finalize
#
# IDEMPOTENT. Every stage reuses its output when it is current and skips:
# the download (DONE stamp), SGDP prep and 1d (recipe stamps), 1a (run only
# when its stamp's fixref mode differs), 1c (soft-union newer than every input
# and stamped with the same settings and code), chunks (newer than the
# soft-union, same sample count), ligate (newer than every chunk), 1e and 4a
# (newer than their inputs). So the same command can be re-submitted after any
# failure: finished work costs seconds, and only what is out of date runs.
# Nothing here deletes data. NOTE: sbatch snapshots each script when it is
# QUEUED, so a stage runs the code as it was at submission.
#
# Usage, from the repo checkout:
#   DRY_RUN=1 bash submit_panel_rebuild.sh     # print every sbatch line
#   bash submit_panel_rebuild.sh               # submit
#
# Knobs (environment):
#   RERUN_MXB=auto|0|1 auto (default): re-lift MXB (1a) only when
#                      mxb_lifted.recipe does not record the current
#                      FIXREF_MODE -- the July/November lift has no stamp (it
#                      used `-m flip -d`), so the first run re-lifts.
#   FIXREF_MODE        flip-all (default) or flip, for EVERY prep (1a, SGDP,
#                      1d, 1c's HGDP+1KG step): `-m flip -d` discards all A/T
#                      and C/G SNPs even when REF matches; flip-all keeps them.
#                      1c refuses to merge cohorts stamped with different modes.
#   SGDP_PREP=<script> default 1b_sgdp_raw_prep_clm.sh (the raw per-sample WGS
#                      calls). If the download is not finished, the prep waits
#                      on the queued sgdp_fetch job (or SGDP_FETCH_JOB=<id>).
#                      1b_sgdp_liftover_clm.sh (phased_data2021, GLIMPSE-imputed)
#                      is REFUSED unless ALLOW_IMPUTED_SGDP=1. none skips SGDP.
#   CHAIN_PUBLICATION  1 (default): after 4a, submit the publication stages 4-6
#                      (submit_publication_rerun.sh START_STAGE=4) against it.
#                      Refused until the previous simulation/RFMix/FLARE outputs
#                      are set aside (scripts/reset-publication-outputs.sh): those
#                      stages reuse whatever they find, whatever its age.
#   HGDP1KG_PATTERN    the HGDP+1KG source; default = the UNFILTERED Koenig et
#                      al. release (original_data_downloads/phased_haplotypes_v2,
#                      4,091 samples, SNVs + indels, no MAF cut) so 1c applies
#                      the 0.005 MAF per super-pop. filter1 is the same release
#                      cut at a global MAF 0.005. The 60 extra samples of the
#                      u235158 archive are Koenig's QC removals and are in
#                      panel_exclude_ids.txt.
#   FILL_REF_COHORTS   cohorts whose file is a complete variant callset, so a
#                      record absent from it is 0/0 for them, checked against
#                      the other cohorts (1c 3a). Default "HGDP1KG MXB SGDP
#                      ALMARRI" with the raw SGDP prep, "HGDP1KG MXB ALMARRI"
#                      otherwise (the imputed release is never filled -- 1c
#                      checks the prep's stamp). "" turns the fill off.
#   FILL_GUARD / FILL_GUARD_MIN_N / FILL_GUARD_MIN_AC / FILL_GUARD_MAX_P
#                      the fill check (default 1 / 20 / 5 / 1e-3).
#   SOFTUNION_GROUPS   whose samples define a super-pop's MAF for the soft-union:
#                      superpop (default; every cohort, pooled by super-pop) or
#                      sample_groups (the July/November rule: HGDP+1KG and MXB
#                      only, so a variant only SGDP or Almarri carry is lost).
#   ADMIX_GROUP_MAXMISS / ADMIX_TOTAL_MAXMISS / ADMIX_MISS_MIN_N
#                      the ADMIXTURE site rule (default 0.15 / 0.07 / 5).
#   ADMIX_INPUT / KING_INPUT   phased (default) or softunion: what 4a / 1e read.
#                      softunion = pre-phasing genotypes, gaps left missing.
#   ADMIX_SITE_FILTER / KING_SITE_FILTER   1 (default) restricts to the site list.
#   HEAVY_CHUNKS       chunk indices given their own array (default 209,302:
#                      the MHC and the chr9 pericentromere, which ran out of
#                      time at the main array's limit in November).
#   THROTTLE / CHUNK_THROTTLE   per-array concurrency (default 15 / 60).
# =============================================================================
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(pwd)}"
cd "$REPO_DIR"
[[ -f 4a_build_homogeneity_panel_kN_clm.sh && -f 1c_phasing-jointcall_clm.sh ]] || {
    echo "ERROR: run from the LAI-sims-accuracy checkout (REPO_DIR=${REPO_DIR})" >&2; exit 2; }

DRY_RUN="${DRY_RUN:-}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
RERUN_MXB="${RERUN_MXB:-auto}"
SGDP_PREP="${SGDP_PREP:-1b_sgdp_raw_prep_clm.sh}"
ALLOW_IMPUTED_SGDP="${ALLOW_IMPUTED_SGDP:-0}"
CHAIN_PUBLICATION="${CHAIN_PUBLICATION:-1}"
# Exported: every prep stage submitted below inherits it (sbatch passes the
# submitting environment), so all cohorts get the same fixref mode.
export FIXREF_MODE="${FIXREF_MODE:-flip-all}"
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chr__CHR__.filtered.SNV_INDEL.phased.shapeit5.bcf}"
# ${VAR-default}, not ${VAR:-default}: FILL_REF_COHORTS="" must mean "no fill",
# not "the default". SGDP is in the default only with the raw per-sample prep.
if [[ "$SGDP_PREP" == "1b_sgdp_raw_prep_clm.sh" ]]; then
    FILL_REF_COHORTS="${FILL_REF_COHORTS-HGDP1KG MXB SGDP ALMARRI}"
else
    FILL_REF_COHORTS="${FILL_REF_COHORTS-HGDP1KG MXB ALMARRI}"
fi
FILL_GUARD="${FILL_GUARD:-1}"
FILL_GUARD_MIN_N="${FILL_GUARD_MIN_N:-20}"
FILL_GUARD_MIN_AC="${FILL_GUARD_MIN_AC:-5}"
FILL_GUARD_MAX_P="${FILL_GUARD_MAX_P:-1e-3}"
SOFTUNION_GROUPS="${SOFTUNION_GROUPS:-superpop}"
ADMIX_GROUP_MAXMISS="${ADMIX_GROUP_MAXMISS:-0.15}"
ADMIX_TOTAL_MAXMISS="${ADMIX_TOTAL_MAXMISS:-0.07}"
ADMIX_MISS_MIN_N="${ADMIX_MISS_MIN_N:-5}"
HEAVY_CHUNKS="${HEAVY_CHUNKS:-209,302}"
THROTTLE="${THROTTLE:-15}"
CHUNK_THROTTLE="${CHUNK_THROTTLE:-60}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"

die () { echo "[rebuild] FATAL: $*" >&2; exit 2; }

# ---- prechecks: everything a stage would discover hours in ----------------
[[ -s "$CHUNKS_TSV" ]] || die "missing ${CHUNKS_TSV} (build with resources/build_b38_chunks.py)"
N_CHUNKS=$(( $(wc -l < "$CHUNKS_TSV") - 1 ))
[[ "$N_CHUNKS" -gt 0 ]] || die "${CHUNKS_TSV} has no chunk rows"
_probe="${HGDP1KG_PATTERN//__CHR__/22}"
[[ -s "$_probe" ]] || die "HGDP1KG_PATTERN does not resolve for chr22: ${_probe}"
[[ -s reference_ids/panel_exclude_ids.txt ]] || die "missing reference_ids/panel_exclude_ids.txt"
[[ "$SGDP_PREP" == "none" || -f "$SGDP_PREP" ]] || die "SGDP_PREP=${SGDP_PREP} not found"
# The imputed SGDP release is what put a whole-cohort gap at 39.6% of chr22's
# soft-union in November; the rebuild is meant to use the raw WGS calls.
if [[ "$SGDP_PREP" == "1b_sgdp_liftover_clm.sh" ]]; then
    [[ "$ALLOW_IMPUTED_SGDP" == "1" ]] || die "SGDP_PREP is the GLIMPSE-imputed phased_data2021 prep." \
        "Point SGDP_PREP at the raw per-sample prep, or set ALLOW_IMPUTED_SGDP=1 to use it anyway."
    [[ " ${FILL_REF_COHORTS} " != *" SGDP "* ]] || die "FILL_REF_COHORTS includes SGDP, but SGDP_PREP" \
        "is the imputed release: a site missing from it is one 1KG phase 3 lacks, not a 0/0."
fi
case "$SOFTUNION_GROUPS" in superpop|sample_groups) ;; *) die "SOFTUNION_GROUPS must be superpop or sample_groups";; esac
case "$FIXREF_MODE" in flip|flip-all) ;; *) die "FIXREF_MODE must be flip or flip-all";; esac
for s in scripts/build-missingness-groups.py scripts/admix-site-filter.py scripts/sbatch-when-ready.sh \
         scripts/fill-contradictions.py scripts/softunion-sites.py scripts/site-provenance.py \
         scripts/sgdp-raw-sample-map.py submit_publication_rerun.sh scripts/reset-publication-outputs.sh; do
    [[ -f "$s" ]] || die "missing $s"
done

# MXB is re-lifted only when its stamp does not record this FIXREF_MODE. The
# lift from before the stamps existed used `-m flip -d` (no A/T or C/G SNPs).
MXB_STAMP="${PANEL_DIR}/mxb_lifted.recipe"
MXB_MODE="flip (no stamp)"
if [[ -s "$MXB_STAMP" ]]; then
    MXB_MODE=$(awk '{for (i = 1; i <= NF; i++) if ($i ~ /^fixref=/) {sub(/^fixref=/, "", $i); print $i; exit}}' "$MXB_STAMP")
fi
case "$RERUN_MXB" in
    auto) if [[ "$MXB_MODE" == "$FIXREF_MODE" ]]; then RERUN_MXB=0; else RERUN_MXB=1; fi ;;
    0)    [[ "$MXB_MODE" == "$FIXREF_MODE" ]] || die "MXB was lifted with fixref '${MXB_MODE}', not" \
              "'${FIXREF_MODE}', and 1c would refuse it. Use RERUN_MXB=auto or 1." ;;
    1)    ;;
    *)    die "RERUN_MXB must be auto, 0 or 1" ;;
esac

# The raw SGDP prep needs the verified download. Not finished yet is fine as
# long as the download job is queued or running: the prep waits on it.
SGDP_DEP=()
if [[ "$SGDP_PREP" == "1b_sgdp_raw_prep_clm.sh" ]]; then
    if [[ -f "${SGDP_RAW_DIR}/vcfs/DONE" ]]; then
        SGDP_STATE="downloaded and verified (${SGDP_RAW_DIR}/vcfs/DONE)"
    else
        FETCH_JID="${SGDP_FETCH_JOB:-}"
        if [[ -z "$FETCH_JID" ]] && command -v squeue > /dev/null; then
            FETCH_JID=$(squeue -h -u "${USER:-$(id -un)}" -n sgdp_fetch -o %i | awk 'NR == 1 {print $1}')
        fi
        [[ -n "$FETCH_JID" ]] || die "the SGDP download is not finished (no ${SGDP_RAW_DIR}/vcfs/DONE) and no" \
            "sgdp_fetch job is queued. Submit it first: sbatch scripts/fetch-sgdp-vcf-variants.sh"
        [[ "$FETCH_JID" =~ ^[0-9]+$ ]] || die "could not read the sgdp_fetch job id (got '${FETCH_JID}');" \
            "pass it as SGDP_FETCH_JOB=<id>"
        SGDP_DEP=(--dependency="afterok:${FETCH_JID}")
        SGDP_STATE="download job ${FETCH_JID} still running -- the SGDP prep waits on it"
    fi
fi
if [[ ! -s reference_ids/almarri_manifest.tsv ]]; then
    echo "[rebuild] WARNING: no reference_ids/almarri_manifest.tsv -- 4a will label Almarri" >&2
    echo "          MEN from the ID pattern, without the country. Build it first:" >&2
    echo "            bash scripts/build-almarri-manifest.sh <almarri_prepped.chrN.bcf> > reference_ids/almarri_manifest.tsv" >&2
fi
# The publication stages reuse any output they find, whatever its age (see
# scripts/reset-publication-outputs.sh). Chained to a rebuilt panel, they must
# start from trees with nothing older than the rebuild epoch, or the previous
# pilot's results are skipped over as done and scored as new.
if [[ "$CHAIN_PUBLICATION" == "1" ]]; then
    echo "[rebuild] publication outputs:" >&2
    PROJECT_ROOT="$PROJECT_ROOT" bash scripts/reset-publication-outputs.sh --check >&2 \
        || die "publication outputs from before the rebuild are still in place (above). Set them aside with" \
               "'bash scripts/reset-publication-outputs.sh' (moves, never deletes), or CHAIN_PUBLICATION=0."
fi

# Main chunk array = every index except the heavy ones.
IFS=',' read -r -a _heavy <<< "$HEAVY_CHUNKS"
MAIN_SPEC=$(python3 -I - "$N_CHUNKS" "$HEAVY_CHUNKS" <<'EOF'
import sys
n = int(sys.argv[1]); heavy = {int(x) for x in sys.argv[2].split(",") if x}
keep = [i for i in range(1, n + 1) if i not in heavy]
runs, start = [], keep[0]
for a, b in zip(keep, keep[1:] + [None]):
    if b != a + 1:
        runs.append(f"{start}-{a}" if a != start else f"{start}"); start = b
print(",".join(runs))
EOF
)

submit () {   # submit <label> <sbatch args...> <script>
    local label="$1"; shift
    echo >&2
    echo ">>> [${label}] sbatch $*" >&2
    if [[ -n "$DRY_RUN" ]]; then echo "DRYRUN_${label}"; return 0; fi
    local jid
    jid=$(scripts/sbatch-when-ready.sh -- "$@") || die "${label}: submission failed"
    [[ "$jid" =~ ^[0-9]+$ ]] || die "${label}: no job id returned (got '${jid}')"
    echo "    ${label} = ${jid}" >&2
    echo "$jid"
}

echo "[rebuild] chunks: ${N_CHUNKS} (main ${MAIN_SPEC}; heavy ${HEAVY_CHUNKS})" >&2
echo "[rebuild] HGDP+1KG source: ${HGDP1KG_PATTERN}" >&2
echo "[rebuild] fixref mode (every prep): ${FIXREF_MODE}; MXB stamp: ${MXB_MODE} -> re-lift MXB: ${RERUN_MXB}" >&2
echo "[rebuild] SGDP prep: ${SGDP_PREP}${SGDP_STATE:+; ${SGDP_STATE}}" >&2
echo "[rebuild] absent-as-0/0 cohorts: ${FILL_REF_COHORTS:-(none)}; fill check ${FILL_GUARD}" \
     "(groups >= ${FILL_GUARD_MIN_N}, >= ${FILL_GUARD_MIN_AC} copies, P < ${FILL_GUARD_MAX_P})" >&2
echo "[rebuild] soft-union MAF groups: ${SOFTUNION_GROUPS}" >&2
echo "[rebuild] ADMIXTURE site rule: groups <= ${ADMIX_GROUP_MAXMISS} (n >= ${ADMIX_MISS_MIN_N}), overall <= ${ADMIX_TOTAL_MAXMISS}" >&2
echo "[rebuild] 4a reads: ${ADMIX_INPUT:-phased} (site filter ${ADMIX_SITE_FILTER:-1}); 1e reads: ${KING_INPUT:-phased} (site filter ${KING_SITE_FILTER:-1})" >&2

PREP_DEPS=()
if [[ "$RERUN_MXB" == "1" ]]; then
    PREP_DEPS+=("$(submit 1a_mxb 1a_prep-mxb-liftover_clm.sh)")
fi
if [[ "$SGDP_PREP" != "none" ]]; then
    PREP_DEPS+=("$(submit sgdp_prep --array=1-22%${THROTTLE} ${SGDP_DEP[@]+"${SGDP_DEP[@]}"} "$SGDP_PREP")")
fi
PREP_DEPS+=("$(submit 1d_almarri --array=1-22%${THROTTLE} 1d_almarri_prep_clm.sh)")
DEP_1C=$(IFS=:; echo "${PREP_DEPS[*]}")

# 1c's settings travel in the environment of THIS submission only (sbatch
# exports the submitting environment by default). Not via --export=A=x,B=y:
# sbatch splits that list on commas, so REQUIRE_COHORTS=SGDP,ALMARRI would
# arrive as REQUIRE_COHORTS=SGDP. The $( ) is a subshell, so nothing leaks into
# the other stages' submissions.
J_1C=$(export PHASE_MODE=chunked REQUIRE_COHORTS="SGDP ALMARRI" \
           HGDP1KG_PATTERN ADMIX_GROUP_MAXMISS ADMIX_TOTAL_MAXMISS ADMIX_MISS_MIN_N \
           FILL_REF_COHORTS FILL_GUARD FILL_GUARD_MIN_N FILL_GUARD_MIN_AC FILL_GUARD_MAX_P \
           SOFTUNION_GROUPS
       submit 1c_merge --array=1-22%${THROTTLE} --dependency=afterok:${DEP_1C} \
           1c_phasing-jointcall_clm.sh)
J_CH=$(submit chunks --array=${MAIN_SPEC}%${CHUNK_THROTTLE} --dependency=afterok:${J_1C} \
    1b_phase_chunks_clm.sh)
J_CHH=$(submit chunks_heavy --array=${HEAVY_CHUNKS} --time=48:00:00 --cpus-per-task=16 --mem=32G \
    --dependency=afterok:${J_1C} 1b_phase_chunks_clm.sh)
J_LG=$(submit ligate --array=1-22%${THROTTLE} --dependency=afterok:${J_CH}:${J_CHH} 1b_phase_ligate_clm.sh)
J_1E=$(submit 1e_king --dependency=afterok:${J_LG} 1e_king_dedup_clm.sh)
J_4A=$(submit 4a_k8 --dependency=afterok:${J_1E} 4a_build_homogeneity_panel_kN_clm.sh)

echo "[rebuild] 4a = ${J_4A}" >&2
if [[ "$CHAIN_PUBLICATION" == "1" ]]; then
    # The publication stages read 4a's outputs; they wait on 4a (afterok).
    echo "[rebuild] chaining the publication stages 4-6 (2b simulation, RFMix/FLARE arms," \
         "finalize) to 4a = ${J_4A}" >&2
    START_STAGE=4 AFTER_JOB="$J_4A" DRY_RUN="$DRY_RUN" bash submit_publication_rerun.sh
else
    cat >&2 <<EOF
[rebuild] Publication stages not chained (CHAIN_PUBLICATION=${CHAIN_PUBLICATION}). To chain them:
    START_STAGE=4 AFTER_JOB=${J_4A} bash submit_publication_rerun.sh
EOF
fi
if [[ -n "$DRY_RUN" ]]; then echo "[rebuild] DRY RUN -- nothing was submitted." >&2; fi
echo "[rebuild] Watch: bash scripts/pipeline-status.sh" >&2
