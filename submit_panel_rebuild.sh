#!/usr/bin/env bash
# =============================================================================
# submit_panel_rebuild.sh -- rebuild the reference panel from the cohort preps
# through 4a, as ONE afterok chain:
#
#   [1a MXB liftover]  [SGDP prep]  1d Almarri          (independent, in parallel)
#        \                 |          /
#         1c  merge + QC + soft-union + ADMIXTURE site lists   (PHASE_MODE=chunked)
#                          |
#         1b_phase_chunks  (main array)  +  heavy chunks (own array, longer limit)
#                          |
#         1b_phase_ligate  -> 1e KING (called sites + metadata duplicates)
#                          |
#         4a K=8  (site-filtered markers, unsupervised, data-first discovery,
#                  supervised, homogeneity)
#
# Every stage now reuses its output only when that output is newer than every
# input (1d also needs its recipe stamp and a file that reads back), so nothing
# from the November panel can be carried into this one: the new soft-union is
# newer than every old chunk, the new chunks are newer than every old ligated
# file, and so on down the chain. Nothing here deletes or moves data.
#
# The stages after 4a (2b simulation, the RFMix/FLARE arms, finalize) are the
# publication orchestrator's stages 4-6; this prints the command that chains
# them to the 4a job it submitted.
#
# Usage, from the repo checkout:
#   DRY_RUN=1 bash submit_panel_rebuild.sh     # print every sbatch line
#   bash submit_panel_rebuild.sh               # submit
#
# Knobs (environment):
#   RERUN_MXB=0|1      1a has not changed and its output reads cleanly, so it
#                      is off by default. 1 re-lifts MXB too (adds hours ahead
#                      of 1c).
#   SGDP_PREP=<script> the SGDP prep to run; default 1b_sgdp_liftover_clm.sh
#                      (phased_data2021 -- it no-ops on chromosomes already
#                      prepped). Set to the per-sample prep once it exists.
#                      SGDP_PREP=none skips the stage.
#   HGDP1KG_PATTERN    the HGDP+1KG source; default = phased_haplotypes_v2_filter1,
#                      the source July used (4,091 samples). The u235158 archive
#                      November used adds 60 samples that are Koenig et al.'s QC
#                      removals; they are now in panel_exclude_ids.txt either way.
#   FILL_REF_COHORTS   cohorts whose file is a complete variant callset, so a
#                      site absent from it is 0/0 for them (see 1c). Default
#                      empty; "MXB ALMARRI" for the current sources.
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
RERUN_MXB="${RERUN_MXB:-0}"
SGDP_PREP="${SGDP_PREP:-1b_sgdp_liftover_clm.sh}"
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/hgdp1kgp_chr__CHR__.shapeit5_phased.filter1_SNP_maf005.vcf.gz}"
FILL_REF_COHORTS="${FILL_REF_COHORTS:-}"
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
for s in scripts/build-missingness-groups.py scripts/admix-site-filter.py scripts/sbatch-when-ready.sh; do
    [[ -f "$s" ]] || die "missing $s"
done
if [[ ! -s reference_ids/almarri_manifest.tsv ]]; then
    echo "[rebuild] WARNING: no reference_ids/almarri_manifest.tsv -- 4a will label Almarri" >&2
    echo "          MEN from the ID pattern, without the country. Build it first:" >&2
    echo "            bash scripts/build-almarri-manifest.sh > reference_ids/almarri_manifest.tsv" >&2
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
echo "[rebuild] SGDP prep: ${SGDP_PREP}; MXB liftover rerun: ${RERUN_MXB}" >&2
echo "[rebuild] absent-as-0/0 cohorts: ${FILL_REF_COHORTS:-(none)}" >&2
echo "[rebuild] ADMIXTURE site rule: groups <= ${ADMIX_GROUP_MAXMISS} (n >= ${ADMIX_MISS_MIN_N}), overall <= ${ADMIX_TOTAL_MAXMISS}" >&2
echo "[rebuild] 4a reads: ${ADMIX_INPUT:-phased} (site filter ${ADMIX_SITE_FILTER:-1}); 1e reads: ${KING_INPUT:-phased} (site filter ${KING_SITE_FILTER:-1})" >&2

PREP_DEPS=()
if [[ "$RERUN_MXB" == "1" ]]; then
    PREP_DEPS+=("$(submit 1a_mxb 1a_prep-mxb-liftover_clm.sh)")
fi
if [[ "$SGDP_PREP" != "none" ]]; then
    PREP_DEPS+=("$(submit sgdp_prep --array=1-22%${THROTTLE} "$SGDP_PREP")")
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
           FILL_REF_COHORTS
       submit 1c_merge --array=1-22%${THROTTLE} --dependency=afterok:${DEP_1C} \
           1c_phasing-jointcall_clm.sh)
J_CH=$(submit chunks --array=${MAIN_SPEC}%${CHUNK_THROTTLE} --dependency=afterok:${J_1C} \
    1b_phase_chunks_clm.sh)
J_CHH=$(submit chunks_heavy --array=${HEAVY_CHUNKS} --time=48:00:00 --cpus-per-task=16 --mem=32G \
    --dependency=afterok:${J_1C} 1b_phase_chunks_clm.sh)
J_LG=$(submit ligate --array=1-22%${THROTTLE} --dependency=afterok:${J_CH}:${J_CHH} 1b_phase_ligate_clm.sh)
J_1E=$(submit 1e_king --dependency=afterok:${J_LG} 1e_king_dedup_clm.sh)
J_4A=$(submit 4a_k8 --dependency=afterok:${J_1E} 4a_build_homogeneity_panel_kN_clm.sh)

[[ -n "$DRY_RUN" ]] && echo "[rebuild] DRY RUN -- nothing was submitted." >&2
cat >&2 <<EOF

[rebuild] 4a = ${J_4A}
[rebuild] Chain the publication stages (2b simulation, RFMix/FLARE arms, finalize)
[rebuild] to it -- they read 4a's outputs:
    START_STAGE=4 AFTER_JOB=${J_4A} bash submit_publication_rerun.sh
[rebuild] Watch: bash scripts/pipeline-status.sh
EOF
