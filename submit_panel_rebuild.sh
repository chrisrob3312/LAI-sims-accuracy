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
#                      (phased_data2021, GLIMPSE-imputed against 1KG phase 3).
#                      The rebuild is meant to use the raw per-sample WGS calls
#                      (vcf_variants) instead, so this default is REFUSED
#                      unless ALLOW_IMPUTED_SGDP=1; set SGDP_PREP to the raw
#                      prep once it exists. SGDP_PREP=none skips the stage.
#   HGDP1KG_PATTERN    the HGDP+1KG source; default = the UNFILTERED Koenig et
#                      al. release (original_data_downloads/phased_haplotypes_v2,
#                      4,091 samples, SNVs + indels, no MAF cut) so 1c applies
#                      the 0.005 MAF per super-pop. filter1 is the same release
#                      cut at a global MAF 0.005. The 60 extra samples of the
#                      u235158 archive are Koenig's QC removals and are in
#                      panel_exclude_ids.txt.
#   FILL_REF_COHORTS   cohorts whose file is a complete variant callset, so a
#                      record absent from it is 0/0 for them, checked against
#                      the other cohorts (1c 3a). Default "HGDP1KG MXB ALMARRI";
#                      add SGDP only with the raw per-sample prep (never the
#                      imputed release -- refused). "" turns the fill off.
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
RERUN_MXB="${RERUN_MXB:-0}"
SGDP_PREP="${SGDP_PREP:-1b_sgdp_liftover_clm.sh}"
ALLOW_IMPUTED_SGDP="${ALLOW_IMPUTED_SGDP:-0}"
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chr__CHR__.filtered.SNV_INDEL.phased.shapeit5.bcf}"
# ${VAR-default}, not ${VAR:-default}: FILL_REF_COHORTS="" must mean "no fill",
# not "the default".
FILL_REF_COHORTS="${FILL_REF_COHORTS-HGDP1KG MXB ALMARRI}"
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
for s in scripts/build-missingness-groups.py scripts/admix-site-filter.py scripts/sbatch-when-ready.sh \
         scripts/fill-contradictions.py scripts/softunion-sites.py scripts/site-provenance.py; do
    [[ -f "$s" ]] || die "missing $s"
done
if [[ ! -s reference_ids/almarri_manifest.tsv ]]; then
    echo "[rebuild] WARNING: no reference_ids/almarri_manifest.tsv -- 4a will label Almarri" >&2
    echo "          MEN from the ID pattern, without the country. Build it first:" >&2
    echo "            bash scripts/build-almarri-manifest.sh <almarri_prepped.chrN.bcf> > reference_ids/almarri_manifest.tsv" >&2
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

[[ -n "$DRY_RUN" ]] && echo "[rebuild] DRY RUN -- nothing was submitted." >&2
cat >&2 <<EOF

[rebuild] 4a = ${J_4A}
[rebuild] Chain the publication stages (2b simulation, RFMix/FLARE arms, finalize)
[rebuild] to it -- they read 4a's outputs:
    START_STAGE=4 AFTER_JOB=${J_4A} bash submit_publication_rerun.sh
[rebuild] Watch: bash scripts/pipeline-status.sh
EOF
