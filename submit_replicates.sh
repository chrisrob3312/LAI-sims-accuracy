#!/usr/bin/env bash
# =============================================================================
# submit_replicates.sh
#
# Run the simulation -> RFMix -> accuracy chain N times as INDEPENDENT
# replicates, then combine them for across-replicate error bars.
#
# WHY: every within-run uncertainty estimator (including the block bootstrap in
# accuracy_v2.R) holds the donor pool fixed -- it resamples inside ONE
# realisation. Replicates re-draw both the recombination breakpoints and the
# sampling of donor haplotypes into mosaics, so the spread across replicates is
# the only estimate that captures the full variance. Expect the
# across-replicate interval to be WIDER than the block bootstrap; if it is not,
# the replicates are not independent.
#
# No edits to 2b/3b/accuracy are needed: every path they use is already
# env-overridable, so a replicate is just a different set of directories.
#
# Usage:
#   DRY_RUN=1 bash submit_replicates.sh        # print the chain, submit nothing
#   bash submit_replicates.sh                  # fire 5 replicates
#   N_REPS=10 bash submit_replicates.sh
#
# After all replicates finish:
#   Rscript scripts/combine-replicates.R \
#     $PROJECT_ROOT/04_accuracy/rep{1..5} --out $PROJECT_ROOT/04_accuracy
# =============================================================================
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
N_REPS="${N_REPS:-5}"
THROTTLE="${THROTTLE:-15}"
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
DRY_RUN="${DRY_RUN:-}"

cd "$REPO_DIR"

submit() {
    local label="$1"; shift
    local args="$1"; shift
    local script="$1"
    local line=(sbatch --parsable $args "$script")
    echo >&2
    echo ">>> [${label}]  ${line[*]}" >&2
    if [[ -n "$DRY_RUN" ]]; then
        echo "        (dry run)" >&2; echo "PLACEHOLDER_${label}"
    else
        local jid; jid=$("${line[@]}")
        echo "        job ${jid}" >&2; echo "$jid"
    fi
}

echo "=========================================================="
echo "replicate driver"
echo "  N_REPS       = $N_REPS"
echo "  ADMIX_POP    = $ADMIX_POP   GEN = $GEN"
echo "  PROJECT_ROOT = $PROJECT_ROOT"
echo "  DRY_RUN      = ${DRY_RUN:-off}"
echo "=========================================================="
cat <<'EOF'

  SEED WARNING, read before trusting the error bars:
  admix-simu's simu-mix.pl is stochastic, so separate runs differ. But if it
  accepts an explicit seed flag, pass a DIFFERENT one per replicate and record
  it -- otherwise the replicates are not reproducible and the run cannot be
  regenerated for a reviewer. Check with:
      perl ${ADMIXSIMU_DIR}/simu-mix.pl 2>&1 | head -30
  and if a seed option exists, add it to SIMU_EXTRA below.

EOF

REP_ACC_DIRS=()
for REP in $(seq 1 "$N_REPS"); do
    echo
    echo "============ REPLICATE ${REP}/${N_REPS} ============"

    SIM_OUT="${PROJECT_ROOT}/02_simulations/rep${REP}"
    SIM_DIR="${SIM_OUT}/${ADMIX_POP}/gen${GEN}"
    RFMIX_OUT="${PROJECT_ROOT}/03_rfmix/rep${REP}"
    ACC_OUT="${PROJECT_ROOT}/04_accuracy/rep${REP}"
    REP_ACC_DIRS+=("$ACC_OUT")

    # Per-replicate seed, if simu-mix.pl supports one. Deterministic in REP so
    # the set is reproducible as a whole.
    SIMU_SEED=$(( 20261015 + REP * 7919 ))
    SIMU_EXTRA="${SIMU_EXTRA:-}"   # e.g. "-seed ${SIMU_SEED}" once verified

    J_SIM=$(submit "rep${REP}_simulate" \
        "--array=1-22%${THROTTLE} --export=ALL,OUTDIR=${SIM_OUT},ADMIX_POP=${ADMIX_POP},GEN=${GEN},SIMU_SEED=${SIMU_SEED},SIMU_EXTRA=${SIMU_EXTRA}" \
        ./2b_simulation_clm.sh)

    J_RF=$(submit "rep${REP}_rfmix" \
        "--array=1-198%${THROTTLE} --mem=64G --dependency=afterok:${J_SIM} --export=ALL,SIM_DIR=${SIM_DIR},RFMIX_OUTDIR=${RFMIX_OUT},ADMIX_POP=${ADMIX_POP},GEN=${GEN}" \
        ./3b_wgs-rfmix-jointcall_clm.sh)

    echo "  replicate ${REP}: sim=${J_SIM} rfmix=${J_RF}"
    echo "  accuracy dir: ${ACC_OUT}"
done

echo
echo "=========================================================="
echo "queued ${N_REPS} replicates."
echo
echo "When every replicate's RFMix sweep has finished, score each one:"
for REP in $(seq 1 "$N_REPS"); do
  echo "  WORKDIR=${PROJECT_ROOT}/03_rfmix/rep${REP}/${ADMIX_POP}/gen${GEN} \\"
  echo "    OUTDIR=${PROJECT_ROOT}/04_accuracy/rep${REP} \\"
  echo "    SIM_DIR=${PROJECT_ROOT}/02_simulations/rep${REP}/${ADMIX_POP}/gen${GEN} \\"
  echo "    bash run_accuracy_pilot.sh"
done
echo
echo "Then combine for across-replicate error bars:"
echo "  Rscript scripts/combine-replicates.R \\"
echo "    ${REP_ACC_DIRS[*]} \\"
echo "    --out ${PROJECT_ROOT}/04_accuracy"
echo "=========================================================="
