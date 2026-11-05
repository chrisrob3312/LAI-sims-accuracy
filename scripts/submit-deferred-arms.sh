#!/usr/bin/env bash
#SBATCH --job-name=defer_submit
#SBATCH --partition=atkinson,mhgcp
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/defer_submit_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/defer_submit_%j.err
# ============================================================================
# submit-deferred-arms.sh
#
# Submit the two extra FLARE arms LATER, from inside a tiny job that waits
# until there is room in the queue.
#
# THE PROBLEM THIS SOLVES. This cluster caps submitted jobs per user, counting
# array ELEMENTS. With the 594-chunk phasing array queued, the main chain uses
# roughly 1,345 of them, and a 198-task arm will not fit: sbatch prints
# "Slurm temporarily unable to accept job, sleeping and retrying", gives up,
# and returns nothing.
#
# The arms are not urgent -- they depend on stage 4, which is days away. They
# only have to be submitted BEFORE stage 4 finishes. So this job is chained to
# LIGATE: by the time ligate completes, the 594 chunks have drained and the
# space exists. It then submits the arms and exits.
#
# Submitting from a job rather than by hand matters because the window opens
# at an unpredictable hour and a missed window has already cost one arm.
#
# Usage -- chain it to ligate, and tell it what stage 4 is:
#   STAGE4_JID=3742133 sbatch --dependency=afterany:3742130 \
#       scripts/submit-deferred-arms.sh
#
# afterany, not afterok: if ligate fails the arms are pointless, but so is
# leaving this job queued forever holding a slot.
# ============================================================================
set -uo pipefail

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
STAGE4_JID="${STAGE4_JID:?set STAGE4_JID to the 2b_simulation job id the arms must wait for}"
THROTTLE="${THROTTLE:-30}"
CHIP_NAME="${CHIP_NAME:-gsa_jessica}"
CHIP_SNPS_TPL="${CHIP_SNPS_TPL:-${REPO_DIR}/reference_ids/chip_${CHIP_NAME}/chr%d.snps}"
FINALIZE_JID="${FINALIZE_JID:-}"   # set to re-chain accuracy onto all arms
MAX_WAIT_MIN="${MAX_WAIT_MIN:-90}"

cd "$REPO_DIR"
echo "[defer] repo=${REPO_DIR}"
echo "[defer] arms will depend on stage 4 = ${STAGE4_JID}"

# try_submit <label> <sbatch args...> -> prints the job id, or empty
try_submit () {
    local label="$1"; shift
    local jid="" waited=0
    while :; do
        jid=$(sbatch --parsable "$@" 2>/dev/null || true)
        jid="${jid%%;*}"
        if [[ "$jid" =~ ^[0-9]+$ ]]; then
            echo "[defer] ${label} submitted as ${jid}"
            printf '%s' "$jid"
            return 0
        fi
        if [[ "$waited" -ge "$MAX_WAIT_MIN" ]]; then
            echo "[defer] ${label}: still no room after ${waited} min -- giving up" >&2
            return 1
        fi
        echo "[defer] ${label}: queue full, retrying in 60s (waited ${waited} min)" >&2
        sleep 60
        waited=$((waited + 1))
    done
}

echo "[defer] queued array elements right now: $(squeue -u "$USER" -h -t all -r 2>/dev/null | wc -l)"

# The 44-task homog arm goes first: it is the smaller of the two, so if only
# one will fit, the one that fits should be attempted first rather than
# blocked behind the larger.
J_HOMOG=$(try_submit "FLARE homog" \
    --array="1-44%${THROTTLE}" --mem=32G \
    --dependency="afterok:${STAGE4_JID}" \
    --export=ALL,PANELS_TO_RUN=NAT_HOMOG \
    ./3b_flare_wgs-jointcall_clm.sh) || J_HOMOG=""

J_CHIP=$(try_submit "FLARE chip" \
    --array="1-198%${THROTTLE}" --mem=32G \
    --dependency="afterok:${STAGE4_JID}" \
    --export="ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}" \
    ./3b_flare_wgs-jointcall_clm.sh) || J_CHIP=""

echo
echo "[defer] FLARE homog = ${J_HOMOG:-FAILED}"
echo "[defer] FLARE chip  = ${J_CHIP:-FAILED}"

# ---------------------------------------------------------------------------
# Re-chain accuracy so one table carries every arm.
#
# The existing finalize job is cancelled only once BOTH new arms exist and the
# replacement has been submitted. Cancelling first and failing to resubmit
# would leave the chain with no accuracy step at all, which is a worse outcome
# than accuracy running on five arms instead of seven.
# ---------------------------------------------------------------------------
if [[ -n "$FINALIZE_JID" && -n "$J_HOMOG" && -n "$J_CHIP" ]]; then
    DEPS="${ARM_JIDS:-}"
    if [[ -z "$DEPS" ]]; then
        echo "[defer] ARM_JIDS not set -- cannot rebuild the accuracy dependency" >&2
        echo "[defer] leaving ${FINALIZE_JID} alone; it will score the arms it knows about" >&2
    else
        NEW=$(try_submit "accuracy (re-chained)" \
            --dependency="afterok:${DEPS}:${J_HOMOG}:${J_CHIP}" \
            ./finalize_pilot.sbatch) || NEW=""
        if [[ -n "$NEW" ]]; then
            scancel "$FINALIZE_JID" && \
              echo "[defer] cancelled old accuracy ${FINALIZE_JID}, replaced by ${NEW}"
        else
            echo "[defer] could not submit the replacement; ${FINALIZE_JID} LEFT IN PLACE" >&2
        fi
    fi
fi

echo "[defer] done."
