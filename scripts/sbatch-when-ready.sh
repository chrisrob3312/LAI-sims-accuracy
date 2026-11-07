#!/usr/bin/env bash
# ============================================================================
# sbatch-when-ready.sh   --   submit when the queue has room, not before
#
# Wraps sbatch. Works out how many array ELEMENTS the job needs, waits until
# the per-user submit cap has that much headroom, then submits -- retrying
# through transient controller backpressure on the way.
#
# TWO DIFFERENT FAILURES, WHICH LOOK THE SAME AND ARE NOT.
#
#   The cap. MaxSubmitJobsPU counts array elements, not jobs, so a 198-task
#   arm needs 198 slots. Over the cap, sbatch refuses and will keep refusing
#   until something finishes. Waiting is the only fix.
#
#   Backpressure. "Slurm temporarily unable to accept job, sleeping and
#   retrying" is slurmctld being busy -- typically while ingesting a large
#   array someone just submitted. It clears on its own in seconds to minutes
#   and has nothing to do with the cap.
#
# Both print much the same thing and both end with sbatch returning nothing.
# Reading one as the other sends you to the wrong fix: this script reports
# which it is each time it waits, so the distinction is visible.
#
# WHY THIS EXISTS. Three times now a trailing arm has failed to submit and the
# chain carried on without it -- once losing a FLARE arm entirely, noticed
# only because the accuracy job's dependency list was one id short. The window
# opens at an unpredictable hour; a human watching for it is the unreliable
# part.
#
# Usage -- everything after -- is passed to sbatch verbatim:
#   scripts/sbatch-when-ready.sh -- --array=1-198%30 --mem=32G \
#       --dependency=afterok:3742133 3b_flare_wgs-jointcall_clm.sh
#
#   # capture the id, as with sbatch --parsable
#   JID=$(scripts/sbatch-when-ready.sh -- --array=1-44 job.sh)
#
#   # run it unattended for a long wait
#   nohup scripts/sbatch-when-ready.sh --max-wait 720 -- ... &> defer.log &
#
# The job id is the ONLY thing on stdout; progress goes to stderr, so command
# substitution captures the id cleanly.
# ============================================================================
set -uo pipefail

CAP="${CAP:-}"            # default: read from sacctmgr
MARGIN="${MARGIN:-20}"    # slots left spare, so this never fills the cap exactly
MAX_WAIT_MIN="${MAX_WAIT_MIN:-240}"
POLL_SEC="${POLL_SEC:-60}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cap)      CAP="$2"; shift 2 ;;
        --margin)   MARGIN="$2"; shift 2 ;;
        --max-wait) MAX_WAIT_MIN="$2"; shift 2 ;;
        --poll)     POLL_SEC="$2"; shift 2 ;;
        --)         shift; break ;;
        -h|--help)  sed -n '2,50p' "$0"; exit 0 ;;
        *)          echo "unknown option $1 (did you forget --?)" >&2; exit 2 ;;
    esac
done
[[ $# -gt 0 ]] || { echo "ERROR: nothing to submit; put sbatch args after --" >&2; exit 2; }

say () { echo "[when-ready] $*" >&2; }

# --- how many array elements does this submission need? --------------------
# --array=1-198%30 -> 198 ; --array=1,3,5 -> 3 ; --array=5 -> 1 ; absent -> 1
# The %N throttle limits CONCURRENCY, not how many elements are submitted, so
# it must not be counted: a 198%30 array occupies 198 slots, not 30.
array_elements () {
    local spec=""
    for a in "$@"; do
        case "$a" in
            --array=*) spec="${a#--array=}" ;;
            -a)        spec="__NEXT__" ;;
            *)         [[ "$spec" == "__NEXT__" ]] && spec="$a" ;;
        esac
    done
    [[ -n "$spec" && "$spec" != "__NEXT__" ]] || { echo 1; return; }
    spec="${spec%%\%*}"                      # drop the %throttle
    local total=0 part lo hi step
    IFS=',' read -ra parts <<< "$spec"
    for part in "${parts[@]}"; do
        if [[ "$part" == *-* ]]; then
            step="${part##*:}"; [[ "$step" == "$part" ]] && step=1
            part="${part%%:*}"
            lo="${part%%-*}"; hi="${part##*-}"
            [[ "$lo" =~ ^[0-9]+$ && "$hi" =~ ^[0-9]+$ ]] || continue
            total=$(( total + (hi - lo) / step + 1 ))
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            total=$(( total + 1 ))
        fi
    done
    [[ "$total" -gt 0 ]] && echo "$total" || echo 1
}

NEED=$(array_elements "$@")

if [[ -z "$CAP" ]]; then
    CAP=$(sacctmgr -n show qos format=MaxSubmitJobsPU 2>/dev/null \
          | awk 'NF && $1 ~ /^[0-9]+$/ {print $1; exit}')
    CAP="${CAP:-2000}"
fi
say "needs ${NEED} array element(s); cap ${CAP}, margin ${MARGIN}"

# $USER is not guaranteed to be exported -- under `set -u` an unset one aborts
# the arithmetic below, and worse, a squeue that errors for ANY reason returns
# nothing, `wc -l` reports 0, and the script concludes the queue is empty and
# submits blind. "I could not measure the queue" must not read as "the queue
# is empty".
WHOAMI="${USER:-$(id -un 2>/dev/null || echo "")}"
[[ -n "$WHOAMI" ]] || { say "cannot determine the current user; refusing to guess at queue depth"; exit 2; }

waited=0
while :; do
    if ! used=$(squeue -u "$WHOAMI" -h -t all -r 2>/dev/null | wc -l) || [[ -z "$used" ]]; then
        say "squeue failed -- treating the queue as FULL rather than empty, and waiting"
        used="$CAP"
    fi
    free=$(( CAP - used ))
    if [[ "$free" -ge $(( NEED + MARGIN )) ]]; then
        jid=$(sbatch --parsable "$@" 2>&1); rc=$?
        # --parsable gives "id" or "id;cluster"; anything else is an error
        # message that happens to have been printed on stdout.
        id="${jid%%;*}"; id="${id%%$'\n'*}"
        if [[ $rc -eq 0 && "$id" =~ ^[0-9]+$ ]]; then
            say "submitted as ${id} (queue was ${used}/${CAP})"
            printf '%s\n' "$id"
            exit 0
        fi
        # Room existed, so this is backpressure, not the cap. Say so: the two
        # need different responses and only one of them is "wait for jobs to
        # finish".
        say "room existed (${used}/${CAP}) but sbatch refused -- controller backpressure, not the cap"
        say "  ${jid//$'\n'/ | }"
    else
        say "queue ${used}/${CAP}, ${free} free, need $(( NEED + MARGIN )) -- waiting"
    fi
    if [[ "$waited" -ge "$MAX_WAIT_MIN" ]]; then
        say "gave up after ${waited} min. NOTHING WAS SUBMITTED."
        exit 1
    fi
    sleep "$POLL_SEC"
    waited=$(( waited + (POLL_SEC + 59) / 60 ))
done
