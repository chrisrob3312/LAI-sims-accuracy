#!/usr/bin/env bash
# =============================================================================
# failure-report.sh -- which pipeline tasks failed, how, and their own logs.
# READ-ONLY on the pipeline (writes only its report file).
#
# For every stage of the rebuild + publication chain:
#   1. the current run's task states (latest submission of each job name;
#      earlier submissions listed as superseded, not counted);
#   2. every task that did not complete: state, exit code, run time against
#      its limit, MaxRSS against the memory requested, node, .out/.err paths;
#   3. the failures grouped by stage, state and exit code, with their nodes --
#      "all on one node" or "all out of memory" shows at once;
#   4. the FULL .err and .out of the first failure of each kind;
#   5. what is still queued and why (a DependencyNeverSatisfied job never runs).
# Logs are found by job id: every stage names them <prefix>_<jobid>.out/.err.
#
#   bash scripts/failure-report.sh                    # since the rebuild epoch
#   SINCE=2026-11-13 bash scripts/failure-report.sh   # since a date
#   CAP=0 bash scripts/failure-report.sh              # never cut a log on screen
# The whole report, logs uncut, goes to ${LOGDIR}/failure_report_<date>.txt.
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
CAP="${CAP:-150}"
U="${USER:-$(id -un)}"
EPOCH="${PROJECT_ROOT}/rebuild_epoch"
if [[ -z "${SINCE:-}" ]]; then
    if [[ -f "$EPOCH" ]]; then SINCE=$(date -r "$EPOCH" +%Y-%m-%d); else SINCE=$(date -d '-7 days' +%Y-%m-%d); fi
fi
PY="${REPO_DIR}/scripts/failure-report.py"
[[ -f "$PY" ]] || { echo "ERROR: ${PY} not found (run from the repo, or set REPO_DIR)" >&2; exit 2; }
WORK=$(mktemp -d "${TMPDIR:-/tmp}/failure_report.XXXXXX") || exit 2
trap 'rm -rf "$WORK"' EXIT

sacct -n -P -u "$U" -S "$SINCE" \
    -o JobID,JobIDRaw,JobName%60,State,ExitCode,Elapsed,Timelimit,MaxRSS,ReqMem,NodeList,Submit \
    > "${WORK}/sacct.txt" || { echo "ERROR: sacct failed (stderr above)" >&2; exit 1; }
squeue -h -u "$U" -r -o '%i|%j|%T|%r' > "${WORK}/squeue.txt" || { echo "ERROR: squeue failed (stderr above)" >&2; exit 1; }
echo "sacct rows since ${SINCE}: $(wc -l < "${WORK}/sacct.txt"); queued/running tasks: $(wc -l < "${WORK}/squeue.txt")"

python3 -I "$PY" --sacct "${WORK}/sacct.txt" --squeue "${WORK}/squeue.txt" --logdir "$LOGDIR" \
    --since "$SINCE" --cap "$CAP" --report "${LOGDIR}/failure_report_$(date +%Y%m%d_%H%M).txt"
