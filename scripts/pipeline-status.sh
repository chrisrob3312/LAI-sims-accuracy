#!/usr/bin/env bash
# ============================================================================
# pipeline-status.sh   --   where is the publication rerun, right now
#
# One read-only view of the whole chain. It exists because "is this broken or
# just not started yet?" is the question that keeps coming up, and an empty
# grep answers it ambiguously: a check that has not run and a check that is
# not firing look identical from the output alone.
#
# So every checkpoint here reports THREE states, never two:
#   not started    the job has not run, so there is nothing to see
#   ok             the job ran and the check fired as expected
#   PROBLEM        the job ran and the check did not fire
#
# Usage:
#   bash scripts/pipeline-status.sh
#   watch -n 300 'bash scripts/pipeline-status.sh'
# ============================================================================
set -uo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
LOG="${LOG:-${PROJECT_ROOT}/logs}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
U="${USER}"

hr () { printf '%s\n' "------------------------------------------------------------------"; }

# state of a job id: RUNNING/PENDING/COMPLETED/FAILED/absent, plus array progress
jstate () {
    local jid="$1"
    [[ -n "$jid" ]] || { printf 'unset'; return; }
    local q
    q=$(squeue -j "$jid" -h -o '%T' 2>/dev/null | sort -u | tr '\n' '/' | sed 's#/$##')
    if [[ -n "$q" ]]; then
        local run pend
        run=$(squeue -j "$jid" -h -t R -r 2>/dev/null | wc -l)
        pend=$(squeue -j "$jid" -h -t PD -r 2>/dev/null | wc -l)
        printf '%s (%s running, %s pending)' "$q" "$run" "$pend"
        return
    fi
    # Not in the queue: finished, or never submitted.
    local acc
    acc=$(sacct -j "$jid" -n -X -o State 2>/dev/null | head -1 | tr -d ' ')
    printf '%s' "${acc:-absent}"
}

echo "=================================================================="
echo " PUBLICATION RERUN STATUS   $(date '+%Y-%m-%d %H:%M')"
echo "=================================================================="
echo
echo " queue: $(squeue -u "$U" -h -t all -r 2>/dev/null | wc -l) array elements"\
     "($(squeue -u "$U" -h -t R -r 2>/dev/null | wc -l) running)"
echo

# --------------------------------------------------------------------------
# Stage table. Job ids are discovered by NAME, so this keeps working across
# the resubmissions this chain has already been through.
# --------------------------------------------------------------------------
hr
printf ' %-22s %-10s %-26s %s\n' "STAGE" "JOBID" "STATE" "TASKS"
hr
# Job ids come from sacct, NOT just squeue: a finished job leaves the queue,
# and "not in queue" is the same ambiguity this script exists to remove --
# it means completed, failed, or never submitted, and those need different
# responses. The MOST RECENT id per name is the live one, which also skips
# the cancelled batches whose delayed notifications keep arriving.
# sacct AND squeue, merged.
#
# sacct omits PENDING jobs: they have no start time, so a queued job is simply
# absent and the table said "never submitted" for eleven jobs that were sitting
# in the queue. Acting on that would mean resubmitting the whole chain and
# duplicating it -- the exact failure this script is supposed to prevent.
SINCE="${SINCE:-$(date -d '3 days ago' +%Y-%m-%d 2>/dev/null || date +%Y-%m-%d)}"
ALL=$(
  sacct -u "$U" -S "$SINCE" -n -X -P -o JobID,JobName,State 2>/dev/null
  squeue -u "$U" -h -o '%i|%j|%T' 2>/dev/null | sed 's/_\[[^]]*\]//'
)

stage_row () {
    local label="$1" name="$2" exact="${3:-}"
    local jid
    jid=$(printf '%s\n' "$ALL" | awk -F'|' -v n="$name" -v ex="$exact" '
        { id=$1; sub(/_.*/, "", id)
          if (ex == "1") { if ($2 == n) print id }
          else          { if (index($2, n) == 1) print id } }' \
        | sort -un | tail -1)
    if [[ -z "$jid" ]]; then
        printf ' %-22s %-10s %-26s %s\n' "$label" "-" "never submitted" ""
        return
    fi
    local counts total
    counts=$(printf '%s\n' "$ALL" | awk -F'|' -v j="$jid" '
        { id=$1; sub(/_.*/, "", id); if (id == j) { st=$3; sub(/ .*/, "", st); c[st]++ } }
        END { for (k in c) printf "%s:%d ", k, c[k] }')
    total=$(printf '%s\n' "$ALL" | awk -F'|' -v j="$jid" '
        { id=$1; sub(/_.*/, "", id); if (id == j) n++ } END { print n+0 }')
    local flag=""
    [[ "$counts" == *FAILED* || "$counts" == *CANCELLED* || "$counts" == *TIMEOUT* ]] && flag=" <<<"
    printf ' %-22s %-10s %-26s %s%s\n' "$label" "$jid" "${counts% }" "$total" "$flag"
}

stage_row "1c prep (chunked)"   "1c_phase"
stage_row "phase chunks"        "1b_chunk"
stage_row "ligate"              "1b_ligat"
stage_row "KING dedup"          "1e_king"
stage_row "4a ADMIXTURE K=8"    "4a_homog"
stage_row "2b simulate"         "2b_simu"
stage_row "RFMix WGS"           "3b_rfmix"
stage_row "RFMix WGS homog"     "3b_homog"
stage_row "RFMix chip"          "3b_chip" 1
stage_row "RFMix chip homog"    "3b_chip_homog" 1
stage_row "FLARE (3 arms)"      "3b_flare"
stage_row "accuracy"            "finalize"
hr
echo " <<< marks a job with failed or cancelled tasks."
echo " FLARE runs as three separate arms; only the newest is shown here."
echo

# --------------------------------------------------------------------------
# Checkpoints. Each distinguishes "not started" from "ran and failed to fire".
# --------------------------------------------------------------------------
echo " CHECKPOINTS"
hr

# 1. the stale-chunk check.
#
# SCOPED TO THE CURRENT JOB ID. Logs are named 1b_chunk_<task>_<jobid>.out and
# they accumulate: June's run left hundreds behind. Globbing them all made
# this report "483 skipped, 0 rebuilt -- June panel may be in use" for an
# array that had not started, which is precisely the false alarm this script
# exists to prevent. A check that cries wolf about old logs is worse than no
# check, because the correct response to a real one is to stop the pipeline.
CHUNK_JID=$(squeue -u "$U" -h -o '%i %j' 2>/dev/null \
            | awk '$2 ~ /^1b_chunk/ {print $1}' | sed 's/_\[.*//; s/_[0-9]*$//' | head -1)
if [[ -z "$CHUNK_JID" ]]; then
    CHUNK_JID=$(sacct -n -X -o JobID,JobName%20 -S today 2>/dev/null \
                | awk '$2 ~ /^1b_chunk/ {print $1}' | tail -1)
fi
if [[ -z "$CHUNK_JID" ]]; then
    printf ' %-34s %s\n' "stale-chunk rebuild" "no chunk array found for this run"
else
    # Array task logs use %j, the TASK's own job id, NOT the array master's,
    # so a glob on the master id matches nothing -- which read as "no logs
    # yet" while 25 chunks had completed. Scope by modification time against
    # the array's submit time instead.
    SUB=$(sacct -j "$CHUNK_JID" -n -X -o Submit 2>/dev/null | head -1)
    shopt -s nullglob
    if [[ -n "$SUB" ]]; then
        cur_logs=($(find "$LOG" -name '1b_chunk_*.out' -newermt "$SUB" 2>/dev/null))
    else
        cur_logs=($(find "$LOG" -name '1b_chunk_*.out' -mtime -3 2>/dev/null))
    fi
    shopt -u nullglob
    if [[ ${#cur_logs[@]} -eq 0 ]]; then
        printf ' %-34s %s\n' "stale-chunk rebuild" \
          "not started (job ${CHUNK_JID} has written no logs yet)"
    else
        n_stale=$(grep -hc 'STALE' "${cur_logs[@]}" 2>/dev/null | paste -sd+ | bc 2>/dev/null)
        n_skip=$(grep -hc 'already done' "${cur_logs[@]}" 2>/dev/null | paste -sd+ | bc 2>/dev/null)
        # Slurm BUFFERS a running task's stdout: its log exists but is empty
        # until the task ends. So "logs exist, counts are zero" means the tasks
        # are still running, not that the check failed to fire -- and saying
        # PROBLEM there would be the same false alarm as reading June's logs.
        if [[ "${n_stale:-0}" -eq 0 && "${n_skip:-0}" -eq 0 ]]; then
            printf ' %-34s %s\n' "stale-chunk rebuild" \
              "running (${#cur_logs[@]} logs, none flushed yet -- Slurm buffers stdout)"
        elif [[ "${n_skip:-0}" -gt 0 && "${n_stale:-0}" -eq 0 ]]; then
            printf ' %-34s %s\n' "stale-chunk rebuild" \
              "*** PROBLEM: ${n_skip} skipped, 0 rebuilt in job ${CHUNK_JID}"
        else
            printf ' %-34s %s\n' "stale-chunk rebuild" \
              "ok (${n_stale:-0} rebuilt, ${n_skip:-0} skipped, ${#cur_logs[@]} logs, job ${CHUNK_JID})"
        fi
    fi
fi

# 2. panel sample count -- the single number that says which panel this is
newest=$(ls -t "$PANEL_DIR"/merged_chr*.shapeit5_phased.softunion_maf005.bcf 2>/dev/null | head -1)
if [[ -z "$newest" ]]; then
    printf ' %-34s %s\n' "phased panel" "not started (no phased BCF)"
else
    n_s=$(bcftools query -l "$newest" 2>/dev/null | wc -l)
    age=$(date -r "$newest" '+%Y-%m-%d %H:%M' 2>/dev/null)
    # A small panel is only a PROBLEM once this run should have replaced it.
    # Until ligate finishes, the old files are simply still there, which is
    # expected and not worth an alarm.
    lig_done=$(squeue -u "$U" -h -o '%j' 2>/dev/null | grep -c '^1b_ligat' || true)
    if [[ "$n_s" -gt 3448 ]]; then
        printf ' %-34s %s\n' "phased panel" "ok (${n_s} samples, $age)"
    elif [[ "${lig_done:-0}" -gt 0 ]]; then
        printf ' %-34s %s\n' "phased panel" \
          "old build present (${n_s} samples, $age) -- ligate has not run yet, expected"
    else
        printf ' %-34s %s\n' "phased panel" \
          "*** PROBLEM: ${n_s} samples ($age) and ligate is done -- old panel was not replaced"
    fi
fi

# 3. KING
if [[ -s "${HOMOG_DIR}/king/duplicate_pairs.tsv" ]]; then
    nd=$(awk 'NR>1' "${HOMOG_DIR}/king/duplicate_pairs.tsv" | wc -l)
    ndrop=$(awk 'NF && $1!~/^#/' reference_ids/panel_duplicate_drop_ids.txt 2>/dev/null | wc -l)
    printf ' %-34s %s\n' "KING de-duplication" "ok (${nd} pairs, ${ndrop} dropped)"
    printf ' %-34s %s\n' "  vs metadata cross-check" \
      "36 pairs are named by alias; KING should find at least these"
else
    printf ' %-34s %s\n' "KING de-duplication" "not started"
fi

# 4. ADMIXTURE / the file everything downstream reads
if [[ -s "${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv" ]]; then
    tot=$(awk 'NR>1' "${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv" | wc -l)
    inc=$(awk -F'\t' 'NR>1 && $NF==1' "${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv" | wc -l)
    printf ' %-34s %s\n' "K=8 panel table" "ok (${tot} samples, ${inc} homogeneous)"
    printf ' %-34s %s\n' "  -> benchmark split" "READY to run (see below)"
else
    printf ' %-34s %s\n' "K=8 panel table" "not started"
fi

# 5. homog list freshness -- the arm that silently used a K=5 list
if [[ -s reference_ids/amr_homog_rfmix.k8.txt ]]; then
    printf ' %-34s %s\n' "homog AMR list (K=8)" \
      "ok ($(awk 'NF' reference_ids/amr_homog_rfmix.k8.txt | wc -l) ids)"
elif [[ -s reference_ids/amr_homog_rfmix.txt ]]; then
    printf ' %-34s %s\n' "homog AMR list" \
      "still K=5 ($(awk 'NF' reference_ids/amr_homog_rfmix.txt | wc -l) ids) -- 4a has not run"
fi
hr
echo
if [[ -s "${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv" ]]; then
cat <<EOF
 NEXT: stage 3 is done, so the benchmark split can run now (seconds):
   python3 -I scripts/split-homog-donors-reference.py \\
     --panel ${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv \\
     --out reference_ids --suffix a  --mode lai
   python3 -I scripts/split-homog-donors-reference.py \\
     --panel ${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv \\
     --out reference_ids --suffix b2 --mode benchmark
   python3 -I scripts/cohort-composition-report.py \\
     --panel ${HOMOG_DIR}/lai_ref_panel_samples.k8.tsv \\
     --refs reference_ids --suffix b2 --out results/cohort_composition.txt
EOF
fi
