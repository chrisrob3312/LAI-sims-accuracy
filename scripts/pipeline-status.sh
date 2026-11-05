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
printf ' %-22s %-10s %s\n' "STAGE" "JOBID" "STATE"
hr
for spec in \
    "1c prep (chunked):1c_phase" \
    "phase chunks:1b_chunk" \
    "ligate:1b_ligat" \
    "KING dedup:1e_king" \
    "4a ADMIXTURE K=8:4a_homog" \
    "2b simulate:2b_simu" \
    "RFMix WGS:3b_rfmix" \
    "RFMix WGS homog:3b_homog" \
    "RFMix chip:3b_chip" \
    "RFMix chip homog:3b_chip_" \
    "FLARE:3b_flare" \
    "deferred submit:defer_su" \
    "accuracy:finalize" ; do
    label="${spec%%:*}"; name="${spec##*:}"
    ids=$(squeue -u "$U" -h -o '%i %j' 2>/dev/null | awk -v n="$name" '$2 ~ "^"n {print $1}' \
          | sed 's/_\[.*//; s/_[0-9]*$//' | sort -u | tr '\n' ',' | sed 's/,$//')
    if [[ -z "$ids" ]]; then
        printf ' %-22s %-10s %s\n' "$label" "-" "not in queue"
    else
        for id in ${ids//,/ }; do
            printf ' %-22s %-10s %s\n' "$label" "$id" "$(jstate "$id")"
        done
    fi
done
hr
echo

# --------------------------------------------------------------------------
# Checkpoints. Each distinguishes "not started" from "ran and failed to fire".
# --------------------------------------------------------------------------
echo " CHECKPOINTS"
hr

# 1. the stale-chunk check
n_logs=$(ls "$LOG"/1b_chunk_*.out 2>/dev/null | wc -l)
if [[ "$n_logs" -eq 0 ]]; then
    printf ' %-34s %s\n' "stale-chunk rebuild" "not started (no chunk logs yet)"
else
    n_stale=$(grep -hc 'STALE' "$LOG"/1b_chunk_*.out 2>/dev/null | paste -sd+ | bc 2>/dev/null)
    n_skip=$(grep -hc 'already done, skipping' "$LOG"/1b_chunk_*.out 2>/dev/null | paste -sd+ | bc 2>/dev/null)
    if [[ "${n_skip:-0}" -gt 0 && "${n_stale:-0}" -eq 0 ]]; then
        printf ' %-34s %s\n' "stale-chunk rebuild" \
          "*** PROBLEM: ${n_skip} skipped, 0 rebuilt -- June panel may be in use"
    else
        printf ' %-34s %s\n' "stale-chunk rebuild" \
          "ok (${n_stale:-0} rebuilt, ${n_skip:-0} legitimately skipped, ${n_logs} logs)"
    fi
fi

# 2. panel sample count -- the single number that says which panel this is
newest=$(ls -t "$PANEL_DIR"/merged_chr*.shapeit5_phased.softunion_maf005.bcf 2>/dev/null | head -1)
if [[ -z "$newest" ]]; then
    printf ' %-34s %s\n' "phased panel" "not started (no phased BCF)"
else
    n_s=$(bcftools query -l "$newest" 2>/dev/null | wc -l)
    age=$(date -r "$newest" '+%Y-%m-%d %H:%M' 2>/dev/null)
    if [[ "$n_s" -le 3448 ]]; then
        printf ' %-34s %s\n' "phased panel" \
          "*** CHECK: ${n_s} samples ($age) -- 3448 is the OLD June panel"
    else
        printf ' %-34s %s\n' "phased panel" "ok (${n_s} samples, $age)"
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
