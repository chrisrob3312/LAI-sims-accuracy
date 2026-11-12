#!/usr/bin/env bash
# =============================================================================
# reset-publication-outputs.sh -- set the previous simulation / RFMix / FLARE
# outputs aside before the publication stages run on a rebuilt panel.
#
# WHY. The publication stages reuse what they find:
#   - every 3b RFMix arm skips a task whose .done marker and .Lat3 exist, and
#     back-fills .done when an old .Lat3 matches its .map line count;
#   - the FLARE arm skips any task whose .anc.vcf.gz exists;
#   - 2b regenerates, but leaves behind files it no longer writes.
# None of them looks at the AGE of what it reuses. After the panel is rebuilt,
# results from the July/November panel would be skipped over as "done" and then
# scored by finalize as if they were new. So these trees must start empty.
#
# WHAT. Each existing tree is MOVED (never deleted) to
# <tree>.pre_rebuild_<YYYYmmdd_HHMM>, then ${PROJECT_ROOT}/rebuild_epoch is
# written. Anything written after the epoch belongs to the rebuild, so
# re-submitting the chain after a failure still reuses the NEW outputs.
#
#   bash scripts/reset-publication-outputs.sh            # move + write the epoch
#   DRY_RUN=1 bash scripts/reset-publication-outputs.sh  # print what it would move
#   bash scripts/reset-publication-outputs.sh --check    # READ-ONLY: exit 0 only if the
#       epoch exists and no tree holds a file older than it (what the driver runs)
#
# Refuses to move anything while a publication-stage job is queued or running.
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
CHIP_NAME="${CHIP_NAME:-gsa_jessica}"
EPOCH="${PROJECT_ROOT}/rebuild_epoch"
TREES=(
    "${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}"
    "${PROJECT_ROOT}/03_rfmix/${ADMIX_POP}/gen${GEN}"
    "${PROJECT_ROOT}/03_rfmix/chip_${CHIP_NAME}/${ADMIX_POP}/gen${GEN}"
    "${PROJECT_ROOT}/03_flare/${ADMIX_POP}/gen${GEN}"
)
PUB_JOBS='^(2b_simu|3b_rfmix|3b_homog|3b_chip|3b_chip_homog|3b_flare|finalize_pilot)$'

if [[ "${1:-}" == "--check" ]]; then
    if [[ ! -f "$EPOCH" ]]; then
        echo "[reset-check] no ${EPOCH}: the previous outputs were never set aside."
        echo "[reset-check] run: bash scripts/reset-publication-outputs.sh"
        exit 1
    fi
    echo "[reset-check] epoch: $(cat "$EPOCH")"
    bad=0
    for t in "${TREES[@]}"; do
        if [[ -d "$t" ]]; then
            old=$(find "$t" -type f ! -newer "$EPOCH" | wc -l)
            new=$(find "$t" -type f -newer "$EPOCH" | wc -l)
            echo "[reset-check] ${t}: ${new} files from this rebuild, ${old} from before it"
            [[ "$old" -eq 0 ]] || bad=1
        else
            echo "[reset-check] ${t}: absent (the stage creates it)"
        fi
    done
    if [[ "$bad" -ne 0 ]]; then
        echo "[reset-check] FAIL: files older than the epoch would be reused as if new."
        echo "[reset-check] run: bash scripts/reset-publication-outputs.sh"
        exit 1
    fi
    echo "[reset-check] OK"
    exit 0
fi

U="${USER:-$(id -un)}"
if command -v squeue > /dev/null; then
    RUNNING=$(squeue -h -u "$U" -o '%i %j %T' | awk -v re="$PUB_JOBS" '$2 ~ re')
    if [[ -n "$RUNNING" ]]; then
        echo "[reset] ERROR: publication-stage jobs are queued or running; moving their" >&2
        echo "        directories now would break them. Wait or scancel them first:" >&2
        echo "$RUNNING" >&2
        exit 1
    fi
fi
STAMP=$(date +%Y%m%d_%H%M)
for t in "${TREES[@]}"; do
    if [[ -d "$t" ]]; then
        n=$(find "$t" -type f | wc -l)
        echo "[reset] ${t} (${n} files) -> ${t}.pre_rebuild_${STAMP}"
        if [[ -z "${DRY_RUN:-}" ]]; then mv "$t" "${t}.pre_rebuild_${STAMP}" || exit 1; fi
    else
        echo "[reset] ${t}: absent, nothing to move"
    fi
done
if [[ -n "${DRY_RUN:-}" ]]; then
    echo "[reset] DRY RUN -- nothing moved, no epoch written"
    exit 0
fi
# (cd ...) rather than git -C: the cluster's git (1.8) has no -C.
echo "rebuild_epoch ${STAMP} commit=$(cd "$REPO_DIR" && git rev-parse --short HEAD)" > "$EPOCH"
echo "[reset] wrote ${EPOCH}: $(cat "$EPOCH")"
