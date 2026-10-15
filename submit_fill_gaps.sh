#!/usr/bin/env bash
# Fire off gap-fill sbatch arrays to get every panel × chr combo done for
# both WGS and chip pilots. Priority order matches what unblocks the
# current committee-deck asterisks / blank chip cells fastest:
#
#   1) WGS HOMOG panel   over chr 1-19       (removes chr20-22-only footnote)
#   2) WGS PEL_EAS panel over chr 1-19       (same)
#   3) Chip main array   over chr 1-19       (fills chip density for HGDP/
#                                             HGDPMXB/HGDPMXB_FULL/PEL/PEL_EAS)
#   4) Chip HOMOG fork   over chr 1-19       (HOMOG on chip too)
#
# All submissions are stackable — later stages wait on earlier stages with
# --dependency=afterany so the SLURM throttle stays sane. Existing .done
# markers are respected (each 3b_*.sh checks for its own .done before running,
# so re-submitting is idempotent).
#
# Usage:
#   bash submit_fill_gaps.sh
#   DRY_RUN=1 bash submit_fill_gaps.sh    # print the sbatch lines only
#
# Env overrides:
#   MAX_CHR    (default 19)   — how far up to fill (chr20-22 already done)
#   THROTTLE   (default 12)   — %throttle for each sub-array
set -euo pipefail

MAX_CHR="${MAX_CHR:-19}"
THROTTLE="${THROTTLE:-15}"
DRY_RUN="${DRY_RUN:-}"

# Chip identity: the chip fork scripts default CHIP_NAME=wgs, which puts
# outputs in chip_wgs/ AND runs at full-WGS density (no chip SNP filter).
# We must export both CHIP_NAME and CHIP_SNPS_TPL so outputs land in
# chip_gsa_jessica/ and the admixed haps get chip-filtered.
CHIP_NAME_DEFAULT="gsa_jessica"
CHIP_SNPS_TPL_DEFAULT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/reference_ids/chip_gsa_jessica/chr%d.snps"
CHIP_NAME="${CHIP_NAME:-$CHIP_NAME_DEFAULT}"
CHIP_SNPS_TPL="${CHIP_SNPS_TPL:-$CHIP_SNPS_TPL_DEFAULT}"

# Combos per panel-set (must match the combo grid in each script):
#   HOMOG fork:       NAT + NATMXB × HOMOG only          => 2 combos
#   PEL_EAS-only run: NAT + NATMXB × PEL_EAS only        => 2 combos
#   chip main:        5 panels × NAT + 4 × NATMXB        => 9 combos
#   chip HOMOG fork:  NAT + NATMXB × HOMOG               => 2 combos
n_task_2() { echo "$((MAX_CHR * 2))"; }
n_task_9() { echo "$((MAX_CHR * 9))"; }

submit() {
    # returns the job ID on stdout; logs to stderr so $(submit …) captures the ID cleanly
    local label="$1"; shift
    local n="$1"; shift
    local script="$1"; shift
    local extra="$*"
    local line=(
        sbatch
        --parsable
        --array="1-${n}%${THROTTLE}"
        $extra
        "$script"
    )
    echo                         >&2
    echo ">>> [$label]  ${line[*]}" >&2
    if [[ -n "$DRY_RUN" ]]; then
        echo "        (dry run — not submitted)" >&2
        echo "PLACEHOLDER_$label"
    else
        local jid
        jid=$("${line[@]}")
        echo "        submitted as job ${jid}" >&2
        echo "$jid"
    fi
}

echo "=========================================================="
echo "gap-fill submitter"
echo "  MAX_CHR       = $MAX_CHR   (chr20-22 already have full coverage)"
echo "  THROTTLE      = %${THROTTLE}"
echo "  CHIP_NAME     = $CHIP_NAME  (must NOT default to 'wgs')"
echo "  CHIP_SNPS_TPL = $CHIP_SNPS_TPL"
echo "  DRY_RUN       = ${DRY_RUN:-off}"
echo "=========================================================="

# ---- 1) WGS HOMOG over chr 1-MAX_CHR -----------------------------------------
J1=$(submit "wgs-homog" "$(n_task_2)" ./3b_homog_wgs-rfmix-jointcall_clm.sh)

# ---- 2) WGS PEL_EAS over chr 1-MAX_CHR ---------------------------------------
# Reuse the main WGS script but constrain to just the PEL_EAS panel.
# NAT_PEL_EAS only; TRACKS_TO_RUN defaults to "NAT NATMXB" inside the script,
# so we don't export it here (avoids sbatch's space-in-value pitfall).
J2=$(submit "wgs-pel-eas" "$(n_task_2)" ./3b_wgs-rfmix-jointcall_clm.sh \
     --dependency=afterany:${J1} \
     --export=ALL,PANELS_TO_RUN=NAT_PEL_EAS)

# ---- 3) chip main over chr 1-MAX_CHR -----------------------------------------
# Must export CHIP_NAME + CHIP_SNPS_TPL, otherwise defaults land in chip_wgs/
# and run at full-WGS density (not chip-filtered).
J3=$(submit "chip-main" "$(n_task_9)" ./3b_chip_wgs-rfmix-jointcall_clm.sh \
     --dependency=afterany:${J2} \
     --export="ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}")

# ---- 4) chip HOMOG over chr 1-MAX_CHR ----------------------------------------
J4=$(submit "chip-homog" "$(n_task_2)" ./3b_chip_homog_wgs-rfmix-jointcall_clm.sh \
     --dependency=afterany:${J3} \
     --export="ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}")

echo
echo "=========================================================="
echo "queued (with afterany chain):"
echo "  1) $J1  wgs-homog"
echo "  2) $J2  wgs-pel-eas  (after $J1)"
echo "  3) $J3  chip-main    (after $J2)"
echo "  4) $J4  chip-homog   (after $J3)"
echo
echo "Watch progress with:"
echo "  squeue -u \$USER --states=PD,R"
echo
echo "When each finishes, refresh the deck:"
echo "  bash run_accuracy_pilot_all.sh"
echo "  git add results/ slides/"
echo "  git commit -m 'accuracy pilot: fill chr1-${MAX_CHR} for HOMOG/PEL_EAS + chip'"
echo "  git push"
echo "=========================================================="
