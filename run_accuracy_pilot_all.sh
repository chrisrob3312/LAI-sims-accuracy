#!/usr/bin/env bash
# Run BOTH the WGS pilot and the chip-GSA pilot, using the same CHRS list,
# then regenerate every slide-facing artifact (figures, publication tables,
# committee pptx) so it's ready to `git commit -am ...`.
#
# Default CHRS = the 11 chromosomes we know have complete sweep coverage
# (chr 1-8 + 20-22 = 121 combos scored, roughly 4x the current chr20-22 pilot).
# Override at call time if the sweep has more done.
#
# Usage:
#   bash run_accuracy_pilot_all.sh
#   CHRS=1,2,3,4,5,6,7,8,9,10,11,20,21,22 bash run_accuracy_pilot_all.sh
#
# Env overrides:
#   CHRS, PROJECT_ROOT, ADMIX_POP, GEN, WGS_WD, CHIP_WD
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
WGS_WD="${WGS_WD:-${PROJECT_ROOT}/03_rfmix/${ADMIX_POP}/gen${GEN}}"
CHIP_WD="${CHIP_WD:-${PROJECT_ROOT}/03_rfmix/chip_gsa_jessica/${ADMIX_POP}/gen${GEN}}"
CHRS="${CHRS:-1,2,3,4,5,6,7,8,20,21,22}"

echo "=================================================================="
echo "[accuracy-pilot-all]"
echo "  REPO_DIR    : $REPO_DIR"
echo "  ADMIX_POP   : $ADMIX_POP"
echo "  GEN         : $GEN"
echo "  CHRS        : $CHRS"
echo "  WGS  WORKDIR: $WGS_WD"
echo "  CHIP WORKDIR: $CHIP_WD"
echo "=================================================================="

# ---------- 1. WGS pilot ----------
echo
echo ">>> [1/3] WGS pilot"
CHRS="$CHRS" \
PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
WORKDIR="$WGS_WD" OUTDIR="${PROJECT_ROOT}/04_accuracy" \
bash "$REPO_DIR/run_accuracy_pilot.sh"

# ---------- 2. Chip pilot (writes to /tmp then copies) ----------
echo
echo ">>> [2/3] Chip-GSA pilot"
CHIP_OUT=/tmp/accuracy_chip_gsa_expanded
CHIP_REVIEW="$REPO_DIR/results/accuracy_pilot_chip_gsa"

CHRS="$CHRS" \
PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
WORKDIR="$CHIP_WD" OUTDIR="$CHIP_OUT" \
bash "$REPO_DIR/run_accuracy_pilot.sh"

# copy chip artifacts into the chip review dir (run_accuracy_pilot.sh writes
# into results/accuracy_pilot; we redirect chip outputs to a separate slot)
rm -rf "$CHIP_REVIEW"
mkdir -p "$CHIP_REVIEW/plots"
cp "$CHIP_OUT/"*.tsv                       "$CHIP_REVIEW/" || true
cp "$CHIP_OUT/plots/"*.pdf                 "$CHIP_REVIEW/plots/" || true

# also replace the top-level review README so it labels chip vs wgs
n_chip=$(ls "$CHIP_WD"/*.gen${GEN}_chr*.done 2>/dev/null | wc -l)
n_wgs=$(ls "$WGS_WD"/*.gen${GEN}_chr*.done 2>/dev/null | wc -l)
cat > "$CHIP_REVIEW/README.md" <<EOF
# accuracy_v2.R chip-GSA pilot run

Scored against **${n_chip}** completed combos in \`chip_gsa_jessica\` as of $(date -u +%Y-%m-%dT%H:%M:%SZ).
CHRS = ${CHRS}

## Contents
- accuracy_summary.tsv, accuracy_per_chr.tsv, panel_vs_baseline.tsv, accuracy_long.tsv
- plots/*.pdf (six PDFs: boxplot / weighted bars / per-chr / heatmap / MXB scatter / delta forest)
EOF

# ---------- 3. Regenerate slide artifacts ----------
echo
echo ">>> [3/3] Regenerating figures, tables, and pptx"
python3 "$REPO_DIR/slides/make_figures.py"
python3 "$REPO_DIR/slides/make_metrics_table.py"
python3 "$REPO_DIR/slides/build_committee_deck.py"

echo
echo "=================================================================="
echo "[accuracy-pilot-all] DONE"
echo "  WGS  combos scored: (see above; source = $WGS_WD, ${n_wgs} .done markers)"
echo "  CHIP combos scored: ${n_chip}"
echo
echo "Ready to commit + push:"
echo "  git add results/ slides/"
echo "  git commit -m 'accuracy pilot expanded to CHRS=$CHRS'"
echo "  git push"
echo "=================================================================="
