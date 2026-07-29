#!/usr/bin/env bash
# Test-run harness for accuracy_v2.R against whatever 3b combos have finished
# so far. Writes tables + PDF plots into results/accuracy_pilot/ inside the
# repo so they can be committed and reviewed on GitHub.
#
# Behaviour:
#   - Reads only combos with a .done marker + non-empty .Lat3 (partial runs OK)
#   - Copies the small tables/PDFs into results/accuracy_pilot/
#   - Leaves accuracy_long.tsv (per-hap, ~200 MB at full scale) OUT of the
#     repo copy; keep it under 04_accuracy for downstream local use.
#
# Usage:
#   bash run_accuracy_pilot.sh                          # defaults below
#   ADMIX_POP=Brasa GEN=12 bash run_accuracy_pilot.sh
#
# Env overrides:
#   PROJECT_ROOT, ADMIX_POP, GEN, RFMIX_ROOT, WORKDIR, OUTDIR, PANELS_TO_RUN
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
WORKDIR="${WORKDIR:-${PROJECT_ROOT}/03_rfmix/${ADMIX_POP}/gen${GEN}}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/04_accuracy}"
REVIEW_DIR="${REPO_DIR}/results/accuracy_pilot"

mkdir -p "$OUTDIR" "$REVIEW_DIR/plots"

# Quick sanity: how many combos will be scored?
n_done=$(ls "$WORKDIR"/*.gen${GEN}_chr*.done 2>/dev/null | wc -l)
echo "[accuracy-pilot] scoring against ${n_done} completed combos in $WORKDIR"

# Run the R analysis
env \
    PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
    WORKDIR="$WORKDIR" OUTDIR="$OUTDIR" \
    Rscript "$REPO_DIR/accuracy_v2.R"

# Stage the review artifacts (small tables + PDFs)
cp "$OUTDIR/accuracy_summary.tsv"    "$REVIEW_DIR/" || true
cp "$OUTDIR/accuracy_per_chr.tsv"    "$REVIEW_DIR/" || true
cp "$OUTDIR/panel_vs_baseline.tsv"   "$REVIEW_DIR/" || true
cp "$OUTDIR/plots/"*.pdf             "$REVIEW_DIR/plots/" || true

# Coverage summary for the review dir README
{
    echo "# accuracy_v2.R pilot run"
    echo
    echo "Ran against **${n_done}/198** completed combos as of $(date -u +%Y-%m-%dT%H:%M:%SZ)."
    echo
    echo "## Contents"
    echo "- accuracy_summary.tsv    -- weighted concordance per (track, panel, ancestry)"
    echo "- accuracy_per_chr.tsv    -- per-chr concordance, feeds the line plot"
    echo "- panel_vs_baseline.tsv   -- paired Wilcoxon vs NAT_HGDP baseline"
    echo "- plots/*.pdf             -- six PDFs (boxplot / weighted bars / per-chr / heatmap / MXB scatter / delta forest)"
    echo
    echo "Full per-hap table lives at \`04_accuracy/accuracy_long.tsv\` on the cluster."
} > "$REVIEW_DIR/README.md"

echo "[accuracy-pilot] review artifacts under $REVIEW_DIR"
echo "[accuracy-pilot] commit + push with:"
echo "  git add results/accuracy_pilot"
echo "  git commit -m 'accuracy pilot: ${n_done}/198 combos'"
echo "  git push -u origin claude/setup-shapeit5-conda-xlcdl"
