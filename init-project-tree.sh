#!/usr/bin/env bash
# ============================================================================
# init-project-tree.sh
#
# Creates the project directory layout under $PROJECT_ROOT and runs the
# one-shot helpers (make_sample_groups.sh, build-pel-panels.sh,
# build-panel-keep-files.sh) so the SLURM stages 1a/1b/2b/3b have all
# their inputs in place.
#
# Usage:
#   ./init-project-tree.sh
#
# Override PROJECT_ROOT via env var if you want a different location.
# ============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
META="${META:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
MXB_POPINFO="${MXB_POPINFO:-reference_ids/MXB50genomes_popinfo.tsv}"

[[ -s "$META"        ]] || { echo "ERROR: missing $META"; exit 1; }
[[ -s "$MXB_POPINFO" ]] || { echo "ERROR: missing $MXB_POPINFO (run from repo root)"; exit 1; }

echo "[$(date +%T)] Creating project tree at $PROJECT_ROOT"
mkdir -p \
    "${PROJECT_ROOT}/01_merged_phased_panel" \
    "${PROJECT_ROOT}/02_simulations" \
    "${PROJECT_ROOT}/03_rfmix" \
    "${PROJECT_ROOT}/04_accuracy" \
    "${PROJECT_ROOT}/panel_keep_files" \
    "${PROJECT_ROOT}/logs"

echo "[$(date +%T)] Building sample_groups.tsv"
./make_sample_groups.sh "$META" "$MXB_POPINFO" > "${PROJECT_ROOT}/sample_groups.tsv"
wc -l "${PROJECT_ROOT}/sample_groups.tsv"

echo "[$(date +%T)] Building pel_rfmix.txt and pel_eas_rfmix.txt"
./build-pel-panels.sh "$META"

echo "[$(date +%T)] Building panel keep-files"
PROJECT_ROOT="$PROJECT_ROOT" ./build-panel-keep-files.sh

echo "[$(date +%T)] # Project tree ready at $PROJECT_ROOT"
echo ""
echo "Layout:"
echo "  ${PROJECT_ROOT}/"
echo "    01_merged_phased_panel/   (1a + 1b outputs)"
echo "    02_simulations/           (2b outputs; place \${ADMIX_POP}.dat / .sample.txt here)"
echo "    03_rfmix/                 (3b outputs)"
echo "    04_accuracy/              (accuracy.R outputs)"
echo "    panel_keep_files/         (built above)"
echo "    logs/                     (SLURM stdout/stderr)"
echo "    sample_groups.tsv         (built above)"
echo ""
echo "Next: sbatch 1a_prep-mxb-liftover_clm.sh"
