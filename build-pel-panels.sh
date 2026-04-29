#!/usr/bin/env bash
# ============================================================================
# build-pel-panels.sh
#
# Derive RFMix-reference sample lists for legacy panels 2 (PEL only) and
# 3 (PEL + EAS) from gnomad_meta_updated.tsv.
#
# Outputs (in $REFS):
#   pel_rfmix.txt        all 1KG PEL samples (~85)
#   pel_eas_rfmix.txt    1KG PEL + 1KG EAS samples (CHB, CHS, JPT, CDX, KHV)
#
# Usage:
#   ./build-pel-panels.sh <gnomad_meta_updated.tsv>
#
# 3b_wgs-rfmix-jointcall_clm.sh picks these up automatically once present.
# ============================================================================

set -euo pipefail

META="${1:-}"
REFS="${REFS:-reference_ids}"

if [[ -z "$META" || ! -s "$META" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv>" >&2
    exit 1
fi

mkdir -p "$REFS"

# 1KG metadata format: <sample> <super_pop> <sub_pop>
# PEL only
awk -F'\t' '$2 == "AMR" && $3 == "PEL" {print $1}' "$META" > "${REFS}/pel_rfmix.txt"

# PEL + EAS (the legacy "PEL+EAS" panel mixes Latino-NAT proxy with EAS as
# an additional reference for distinguishing NAT from EAS ancestry).
{
    awk -F'\t' '$2 == "AMR" && $3 == "PEL" {print $1}' "$META"
    awk -F'\t' '$2 == "EAS" {print $1}' "$META"
} > "${REFS}/pel_eas_rfmix.txt"

echo "PEL only        : $(wc -l < ${REFS}/pel_rfmix.txt) IDs"
echo "PEL + EAS       : $(wc -l < ${REFS}/pel_eas_rfmix.txt) IDs"
echo "Outputs:"
echo "  ${REFS}/pel_rfmix.txt"
echo "  ${REFS}/pel_eas_rfmix.txt"
echo ""
echo "Now 3b_wgs-rfmix-jointcall_clm.sh will run all 4 panels by default."
echo "Override panel selection at submit time, e.g.:"
echo "  PANELS_TO_RUN='NAT_HGDP NAT_PEL NAT_PEL_EAS NAT_HGDPMXB' \\"
echo "    ADMIX_POP=Brasa GEN=12 sbatch 3b_wgs-rfmix-jointcall_clm.sh"
