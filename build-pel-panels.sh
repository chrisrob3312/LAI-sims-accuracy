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
#
# The gnomAD meta file is wide (>150 cols); we resolve the relevant column
# indices by header name to be robust to column-order changes:
#   project_meta.sample_id        -> sample ID
#   hgdp_tgp_meta.Population      -> 1KG sub-pop (PEL, CEU, IBS, YRI, ...)
#                                     or HGDP pop name (Brahui, Yoruba, ...)
#   hgdp_tgp_meta.Genetic.region  -> super-pop (AFR/AMR/EUR/EAS/SAS/CSA/OCE/MID)
# ============================================================================

set -euo pipefail

META="${1:-}"
REFS="${REFS:-reference_ids}"

if [[ -z "$META" || ! -s "$META" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv>" >&2
    exit 1
fi

mkdir -p "$REFS"

SAMPLE_COL=$(head -1 "$META" | tr '\t' '\n' | awk '$0=="project_meta.sample_id"{print NR; exit}')
POP_COL=$(   head -1 "$META" | tr '\t' '\n' | awk '$0=="hgdp_tgp_meta.Population"{print NR; exit}')
REGION_COL=$(head -1 "$META" | tr '\t' '\n' | awk '$0=="hgdp_tgp_meta.Genetic.region"{print NR; exit}')

[[ -z "$SAMPLE_COL" ]] && { echo "ERROR: column 'project_meta.sample_id' not found"; exit 1; }
[[ -z "$POP_COL"    ]] && { echo "ERROR: column 'hgdp_tgp_meta.Population' not found"; exit 1; }
[[ -z "$REGION_COL" ]] && { echo "ERROR: column 'hgdp_tgp_meta.Genetic.region' not found"; exit 1; }

echo "Using sample col $SAMPLE_COL, pop col $POP_COL, region col $REGION_COL" >&2

# PEL only -- match on the 1KG sub-population label
awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" '
NR > 1 && $pc == "PEL" { print $sc }' "$META" > "${REFS}/pel_rfmix.txt"

# PEL + EAS  (legacy panel 3 mixes Latino-NAT proxy with EAS samples)
{
    awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" '
        NR > 1 && $pc == "PEL" { print $sc }' "$META"
    awk -F'\t' -v sc="$SAMPLE_COL" -v rc="$REGION_COL" '
        NR > 1 && $rc == "EAS" { print $sc }' "$META"
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
