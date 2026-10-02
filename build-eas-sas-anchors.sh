#!/usr/bin/env bash
# ============================================================================
# build_eas_sas_anchors.sh
#
# Derive reference_ids/eas_rfmix.txt and reference_ids/sas_rfmix.txt for
# supervised ADMIXTURE K=5 anchoring, matching the recipe in
# build-pel-panels.sh. Filters gnomad_meta_updated.tsv by:
#
#   EAS = release==true AND Population in {CHB, CHS, JPT, CDX, KHV}
#   SAS = release==true AND Population in {GIH, PJL, BEB, STU, ITU}
#         (gnomAD's Genetic.region field lumps 1KG-SAS into CSA, so we
#         restrict to the canonical 1KG-SAS sub-pops by name.)
#
# Usage:
#   ./build_eas_sas_anchors.sh <gnomad_meta_updated.tsv>
# ============================================================================

set -euo pipefail
META="${1:-}"
REFS="${REFS:-reference_ids}"

if [[ -z "$META" || ! -s "$META" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv>" >&2
    exit 1
fi
mkdir -p "$REFS"

header=$(head -1 "$META")
SAMPLE_COL=$(awk -F'\t' -v t="project_meta.sample_id"   '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")
POP_COL=$(   awk -F'\t' -v t="hgdp_tgp_meta.Population" '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")
REL_COL=$(   awk -F'\t' -v t="release"                  '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")

[[ -z "$SAMPLE_COL" ]] && { echo "ERROR: 'project_meta.sample_id' column not found"; exit 1; }
[[ -z "$POP_COL"    ]] && { echo "ERROR: 'hgdp_tgp_meta.Population' column not found"; exit 1; }
[[ -z "$REL_COL"    ]] && { echo "ERROR: 'release' column not found"; exit 1; }

extract_anchors () {
    local pops_regex="$1"
    awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" -v rc="$REL_COL" -v re="$pops_regex" '
        function is_true(v) { v=tolower(v); return v=="true" || v=="t" || v=="1" }
        NR > 1 && $pc ~ re && is_true($rc) { print $sc }' "$META"
}

extract_anchors '^(CHB|CHS|JPT|CDX|KHV)$' > "${REFS}/eas_rfmix.txt"
extract_anchors '^(GIH|PJL|BEB|STU|ITU)$' > "${REFS}/sas_rfmix.txt"

echo "EAS anchors: $(wc -l < "${REFS}/eas_rfmix.txt") IDs   (CHB+CHS+JPT+CDX+KHV, release-set)"
echo "SAS anchors: $(wc -l < "${REFS}/sas_rfmix.txt") IDs   (GIH+PJL+BEB+STU+ITU, release-set)"
