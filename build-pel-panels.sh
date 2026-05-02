#!/usr/bin/env bash
# ============================================================================
# build-pel-panels.sh
#
# Derive RFMix-reference sample lists for legacy panels 2 (PEL) and 3
# (PEL + 1KG-EAS) from gnomad_meta_updated.tsv. Matches the published
# Honorato-Mauer et al. 2025 (AJHG, PMID 39753130) panel composition:
#   - Panel 2 = 1KG PEL, unrelated only           (paper: n = 85)
#   - Panel 3 = Panel 2 + 1KG EAS, unrelated only (paper: n = 85 + 505 = 590)
#
# "Unrelated" filter strategy:
#   Primary:  gnomAD's `release` flag (true means included in the v3 release
#             set, which is the kinship-pruned, QC-passed unrelated subset).
#             This matches the paper's "We filtered to keep only unrelated
#             individuals" criterion.
#   Fallback: if the `release` column isn't found, fall back to filtering
#             via $OUTLIERS (related_outliers.txt).
#
# For EAS, restrict to 1KG sub-pops (CHB/CHS/JPT/CDX/KHV) -- excludes HGDP
# East_Asia samples (Han, Japanese, Cambodian, Dai, ...) which share
# Genetic.region=EAS but are not what the paper used.
#
# Outputs (in $REFS):
#   pel_rfmix.txt        1KG PEL unrelated         (~85)
#   pel_eas_rfmix.txt    1KG PEL + 1KG EAS, unrel  (~590)
#
# Usage:
#   ./build-pel-panels.sh <gnomad_meta_updated.tsv>
# ============================================================================

set -euo pipefail

META="${1:-}"
REFS="${REFS:-reference_ids}"
OUTLIERS="${OUTLIERS:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/sample_map_files/related_outliers.txt}"

if [[ -z "$META" || ! -s "$META" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv>" >&2
    exit 1
fi

mkdir -p "$REFS"

# Read header once (avoid SIGPIPE under pipefail).
header=$(head -1 "$META")
SAMPLE_COL=$(  awk -F'\t' -v t="project_meta.sample_id"   '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")
POP_COL=$(     awk -F'\t' -v t="hgdp_tgp_meta.Population" '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")
RELEASE_COL=$( awk -F'\t' -v t="release"                 '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")

[[ -z "$SAMPLE_COL" ]] && { echo "ERROR: column 'project_meta.sample_id' not found"; exit 1; }
[[ -z "$POP_COL"    ]] && { echo "ERROR: column 'hgdp_tgp_meta.Population' not found"; exit 1; }

echo "Using sample col $SAMPLE_COL, pop col $POP_COL" >&2

# Decide which unrelated-filter to use.
if [[ -n "$RELEASE_COL" ]]; then
    UNREL_MODE="release_col"
    echo "Unrelated filter: gnomAD 'release' column ($RELEASE_COL) == true" >&2
    # Helper: emit IDs matching a pop regex AND release=true.
    extract_unrelated () {
        local pops_regex="$1"
        awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" -v rc="$RELEASE_COL" -v re="$pops_regex" '
            NR > 1 && $pc ~ re && tolower($rc) == "true" { print $sc }' "$META"
    }
else
    UNREL_MODE="related_outliers_file"
    [[ -s "$OUTLIERS" ]] || { echo "ERROR: 'release' col missing AND no $OUTLIERS"; exit 1; }
    echo "Unrelated filter: $OUTLIERS ($(wc -l < "$OUTLIERS") IDs)" >&2
    extract_unrelated () {
        local pops_regex="$1"
        awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" -v re="$pops_regex" '
            NR > 1 && $pc ~ re { print $sc }' "$META" \
            | grep -vxFf "$OUTLIERS" || true
    }
fi

# Panel 2: PEL only (1KG, unrelated)
extract_unrelated '^PEL$' > "${REFS}/pel_rfmix.txt"

# Panel 3: PEL + 1KG EAS (CHB, CHS, JPT, CDX, KHV), unrelated.
{
    extract_unrelated '^PEL$'
    extract_unrelated '^(CHB|CHS|JPT|CDX|KHV)$'
} > "${REFS}/pel_eas_rfmix.txt"

echo "Panel 2 (PEL,         unrel) : $(wc -l < ${REFS}/pel_rfmix.txt) IDs   [paper: 85]"
echo "Panel 3 (PEL + 1KG EAS, unrel) : $(wc -l < ${REFS}/pel_eas_rfmix.txt) IDs   [paper: 590]"
echo ""
echo "Filter mode: $UNREL_MODE"
echo ""
echo "Now 3b_wgs-rfmix-jointcall_clm.sh will run all 4 panels by default."
echo "Override panel selection at submit time, e.g.:"
echo "  PANELS_TO_RUN='NAT_HGDP NAT_PEL NAT_PEL_EAS NAT_HGDPMXB' \\"
echo "    ADMIX_POP=Brasa GEN=12 sbatch 3b_wgs-rfmix-jointcall_clm.sh"
