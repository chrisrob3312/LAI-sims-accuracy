#!/usr/bin/env bash
# build-almarri-manifest.sh -- panel sample name -> population, Almarri 2021
#
# 4a resolves a sample's super-pop through MXB, then SGDP, then the gnomAD
# HGDP+1KG metadata. Almarri is in NONE of those; without this manifest 4a can
# still tell an Almarri genome from its ID (-> MEN) but not its population, so
# every one is labelled ALMARRI_UNKNOWN_POP in the discovery tables.
#
# KEY = THE PANEL SAMPLE NAME, matched to appg.meta by the APPG NUMBER.
#   appg.meta columns: long_ID  short_ID  pop_ID  population  annotation
#     longranger222_wgs_26958_APPG7555879_GRCh38_gatk  APPG7555879  UAE19  UAE  EmiratiD
# The panel carries the sample names of the per-sample VCF headers 1d merged.
# Those are USUALLY appg.meta's long_ID, but not always -- one panel sample is
#     longranger222_wgs_27426_APPG7555931_GRCh38_gatk
# where appg.meta has
#     longranger222_wgs_APPG7555931_GRCh38_gatk
# so keying on long_ID left it unlabelled. The APPG number (short_ID) is what
# both agree on: each panel name is matched by the one APPG<digits> token it
# contains, and the PANEL name is written as the key, which is what 4a looks up.
#
# Usage, on the cluster, inside the shapeit5 env (bcftools for a BCF argument):
#   bash scripts/build-almarri-manifest.sh <PANEL> > reference_ids/almarri_manifest.tsv
# PANEL = any readable Almarri prepped BCF (its sample list is the panel's
# Almarri sample names), or a text file of those names, one per line.
# Exit status 3 = some panel name could not be matched (listed on stderr; the
# rows that did match are still written).
set -euo pipefail
META="${ALMARRI_APPG_META:-$HOME/almarri_vcfs/meta/appg.meta}"
PANEL="${1:-}"
[[ -n "$PANEL" ]] || { echo "usage: $0 <almarri_prepped.chrN.bcf | panel_ids.txt> > almarri_manifest.tsv" >&2; exit 2; }
[[ -s "$PANEL" ]] || { echo "ERROR: ${PANEL} is missing or empty" >&2; exit 2; }
[[ -s "$META" ]] || { echo "ERROR: no Almarri metadata at ${META}" >&2
                      echo "       set ALMARRI_APPG_META to its path" >&2; exit 1; }

# Check the layout instead of assuming it: column 1 must look like a long_ID
# and column 2 like a short_ID on every data row.
BAD=$(awk -F'[ \t]+' 'NR>1 && NF && !($1 ~ /APPG[0-9]+/ && $2 ~ /^APPG[0-9]+$/) {n++} END {print n+0}' "$META")
if [[ "$BAD" -gt 0 ]]; then
    echo "ERROR: ${BAD} data row(s) of ${META} do not have long_ID/short_ID in" >&2
    echo "       columns 1/2 -- the layout is not what this script expects:" >&2
    awk 'NR<=3' "$META" >&2
    exit 1
fi

IDS=$(mktemp)
trap 'rm -f "$IDS"' EXIT
case "$PANEL" in
    *.bcf|*.vcf.gz|*.vcf) bcftools query -l "$PANEL" > "$IDS" ;;
    *)                    awk 'NF {print $1}' "$PANEL" > "$IDS" ;;
esac
[[ -s "$IDS" ]] || { echo "ERROR: no sample names read from ${PANEL}" >&2; exit 2; }

printf 'sample\tpop\tpop_id\tshort_id\tannotation\n'
set +e
awk -F'[ \t]+' -v OFS='\t' '
    FNR == NR { if (FNR > 1 && NF >= 4) { long[$2] = $1; pop[$2] = $4; pid[$2] = $3
                                          ann[$2] = (NF >= 5 ? $5 : "") }
                next }
    {
        s = $1; n = 0; t = s; tok = ""
        while (match(t, /APPG[0-9]+/)) { n++; tok = substr(t, RSTART, RLENGTH); t = substr(t, RSTART + RLENGTH) }
        if (n != 1)          { printf "[almarri-manifest] UNMATCHED %s: %d APPG tokens\n", s, n > "/dev/stderr"; bad++; next }
        if (!(tok in pop))   { printf "[almarri-manifest] UNMATCHED %s: %s is not in appg.meta\n", s, tok > "/dev/stderr"; bad++; next }
        if (tok in seen)     { printf "[almarri-manifest] DUPLICATE %s: %s already matched %s\n", s, tok, seen[tok] > "/dev/stderr"; bad++; next }
        seen[tok] = s
        if (long[tok] != s)  { printf "[almarri-manifest] note: panel name %s differs from appg.meta long_ID %s (matched on %s)\n", s, long[tok], tok > "/dev/stderr"; diff++ }
        print s, pop[tok], pid[tok], tok, ann[tok]; ok++
    }
    END {
        for (k in pop) if (!(k in seen)) { notin++; if (notin <= 10) printf "[almarri-manifest] in appg.meta, not in the panel: %s (%s)\n", k, pop[k] > "/dev/stderr" }
        printf "[almarri-manifest] %d panel samples labelled (%d under a name that differs from appg.meta), %d unmatched; %d appg.meta rows not in the panel\n", ok, diff, bad, notin > "/dev/stderr"
        exit (bad > 0 ? 3 : 0)
    }' "$META" "$IDS" | sort -t $'\t' -k1,1
rc=${PIPESTATUS[0]}
set -e
exit "$rc"
