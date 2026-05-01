#!/usr/bin/env bash
# Emit a 2-column sample -> superpop TSV for `bcftools +fill-tags -S`,
# mapping the gnomad_meta_updated.tsv `hgdp_tgp_meta.Genetic.region` field
# to 8 superpops:
#   AFR, AMR, EUR, EAS, SAS, CSA, OCE, MEN
# (MEN = MID/MENA in the source file -- renamed for our convention.)
# Then append MXB sample IDs as AMR.
#
# Usage:
#   make_sample_groups.sh gnomad_meta_updated.tsv reference_ids/MXB50genomes_popinfo.tsv > sample_groups.tsv
#
# The gnomAD HGDP+1KG metadata is wide (>150 columns); we look up the column
# index of the desired field by name in the header, so this script doesn't
# break if column order changes.

set -euo pipefail

META="${1:-}"
MXB="${2:-}"

if [[ -z "$META" || -z "$MXB" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv> <MXB50genomes_popinfo.tsv> > sample_groups.tsv" >&2
    exit 1
fi

# Resolve column indices (1-based) by header name.
# Read the header into a variable first to avoid SIGPIPE racing between
# `head` and `awk ... exit` under `set -euo pipefail`.
header=$(head -1 "$META")
SAMPLE_COL=$(awk -F'\t' -v t="project_meta.sample_id"       '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")
POP_COL=$(   awk -F'\t' -v t="hgdp_tgp_meta.Genetic.region" '{for(i=1;i<=NF;i++) if($i==t){print i; exit}}' <<< "$header")

if [[ -z "$SAMPLE_COL" ]]; then
    echo "ERROR: column 'project_meta.sample_id' not found in $META header" >&2
    exit 1
fi
if [[ -z "$POP_COL" ]]; then
    echo "ERROR: column 'hgdp_tgp_meta.Genetic.region' not found in $META header" >&2
    exit 1
fi

echo "Using sample column $SAMPLE_COL ('project_meta.sample_id')" >&2
echo "Using pop    column $POP_COL ('hgdp_tgp_meta.Genetic.region')" >&2

awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" '
NR == 1 { next }                # skip header
{
    pop = $pc
    if      (pop == "AFR") super = "AFR"
    else if (pop == "AMR") super = "AMR"
    else if (pop == "EUR") super = "EUR"
    else if (pop == "EAS") super = "EAS"
    else if (pop == "SAS") super = "SAS"
    else if (pop == "CSA") super = "CSA"
    else if (pop == "OCE") super = "OCE"
    else if (pop == "MID") super = "MEN"
    else { printf("WARN: unmapped pop \"%s\" for sample %s\n", pop, $sc) > "/dev/stderr"; next }
    print $sc "\t" super
}' "$META"

awk -F'\t' 'NR > 1 { print $1 "\tAMR" }' "$MXB"
