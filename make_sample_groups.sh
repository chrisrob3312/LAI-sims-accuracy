#!/usr/bin/env bash
# Emit a 2-column sample -> superpop TSV for `bcftools +fill-tags -S`,
# mapping gnomad_meta_updated.tsv labels to 8 superpops:
#   AFR, AMR, EUR, EAS, SAS, CSA, OCE, MEN
# and appending MXB sample IDs as AMR.
#
# Usage:
#   make_sample_groups.sh gnomad_meta_updated.tsv reference_ids/MXB50genomes_popinfo.tsv > sample_groups.tsv
#
# Input formats:
#   gnomad_meta_updated.tsv: <sample>\t<pop>\t<subpop>
#     1KG rows use AFR/AMR/EUR/EAS/SAS in column 2.
#     HGDP rows use Africa/America/Europe/East_Asia/Central_South_Asia/Oceania/Middle_East.
#   MXB popinfo: header row + Sample_ID in column 1 (rest ignored).
#
# Unmapped populations print a WARN to stderr and are skipped.

set -euo pipefail

META="${1:-}"
MXB="${2:-}"

if [[ -z "$META" || -z "$MXB" ]]; then
    echo "usage: $0 <gnomad_meta_updated.tsv> <MXB50genomes_popinfo.tsv> > sample_groups.tsv" >&2
    exit 1
fi

awk -F'\t' '
{
    pop = $2
    if      (pop == "AFR" || pop == "Africa")            super = "AFR"
    else if (pop == "AMR" || pop == "America")           super = "AMR"
    else if (pop == "EUR" || pop == "Europe")            super = "EUR"
    else if (pop == "EAS" || pop == "East_Asia")         super = "EAS"
    else if (pop == "SAS")                               super = "SAS"
    else if (pop == "Central_South_Asia")                super = "CSA"
    else if (pop == "Oceania")                           super = "OCE"
    else if (pop == "Middle_East")                       super = "MEN"
    else { printf("WARN: unmapped pop \"%s\" for sample %s\n", pop, $1) > "/dev/stderr"; next }
    print $1 "\t" super
}' "$META"

awk -F'\t' 'NR > 1 { print $1 "\tAMR" }' "$MXB"
