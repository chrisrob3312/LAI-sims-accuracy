#!/usr/bin/env bash
# build-almarri-manifest.sh -- sample -> population for the Almarri 2021 genomes
#
# 4a resolves a sample's super-pop through MXB, then SGDP, then the gnomAD
# HGDP+1KG metadata. Almarri is in NONE of those, so without this manifest
# every Almarri sample lands on super_pop=NA and is hard-excluded by
#   elif sp == "NA": included = 0
# -- all 137, the largest single MEN contribution, prepped by 1d and merged by
# 1c only to be dropped at the homogeneity filter in silence. 4a also reads
# this file to put the Almarri genomes in the MEN anchor draw pool.
#
# KEY = THE PANEL SAMPLE NAME. appg.meta has five columns:
#   long_ID  short_ID  pop_ID  population  annotation
#   longranger222_wgs_26958_APPG7555879_GRCh38_gatk  APPG7555879  UAE19  UAE  EmiratiD
# The phased panel carries long_ID (the sample name in the longranger VCF
# headers that 1d merges). This script used to key on pop_ID (UAE19), the name
# in the distributed .fam, which matches NO panel sample -- so the manifest
# existed and still labelled nothing.
#
# Usage, on the cluster:
#   bash scripts/build-almarri-manifest.sh > reference_ids/almarri_manifest.tsv
#   wc -l reference_ids/almarri_manifest.tsv     # expect 137 + header
# Then confirm every key is a panel sample (see the check printed on stderr).
set -euo pipefail
META="${ALMARRI_APPG_META:-$HOME/almarri_vcfs/meta/appg.meta}"
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

printf 'sample\tpop\tpop_id\tshort_id\tannotation\n'
awk -F'[ \t]+' 'NR>1 && NF>=4 {print $1"\t"$4"\t"$3"\t"$2"\t"(NF>=5 ? $5 : "")}' "$META" \
  | sort -u

N=$(awk 'NR>1 && NF>=4' "$META" | wc -l)
echo "[almarri-manifest] ${N} rows written from ${META} (keyed on long_ID)" >&2
echo "[almarri-manifest] verify against the panel, e.g.:" >&2
echo "  comm -3 <(awk -F'\\t' 'NR>1{print \$1}' reference_ids/almarri_manifest.tsv | sort) \\" >&2
echo "          <(bcftools query -l <PROJECT_ROOT>/01_merged_phased_panel/almarri_prepped.chr22.bcf | sort)" >&2
echo "  (no output = every panel Almarri sample is labelled and nothing extra)" >&2
