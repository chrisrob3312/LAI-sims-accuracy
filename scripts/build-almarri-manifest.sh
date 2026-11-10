#!/usr/bin/env bash
# build-almarri-manifest.sh -- sample -> population for the Almarri 2021 genomes
#
# 4a resolves a sample's super-pop through MXB, then SGDP, then the gnomAD
# HGDP+1KG metadata. Almarri is in NONE of those, so without this manifest
# every Almarri sample lands on super_pop=NA and is hard-excluded by
#   elif sp == "NA": included = 0
# -- all 137, the largest single MEN contribution, prepped by 1d and merged by
# 1c only to be dropped at the homogeneity filter in silence.
#
# 4a falls back to the ID pattern <COUNTRY><number> when this file is absent,
# which recovers the samples but gives coarser labels. This gives the real ones.
#
# Usage, on the cluster:
#   bash scripts/build-almarri-manifest.sh > reference_ids/almarri_manifest.tsv
#   wc -l reference_ids/almarri_manifest.tsv     # expect 137 + header
set -euo pipefail
META="${ALMARRI_APPG_META:-$HOME/almarri_vcfs/meta/appg.meta}"
[[ -s "$META" ]] || { echo "ERROR: no Almarri metadata at ${META}" >&2
                      echo "       set ALMARRI_APPG_META to its path" >&2; exit 1; }
printf 'sample\tpop\n'
# Column 3 is the sample id (IRAQ1, UAE57, ...), column 4 the population
# (IRAQ, UAE, ...). Both are checked rather than assumed: a layout change
# would otherwise write a manifest of empty labels.
awk -F'[ \t]+' 'NR>1 && NF>=4 && $3 != "" && $4 != "" {print $3"\t"$4}' "$META" \
  | sort -u
