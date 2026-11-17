#!/usr/bin/env bash
# =============================================================================
# admix-site-list.sh -- the ADMIXTURE site list for one chromosome, from the
# published PRE-phasing soft-union (01_merged_phased_panel/merged_chrN.softunion.bcf).
#
#   bash scripts/admix-site-list.sh <chr> <out.bed>     # site list + <out>.summary.tsv
#   bash scripts/admix-site-list.sh --groups <chr> <out.tsv>   # only the groups
#
# THE RULE (17 Nov 2026): every super-pop, its cohorts pooled, <= 7% missing
# (ADMIX_SUPERPOP_MAXMISS) and every cohort, its super-pops pooled, <= 15%
# (ADMIX_COHORT_MAXMISS); margins under ADMIX_MISS_MIN_N (5) samples are not
# judged alone. scripts/admix-site-margins.py applies it.
#
# THE GROUPS are rebuilt exactly as 1c built them (scripts/admix-rule-compare.sh
# reproduces 1c's own list byte for byte with them): cohort = the file a sample
# came from, with 1c's ID de-duplication order -- HGDP+1KG and MXB, then SGDP,
# then Almarri -- and the HGDP+1KG source read from the soft-union's recipe
# stamp; super-pop from sample_groups.tsv (CSA->SAS, MID->MEN), the SGDP
# manifest, MXB = AMR, Almarri = MEN (scripts/build-missingness-groups.py).
#
# Why here and not in 1c: 1c measured a different rule (cohort x super-pop
# groups <= 15%, panel <= 7%), and changing 1c re-phases the panel. The site
# list only decides which markers ADMIXTURE reads, so 4a computes it itself.
# =============================================================================
set -uo pipefail
MODE=sites
if [[ "${1:-}" == "--groups" ]]; then MODE=groups; shift; fi
CHR="${1:?chromosome}"; OUT="${2:?output path}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
SGDP_MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
EXCLUDE_IDS="${EXCLUDE_IDS:-${REPO_DIR}/reference_ids/panel_exclude_ids.txt}"
SP_MAX="${ADMIX_SUPERPOP_MAXMISS:-0.07}"
CO_MAX="${ADMIX_COHORT_MAXMISS:-0.15}"
MIN_N="${ADMIX_MISS_MIN_N:-5}"
GROUPER="${REPO_DIR}/scripts/build-missingness-groups.py"
FILTER="${REPO_DIR}/scripts/admix-site-margins.py"

die () { echo "ERROR [chr${CHR}]: $*" >&2; exit 1; }
for f in "$GROUPER" "$FILTER"; do [[ -f "$f" ]] || die "$f not found (set REPO_DIR to the repo)"; done
for f in "$SAMPLE_GROUPS" "$SGDP_MANIFEST"; do [[ -s "$f" ]] || die "missing $f"; done
SU="${PHASED_DIR}/merged_chr${CHR}.softunion.bcf"
[[ -s "$SU" ]] || die "missing $SU (1c publishes it in PHASE_MODE=chunked)"
[[ -s "${SU}.recipe" ]] || die "missing ${SU}.recipe (1c's stamp names the HGDP+1KG source)"
SRC=$(sed -n 's/.* src=\([^ ]*\).*/\1/p' "${SU}.recipe")
HG="${SRC//__CHR__/$CHR}"
[[ -s "$HG" ]] || die "cannot read the HGDP+1KG source named in ${SU}.recipe: '${HG}'"

WORK=$(mktemp -d "${OUT}.work.XXXXXX") || die "cannot make a work folder beside $OUT"
trap 'rm -rf "$WORK"' EXIT

# Cohort membership in 1c's order (1c step 2: base = HGDP+1KG + MXB; each added
# cohort loses any ID already present or excluded by name).
bcftools query -l "$SU" > "${WORK}/samples.txt" || die "cannot read samples of $SU"
sort -u "${WORK}/samples.txt" > "${WORK}/samples.sorted"
if [[ -s "$EXCLUDE_IDS" ]]; then
    awk 'NF && $1 !~ /^#/ {print $1}' "$EXCLUDE_IDS" | sort -u > "${WORK}/excl.txt"
else
    echo "NOTE [chr${CHR}]: no ${EXCLUDE_IDS}; no named exclusions applied" >&2; : > "${WORK}/excl.txt"
fi
MXB_F="${PHASED_DIR}/mxb_lifted.chr${CHR}.bcf"
[[ -s "$MXB_F" ]] || die "missing $MXB_F"
bcftools query -l "$MXB_F" | sort -u > "${WORK}/mxb.txt" || die "cannot read samples of $MXB_F"
{ bcftools query -l "$HG"; cat "${WORK}/mxb.txt"; } | sort -u > "${WORK}/seen.txt" || die "cannot read samples of $HG"
COH=(--cohort "MXB=${WORK}/mxb.txt")
for c in sgdp almarri; do
    F="${PHASED_DIR}/${c}_prepped.chr${CHR}.bcf"
    [[ -s "$F" ]] || { echo "NOTE [chr${CHR}]: no ${F}; ${c} not in this panel" >&2; continue; }
    bcftools query -l "$F" | sort -u | comm -23 - "${WORK}/seen.txt" | comm -23 - "${WORK}/excl.txt" \
        | comm -12 - "${WORK}/samples.sorted" > "${WORK}/${c}.txt"
    sort -u "${WORK}/seen.txt" "${WORK}/${c}.txt" -o "${WORK}/seen.txt"
    COH+=(--cohort "$(echo "$c" | tr a-z A-Z)=${WORK}/${c}.txt")
done
python3 -I "$GROUPER" --samples "${WORK}/samples.txt" "${COH[@]}" --default-cohort HGDP1KG \
    --groups "$SAMPLE_GROUPS" --sgdp-manifest "$SGDP_MANIFEST" > "${WORK}/groups.tsv" \
    || die "build-missingness-groups.py failed"
if [[ "$MODE" == groups ]]; then
    cp "${WORK}/groups.tsv" "$OUT"; exit 0
fi

FMT=$(python3 -I "$FILTER" --groups "${WORK}/groups.tsv" --print-format) || die "no query format"
bcftools +fill-tags "$SU" -Ou -- -S "${WORK}/groups.tsv" -t AN \
  | bcftools query -f "$FMT" \
  | python3 -I "$FILTER" --groups "${WORK}/groups.tsv" \
        --max-superpop-miss "$SP_MAX" --max-cohort-miss "$CO_MAX" --min-n "$MIN_N" \
        --out "${OUT}.partial" --summary "${OUT}.summary.tsv.partial" \
  || die "fill-tags | query | admix-site-margins.py failed"
[[ -s "${OUT}.partial" ]] || die "no site passed the rule"
cp "${WORK}/groups.tsv" "${OUT}.groups.tsv"
mv "${OUT}.summary.tsv.partial" "${OUT}.summary.tsv"
mv "${OUT}.partial" "$OUT"
