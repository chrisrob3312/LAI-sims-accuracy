#!/usr/bin/env bash
# =============================================================================
# admix-rule-compare.sh -- which soft-union sites ADMIXTURE may use, under
# four versions of the missingness rule, measured the way 1c measures it.
# READ-ONLY on the pipeline: writes only into $OUT (a folder of its own).
#
#   coded       every cohort x super-pop group <= 15% missing, panel <= 7%
#               (1c since 12 Nov; must reproduce 1c's published site list --
#               checked, so a wrong reconstruction of the groups shows)
#   sp07_tot15  every super-pop <= 7% missing, panel <= 15%  (16 Nov wording)
#   grp07       every cohort x super-pop group <= 7%, panel <= 7%
#   sp15_tot07  every super-pop <= 15%, panel <= 7%
# Groups below ADMIX_MISS_MIN_N (5) samples are not judged on their own but
# count toward the panel. Missingness is 1 - AN / (2n) on the PRE-phasing
# soft-union (01_merged_phased_panel/merged_chrN.softunion.bcf), exactly as
# 1c step 5b does it (scripts/admix-site-filter.py).
#
# The groups are rebuilt as 1c built them: cohort = the file a sample came
# from, with 1c's ID de-duplication order (HGDP+1KG and MXB, then SGDP, then
# Almarri); the HGDP+1KG source is read from the soft-union's recipe stamp.
#
#   bash scripts/admix-rule-compare.sh                      # chr22, all four rules
#   CHRS="$(seq 1 22)" RULES=sp07_tot15 bash scripts/admix-rule-compare.sh
#       -> $OUT/sp07_tot15/merged_chrN.admix_sites.bed for all 22, which 4a
#          reads with ADMIX_SITES_TPL=$OUT/sp07_tot15/merged_chr%s.admix_sites.bed
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
SGDP_MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
EXCLUDE_IDS="${EXCLUDE_IDS:-${REPO_DIR}/reference_ids/panel_exclude_ids.txt}"
OUT="${OUT:-${PROJECT_ROOT}/04_homogeneity_panel/admix_rule_compare}"
CHRS="${CHRS:-22}"
RULES="${RULES:-coded sp07_tot15 grp07 sp15_tot07}"
MIN_N="${ADMIX_MISS_MIN_N:-5}"
FILTER="${REPO_DIR}/scripts/admix-site-filter.py"
GROUPER="${REPO_DIR}/scripts/build-missingness-groups.py"
for f in "$FILTER" "$GROUPER"; do [[ -f "$f" ]] || { echo "ERROR: $f not found (run from the repo, or set REPO_DIR)" >&2; exit 2; }; done
for f in "$SAMPLE_GROUPS" "$SGDP_MANIFEST"; do [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 2; }; done

rule_spec () {   # -> "<groups: co|sp> <max group> <max total>"
    case "$1" in
        coded)      echo "co 0.15 0.07" ;;
        sp07_tot15) echo "sp 0.07 0.15" ;;
        grp07)      echo "co 0.07 0.07" ;;
        sp15_tot07) echo "sp 0.15 0.07" ;;
        *) echo "ERROR: unknown rule '$1'" >&2; return 1 ;;
    esac
}
for r in $RULES; do rule_spec "$r" > /dev/null || exit 2; mkdir -p "${OUT}/${r}"; done
WORK="${OUT}/work"; mkdir -p "$WORK"

for CHR in $CHRS; do
    SU="${PHASED_DIR}/merged_chr${CHR}.softunion.bcf"
    PUB_BED="${PHASED_DIR}/merged_chr${CHR}.admix_sites.bed"
    [[ -s "$SU" ]] || { echo "ERROR: missing $SU (1c publishes it in PHASE_MODE=chunked)" >&2; exit 1; }
    SRC=$(sed -n 's/.* src=\([^ ]*\).*/\1/p' "${SU}.recipe")
    HG="${SRC//__CHR__/$CHR}"
    echo "=================== chr${CHR}"
    echo "soft-union: $SU"
    echo "HGDP+1KG source (from ${SU}.recipe): ${HG:-NOT FOUND}"
    [[ -s "$HG" ]] || { echo "ERROR: cannot read the HGDP+1KG source named in ${SU}.recipe" >&2; exit 1; }

    # Cohort membership in 1c's order.
    bcftools query -l "$SU" > "${WORK}/samples.txt" || exit 1
    sort -u "${WORK}/samples.txt" > "${WORK}/samples.sorted"
    if [[ -s "$EXCLUDE_IDS" ]]; then
        awk 'NF && $1 !~ /^#/ {print $1}' "$EXCLUDE_IDS" | sort -u > "${WORK}/excl.txt"
    else
        echo "  NOTE: no ${EXCLUDE_IDS}; no named exclusions applied"; : > "${WORK}/excl.txt"
    fi
    { bcftools query -l "$HG"; bcftools query -l "${PHASED_DIR}/mxb_lifted.chr${CHR}.bcf"; } | sort -u > "${WORK}/seen.txt" || exit 1
    bcftools query -l "${PHASED_DIR}/mxb_lifted.chr${CHR}.bcf" | sort -u > "${WORK}/mxb.txt" || exit 1
    COH=(--cohort "MXB=${WORK}/mxb.txt")
    for c in sgdp almarri; do
        F="${PHASED_DIR}/${c}_prepped.chr${CHR}.bcf"
        [[ -s "$F" ]] || { echo "  (no ${F}: ${c} not in this panel)"; continue; }
        bcftools query -l "$F" | sort -u | comm -23 - "${WORK}/seen.txt" | comm -23 - "${WORK}/excl.txt" \
            | comm -12 - "${WORK}/samples.sorted" > "${WORK}/${c}.txt"
        sort -u "${WORK}/seen.txt" "${WORK}/${c}.txt" -o "${WORK}/seen.txt"
        COH+=(--cohort "$(echo "$c" | tr a-z A-Z)=${WORK}/${c}.txt")
        echo "  ${c}: $(wc -l < "${WORK}/${c}.txt") samples in the soft-union"
    done
    python3 -I "$GROUPER" --samples "${WORK}/samples.txt" "${COH[@]}" --default-cohort HGDP1KG \
        --groups "$SAMPLE_GROUPS" --sgdp-manifest "$SGDP_MANIFEST" > "${WORK}/groups.co.tsv" || exit 1
    # Super-pop groups: the same samples, cohort prefix removed (NA kept so it
    # still counts toward the panel).
    awk -F'\t' '{g = $2; sub(/^[^_]*_/, "", g); print $1 "\tSP" g}' "${WORK}/groups.co.tsv" > "${WORK}/groups.sp.tsv"
    echo "  super-pop groups: $(cut -f2 "${WORK}/groups.sp.tsv" | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"

    for G in co sp; do
        need=0
        for r in $RULES; do read -r g _ _ <<<"$(rule_spec "$r")"; [[ "$g" == "$G" ]] && need=1; done
        [[ "$need" == 1 ]] || continue
        FMT=$(python3 -I "$FILTER" --groups "${WORK}/groups.${G}.tsv" --print-format) || exit 1
        bcftools +fill-tags "$SU" -Ou -- -S "${WORK}/groups.${G}.tsv" -t AN \
            | bcftools query -f "$FMT" | gzip -1 > "${WORK}/an.${G}.tsv.gz" || exit 1
    done

    printf '  %-11s %-20s %6s %6s %12s %7s  %s\n' rule groups grp tot positions pct "vs 1c's published list ($(wc -l < "$PUB_BED") positions)"
    for r in $RULES; do
        read -r G GM TM <<<"$(rule_spec "$r")"
        BED="${OUT}/${r}/merged_chr${CHR}.admix_sites.bed"
        zcat "${WORK}/an.${G}.tsv.gz" | python3 -I "$FILTER" --groups "${WORK}/groups.${G}.tsv" \
            --max-group-miss "$GM" --max-total-miss "$TM" --min-group-n "$MIN_N" \
            --out "$BED" --summary "${OUT}/${r}/merged_chr${CHR}.admix_sites.summary.tsv" \
            2> "${OUT}/${r}/chr${CHR}.log" || { cat "${OUT}/${r}/chr${CHR}.log" >&2; exit 1; }
        TOTPOS=$(sed -n 's/.* positions=\([0-9]*\) .*/\1/p' "${OUT}/${r}/merged_chr${CHR}.admix_sites.summary.tsv")
        N=$(wc -l < "$BED")
        BOTH=$(comm -12 <(sort "$BED") <(sort "$PUB_BED") | wc -l)
        CMP="same $BOTH; only here $((N - BOTH)); only in 1c's $(( $(wc -l < "$PUB_BED") - BOTH ))"
        [[ "$r" == coded ]] && { if cmp -s "$BED" "$PUB_BED"; then CMP="IDENTICAL (reconstruction checks out)"; else CMP="DIFFERS: ${CMP} -- the groups were not rebuilt as 1c built them"; fi; }
        printf '  %-11s %-20s %6s %6s %12s %6.1f%%  %s\n' "$r" "$([[ $G == co ]] && echo 'cohort x super-pop' || echo 'super-pop')" \
            "$GM" "$TM" "$N" "$(awk -v n="$N" -v t="$TOTPOS" 'BEGIN {print (t ? 100 * n / t : 0)}')" "$CMP"
    done
    echo "  each rule's own report (groups over the limit, and how many sites for no other reason):"
    for r in $RULES; do sed "s/^/    [${r}] /" "${OUT}/${r}/chr${CHR}.log"; done
    rm -f "${WORK}"/an.*.tsv.gz
done
echo "site lists and per-group summaries -> ${OUT}/<rule>/"
