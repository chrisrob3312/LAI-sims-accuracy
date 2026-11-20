#!/usr/bin/env bash
# =============================================================================
# admix-rule-compare.sh -- which soft-union sites ADMIXTURE may use, under
# several versions of the missingness rule, measured the way 1c measures it.
# READ-ONLY on the pipeline: writes only into $OUT (a folder of its own).
#
#   sp07_coh15  every super-pop <= 7% missing, every cohort overall <= 15%
#               (THE RULE from 17 Nov; 4a computes it with scripts/admix-site-list.sh)
#   coded       every cohort x super-pop group <= 15% missing, panel <= 7%
#               (1c's rule from 12 Nov; must reproduce 1c's published site list --
#               checked, so a wrong reconstruction of the groups shows)
#   sp07_tot15  every super-pop <= 7% missing, panel <= 15%
#   grp07       every cohort x super-pop group <= 7%, panel <= 7%
#   sp15_tot07  every super-pop <= 15%, panel <= 7%
# Groups or margins below ADMIX_MISS_MIN_N (5) samples are not judged alone.
# Missingness is 1 - AN / (2n) on the PRE-phasing soft-union
# (01_merged_phased_panel/merged_chrN.softunion.bcf).
#
# The groups come from scripts/admix-site-list.sh --groups: cohort = the file a
# sample came from, in 1c's ID de-duplication order; super-pop from the
# metadata, exactly as 1c built them.
#
#   bash scripts/admix-rule-compare.sh                      # chr22, every rule
#   CHRS="$(seq 1 22)" RULES=coded bash scripts/admix-rule-compare.sh
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
export PROJECT_ROOT PHASED_DIR REPO_DIR
OUT="${OUT:-${PROJECT_ROOT}/04_homogeneity_panel/admix_rule_compare}"
CHRS="${CHRS:-22}"
RULES="${RULES:-sp07_coh15 coded sp07_tot15 grp07 sp15_tot07}"
MIN_N="${ADMIX_MISS_MIN_N:-5}"
FILTER="${REPO_DIR}/scripts/admix-site-filter.py"
MARGINS="${REPO_DIR}/scripts/admix-site-margins.py"
LISTER="${REPO_DIR}/scripts/admix-site-list.sh"
for f in "$FILTER" "$MARGINS" "$LISTER"; do [[ -f "$f" ]] || { echo "ERROR: $f not found (run from the repo, or set REPO_DIR)" >&2; exit 2; }; done

rule_spec () {   # -> "<filter: grp|mg> <groups: co|sp> <cap 1> <cap 2>"
    case "$1" in
        sp07_coh15) echo "mg co 0.07 0.15" ;;      # super-pop margin, cohort margin
        coded)      echo "grp co 0.15 0.07" ;;     # each group, panel
        sp07_tot15) echo "grp sp 0.07 0.15" ;;
        grp07)      echo "grp co 0.07 0.07" ;;
        sp15_tot07) echo "grp sp 0.15 0.07" ;;
        *) echo "ERROR: unknown rule '$1'" >&2; return 1 ;;
    esac
}
describe () {
    case "$1" in
        sp07_coh15) echo "super-pop <= 7%, cohort <= 15%" ;;
        coded)      echo "cohort x super-pop <= 15%, panel <= 7%" ;;
        sp07_tot15) echo "super-pop <= 7%, panel <= 15%" ;;
        grp07)      echo "cohort x super-pop <= 7%, panel <= 7%" ;;
        sp15_tot07) echo "super-pop <= 15%, panel <= 7%" ;;
    esac
}
for r in $RULES; do rule_spec "$r" > /dev/null || exit 2; mkdir -p "${OUT}/${r}"; done
WORK="${OUT}/work"; mkdir -p "$WORK"

for CHR in $CHRS; do
    SU="${PHASED_DIR}/merged_chr${CHR}.softunion.bcf"
    PUB_BED="${PHASED_DIR}/merged_chr${CHR}.admix_sites.bed"
    echo "=================== chr${CHR}"
    echo "soft-union: $SU"
    bash "$LISTER" --groups "$CHR" "${WORK}/groups.co.tsv" || exit 1
    # Super-pop-only groups: the same samples, cohort prefix removed (NA kept so
    # it still counts toward the panel).
    awk -F'\t' '{g = $2; sub(/^[^_]*_/, "", g); print $1 "\tSP" g}' "${WORK}/groups.co.tsv" > "${WORK}/groups.sp.tsv"
    echo "  super-pop sizes: $(cut -f2 "${WORK}/groups.sp.tsv" | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"
    echo "  cohort sizes:    $(cut -f2 "${WORK}/groups.co.tsv" | sed 's/_.*//' | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"

    for G in co sp; do
        need=0
        for r in $RULES; do read -r _ g _ _ <<<"$(rule_spec "$r")"; [[ "$g" == "$G" ]] && need=1; done
        [[ "$need" == 1 ]] || continue
        FMT=$(python3 -I "$FILTER" --groups "${WORK}/groups.${G}.tsv" --print-format) || exit 1
        bcftools +fill-tags "$SU" -Ou -- -S "${WORK}/groups.${G}.tsv" -t AN \
            | bcftools query -f "$FMT" | gzip -1 > "${WORK}/an.${G}.tsv.gz" || exit 1
    done

    printf '  %-11s %-40s %12s %7s  %s\n' rule "rule" positions pct "vs 1c's published list ($(wc -l < "$PUB_BED") positions)"
    for r in $RULES; do
        read -r KIND G C1 C2 <<<"$(rule_spec "$r")"
        BED="${OUT}/${r}/merged_chr${CHR}.admix_sites.bed"
        SUM="${OUT}/${r}/merged_chr${CHR}.admix_sites.summary.tsv"
        if [[ "$KIND" == mg ]]; then
            zcat "${WORK}/an.${G}.tsv.gz" | python3 -I "$MARGINS" --groups "${WORK}/groups.${G}.tsv" \
                --max-superpop-miss "$C1" --max-cohort-miss "$C2" --min-n "$MIN_N" \
                --out "$BED" --summary "$SUM" 2> "${OUT}/${r}/chr${CHR}.log" || { cat "${OUT}/${r}/chr${CHR}.log" >&2; exit 1; }
        else
            zcat "${WORK}/an.${G}.tsv.gz" | python3 -I "$FILTER" --groups "${WORK}/groups.${G}.tsv" \
                --max-group-miss "$C1" --max-total-miss "$C2" --min-group-n "$MIN_N" \
                --out "$BED" --summary "$SUM" 2> "${OUT}/${r}/chr${CHR}.log" || { cat "${OUT}/${r}/chr${CHR}.log" >&2; exit 1; }
        fi
        TOTPOS=$(sed -n 's/.* positions=\([0-9]*\) .*/\1/p' "$SUM")
        N=$(wc -l < "$BED")
        BOTH=$(comm -12 <(sort "$BED") <(sort "$PUB_BED") | wc -l)
        CMP="same $BOTH; only here $((N - BOTH)); only in 1c's $(( $(wc -l < "$PUB_BED") - BOTH ))"
        [[ "$r" == coded ]] && { if cmp -s "$BED" "$PUB_BED"; then CMP="IDENTICAL (reconstruction checks out)"; else CMP="DIFFERS: ${CMP} -- the groups were not rebuilt as 1c built them"; fi; }
        printf '  %-11s %-40s %12s %6.1f%%  %s\n' "$r" "$(describe "$r")" "$N" \
            "$(awk -v n="$N" -v t="$TOTPOS" 'BEGIN {print (t ? 100 * n / t : 0)}')" "$CMP"
    done
    echo "  each rule's own report (groups or margins over the limit, and how many sites for no other reason):"
    for r in $RULES; do sed "s/^/    [${r}] /" "${OUT}/${r}/chr${CHR}.log"; done
    rm -f "${WORK}"/an.*.tsv.gz
done
echo "site lists and summaries -> ${OUT}/<rule>/"
