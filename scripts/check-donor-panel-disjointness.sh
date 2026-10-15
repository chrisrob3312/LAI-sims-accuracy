#!/usr/bin/env bash
# ============================================================================
# check-donor-panel-disjointness.sh
#
# Assert that no simulation DONOR is a member of a third-party IMPUTATION
# reference panel. Run this before generating any benchmarking data for the
# genotyping_pipeline project (Project B).
#
# WHY THIS MATTERS, and why it is not the same as the donor/reference check
# already in 2b_simulation_clm.sh:
#   * 2b checks donors against OUR OWN LAI reference panel. We control that
#     panel, so the check is a true set intersection and a hard gate.
#   * This script checks donors against an EXTERNAL imputation panel that we
#     never see. If a donor is inside TOPMed / HRC / 1000G / All of Us, then
#     imputation reconstructs a haplotype the panel already holds, accuracy
#     looks near-perfect, and NO downstream QC step can detect it. The
#     reference panel is server-side; merging and relatedness filtering on our
#     end cannot remove a sample that lives in the panel.
#
# WHAT IS AND IS NOT VERIFIABLE (verified 2026-10):
#   1000 Genomes  -- VERIFIABLE. Public per-sample manifest. Checked here.
#   HRC r1.1      -- NOT verifiable (no public manifest), BUT all of 1000G
#                    Phase 3 is embedded in HRC (McCarthy et al. 2016:
#                    "the 1000 Genomes Phase 3 cohort ... is also included").
#                    So excluding 1KG donors is the necessary condition for
#                    HRC to be a fair test, and this check covers it.
#   TOPMed r2/r3  -- NOT verifiable. No public participant manifest. HGDP is
#                    absent from r2 on strong evidence (it appears in Taliun
#                    et al. 2021 only as an ancestry PCA reference, the
#                    inverse of panel membership); r3 added ~36k samples from
#                    an undocumented freeze and cannot be checked at all.
#   All of Us     -- NOT verifiable, structurally. No participant ID space to
#                    diff against, and the Data User Code of Conduct forbids
#                    the query as re-identification.
#
# Because three of the four cannot be verified, this gate is necessary but NOT
# sufficient. The complement is the donor-vs-holdout probe: build some
# homogeneous individuals from donor haplotypes and some from a disjoint
# same-population holdout. If panel leakage exists, donor-derived individuals
# impute measurably better, and that difference IS the leakage estimate.
#
# Usage:
#   scripts/check-donor-panel-disjointness.sh reference_ids/*_simulation_b2.txt
#   KG_PANEL=/path/to/cached.panel scripts/check-donor-panel-disjointness.sh LIST...
#
# Exit codes: 0 disjoint; 2 overlap found; 3 could not obtain the manifest.
# ============================================================================
set -euo pipefail

KG_URL="${KG_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel}"
CACHE_DIR="${CACHE_DIR:-${TMPDIR:-/tmp}}"
KG_PANEL="${KG_PANEL:-${CACHE_DIR}/integrated_call_samples_v3.20130502.ALL.panel}"

if [[ $# -eq 0 ]]; then
    echo "usage: $0 <donor_list.txt> [more lists...]" >&2
    exit 1
fi

# --- obtain the 1000G manifest -------------------------------------------
if [[ ! -s "$KG_PANEL" ]]; then
    echo "[disjoint] fetching 1000G sample manifest -> $KG_PANEL"
    if ! curl -sSLf --retry 3 --retry-delay 5 -o "$KG_PANEL" "$KG_URL" 2>/dev/null; then
        rm -f "$KG_PANEL"
        cat >&2 <<EOF
ERROR: could not download the 1000G manifest from
         $KG_URL
       Some networks block ftp.1000genomes.ebi.ac.uk. Fetch it by hand on a
       machine that can reach it and pass the path:
         KG_PANEL=/path/to/integrated_call_samples_v3.20130502.ALL.panel $0 ...
       Mirror: https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/
EOF
        exit 3
    fi
else
    echo "[disjoint] using cached manifest $KG_PANEL"
fi

# The manifest is whitespace-separated: sample pop super_pop gender (with a
# header row). Take column 1, skipping the header.
KG_IDS="${CACHE_DIR}/.kg_ids.$$"
trap 'rm -f "$KG_IDS" "$DONORS"' EXIT
awk 'NR>1 && NF {print $1}' "$KG_PANEL" | sort -u > "$KG_IDS"
N_KG=$(wc -l < "$KG_IDS")
if [[ "$N_KG" -lt 2000 ]]; then
    echo "ERROR: manifest parsed to only ${N_KG} ids (expected 2504). Refusing to" >&2
    echo "       certify disjointness against a manifest that did not parse." >&2
    exit 3
fi
echo "[disjoint] 1000G manifest: ${N_KG} sample ids"

# --- collect the donor ids -----------------------------------------------
DONORS="${CACHE_DIR}/.donors.$$"
: > "$DONORS"
for f in "$@"; do
    if [[ ! -s "$f" ]]; then
        echo "[disjoint] WARNING: $f missing or empty -- skipped" >&2
        continue
    fi
    n=$(awk 'NF' "$f" | wc -l)
    echo "[disjoint]   $(basename "$f"): ${n} ids"
    awk 'NF {print $1}' "$f" >> "$DONORS"
done
sort -u -o "$DONORS" "$DONORS"
N_DONOR=$(wc -l < "$DONORS")
echo "[disjoint] donors to check: ${N_DONOR}"
[[ "$N_DONOR" -gt 0 ]] || { echo "ERROR: no donor ids collected"; exit 1; }

# --- the check ------------------------------------------------------------
OVERLAP=$(comm -12 "$DONORS" "$KG_IDS" | wc -l)
if [[ "$OVERLAP" -gt 0 ]]; then
    echo
    echo "FAIL: ${OVERLAP} donor(s) are 1000 Genomes samples."
    echo "      1000G is the Michigan 1000g-phase3-v5 panel AND is embedded"
    echo "      whole inside HRC r1.1, so these donors would be imputed from"
    echo "      their own haplotypes. Offending ids:"
    comm -12 "$DONORS" "$KG_IDS" | head -40 | sed 's/^/        /'
    [[ "$OVERLAP" -gt 40 ]] && echo "        ... and $((OVERLAP - 40)) more"
    echo
    echo "      Re-draw these donors from HGDP, SGDP or MXB, which are absent"
    echo "      from 1000G and HRC."
    exit 2
fi

echo "[disjoint] PASS: 0 of ${N_DONOR} donors appear in the 1000G manifest."
echo "[disjoint] NOTE: TOPMed, HRC and All of Us cannot be checked mechanically."
echo "[disjoint]       This gate is necessary, not sufficient -- pair it with the"
echo "[disjoint]       donor-vs-holdout leakage probe."
