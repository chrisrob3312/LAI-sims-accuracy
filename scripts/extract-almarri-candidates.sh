#!/usr/bin/env bash
# ============================================================================
# extract-almarri-candidates.sh
#
# Emit per-super-pop CANDIDATE lists from the HGDP-1kGP-Almarri .fam
# (AlAli / Almarri / Atkinson, "An Expanded Reference for Middle Eastern
# Populations Improves GWAS"), so they join the candidate pools that
# 4a_build_homogeneity_panel_kN_clm.sh unions by glob.
#
# THESE ARE CANDIDATES, NOT ANCHORS. Being in this list confers nothing. The
# pool does two jobs -- it names the unsupervised clusters by mean-Q, and it is
# the draw pool for anchors -- but which samples become ANCHORS is decided by
# the top-N >= ANCHOR_MIN_Q cut in stage A, and which enter the HOMOGENEOUS
# panel by the supervised Q >= HOMOG_THRESHOLD cut in stage C. Nothing here
# pre-commits either.
#
# Why this source matters: MEN currently rests on SGDP's 19 candidates, which
# is thin for a mean-Q cluster assignment. The expanded reference adds 136 WGS
# samples from seven previously unrepresented regions (Syria, Jordan, Iraq,
# Saudi Arabia, UAE, Yemen, Oman), and its .fam already carries harmonized
# HGDP-1kGP super-pop labels in column 1.
#
# NOTE ON THE GENOTYPES: the .fam accompanies an LD-pruned, MAF>=0.05,
# MHC-excluded PCA panel (833,878 variants of 7,761,449 intersecting). That is
# the right input for ADMIXTURE but NOT a source of reference haplotypes for
# RFMix/FLARE -- those need the unpruned, phased WGS. This script only takes
# sample IDs, so it is unaffected; the WGS ingest is a separate stage.
#
# Usage:
#   scripts/extract-almarri-candidates.sh /storage/atkinson/shared_projects/mid/almarri_hgdp1kgpqc.fam
#   OUTDIR=reference_ids scripts/extract-almarri-candidates.sh <fam>
# ============================================================================
set -euo pipefail

FAM="${1:-}"
OUTDIR="${OUTDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/reference_ids}"
TAG="${TAG:-almarri}"

if [[ -z "$FAM" || ! -s "$FAM" ]]; then
    echo "usage: $0 <almarri_hgdp1kgpqc.fam>" >&2
    exit 1
fi
mkdir -p "$OUTDIR"

# Their super-pop label -> our candidate-file infix. CSA maps to sas because
# 4a's super_pop_of() already collapses CSA -> SAS; revisit here if CSA is ever
# split into South and Central Asian.
declare -A MAP=( [AFR]=afr [AMR]=amr [EAS]=eas [EUR]=eur [MID]=men [OCE]=oceania [CSA]=sas )

echo "[almarri] source: $FAM ($(wc -l < "$FAM") samples)"
echo "[almarri] super-pop labels present:"
awk '{print $1}' "$FAM" | sort | uniq -c | sort -rn | sed 's/^/    /'
echo

MANIFEST="${OUTDIR}/${TAG}_manifest.tsv"
printf 'sample\tpopulation\tsuper_pop\n' > "$MANIFEST"

for sp in "${!MAP[@]}"; do
    infix="${MAP[$sp]}"
    out="${OUTDIR}/${TAG}_${infix}_candidates.txt"
    awk -v s="$sp" '$1==s {print $2}' "$FAM" | sort -u > "$out"
    n=$(wc -l < "$out")
    if [[ "$n" -eq 0 ]]; then
        rm -f "$out"
        echo "[almarri] ${sp}: 0 samples -- no file written"
        continue
    fi
    # population is unknown from the .fam alone (col1 is the super-pop), so it
    # is recorded as the super-pop. Replace with real population labels from
    # the paper's Table S5 if finer resolution is wanted for the [refine] pass.
    awk -v s="$sp" -v OFS='\t' '$1==s {print $2, s, s}' "$FAM" | sort -u >> "$MANIFEST"
    echo "[almarri] ${sp} -> ${n} candidates -> $(basename "$out")"
done

echo
echo "[almarri] manifest -> $MANIFEST ($(($(wc -l < "$MANIFEST") - 1)) rows)"
echo "[almarri] 4a picks these up automatically via its *_men_candidates.txt /"
echo "[almarri] *_oceania_candidates.txt globs. They remain CANDIDATES until"
echo "[almarri] stage A's top-N >= ANCHOR_MIN_Q cut selects anchors."
echo
echo "[almarri] REMINDER: these IDs must match the sample names in the merged"
echo "[almarri] panel. The Almarri WGS is not ingested yet, so until it is,"
echo "[almarri] every one of these will be absent from the .fam and will"
echo "[almarri] contribute nothing -- which is harmless but also pointless."
echo "[almarri] Ingest the WGS first, then re-run 4a."
