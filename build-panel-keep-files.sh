#!/usr/bin/env bash
# ============================================================================
# build-panel-keep-files.sh
#
# Compose the RFMix-reference and SIMU-donor sample-ID lists for the 4-panel
# comparison from the existing reference_ids/* lists. Run from repo root after
# 1b_phasing-jointcall_clm.sh has produced the joint phased panel.
#
# Output (in $OUTDIR, default panel_keep_files/):
#   --- RFMix reference panels ---
#   panel1_HGDPNAT_IBS_YRI.keep         (HGDP-NAT-rfmix + IBS + YRI)
#   panel4_HGDPMXB_IBS_YRI.keep         (HGDP-NAT-rfmix + MXB-rfmix + IBS + YRI)
#   --- SIMU donor sets ---
#   simu_HGDPNAT_donors.keep            (HGDP-NAT-simu)
#   simu_HGDPMXB_donors.keep            (HGDP-NAT-simu + MXB-simu)
#
# Panels 2 (PEL) and 3 (PEL+EAS) are not built here -- they reuse 1KG IDs
# already pulled by the legacy wgs-simulation-rfmix-jointcall.sh. If you want
# them rebuilt against the new merged panel, point that script's --vcf at
# merged_chr${i}.shapeit5_phased.softunion_maf005.rechr.bcf (no other changes
# needed beyond the new 4th panel block).
# ============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REFS="${REFS:-reference_ids}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/panel_keep_files}"
mkdir -p "$OUTDIR"

# Required inputs (new on this branch + existing)
for f in amr_rfmix.txt mxb_rfmix.txt amr_hgdpmxb_rfmix.txt \
         amr_simulation.txt mxb_simu.txt amr_hgdpmxb_simu.txt \
         eur_rfmix.txt afr_rfmix.txt; do
    [[ -s "${REFS}/${f}" ]] || { echo "ERROR: missing ${REFS}/${f}"; exit 1; }
done

# --- Panel 1 (existing) : HGDP-NAT + IBS + YRI ---
cat "${REFS}/amr_rfmix.txt" "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" \
    > "${OUTDIR}/panel1_HGDPNAT_IBS_YRI.keep"

# --- Panel 4 (new) : HGDP-NAT-rfmix + MXB-rfmix + IBS + YRI ---
cat "${REFS}/amr_hgdpmxb_rfmix.txt" "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" \
    > "${OUTDIR}/panel4_HGDPMXB_IBS_YRI.keep"

# --- SIMU donor lists ---
cp "${REFS}/amr_simulation.txt"        "${OUTDIR}/simu_HGDPNAT_donors.keep"
cp "${REFS}/amr_hgdpmxb_simu.txt"      "${OUTDIR}/simu_HGDPMXB_donors.keep"

# Sanity counts
echo "Panel 1 (HGDP-NAT + IBS + YRI):           $(wc -l < ${OUTDIR}/panel1_HGDPNAT_IBS_YRI.keep) IDs"
echo "Panel 4 (HGDP-NAT + MXB + IBS + YRI):     $(wc -l < ${OUTDIR}/panel4_HGDPMXB_IBS_YRI.keep) IDs"
echo "SIMU donors (HGDP-NAT only):              $(wc -l < ${OUTDIR}/simu_HGDPNAT_donors.keep) IDs"
echo "SIMU donors (HGDP-NAT + MXB):             $(wc -l < ${OUTDIR}/simu_HGDPMXB_donors.keep) IDs"
echo "Outputs in ${OUTDIR}/"
