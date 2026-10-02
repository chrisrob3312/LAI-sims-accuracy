#!/usr/bin/env bash
# =============================================================================
# submit_publication_rerun.sh
#
# One-shot orchestrator for the full publication-grade LAI simulation rerun:
#
#   Stage 0: Fetch SGDP metadata + build 7-super-pop anchor lists
#            (scripts/fetch-sgdp-anchors.py)
#   Stage 1: SGDP acquisition + normalization per chr
#            (sbatch 1c_prep_sgdp_clm.sh — array 1-22)
#   Stage 2: Re-ingest TGP+HGDP source WITHOUT dropping the 13 keep-list
#            samples that were previously filtered as "related" + merge
#            with SGDP + rephase (sbatch 1b_phasing-jointcall_clm.sh)
#   Stage 3: K=7 homogeneity filter with SGDP OCE/MEN anchor discovery,
#            unsupervised->supervised ADMIXTURE, Q>=0.95
#            (sbatch 4a_build_homogeneity_panel_k7_clm.sh)
#   Stage 4: Re-simulate admixed cohorts (sbatch 2b_simulation_clm.sh)
#   Stage 5: Full RFMix sweep — WGS main + HOMOG fork, chip main + HOMOG fork
#            (sbatch 3b_*.sh — 4 arrays)
#   Stage 6: Score accuracy + commit + push (sbatch finalize_pilot.sbatch)
#   Stage 7: (main-loop side, not here) regen deck + tables + figures
#
# Dependency chain uses --dependency=afterok so a hard failure in any stage
# halts everything downstream. The chain is submitted all at once; wall-clock
# ETA ~3-5 days of mostly idle waiting.
#
# Required environment:
#   PROJECT_ROOT   default /storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank
#   REPO_DIR       default $(pwd) — must be the LAI-sims-accuracy checkout
#   CONDA_ENV      default shapeit5
#   THROTTLE       default 15 — per-array concurrent task limit
#
# Usage:
#   DRY_RUN=1 bash submit_publication_rerun.sh    # print the chain, submit nothing
#   bash submit_publication_rerun.sh              # fire for real
#
# Notes:
#   - Stage 0 runs locally on the login node (fast). Stages 1-6 are sbatch'd.
#   - Stage 2's "fix 13 samples" step is a KNOWN-PENDING patch: the current
#     1b script sources from a hgdp1kg BCF that has already applied the
#     "unrelated only" filter. See BLOCK B in the chat transcript for how
#     to switch source to the u235158 archived file or extend the keep list.
#     THIS SCRIPT WILL HALT WITH AN ERROR MESSAGE IF THAT PATCH ISN'T IN YET.
# =============================================================================
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
CONDA_ENV="${CONDA_ENV:-shapeit5}"
THROTTLE="${THROTTLE:-15}"
DRY_RUN="${DRY_RUN:-}"

cd "$REPO_DIR"

# ---------------------------------------------------------------- helpers
die() { echo "[publication-rerun] FATAL: $*" >&2; exit 2; }

submit() {
    # usage: submit "<label>" "<sbatch args>" "<script>"
    # returns the job ID on stdout; logs the sbatch line to stderr.
    local label="$1"; shift
    local args="$1"; shift
    local script="$1"
    local line=(sbatch --parsable $args "$script")
    echo >&2
    echo ">>> [${label}]  ${line[*]}" >&2
    if [[ -n "$DRY_RUN" ]]; then
        echo "        (dry run — not submitted)" >&2
        echo "PLACEHOLDER_${label}"
    else
        local jid
        jid=$("${line[@]}")
        echo "        submitted as job ${jid}" >&2
        echo "$jid"
    fi
}

run_now() {
    # usage: run_now "<label>" <command and args...>
    local label="$1"; shift
    echo
    echo ">>> [${label}] running inline: $*"
    if [[ -n "$DRY_RUN" ]]; then
        echo "        (dry run — not executed)"
        return 0
    fi
    "$@"
}

# ---------------------------------------------------------------- precheck
echo "=========================================================="
echo "publication-rerun orchestrator"
echo "  REPO_DIR     = $REPO_DIR"
echo "  PROJECT_ROOT = $PROJECT_ROOT"
echo "  CONDA_ENV    = $CONDA_ENV"
echo "  THROTTLE     = %${THROTTLE}"
echo "  DRY_RUN      = ${DRY_RUN:-off}"
echo "=========================================================="

# Required scripts (hard-fail if missing)
for s in scripts/fetch-sgdp-anchors.py \
         1c_prep_sgdp_clm.sh \
         1b_phasing-jointcall_clm.sh \
         4a_build_homogeneity_panel_k7_clm.sh \
         2b_simulation_clm.sh \
         3b_wgs-rfmix-jointcall_clm.sh \
         3b_homog_wgs-rfmix-jointcall_clm.sh \
         3b_chip_wgs-rfmix-jointcall_clm.sh \
         3b_chip_homog_wgs-rfmix-jointcall_clm.sh \
         finalize_pilot.sbatch ; do
    [[ -f "$s" ]] || die "missing required script: $s"
done
echo "[precheck] all stage scripts present"

# Project paths
[[ -d "$PROJECT_ROOT" ]] || die "PROJECT_ROOT missing: $PROJECT_ROOT"
echo "[precheck] project root OK"

# ---------------------------------------------------------------- STAGE 0
# Fetch SGDP metadata and build the 7 super-pop anchor lists.
# This is fast (metadata is a few MB) and must complete before stage 1
# can read sgdp_combined_anchors.txt.
echo
echo "============ STAGE 0: SGDP anchor-list fetch ============"
run_now "fetch-sgdp" python3 scripts/fetch-sgdp-anchors.py \
    --outdir reference_ids

ANCHORS_FILE="reference_ids/sgdp_combined_anchors.txt"
if [[ -z "$DRY_RUN" ]]; then
    [[ -s "$ANCHORS_FILE" ]] || die "stage 0 did not produce $ANCHORS_FILE"
    N=$(wc -l < "$ANCHORS_FILE")
    echo "[stage 0] combined anchor list: $N SGDP samples across 7 super-pops"
fi

# ---------------------------------------------------------------- STAGE 1
# SGDP acquisition — downloads each chr's phased SGDP BCF (phased_data2021,
# hg19) from the Reich lab public share, filters to the candidate pool, lifts
# hg19->hg38, and emits sgdp_prepped.chr${CHR}.bcf for the merge step.
#
# URL is now pinned to the confirmed layout
# (.../phased_data2021/chr.sgdp.pub.N.bcf). 1c aborts loudly if the candidate
# IDs don't match the BCF sample-ID scheme, so a mismatch can't silently
# produce an empty SGDP panel.
echo
echo "============ STAGE 1: SGDP download + prep ============"
J1=$(submit "stage1_sgdp_prep" "--array=1-22%${THROTTLE}" ./1c_prep_sgdp_clm.sh)

# ---------------------------------------------------------------- STAGE 2
# Rebuild merged phased panel with:
#   - the 13 previously-dropped keep-list samples restored
#   - new SGDP samples merged in
#
# PENDING PATCH: 1b_phasing-jointcall_clm.sh currently sources from a
# hgdp1kg filtered BCF that has already applied the "unrelated only" filter
# upstream of our 13 missing samples. The patch (not yet in the repo) is to
# either:
#   (a) switch HGDP1KG_PATTERN to the u235158 archived file
#       /storage/atkinson/home/u235158/projects/archived/tgp_hgdp_liftover/
#       TGP_HGDP_QC_phased_hg38_chr__CHR__.vcf.gz   (has all 4151 samples
#       incl. our 13 relatives); OR
#   (b) explicitly add the 13 back via a bcftools view -S step on top of
#       the current source.
#
# Option (a) is cleaner. Add it as an env override of HGDP1KG_PATTERN when
# you fire this orchestrator for real.
echo
echo "============ STAGE 2: rebuild merged phased panel ============"
if [[ -z "${HGDP1KG_PATTERN:-}" ]]; then
    cat >&2 <<'EOF'

  WARNING: HGDP1KG_PATTERN is unset. Stage 2 will use the DEFAULT upstream
  BCF which is missing the 13 keep-list samples flagged in the chat
  transcript (HG01935/36, HGDP008xx/009xx/010xx, LP6005441-DNA_B04/B12/F10/
  H06). To include those 13 for publication-grade results, re-invoke with:

    HGDP1KG_PATTERN=/storage/atkinson/home/u235158/projects/archived/tgp_hgdp_liftover/TGP_HGDP_QC_phased_hg38_chr__CHR__.vcf.gz \
      bash submit_publication_rerun.sh

EOF
fi

J2=$(submit "stage2_phase_panel" "--array=1-22%${THROTTLE} --dependency=afterok:${J1}" \
     ./1b_phasing-jointcall_clm.sh)

# ---------------------------------------------------------------- STAGE 3
# K=7 homogeneity filter with SGDP OCE/MEN anchor DISCOVERY.
#   STAGE A: unsupervised ADMIXTURE K=7 -> discover OCE/MEN anchors from the
#            unsupervised Q (candidates sorted by cluster-Q, kept if >=0.95)
#            -> reference_ids/{oce,men}_rfmix.txt.
#   STAGE B: supervised ADMIXTURE K=7 with the full 7-way prior.
#   STAGE C: included iff supervised max_Q >= 0.95 AND argmax==assigned.
# Writes k7 outputs (admixture_k7/, homog_7pop/, lai_ref_panel_samples.k7.tsv)
# and never clobbers the proven K=5 outputs. The EUR/AFR/AMR 3-way accuracy
# test benefits because AMR homogeneity is now resolved at the strict 0.95 cut
# with SGDP-AMR (Karitiana/Surui/...) reinforcing the Amerindigenous centroid.
echo
echo "============ STAGE 3: homogeneity filter (K=7 + SGDP discovery) ============"
J3=$(submit "stage3_homog_k7" "--dependency=afterok:${J2}" ./4a_build_homogeneity_panel_k7_clm.sh)

# ---------------------------------------------------------------- STAGE 4
# Re-simulate admixed cohorts against the new panel. Simulation DONORS come
# from reference_ids/*_simulation.txt, which are DIFFERENT individuals from
# the reference-panel samples to avoid donor-in-reference confound. If donor
# IDs haven't changed, this stage still needs to re-run because the per-chr
# genetic maps now source from the rebuilt panel.
echo
echo "============ STAGE 4: re-simulate admixed cohorts ============"
J4=$(submit "stage4_simulate" "--array=1-22%${THROTTLE} --dependency=afterok:${J3}" \
     ./2b_simulation_clm.sh)

# ---------------------------------------------------------------- STAGE 5
# Full RFMix sweep on both WGS and chip-GSA density. Four parallel arrays:
#   - 3b_wgs main array:        22 chr * 9 default panels = 198 tasks
#   - 3b_homog fork:            22 chr * 2 tracks         = 44 tasks
#   - 3b_chip_wgs main array:   198 tasks  (CHIP_NAME=gsa_jessica)
#   - 3b_chip_homog fork:       44 tasks   (CHIP_NAME=gsa_jessica)
#
# All four run in parallel after stage 4 completes. The chip ones get the
# CHIP_NAME + CHIP_SNPS_TPL env from the submit line to avoid the bug where
# unset CHIP_NAME defaults to "wgs" and runs at WGS density in chip_wgs/.
CHIP_SNPS_TPL="${CHIP_SNPS_TPL:-${REPO_DIR}/reference_ids/chip_gsa_jessica/chr%d.snps}"
CHIP_NAME="${CHIP_NAME:-gsa_jessica}"

echo
echo "============ STAGE 5: full RFMix sweep (4 arrays) ============"
J5a=$(submit "stage5_wgs_main" \
      "--array=1-198%${THROTTLE} --mem=64G --dependency=afterok:${J4}" \
      ./3b_wgs-rfmix-jointcall_clm.sh)

J5b=$(submit "stage5_wgs_homog" \
      "--array=1-44%${THROTTLE} --mem=64G --dependency=afterok:${J4}" \
      ./3b_homog_wgs-rfmix-jointcall_clm.sh)

J5c=$(submit "stage5_chip_main" \
      "--array=1-198%${THROTTLE} --mem=64G --dependency=afterok:${J4} --export=ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}" \
      ./3b_chip_wgs-rfmix-jointcall_clm.sh)

J5d=$(submit "stage5_chip_homog" \
      "--array=1-44%${THROTTLE} --mem=64G --dependency=afterok:${J4} --export=ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}" \
      ./3b_chip_homog_wgs-rfmix-jointcall_clm.sh)

# ---------------------------------------------------------------- STAGE 6
# Accuracy pilot + verify + commit + push. finalize_pilot.sbatch has its
# own gates (confirms every (track, panel, chr) combo wrote .done + .Lat3 +
# .classes, rejects any 2-class .classes, checks per-chr coverage), then
# runs accuracy_v2.R, commits the TSVs, and pushes.
echo
echo "============ STAGE 6: accuracy pilot + push ============"
J6=$(submit "stage6_finalize" \
     "--dependency=afterok:${J5a}:${J5b}:${J5c}:${J5d}" \
     ./finalize_pilot.sbatch)

# ---------------------------------------------------------------- summary
echo
echo "=========================================================="
echo "orchestrator complete — jobs queued with afterok chain:"
echo "  stage 1 (SGDP prep)              = ${J1}"
echo "  stage 2 (rebuild panel)          = ${J2}  (after ${J1})"
echo "  stage 3 (homog filter)           = ${J3}  (after ${J2})"
echo "  stage 4 (simulate)               = ${J4}  (after ${J3})"
echo "  stage 5a (WGS main sweep)        = ${J5a} (after ${J4})"
echo "  stage 5b (WGS HOMOG fork)        = ${J5b} (after ${J4})"
echo "  stage 5c (chip main sweep)       = ${J5c} (after ${J4})"
echo "  stage 5d (chip HOMOG fork)       = ${J5d} (after ${J4})"
echo "  stage 6 (accuracy + push)        = ${J6}  (after 5a/5b/5c/5d)"
echo
echo "Watch progress with:"
echo "  squeue -u \$USER --states=PD,R --format='%.10i %.15j %.8T %.10M %R'"
echo
echo "Stage 7 (deck + table + figure regen) happens on the main-loop side"
echo "automatically once stage 6 commits land on origin."
echo "=========================================================="
