#!/usr/bin/env bash
# =============================================================================
# submit_publication_rerun.sh
#
# One-shot orchestrator for the full publication-grade LAI simulation rerun:
#
#   Stage 0: Fetch SGDP metadata + build 7-super-pop CANDIDATE lists
#            (scripts/fetch-sgdp-anchors.py). Candidates != anchors: anchors
#            are discovered in stage 3 from these.
#   Stage 1: SGDP acquisition + normalization per chr
#            (sbatch 1b_sgdp_liftover_clm.sh — array 1-22)
#   Stage 2: Re-ingest TGP+HGDP source WITHOUT dropping the 60 keep-list
#            samples that were previously filtered as "related" + merge
#            with SGDP + rephase (sbatch 1c_phasing-jointcall_clm.sh)
#   Stage 3: K=8 homogeneity filter with SGDP anchor discovery,
#            unsupervised->supervised ADMIXTURE, Q>=0.93
#            (sbatch 4a_build_homogeneity_panel_kN_clm.sh)
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

# START_STAGE / STOP_STAGE / AFTER_JOB -- submit part of the chain.
#
# Needed because sbatch SNAPSHOTS the batch script at submit time: a stage
# queued now runs the script as it is now, so a stage whose script still needs
# work must not be queued yet. Start the long stages that are ready, fix the
# rest, then queue the remainder against the running job:
#
#   START_STAGE=0 STOP_STAGE=1 bash submit_publication_rerun.sh   # stages 0-1
#   START_STAGE=2 AFTER_JOB=3738999 bash submit_publication_rerun.sh  # 2-6
#
# AFTER_JOB is the job the first submitted stage waits on (afterok). It is
# required whenever START_STAGE > 1, because the stage it would have depended
# on was not submitted in this invocation.
START_STAGE="${START_STAGE:-0}"
STOP_STAGE="${STOP_STAGE:-6}"
AFTER_JOB="${AFTER_JOB:-}"

cd "$REPO_DIR"

# Should stage N be submitted in this invocation?
want () { [[ "$1" -ge "$START_STAGE" && "$1" -le "$STOP_STAGE" ]]; }

# The dependency clause for stage N: afterok on the job the previous stage
# produced, or on AFTER_JOB when the previous stage was skipped.
dep_of () {
    local prev_jid="$1"
    if [[ -n "$prev_jid" ]]; then echo "--dependency=afterok:${prev_jid}"
    elif [[ -n "$AFTER_JOB" ]]; then echo "--dependency=afterok:${AFTER_JOB}"
    else echo ""; fi
}

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

# Stage 0 runs python on the LOGIN NODE, where `module load anaconda3` and an
# activated conda env coexist. If PYTHONHOME/PYTHONPATH point at the system
# anaconda while the binary on PATH is the conda env's, python loads the
# system stdlib against the env's C extensions and dies at import with
#
#   AssertionError: SRE module mismatch
#
# (_sre.MAGIC from one build, re/_compiler.py from another). 4a already
# sanitises this for the batch stages; stage 0 did not, so the orchestrator
# crashed before submitting anything. Resolve the interpreter explicitly and
# strip both variables for the call.
PY="python3"
[[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python3" ]] && PY="${CONDA_PREFIX}/bin/python3"
py_run () { env -u PYTHONHOME -u PYTHONPATH "$PY" "$@"; }

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
         1b_sgdp_liftover_clm.sh \
         1c_phasing-jointcall_clm.sh \
         4a_build_homogeneity_panel_kN_clm.sh \
         2b_simulation_clm.sh \
         3b_wgs-rfmix-jointcall_clm.sh \
         3b_homog_wgs-rfmix-jointcall_clm.sh \
         3b_chip_wgs-rfmix-jointcall_clm.sh \
         3b_chip_homog_wgs-rfmix-jointcall_clm.sh \
         finalize_pilot.sbatch ; do
    [[ -f "$s" ]] || die "missing required script: $s"
done
echo "[precheck] all stage scripts present"

# A partial submit with no AFTER_JOB would queue stage N with NO dependency,
# so it would start immediately against inputs the skipped stages have not
# produced yet -- and fail hours later, or worse, succeed on stale inputs.
if [[ "$START_STAGE" -gt 1 && -z "$AFTER_JOB" ]]; then
    die "START_STAGE=${START_STAGE} needs AFTER_JOB=<jobid of the stage before it>.
       Without it stage ${START_STAGE} would be queued with no dependency and
       start at once, before its inputs exist. Find the job id with:
         squeue -u \$USER --format='%.12i %.20j %.8T'"
fi
[[ "$START_STAGE" -le "$STOP_STAGE" ]] || die "START_STAGE > STOP_STAGE"
echo "[precheck] submitting stages ${START_STAGE}..${STOP_STAGE}${AFTER_JOB:+ after job ${AFTER_JOB}}"

# Project paths
[[ -d "$PROJECT_ROOT" ]] || die "PROJECT_ROOT missing: $PROJECT_ROOT"
echo "[precheck] project root OK"

# ---------------------------------------------------------------- STAGE 0
# Fetch SGDP metadata and build the 7 super-pop CANDIDATE lists (BCF-native
# IDs). Fast (metadata is a few MB); must complete before stage 1 can read
# sgdp_combined_candidates.txt.
echo
echo "============ STAGE 0: SGDP candidate-list fetch ============"
CAND_FILE="reference_ids/sgdp_combined_candidates.txt"
if want 0; then
    echo "[stage 0] interpreter: ${PY}"
    run_now "fetch-sgdp" py_run scripts/fetch-sgdp-anchors.py \
        --outdir reference_ids
    if [[ -z "$DRY_RUN" ]]; then
        [[ -s "$CAND_FILE" ]] || die "stage 0 did not produce $CAND_FILE"
        N=$(wc -l < "$CAND_FILE")
        echo "[stage 0] combined candidate list: $N SGDP samples across 7 super-pops"
    fi
else
    echo "  (skipped: START_STAGE=${START_STAGE})"
    # Stage 1 reads this file, so it must already exist even when stage 0 is
    # skipped. Checking now beats discovering it 22 array tasks later.
    want 1 && [[ ! -s "$CAND_FILE" ]] && die "stage 0 skipped but $CAND_FILE is missing"
fi

# ---------------------------------------------------------------- STAGE 1
# SGDP acquisition — downloads each chr's phased SGDP BCF (phased_data2021,
# hg19) from the Reich lab public share, filters to the candidate pool, lifts
# hg19->hg38, and emits sgdp_prepped.chr${CHR}.bcf for the merge step.
#
# URL is now pinned to the confirmed layout
# (.../phased_data2021/chr.sgdp.pub.N.bcf). 1b aborts loudly if the candidate
# IDs don't match the BCF sample-ID scheme, so a mismatch can't silently
# produce an empty SGDP panel.
echo
echo "============ STAGE 1: SGDP download + prep ============"
# Stage 1 acquires the two external cohorts. They are independent of each
# other, so both go in now and stage 2 waits on BOTH: afterok takes a
# colon-separated list.
J1=""
if want 1; then
    J1=$(submit "stage1_sgdp_prep" "--array=1-22%${THROTTLE} $(dep_of "")" \
         ./1b_sgdp_liftover_clm.sh)
    if [[ "${RUN_ALMARRI:-1}" == "1" ]]; then
        # The per-sample VCFs ship without indexes and 1d region-queries every
        # one of them, so the index pass runs once, inline, before the array.
        run_now "almarri-index" bash ./1d_almarri_prep_clm.sh --index
        J1B=$(submit "stage1_almarri_prep" "--array=1-22%${THROTTLE} $(dep_of "")" \
              ./1d_almarri_prep_clm.sh)
        J1="${J1}:${J1B}"
    else
        echo "  (RUN_ALMARRI=0 -- Almarri cohort not acquired)"
    fi
else
    echo "  (skipped: START_STAGE=${START_STAGE})"
fi

# ---------------------------------------------------------------- STAGE 2
# Rebuild merged phased panel with:
#   - the 60 previously-dropped keep-list samples restored
#   - new SGDP samples merged in
#
# PENDING PATCH: 1c_phasing-jointcall_clm.sh currently sources from a
# hgdp1kg filtered BCF that has already applied the "unrelated only" filter
# upstream of those samples. To restore them, set HGDP1KG_PATTERN to the
# u235158 archived source:
#
#   /storage/atkinson/home/u235158/projects/archived/tgp_hgdp_liftover/
#   TGP_HGDP_QC_phased_hg38_chr__CHR___rechr.vcf.gz
#
# NOTE THE _rechr SUFFIX -- it is part of the real filename. A pattern
# without it resolves to nothing and the precheck below rejects it.
#
# It restores 60 samples, not 13. One of them, CHMI_CHMI3_WGS2, is a
# hydatidiform-mole cell line (a complete-hydatidiform-mole genome, so
# effectively haploid) and must NOT enter an ADMIXTURE panel or a donor
# pool. Exclude it explicitly when using this source.
echo
echo "============ STAGE 2: rebuild merged phased panel ============"
if [[ -z "${HGDP1KG_PATTERN:-}" ]]; then
    cat >&2 <<'EOF'

  NOTE: HGDP1KG_PATTERN is unset. Stage 2 uses the DEFAULT upstream BCF, which
  has already had an "unrelated only" filter applied upstream and is therefore
  missing 60 keep-list samples. That is a valid run -- the panel is simply the
  relatedness-filtered set, and the methods should say so. To restore them:

    HGDP1KG_PATTERN=/storage/atkinson/home/u235158/projects/archived/\
  tgp_hgdp_liftover/TGP_HGDP_QC_phased_hg38_chr__CHR___rechr.vcf.gz

  (the _rechr suffix is part of the real filename), and exclude
  CHMI_CHMI3_WGS2 -- a hydatidiform-mole cell line, not a population sample.

EOF
else
    # VALIDATE that the override actually resolves. Checking only that the
    # variable is SET lets a wrong path through: stage 1 then runs the entire
    # SGDP download and liftover (hours) before stage 2 dies on a missing
    # file. Expand __CHR__ for one chromosome and look for the file now.
    _probe="${HGDP1KG_PATTERN//__CHR__/22}"
    if [[ ! -s "$_probe" ]]; then
        cat >&2 <<EOF

  FATAL: HGDP1KG_PATTERN is set but does not resolve to a readable file.
         pattern : ${HGDP1KG_PATTERN}
         probed  : ${_probe}

         Firing anyway would waste stage 1 (the full SGDP download and
         hg19->hg38 liftover) and then fail at stage 2.

         Either point it at a real per-chr source containing __CHR__, or
         UNSET it to run with the default relatedness-filtered panel:

           unset HGDP1KG_PATTERN
           bash submit_publication_rerun.sh

EOF
        exit 2
    fi
    echo "[precheck] HGDP1KG_PATTERN resolves: $_probe"
fi

J2=""
if want 2; then
    # Pass the panel-composition settings explicitly rather than relying on
    # sbatch's default environment inheritance: these decide WHICH SAMPLES get
    # jointly phased, and a dropped variable means re-phasing everything.
    # submit() word-splits its args on purpose, so nothing in --export may
    # contain a space. REQUIRE_COHORTS="SGDP ALMARRI" would otherwise put
    # ALMARRI in the argument position where sbatch expects the script.
    # 1c accepts the list comma- or space-separated.
    _e2="ALL"
    [[ -n "${HGDP1KG_PATTERN:-}" ]] && _e2="${_e2},HGDP1KG_PATTERN=${HGDP1KG_PATTERN}"
    [[ -n "${REQUIRE_COHORTS:-}" ]] && _e2="${_e2},REQUIRE_COHORTS=${REQUIRE_COHORTS// /,}"
    [[ -n "${EXCLUDE_IDS:-}"     ]] && _e2="${_e2},EXCLUDE_IDS=${EXCLUDE_IDS}"
    case "$_e2" in
        *" "*) die "a value passed to stage 2 contains a space: ${_e2}" ;;
    esac
    J2=$(submit "stage2_phase_panel" \
         "--array=1-22%${THROTTLE} $(dep_of "$J1") --export=${_e2}" \
         ./1c_phasing-jointcall_clm.sh)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

# ---------------------------------------------------------------- STAGE 3
# K=8 homogeneity filter with SGDP anchor DISCOVERY.
#   STAGE A: unsupervised ADMIXTURE K=8 -> discover anchors from the
#            unsupervised Q (candidates sorted by cluster-Q, kept if
#            >= ANCHOR_MIN_Q = 0.93, capped at ANCHOR_CAP = 50).
#   STAGE B: supervised ADMIXTURE K=8 with the full 8-way prior.
#   STAGE C: included iff supervised max_Q >= HOMOG_THRESHOLD (0.93) AND the
#            supervised argmax cluster == the assigned super-pop.
# Writes k8 outputs (admixture_k8/, homog_8pop/, lai_ref_panel_samples.k8.tsv)
# and never clobbers the proven K=5 outputs.
#
# KNOWN GAP as of this revision: stage A discovers anchors for OCE, MEN and
# OTH only; EUR/AFR/AMR/EAS/SAS still use the pre-existing *_rfmix.txt
# exemplar lists rather than being discovered the same way. Uniform discovery
# across all 7 candidate pools is the agreed design and is not in yet.
#
# THIS MATTERS AT SUBMIT TIME, NOT AT RUN TIME. sbatch takes a COPY of the
# batch script when the job is queued, and the job runs that copy -- editing
# 4a_build_homogeneity_panel_kN_clm.sh after this orchestrator has queued
# stage 3 does NOT change what stage 3 runs. Verify what a queued job will
# actually execute with:
#     scontrol write batch_script <jobid> -
# So either fix 4a before submitting stage 3, or submit stage 1 on its own
# (START_STAGE=1) and submit stages 2-6 later with AFTER_JOB=<stage1 jobid>.
echo
echo "============ STAGE 3: homogeneity filter (K=8 + SGDP discovery) ============"
J3=""
if want 3; then
    J3=$(submit "stage3_homog_kN" "$(dep_of "$J2")" \
         ./4a_build_homogeneity_panel_kN_clm.sh)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

# ---------------------------------------------------------------- STAGE 4
# Re-simulate admixed cohorts against the new panel. Simulation DONORS come
# from reference_ids/*_simulation.txt, which are DIFFERENT individuals from
# the reference-panel samples to avoid donor-in-reference confound. If donor
# IDs haven't changed, this stage still needs to re-run because the per-chr
# genetic maps now source from the rebuilt panel.
echo
echo "============ STAGE 4: re-simulate admixed cohorts ============"
J4=""
if want 4; then
    J4=$(submit "stage4_simulate" "--array=1-22%${THROTTLE} $(dep_of "$J3")" \
         ./2b_simulation_clm.sh)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

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
J5a=""; J5b=""; J5c=""; J5d=""
if want 5; then
    _d5="$(dep_of "$J4")"
    J5a=$(submit "stage5_wgs_main" \
          "--array=1-198%${THROTTLE} --mem=64G ${_d5}" \
          ./3b_wgs-rfmix-jointcall_clm.sh)

    J5b=$(submit "stage5_wgs_homog" \
          "--array=1-44%${THROTTLE} --mem=64G ${_d5}" \
          ./3b_homog_wgs-rfmix-jointcall_clm.sh)

    J5c=$(submit "stage5_chip_main" \
          "--array=1-198%${THROTTLE} --mem=64G ${_d5} --export=ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}" \
          ./3b_chip_wgs-rfmix-jointcall_clm.sh)

    J5d=$(submit "stage5_chip_homog" \
          "--array=1-44%${THROTTLE} --mem=64G ${_d5} --export=ALL,CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}" \
          ./3b_chip_homog_wgs-rfmix-jointcall_clm.sh)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

# ---------------------------------------------------------------- STAGE 6
# Accuracy pilot + verify + commit + push. finalize_pilot.sbatch has its
# own gates (confirms every (track, panel, chr) combo wrote .done + .Lat3 +
# .classes, rejects any 2-class .classes, checks per-chr coverage), then
# runs accuracy_v2.R, commits the TSVs, and pushes.
echo
echo "============ STAGE 6: accuracy pilot + push ============"
J6=""
if want 6; then
    if [[ -n "$J5a" ]]; then
        _d6="--dependency=afterok:${J5a}:${J5b}:${J5c}:${J5d}"
    else
        _d6="$(dep_of "")"
    fi
    J6=$(submit "stage6_finalize" "${_d6}" ./finalize_pilot.sbatch)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

# ---------------------------------------------------------------- summary
echo
echo "=========================================================="
echo "orchestrator complete — jobs queued with afterok chain:"
row () { [[ -n "$2" ]] && printf '  %-32s = %s\n' "$1" "$2" || printf '  %-32s   (not submitted)\n' "$1"; }
row "stage 1 (SGDP prep)"        "$J1"
row "stage 2 (rebuild panel)"    "$J2"
row "stage 3 (homog filter)"     "$J3"
row "stage 4 (simulate)"         "$J4"
row "stage 5a (WGS main sweep)"  "$J5a"
row "stage 5b (WGS HOMOG fork)"  "$J5b"
row "stage 5c (chip main sweep)" "$J5c"
row "stage 5d (chip HOMOG fork)" "$J5d"
row "stage 6 (accuracy + push)"  "$J6"
if [[ "$STOP_STAGE" -lt 6 ]]; then
    echo
    echo "PARTIAL SUBMIT. Queue the rest once the stage scripts are ready:"
    echo "  START_STAGE=$((STOP_STAGE+1)) AFTER_JOB=<last job id above> \\"
    echo "    bash submit_publication_rerun.sh"
fi
echo
echo "Watch progress with:"
echo "  squeue -u \$USER --states=PD,R --format='%.10i %.15j %.8T %.10M %R'"
echo
echo "Stage 7 (deck + table + figure regen) happens on the main-loop side"
echo "automatically once stage 6 commits land on origin."
echo "=========================================================="
