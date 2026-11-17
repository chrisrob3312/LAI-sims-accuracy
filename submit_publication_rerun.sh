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
#   Stage 3: two homogeneity filters side by side (sbatch
#            4a_build_homogeneity_panel_kN_clm.sh, HOMOG_SET=global / lai3way):
#              global  K=7, every sample, central high-Q anchors (Q >= 0.94),
#                      homogeneous at Q >= 0.94 -- for the benchmark
#              lai3way K=5 as July (EUR AFR AMR EAS SAS, July's anchors, MEN/
#                      OCE/no-super-pop left out), Q >= 0.94 -- writes the
#                      3-way HOMOG arms' amr_homog_rfmix.txt
#            Both: every super-pop <= 7% missing, every cohort <= 15%; LD prune
#            on the unrelated set; relatives projected.
#   Stages 4-6, once per cohort in COHORTS (default Brasa:12:NAT
#            MexLike:16:NATMXB -- Brazilian-like and Mexican-like):
#   Stage 4: simulate the cohort (sbatch 2b_simulation_clm.sh)
#   Stage 5: RFMix + FLARE, WGS + GSA chip, main panels + HOMOG fork
#            (sbatch 3b_*.sh -- 8 arrays)
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
        return 0
    fi

    # Retry, then DIE. This cluster enforces a QOS submit limit, and a chain
    # of array jobs can hit it partway through: sbatch prints "Slurm
    # temporarily unable to accept job, sleeping and retrying", exhausts its
    # own retries and returns nothing.
    #
    # Returning an empty job id here is worse than failing. The caller stores
    # it, the next stage's --dependency is built from it, and the chain
    # continues around the hole -- which is how a FLARE arm vanished from a
    # submission while every other stage queued normally and the summary
    # looked complete.
    local jid="" attempt out rc errf
    for attempt in 1 2 3 4 5; do
        # The id is on STDOUT ("id" or "id;cluster"). This cluster's job_submit
        # plugin prints an informational "sbatch: slurm_job_submit: ..." line on
        # STDERR for every job: shown here, never parsed. A successful sbatch is
        # never retried -- that is how a stage gets queued twice.
        errf=$(mktemp)
        out=$("${line[@]}" 2> "$errf"); rc=$?
        [[ -s "$errf" ]] && sed 's/^/        /' "$errf" >&2
        rm -f "$errf"
        jid=$(printf '%s\n' "$out" | awk -F';' '$1 ~ /^[0-9]+$/ {id = $1} END {print id}')
        [[ $rc -eq 0 && -n "$jid" ]] && break
        if [[ $rc -eq 0 ]]; then
            die "sbatch exited 0 for [${label}] but printed no job id ('${out}'). A job MAY have been
       queued: check squeue before re-running anything."
        fi
        jid=""
        local wait=$((attempt * 20))
        echo "        sbatch failed (exit ${rc}, attempt ${attempt}/5); retrying in ${wait}s" >&2
        sleep "$wait"
    done
    if [[ -z "$jid" ]]; then
        die "could not submit [${label}] after 5 attempts.

       sbatch returned no job id. The usual cause here is the QOS submit
       limit: too many array jobs queued at once. Jobs already submitted in
       this run are LEFT QUEUED -- note their ids from the lines above, wait
       for some to finish, and resume with
         START_STAGE=<this stage> AFTER_JOB=<the id this stage depends on>
       rather than re-running the whole chain, which would duplicate them."
    fi
    echo "        submitted as job ${jid}" >&2
    echo "$jid"
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
# AFTER_JOB=none means "the inputs are already on disk, start now". It exists
# because the normal case after a staged recovery is that every upstream array
# has ALREADY COMPLETED and there is no job left to depend on -- at which point
# demanding a job id would force inventing one, or waiting on an unrelated job.
# It is spelled out rather than silently allowed, because a stage queued with
# no dependency starts immediately, and if the inputs are not actually there it
# fails 22 times over.
if [[ "$AFTER_JOB" == "none" ]]; then
    AFTER_JOB=""
    echo "[precheck] AFTER_JOB=none -- no dependency; inputs must already exist"
    if want 2; then
        _missing=0
        for _c in $(seq 1 22); do
            for _f in "${PROJECT_ROOT}/01_merged_phased_panel/mxb_lifted.chr${_c}.bcf" \
                      "${PROJECT_ROOT}/01_merged_phased_panel/sgdp_prepped.chr${_c}.bcf"; do
                [[ -s "$_f" ]] || { echo "  MISSING: $_f" >&2; _missing=$((_missing+1)); }
            done
            _a="${PROJECT_ROOT}/01_merged_phased_panel/almarri_prepped.chr${_c}.bcf"
            if [[ -n "${REQUIRE_COHORTS:-}" && "$REQUIRE_COHORTS" == *ALMARRI* ]]; then
                [[ -s "$_a" ]] || { echo "  MISSING: $_a" >&2; _missing=$((_missing+1)); }
            fi
        done
        [[ "$_missing" -eq 0 ]] || die "AFTER_JOB=none but ${_missing} stage-2 input(s) are absent.
       Stage 2 would start at once and fail on every chromosome.
       Wait for the upstream arrays, or pass their job id as AFTER_JOB."
        echo "[precheck] all stage-2 inputs present on disk"
    fi
elif [[ "$START_STAGE" -gt 1 && -z "$AFTER_JOB" ]]; then
    die "START_STAGE=${START_STAGE} needs AFTER_JOB=<jobid of the stage before it>,
       or AFTER_JOB=none if every upstream array has already completed and its
       outputs are on disk (the inputs are then checked before anything is
       queued). Without either, stage ${START_STAGE} would be queued with no
       dependency and start before its inputs exist. Find a job id with:
         squeue -u \$USER --format='%.12i %.20j %.8T'"
fi
[[ "$START_STAGE" -le "$STOP_STAGE" ]] || die "START_STAGE > STOP_STAGE"
echo "[precheck] submitting stages ${START_STAGE}..${STOP_STAGE}${AFTER_JOB:+ after job ${AFTER_JOB}}"

# Re-running this while an earlier submission of the SAME stages is still
# queued would put two copies of each on the same files -- two 4a jobs after
# one 1e, two 2b arrays writing one simulation tree. Refuse (a dry run only
# warns). Only the stages this run would submit are checked, so the job named
# in AFTER_JOB -- e.g. a queued 1e that 4a is to wait on -- is not flagged.
if command -v squeue > /dev/null; then
    _names=" "
    want 2 && _names+="1c_phase "
    if want 3 && { want 2 || [[ "${RUN_KING:-}" == "1" ]]; }; then _names+="1e_king_dedup "; fi
    want 3 && _names+="4a_homog_kN 4a_homog_lai3way "
    want 4 && _names+="2b_simu "
    want 5 && _names+="3b_rfmix 3b_homog 3b_chip 3b_chip_homog 3b_flare "
    want 6 && _names+="finalize_pilot "
    _queued=$(squeue -h -u "${USER:-$(id -un)}" -o '%A %j %T' | awk -v n="$_names" 'index(n, " " $2 " ")' | sort -u)
    if [[ -n "$_queued" ]]; then
        echo "[precheck] stages this run would submit are still queued from an earlier submission:" >&2
        sed 's/^/            /' <<<"$_queued" >&2
        if [[ -n "$DRY_RUN" ]]; then
            echo "[precheck] DRY RUN: a real run REFUSES until they are cancelled" >&2
        else
            die "cancel them first (scancel <id> ...), then re-run"
        fi
    fi
fi

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
# Rebuild the merged phased panel. The panel rebuild is normally submitted by
# submit_panel_rebuild.sh, which chains stages 4-6 of this script after its 4a
# job (START_STAGE=4), so this stage is skipped. HGDP1KG_PATTERN defaults (in
# 1c and in that driver) to the UNFILTERED Koenig et al. phased_haplotypes_v2
# release; the u235158 archive's 60 extra samples are Koenig's QC removals and
# are excluded by name (reference_ids/panel_exclude_ids.txt, incl. CHMI).
echo
echo "============ STAGE 2: rebuild merged phased panel ============"
if [[ -z "${HGDP1KG_PATTERN:-}" ]]; then
    if want 2; then
        echo "  HGDP1KG_PATTERN unset: 1c uses its default, the unfiltered Koenig v2 release." >&2
    fi
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

J2=""; J2B=""
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
    # PHASE_MODE=chunked makes stage 2 stop after the soft-union and hand the
    # phasing to 1b_phase_chunks + 1b_phase_ligate. Stage 3 must then depend on
    # LIGATE, not on this job, so the chain below is not the right shape --
    # submit stages 3+ separately with AFTER_JOB set to the ligate id.
    [[ -n "${PHASE_MODE:-}"      ]] && _e2="${_e2},PHASE_MODE=${PHASE_MODE}"
    if [[ "${PHASE_MODE:-whole}" == "chunked" ]] && want 3; then
        die "PHASE_MODE=chunked with stage 3 in range.

       Stage 2 will only publish the soft-union; the phasing is then
       1b_phase_chunks and 1b_phase_ligate, and stage 3 has to wait for
       LIGATE rather than for stage 2. Submit stage 2 alone:
         PHASE_MODE=chunked START_STAGE=2 STOP_STAGE=2 ... bash $0
       then the chunk and ligate arrays, then
         START_STAGE=3 RUN_KING=1 AFTER_JOB=<ligate id> bash $0"
    fi
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
# Homogeneity filter: TWO 4a jobs, run side by side (17 Nov), one per project.
#   global  (HOMOG_SET=global, K=7) -- every panel sample; 4a's own ADMIXTURE
#           site list (every super-pop <= 7% missing, every cohort <= 15%), LD
#           pruned on the unrelated fit set; unsupervised -> central high-Q
#           anchors (Q >= 0.94) -> supervised; homogeneous at Q >= 0.94.
#           For the genotype-harmonization benchmark (run_benchmark_cohort.sh).
#   lai3way (HOMOG_SET=lai3way, K=5) -- the 3-way LAI simulation's AMR panel,
#           judged as July judged it: EUR AFR AMR EAS SAS with July's anchors,
#           MEN / OCE / no-super-pop samples left out of the run, Q >= 0.94.
#           Publishes reference_ids/amr_homog_rfmix.txt (July's AMR panel + the
#           AMR samples new since July that pass) for the HOMOG arms.
# Stages 4-6 wait on lai3way only; the benchmark waits on global. They share
# 1e's KING table and one concatenated input (built under a lock).
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
echo "============ STAGE 3: homogeneity filter (4a: global K=7 + lai3way K=5) ============"
# ----------------------------------------------------------------- STAGE 2B
# KING de-duplication, BETWEEN the merge and ADMIXTURE.
#
# It must precede stage 3 and cannot follow it: a person present twice is
# counted twice by ADMIXTURE and pulls their population's centroid twice as
# hard, so the Q matrix that the homogeneity cut and the donor/reference split
# both read is wrong before either is computed. Removing the duplicate
# afterwards means redoing both.
# Submitted only when STAGE 2 is in this invocation. Keying it to stage 3
# instead meant that resuming at stage 3 -- the normal case when stage 2 is
# already running and KING is already queued behind it -- submitted a SECOND
# KING that re-did hours of concat and plink work, and chained stage 3 to the
# duplicate rather than to the original. RUN_KING=1 forces it when stage 2 is
# out of range and KING genuinely has not run.
J2B=""
if want 3 && { want 2 || [[ "${RUN_KING:-}" == "1" ]]; }; then
    echo
    echo "============ STAGE 2B: KING de-duplication ============"
    J2B=$(submit "stage2b_king_dedup" "$(dep_of "$J2")" ./1e_king_dedup_clm.sh)
fi

J3=""; J3L=""
if want 3; then
    # Both wait on KING when it ran, so 4a reads a finished drop list rather
    # than whatever happens to be on disk.
    J3=$(submit "stage3_homog_global_k7" "$(dep_of "${J2B:-$J2}") --export=ALL,HOMOG_SET=global" \
         ./4a_build_homogeneity_panel_kN_clm.sh)
    J3L=$(submit "stage3_homog_lai3way_k5" \
          "$(dep_of "${J2B:-$J2}") --job-name=4a_homog_lai3way --export=ALL,HOMOG_SET=lai3way" \
          ./4a_build_homogeneity_panel_kN_clm.sh)
else
    echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
fi

# ------------------------------------------------------------------ STAGE 3B
# The benchmark donor/reference split, IN PARALLEL with stages 4-6.
#
# It depends on the GLOBAL 4a and on nothing after it: the split is a function
# of the K=7 panel alone. Chaining it behind the RFMix sweep and the accuracy
# tables would hold the genotyping-pipeline benchmark hostage to a multi-day
# run it does not need, so it branches off stage 3 beside stage 4.
if want 4 && [[ "${RUN_SPLIT:-1}" == "1" ]]; then
    echo
    echo "============ STAGE 3B: benchmark donor/reference split ============"
    if [[ -n "$J3" ]]; then
        echo "  the global 4a is queued as ${J3}; the benchmark cohort needs its output, so"
        echo "  once it completes run (the split takes seconds; 6b-6e are arrays):"
        echo
        echo "    ./run_benchmark_cohort.sh --dry-run"
        echo "    ./run_benchmark_cohort.sh"
        echo
        echo "  Nothing downstream in THIS chain waits on it."
    fi
fi

# ------------------------------------------------------------- STAGES 4-6
# One simulated cohort per entry of COHORTS, "<ADMIX_POP>:<GEN>:<tracks>",
# tracks joined with '+' (sbatch --export splits on commas):
#   Brasa:12:NAT        Brazilian-like  NAT .15 / EUR .60 / AFR .25, HGDP-NAT donors
#   MexLike:16:NATMXB   Mexican-like    NAT .50 / EUR .45 / AFR .05, HGDP-NAT + MXB donors
# the pairing of the August pilot, whose two cohorts both ran off Brasa.dat
# (so its Mexican-like cohort had Brazilian proportions and 12 generations).
# Each cohort: a 2b array, the eight RFMix/FLARE arrays, a finalize. The
# finalize jobs run one after another (afterany on the previous one), since
# each commits and pushes from the same checkout.
COHORTS="${COHORTS:-Brasa:12:NAT MexLike:16:NATMXB}"
CHIP_SNPS_TPL="${CHIP_SNPS_TPL:-${REPO_DIR}/reference_ids/chip_gsa_jessica/chr%d.snps}"
CHIP_NAME="${CHIP_NAME:-gsa_jessica}"
# The main arms' panels, as the 3b scripts build their default list.
MAIN_PANELS="NAT_HGDP"
[[ -s reference_ids/pel_rfmix.txt     ]] && MAIN_PANELS="${MAIN_PANELS} NAT_PEL"
[[ -s reference_ids/pel_eas_rfmix.txt ]] && MAIN_PANELS="${MAIN_PANELS} NAT_PEL_EAS"
MAIN_PANELS="${MAIN_PANELS} NAT_HGDPMXB NAT_HGDPMXB_FULL"
declare -A C_GEN C_TRACKS C_MAIN C_HOMOG
COHORT_LIST=()
for _c in $COHORTS; do
    IFS=':' read -r _pop _gen _tr <<< "$_c"
    [[ -n "$_pop" && "$_gen" =~ ^[0-9]+$ && -n "$_tr" ]] || die "COHORTS entry '${_c}' is not <ADMIX_POP>:<GEN>:<tracks>"
    _nmain=0; _nhom=0
    for _t in ${_tr//+/ }; do
        case "$_t" in NAT|NATMXB) ;; *) die "COHORTS entry '${_c}': track '${_t}' is not NAT or NATMXB" ;; esac
        _nhom=$((_nhom + 1))
        for _p in $MAIN_PANELS; do
            [[ "$_t" == NATMXB && "$_p" == NAT_HGDPMXB_FULL ]] && continue   # donors would be in the reference
            _nmain=$((_nmain + 1))
        done
    done
    _dat="${PROJECT_ROOT}/02_simulations/${_pop}.dat"
    [[ -s "$_dat" ]] || _dat="${REPO_DIR}/simulation_models/${_pop}.dat"
    [[ -s "$_dat" ]] || die "COHORTS entry '${_c}': no ${PROJECT_ROOT}/02_simulations/${_pop}.dat and no simulation_models/${_pop}.dat"
    C_GEN[$_pop]="$_gen"; C_TRACKS[$_pop]="$_tr"
    C_MAIN[$_pop]=$((22 * _nmain)); C_HOMOG[$_pop]=$((22 * _nhom))
    COHORT_LIST+=("$_pop")
    echo "[precheck] cohort ${_pop}: gen ${_gen}, tracks ${_tr//+/ }; model ${_dat}: $(tr '\n' '|' < "$_dat")"
    echo "[precheck]   arrays: main ${C_MAIN[$_pop]} tasks (${_nmain} track x panel combos x 22), HOMOG ${C_HOMOG[$_pop]}"
done

declare -A J4 J5 J6
_prev6=""
for _pop in "${COHORT_LIST[@]}"; do
    _gen="${C_GEN[$_pop]}"; _tr="${C_TRACKS[$_pop]}"
    _ex="ADMIX_POP=${_pop},GEN=${_gen},TRACKS_TO_RUN=${_tr}"
    _sfx=""; [[ "$_pop" == Brasa ]] || _sfx="_${_pop}"

    # ------------------------------------------------------------ STAGE 4
    # Re-simulate the cohort against the new panel. Simulation DONORS come
    # from reference_ids/*_simulation.txt, which are DIFFERENT individuals from
    # the reference-panel samples to avoid donor-in-reference confound.
    echo
    echo "============ STAGE 4: simulate ${_pop} (gen ${_gen}, ${_tr//+/ }) ============"
    J4[$_pop]=""
    if want 4; then
        J4[$_pop]=$(submit "stage4_simulate_${_pop}" "--array=1-22%${THROTTLE} $(dep_of "$J3L") --export=ALL,${_ex}" \
                    ./2b_simulation_clm.sh)
    else
        echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
    fi

    # ------------------------------------------------------------ STAGE 5
    # RFMix + FLARE on both WGS and chip-GSA density, the same grid for both
    # methods: main panels and the HOMOG fork, at each density. The chip ones
    # get CHIP_NAME + CHIP_SNPS_TPL on the submit line, so an unset CHIP_NAME
    # cannot default to "wgs" and run at WGS density in chip_wgs/.
    echo
    echo "============ STAGE 5: RFMix + FLARE sweep for ${_pop} (8 arrays) ============"
    J5[$_pop]=""
    if want 5; then
        _d5="$(dep_of "${J4[$_pop]}")"
        _m="--array=1-${C_MAIN[$_pop]}%${THROTTLE}"; _h="--array=1-${C_HOMOG[$_pop]}%${THROTTLE}"
        _chip="CHIP_NAME=${CHIP_NAME},CHIP_SNPS_TPL=${CHIP_SNPS_TPL}"
        _j=$(submit "stage5_wgs_main_${_pop}"   "${_m} --mem=64G ${_d5} --export=ALL,${_ex}" ./3b_wgs-rfmix-jointcall_clm.sh)
        _js="$_j"
        _j=$(submit "stage5_wgs_homog_${_pop}"  "${_h} --mem=64G ${_d5} --export=ALL,${_ex}" ./3b_homog_wgs-rfmix-jointcall_clm.sh)
        _js="${_js}:${_j}"
        _j=$(submit "stage5_chip_main_${_pop}"  "${_m} --mem=64G ${_d5} --export=ALL,${_ex},${_chip}" ./3b_chip_wgs-rfmix-jointcall_clm.sh)
        _js="${_js}:${_j}"
        _j=$(submit "stage5_chip_homog_${_pop}" "${_h} --mem=64G ${_d5} --export=ALL,${_ex},${_chip}" ./3b_chip_homog_wgs-rfmix-jointcall_clm.sh)
        _js="${_js}:${_j}"
        # FLARE on the SAME grid as the four RFMix arms, so the only thing
        # differing between the two sweeps is the method. Off with RUN_FLARE=0
        # if the jar is not installed.
        if [[ "${RUN_FLARE:-1}" == "1" && -f ./3b_flare_wgs-jointcall_clm.sh ]]; then
            _j=$(submit "stage5_flare_wgs_${_pop}"        "${_m} --mem=32G ${_d5} --export=ALL,${_ex}" ./3b_flare_wgs-jointcall_clm.sh)
            _js="${_js}:${_j}"
            _j=$(submit "stage5_flare_wgs_homog_${_pop}"  "${_h} --mem=32G ${_d5} --export=ALL,${_ex},PANELS_TO_RUN=NAT_HOMOG" ./3b_flare_wgs-jointcall_clm.sh)
            _js="${_js}:${_j}"
            _j=$(submit "stage5_flare_chip_${_pop}"       "${_m} --mem=32G ${_d5} --export=ALL,${_ex},${_chip}" ./3b_flare_wgs-jointcall_clm.sh)
            _js="${_js}:${_j}"
            _j=$(submit "stage5_flare_chip_homog_${_pop}" "${_h} --mem=32G ${_d5} --export=ALL,${_ex},${_chip},PANELS_TO_RUN=NAT_HOMOG" ./3b_flare_wgs-jointcall_clm.sh)
            _js="${_js}:${_j}"
        else
            echo "  (RUN_FLARE=0 or script absent -- FLARE arms not submitted)"
        fi
        J5[$_pop]="$_js"
    else
        echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
    fi

    # ------------------------------------------------------------ STAGE 6
    # Accuracy + verify + commit + push. finalize_pilot.sbatch has its own
    # gates (every (track, panel, chr) combo wrote .done + .Lat3 + .classes,
    # no 2-class .classes, per-chr coverage), then runs accuracy_v2.R and
    # pushes results/accuracy_pilot*${_sfx}/. FLARE is scored when it was
    # submitted (SCORE_FLARE), so the tables carry RFMix and FLARE side by side.
    echo
    echo "============ STAGE 6: accuracy + push for ${_pop} ============"
    J6[$_pop]=""
    if want 6; then
        if [[ -n "${J5[$_pop]}" ]]; then _d6="--dependency=afterok:${J5[$_pop]}"
        else _d6="$(dep_of "")"; fi
        if [[ -n "$_prev6" ]]; then
            if [[ "$_d6" == --dependency=* ]]; then _d6="${_d6},afterany:${_prev6}"
            else _d6="--dependency=afterany:${_prev6}"; fi
        fi
        _sf=0
        if [[ "${RUN_FLARE:-1}" == "1" ]] && { [[ -z "${J5[$_pop]}" ]] || [[ "${J5[$_pop]}" == *:*:*:*:* ]]; }; then _sf=1; fi
        J6[$_pop]=$(submit "stage6_finalize_${_pop}" "${_d6} --export=ALL,SCORE_FLARE=${_sf},CHIP_NAME=${CHIP_NAME},${_ex}" \
                    ./finalize_pilot.sbatch)
        _prev6="${J6[$_pop]}"
    else
        echo "  (skipped: START_STAGE=${START_STAGE} STOP_STAGE=${STOP_STAGE})"
    fi
done

# ---------------------------------------------------------------- summary
echo
echo "=========================================================="
echo "orchestrator complete — jobs queued with afterok chain:"
row () { [[ -n "$2" ]] && printf '  %-32s = %s\n' "$1" "$2" || printf '  %-32s   (not submitted)\n' "$1"; }
row "stage 1 (SGDP prep)"        "$J1"
row "stage 2 (rebuild panel)"    "$J2"
row "stage 2b (KING dedup)"      "$J2B"
row "stage 3 (4a global K=7)"    "$J3"
row "stage 3 (4a lai3way K=5)"   "$J3L"
for _pop in "${COHORT_LIST[@]}"; do
    row "stage 4 (simulate ${_pop})"          "${J4[$_pop]}"
    row "stage 5 (8 arrays ${_pop})"          "${J5[$_pop]//:/ }"
    row "stage 6 (accuracy + push ${_pop})"   "${J6[$_pop]}"
done
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
