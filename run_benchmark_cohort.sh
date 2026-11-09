#!/usr/bin/env bash
# ============================================================================
# run_benchmark_cohort.sh  --  idempotent driver for the benchmark cohort
#
# Runs the cohort build from the K=8 Q table to a transfer-ready set. Every
# stage checks whether its own output is already present AND current before
# doing anything, so re-running after a failure, a timeout or a fix costs only
# the stages that still need doing. Nothing is deleted without --force.
#
# WHY A DRIVER AND NOT A CHAIN OF sbatch --dependency. The fast stages (split,
# assemble, report, preflight) take seconds and belong in the foreground where
# their output can be read before the slow ones are submitted. The slow stages
# are arrays whose dependencies are genuinely a DAG, not a line: 6b and 6c are
# independent, 6d and 6e both wait on the pair. A hand-written chain gets that
# wrong in a way that only shows up hours later.
#
# THE PREFLIGHT IS THE POINT. ingest_simset.sh refuses a cohort that breaks its
# rules, and it runs at the END -- after simulation, chip files and a transfer.
# Several of those rules are decided at the split and the cell deal, which take
# seconds. scripts/check-simset-contract.py re-states the ones knowable from
# the plan and runs them BEFORE anything is submitted.
#
# Usage:
#   ./run_benchmark_cohort.sh                 # do whatever still needs doing
#   ./run_benchmark_cohort.sh --dry-run       # say what it would do
#   ./run_benchmark_cohort.sh --from simulate # skip ahead
#   ./run_benchmark_cohort.sh --only preflight
#   ./run_benchmark_cohort.sh --force split   # redo one stage
#
# Stages, in order:
#   check      4a landed, vocabulary sane, anchors present
#   split      donor/reference lists, both modes, distinct suffixes
#   assemble   cohort manifests + chip/build assignment
#   report     composition + counts
#   preflight  the intake's rules, before anything is submitted
#   simulate   6b (22 tasks)   + 6c (22 tasks), independent
#   cells      6d (22 tasks)   after simulate
#   kinship    6e (1 job)      after simulate
#   verify     contract preflight again, now with the built cohort
# ============================================================================
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
K="${K:-8}"
PANEL_TSV="${PANEL_TSV:-${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv}"
KING_TSV="${KING_TSV:-${HOMOG_DIR}/king/king.kin0}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"
RESULTS="${RESULTS:-${REPO_DIR}/results}"
LAI_SUFFIX="${LAI_SUFFIX:-lai}"
BENCH_SUFFIX="${BENCH_SUFFIX:-bench}"
# THE CONTRACT'S DESIGN IS one-cell. Doc 20 rule 3: every group split equally
# over the 3 chips and, within a chip, half hg19 and half hg38 -- six cells,
# each person in exactly ONE. Intake check 2.chips enforces it. paired-builds
# is statistically stronger for the liftover comparison but the intake refuses
# it, so it is not the default here.
DESIGN="${DESIGN:-one-cell}"
# The contract's threshold is FLAT 0.93 (Doc 20, donors.tsv and
# real_samples.tsv; plan_simset.R --threshold 0.93). Intake 2.real_samples
# compares sample_type against top_fraction at this one number, so the tiered
# 0.95/0.93 the panel build uses would disagree for EUR/EAS/AFR between them.
FLAT_THRESHOLD="${FLAT_THRESHOLD:-0.93}"
NO_DONOR_SUPERPOPS="${NO_DONOR_SUPERPOPS:-OCE}"
# The contract's reference share: AFR/EUR/EAS/SAS take their reference from
# 1000 Genomes ONLY, NAT adds half the MXB, MEN a third of its pool. 1KG can
# never BE a donor, so spending non-1KG samples on a role 1KG already fills
# wastes the only donor-eligible material there is -- it roughly doubles the
# SAS donor pool, which is the binding constraint on N_eff.
REF_SHARE="${REF_SHARE:-default=0,NAT@MXB=0.5,MEN=0.33}"
# ...but 1000 Genomes carries none of the deeply divergent lineages the other
# sources do (no Biaka/Mbuti/San/Mandenka in AFR, no Basque/Sardinian/Orcadian
# in EUR, none of the HGDP minority populations in EAS). A donor haplotype
# whose ancestry no reference covers measures the panel's gap, not the method.
# 3 per population anchors each lineage and still leaves SAS donors well above
# the even split. 0 is the strict contract reading. Coverage is counted on
# ANCESTRY GROUPS (reference_ids/population_aliases.tsv), not raw labels, so
# HGDP Yoruba is not reported as a gap when 1000 Genomes YRI is in the panel.
REF_MIN_PER_POP="${REF_MIN_PER_POP:-3}"

DRY=0; FROM=""; ONLY=""; FORCE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY=1; shift ;;
        --from)    FROM="$2"; shift 2 ;;
        --only)    ONLY="$2"; shift 2 ;;
        --force)   FORCE="$2"; shift 2 ;;
        -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

STAGES=(check split assemble report preflight simulate cells kinship verify)
want () {   # should this stage run at all?
    local s="$1"
    [[ -n "$ONLY" ]] && { [[ "$ONLY" == "$s" ]] && return 0 || return 1; }
    if [[ -n "$FROM" ]]; then
        local seen=0
        for x in "${STAGES[@]}"; do
            [[ "$x" == "$FROM" ]] && seen=1
            [[ "$x" == "$s" ]] && { [[ "$seen" == 1 ]] && return 0 || return 1; }
        done
    fi
    return 0
}
forced () { [[ "$FORCE" == "$1" ]]; }
say  () { printf '[%s] %s\n' "$(date +%T)" "$*"; }
head_ () { printf '\n==== %-12s %s\n' "$1" "$2"; }
run  () {
    if [[ "$DRY" == 1 ]]; then printf '      WOULD RUN: %s\n' "$*"; return 0; fi
    "$@"
}

mkdir -p "$RESULTS"
FAILED=0

# ---------------------------------------------------------------------------
# check -- is 4a actually done, and is its output usable
# ---------------------------------------------------------------------------
if want check; then
    head_ check "4a output"
    if [[ ! -s "$PANEL_TSV" ]]; then
        say "MISSING: $PANEL_TSV"
        say "4a has not finished. Nothing downstream can run."
        exit 1
    fi
    N_ROWS=$(( $(wc -l < "$PANEL_TSV") - 1 ))
    N_INCL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="included") c=i; next} $c==1' "$PANEL_TSV" | wc -l)
    N_ANCHOR=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="assigned") c=i; next}
                           c && $c!="" && $c!="-" && $c!="NA" && $c!="."' "$PANEL_TSV" | wc -l)
    N_CHMI=$(grep -c CHMI "$PANEL_TSV" || true)
    say "panel table: ${N_ROWS} rows, ${N_INCL} included, ${N_ANCHOR} anchors"
    [[ "$N_CHMI" -eq 0 ]] || { say "WARNING: ${N_CHMI} CHMI row(s) still present"; }
    if [[ "$N_ANCHOR" -eq 0 ]]; then
        say "WARNING: the 'assigned' column is empty -- 4a's .pop prior was not"
        say "         written, so no Q in this table is a supervised estimate."
    fi
    if [[ -s "${HOMOG_DIR}/anchor_summary.k${K}.tsv" ]]; then
        THIN=$(awk -F'\t' 'NR>1 && $5==1 {printf " %s(%d)", $1, $2}' \
               "${HOMOG_DIR}/anchor_summary.k${K}.tsv")
        [[ -n "$THIN" ]] && say "THIN ANCHOR SETS:${THIN} -- see anchor_summary"
    fi
    if [[ ! -s "$KING_TSV" ]]; then
        say "WARNING: no KING table at ${KING_TSV}."
        say "         The split cannot group relatives, and the intake's"
        say "         cross-role rule (kinship >= 0.177 in two roles) fails"
        say "         only AFTER the whole build. Run 1e first."
    fi
fi

# ---------------------------------------------------------------------------
# split -- both modes, DISTINCT suffixes
# ---------------------------------------------------------------------------
if want split; then
    head_ split "donor / reference lists"
    SPLITTER="${REPO_DIR}/scripts/split-homog-donors-reference.py"
    need_split=0
    for sfx in "$LAI_SUFFIX" "$BENCH_SUFFIX"; do
        ls "${REFS}"/*_reference_"${sfx}".txt >/dev/null 2>&1 || need_split=1
    done
    forced split && need_split=1
    if [[ "$need_split" == 0 ]]; then
        say "present for suffixes '${LAI_SUFFIX}' and '${BENCH_SUFFIX}', skipping"
    else
        REL=(); [[ -s "$KING_TSV" ]] && REL=(--relatedness "$KING_TSV")
        run python3 -I "$SPLITTER" --panel "$PANEL_TSV" --out "$REFS" \
            --mode lai --suffix "$LAI_SUFFIX" "${REL[@]}" || FAILED=1
        run python3 -I "$SPLITTER" --panel "$PANEL_TSV" --out "$REFS" \
            --mode benchmark --suffix "$BENCH_SUFFIX" \
            --no-donor-superpops "$NO_DONOR_SUPERPOPS" \
            --ref-share "$REF_SHARE" --ref-min-per-pop "$REF_MIN_PER_POP" \
            --pop-aliases "${REPO_DIR}/reference_ids/population_aliases.tsv" \
            "${REL[@]}" || FAILED=1
    fi
    # The two suffixes MUST differ: both modes write the same filenames, so a
    # shared suffix means the second silently overwrote the first.
    if [[ "$LAI_SUFFIX" == "$BENCH_SUFFIX" ]]; then
        say "FATAL: LAI_SUFFIX and BENCH_SUFFIX are both '${LAI_SUFFIX}'."
        say "       The second split would overwrite the first with nothing"
        say "       in the filename to say which policy produced it."
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# assemble -- cohort manifests
# ---------------------------------------------------------------------------
MAN="${RESULTS}/benchmark_cohort_manifest.tsv"
IND="${RESULTS}/benchmark_cohort_simset_individuals.tsv"
if want assemble; then
    head_ assemble "cohort manifests (design=${DESIGN})"
    if [[ -s "$MAN" && -s "$IND" ]] && ! forced assemble; then
        say "present ($(( $(wc -l < "$MAN") - 1 )) real, $(( $(wc -l < "$IND") - 1 )) simulated), skipping"
    else
        run python3 -I "${REPO_DIR}/scripts/assemble-benchmark-cohort.py" \
            --panel "$PANEL_TSV" --refs "$REFS" --suffix "$BENCH_SUFFIX" \
            --design "$DESIGN" --flat-threshold "$FLAT_THRESHOLD" \
            --out "${RESULTS}/benchmark_cohort" || FAILED=1
    fi
fi

# ---------------------------------------------------------------------------
# report -- composition
# ---------------------------------------------------------------------------
if want report; then
    head_ report "composition"
    COMP="${RESULTS}/cohort_composition.txt"
    if [[ -s "$COMP" ]] && ! forced report; then
        say "present, skipping"
    else
        run python3 -I "${REPO_DIR}/scripts/cohort-composition-report.py" \
            --panel "$PANEL_TSV" --refs "$REFS" --suffix "$BENCH_SUFFIX" \
            --out "$COMP" || FAILED=1
    fi
fi

# ---------------------------------------------------------------------------
# preflight -- the intake's rules, BEFORE anything is submitted
# ---------------------------------------------------------------------------
if want preflight; then
    head_ preflight "ingest_simset.sh rules"
    PF=("${REPO_DIR}/scripts/check-simset-contract.py"
        --manifest "$MAN" --individuals "$IND"
        --refs "$REFS" --suffix "$BENCH_SUFFIX"
        --threshold "$FLAT_THRESHOLD" --panel "$PANEL_TSV")
    [[ -s "$KING_TSV" ]] && PF+=(--king "$KING_TSV")
    if [[ "$DRY" == 1 ]]; then
        printf '      WOULD RUN: python3 -I %s ...\n' "${PF[0]}"
    else
        if ! python3 -I "${PF[@]}"; then
            say ""
            say "PREFLIGHT FAILED -- not submitting anything."
            say "Every rule above is decided at the split or the cell deal,"
            say "which take seconds. The intake runs after simulation, chip"
            say "files and a transfer, so fixing it now is the cheap moment."
            exit 1
        fi
    fi
fi

# ---------------------------------------------------------------------------
# simulate -- 6b and 6c, independent of each other
# ---------------------------------------------------------------------------
SIMSET_DIR="${SIMSET_DIR:-${PROJECT_ROOT}/06_benchmark_simset}"
REAL_DIR="${REAL_DIR:-${PROJECT_ROOT}/07_benchmark_real}"
JOB_6B=""; JOB_6C=""
sub () {   # sub <script> [extra sbatch args...] -> prints the job id
    local script="$1"; shift
    if [[ "$DRY" == 1 ]]; then printf 'DRYRUN'; return 0; fi
    sbatch --parsable "$@" "$script" 2>/dev/null
}
n_done () { ls "$1" 2>/dev/null | grep -c . || true; }

if want simulate; then
    head_ simulate "6b simulate + 6c real genomes"
    N_SIM_CHR=$(ls -d "${SIMSET_DIR}"/chr* 2>/dev/null | wc -l)
    if [[ "$N_SIM_CHR" -ge 22 ]] && ! forced simulate; then
        say "6b: ${N_SIM_CHR} chromosome dirs present, skipping"
    else
        say "6b: ${N_SIM_CHR}/22 chromosomes present -- submitting"
        JOB_6B=$(sub "${REPO_DIR}/6b_simulate_benchmark_cohort_clm.sh")
        say "6b: job ${JOB_6B}"
    fi
    N_REAL=$(n_done "${REAL_DIR}/benchmark_real.chr*.bcf")
    if [[ "$N_REAL" -ge 22 ]] && ! forced simulate; then
        say "6c: ${N_REAL} BCFs present, skipping"
    else
        say "6c: ${N_REAL}/22 BCFs present -- submitting"
        JOB_6C=$(sub "${REPO_DIR}/6c_extract_real_test_genomes_clm.sh")
        say "6c: job ${JOB_6C}"
    fi
fi

# dependency string for the stages that need BOTH
DEP=""
for j in $JOB_6B $JOB_6C; do
    [[ "$j" == "DRYRUN" || -z "$j" ]] && continue
    DEP="${DEP}${DEP:+:}${j}"
done
DEPARG=(); [[ -n "$DEP" ]] && DEPARG=(--dependency="afterok:${DEP}")

# ---------------------------------------------------------------------------
# cells -- 6d, after both halves exist
# ---------------------------------------------------------------------------
if want cells; then
    head_ cells "6d chip x build x chr"
    CELLS_DIR="${CELLS_DIR:-${PROJECT_ROOT}/08_chip_cells}"
    N_CELLS=$(n_done "${CELLS_DIR}/benchmark.*.bcf")
    if [[ "$N_CELLS" -gt 0 ]] && ! forced cells; then
        say "${N_CELLS} cell file(s) present, skipping"
    else
        JOB_6D=$(sub "${REPO_DIR}/6d_write_chip_cells_clm.sh" "${DEPARG[@]}")
        say "6d: job ${JOB_6D}${DEP:+ (after ${DEP})}"
        say "    only chips with a site list produce output; the rest skip"
    fi
fi

# ---------------------------------------------------------------------------
# kinship -- 6e, parallel to 6d
# ---------------------------------------------------------------------------
if want kinship; then
    head_ kinship "6e genome-wide kinship"
    KOUT="${PROJECT_ROOT}/09_cohort_kinship/cohort_king.kin0"
    if [[ -s "$KOUT" ]] && ! forced kinship; then
        say "present, skipping"
    else
        JOB_6E=$(sub "${REPO_DIR}/6e_cohort_kinship_clm.sh" "${DEPARG[@]}")
        say "6e: job ${JOB_6E}${DEP:+ (after ${DEP})}"
    fi
fi

# ---------------------------------------------------------------------------
# verify -- the preflight again, now against the built cohort
# ---------------------------------------------------------------------------
if want verify; then
    head_ verify "post-build checks"
    KOUT="${PROJECT_ROOT}/09_cohort_kinship/cohort_king.kin0"
    if [[ ! -s "$KOUT" ]]; then
        say "6e has not produced a kinship table yet -- re-run this stage"
        say "once it has: ./run_benchmark_cohort.sh --only verify"
    else
        run python3 -I "${REPO_DIR}/scripts/check-simset-contract.py" \
            --manifest "$MAN" --individuals "$IND" --refs "$REFS" \
            --suffix "$BENCH_SUFFIX" --threshold "$FLAT_THRESHOLD" \
            --king "$KOUT" --panel "$PANEL_TSV" || FAILED=1
        run python3 -I "${REPO_DIR}/scripts/cohort-power-summary.py" \
            --kinship "$KOUT" --manifest "$MAN" --individuals "$IND" \
            --out "${RESULTS}/cohort_power_summary.txt" || true
    fi
fi

echo ""
echo "===================================================================="
if [[ "$DRY" == 1 ]]; then
    echo " dry run -- nothing was submitted or written"
elif [[ "$FAILED" == 1 ]]; then
    echo " FINISHED WITH ERRORS -- read the stages above"
else
    echo " ok"
fi
echo " re-run this script at any point; finished stages are skipped."
echo " transfer when the cohort is built:"
echo "   ./6a_transfer_to_mab_clm.sh verify"
echo "   ./6a_transfer_to_mab_clm.sh benchmark"
echo "===================================================================="
[[ "$FAILED" == 1 ]] && exit 1
exit 0
