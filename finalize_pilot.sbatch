#!/usr/bin/env bash
#SBATCH --job-name=finalize_pilot
#SBATCH --partition=atkinson,mhgcp
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/finalize_pilot_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/finalize_pilot_%j.err
#
# Fires after all fix arrays exit ok (--dependency=afterok:...). Verifies no
# broken .classes remain, scores RFMix WGS + chip via accuracy_v2.R and, when
# SCORE_FLARE=1 (the driver sets it when the FLARE arms were submitted), FLARE
# WGS + chip the same way (scripts/flare-to-lat3.py puts FLARE's calls in the
# .Lat3/.classes/.done layout accuracy_v2.R reads), then commits + pushes the
# fresh TSVs so the interactive session can walk away.
#
# Deliberately skips Python figure/table/pptx regen — the cluster's
# shapeit5 env is missing matplotlib. The main-loop side does that once
# the push lands.
set -euo pipefail

REPO_DIR="/storage/atkinson/home/magyar/LAI-sims-accuracy"
PROJECT_ROOT="/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank"
ADMIX_POP="Brasa"
GEN="12"
CHRS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"

WGS_WD="${PROJECT_ROOT}/03_rfmix/${ADMIX_POP}/gen${GEN}"
CHIP_WD="${PROJECT_ROOT}/03_rfmix/chip_gsa_jessica/${ADMIX_POP}/gen${GEN}"
SCORE_FLARE="${SCORE_FLARE:-0}"
CHIP_NAME="${CHIP_NAME:-gsa_jessica}"
FLARE_WD="${PROJECT_ROOT}/03_flare/${ADMIX_POP}/gen${GEN}"
FLARE_CHIP_WD="${PROJECT_ROOT}/03_flare_${CHIP_NAME}/${ADMIX_POP}/gen${GEN}"

# --- conda / R environment ---
# unset any stale python env so shapeit5's R gets the right libraries
unset PYTHONHOME PYTHONPATH
source "$HOME/.bashrc" 2>/dev/null || true
if command -v conda &>/dev/null; then
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate shapeit5 2>/dev/null || echo "[warn] conda activate shapeit5 failed"
fi

echo "==================================================================="
echo "[finalize_pilot] host=$(hostname)  time=$(date -u +%FT%TZ)"
echo "[finalize_pilot] repo=$REPO_DIR"
echo "[finalize_pilot] chrs=$CHRS"
echo "==================================================================="

# ---------- 0. Hard precheck — every required path/tool/input must exist ----------
die() { echo "[finalize_pilot] FATAL: $*" >&2; exit 2; }

echo
echo ">>> [0/4] Precheck"

# 0a. Repo layout
[ -d "$REPO_DIR" ]                              || die "REPO_DIR missing: $REPO_DIR"
[ -f "$REPO_DIR/run_accuracy_pilot.sh" ]        || die "run_accuracy_pilot.sh missing"
[ -f "$REPO_DIR/accuracy_v2.R" ]                || die "accuracy_v2.R missing"
[ -d "$REPO_DIR/results" ]                      || die "results/ missing"
cd "$REPO_DIR"
echo "  repo: ok"

# 0b. Cluster paths
[ -d "$PROJECT_ROOT" ]                          || die "PROJECT_ROOT missing: $PROJECT_ROOT"
[ -d "$WGS_WD" ]                                || die "WGS WORKDIR missing: $WGS_WD"
[ -d "$CHIP_WD" ]                               || die "CHIP WORKDIR missing: $CHIP_WD"
[ -d "${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}" ] \
    || die "SIM_DIR missing: ${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}"
echo "  cluster paths: ok"

# 0c. Tools
# R_LIBS_USER exactly as run_accuracy_pilot.sh sets it, BEFORE the dplyr check:
# /opt/biobuilds-2016.04 otherwise shadows the env's tidyverse with a 2016
# magrittr R 4 refuses to load, and this precheck would fail on a library the
# scoring step itself never uses. Exported, so the scoring step sees the same.
if [[ -z "${R_LIBS_USER:-}" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib/R/library" ]]; then
        export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"
    elif [[ -d "$HOME/.conda/envs/shapeit5/lib/R/library" ]]; then
        export R_LIBS_USER="$HOME/.conda/envs/shapeit5/lib/R/library"
    fi
fi
echo "  R_LIBS_USER=${R_LIBS_USER:-<unset>}"
command -v Rscript >/dev/null                   || die "Rscript not on PATH (conda activate shapeit5?)"
command -v git     >/dev/null                   || die "git not on PATH"
Rscript -e 'suppressMessages(library(dplyr))' 2>/dev/null \
    || die "R env broken: can't load dplyr (magrittr version issue?)"
echo "  tools: ok (Rscript, git, dplyr all loadable)"

# 0d. Every (track, panel, chr) combo we expect to score must have all 3 files
#     (.done marker + .Lat3 + .classes). Fail if any are missing.
PANELS_MAIN="NAT_HGDP NAT_PEL NAT_PEL_EAS NAT_HGDPMXB NAT_HGDPMXB_FULL NAT_HOMOG"
TRACKS="NAT NATMXB"
IFS=',' read -ra CHR_ARR <<< "$CHRS"
missing_wgs=0; missing_chip=0
for chr in ${CHR_ARR[@]+"${CHR_ARR[@]}"}; do
    for track in $TRACKS; do
        for panel in $PANELS_MAIN; do
            # skip the deliberately-excluded combo (NATMXB + HGDPMXB_FULL double-dips donors)
            [ "$track" = "NATMXB" ] && [ "$panel" = "NAT_HGDPMXB_FULL" ] && continue
            for ext in done Lat3 classes; do
                f_wgs="$WGS_WD/${track}.${panel}.gen${GEN}_chr${chr}.${ext}"
                f_chp="$CHIP_WD/${track}.${panel}.gen${GEN}_chr${chr}.${ext}"
                [ -f "$f_wgs" ] || { echo "  WGS missing: $f_wgs"  ; missing_wgs=$((missing_wgs+1)); }
                [ -f "$f_chp" ] || { echo "  CHIP missing: $f_chp"; missing_chip=$((missing_chip+1)); }
            done
        done
    done
done
if [ "$missing_wgs" -gt 0 ] || [ "$missing_chip" -gt 0 ]; then
    die "$missing_wgs WGS + $missing_chip CHIP required files missing — do not score partial output"
fi
echo "  input files: ok (every chr × combo × 3 files present in both WGS and CHIP)"

# 0e. FLARE: the same grid, one .anc.vcf.gz per combo, at both densities.
if [ "$SCORE_FLARE" = "1" ]; then
    missing_flare=0
    for chr in ${CHR_ARR[@]+"${CHR_ARR[@]}"}; do
        for track in $TRACKS; do
            for panel in $PANELS_MAIN; do
                [ "$track" = "NATMXB" ] && [ "$panel" = "NAT_HGDPMXB_FULL" ] && continue
                for wd in "$FLARE_WD" "$FLARE_CHIP_WD"; do
                    f="$wd/${track}.${panel}.gen${GEN}_chr${chr}.anc.vcf.gz"
                    [ -s "$f" ] || { echo "  FLARE missing: $f"; missing_flare=$((missing_flare+1)); }
                done
            done
        done
    done
    [ "$missing_flare" -eq 0 ] || die "$missing_flare FLARE combos missing — do not score partial output"
    [ -f "$REPO_DIR/scripts/flare-to-lat3.py" ] || die "scripts/flare-to-lat3.py missing"
    echo "  FLARE input files: ok (WGS + chip)"
else
    echo "  FLARE: not scored (SCORE_FLARE=${SCORE_FLARE})"
fi

# ---------- 1. Verify .classes are all healthy (fail if any broken) ----------
echo
echo ">>> [1/4] Verifying .classes healthy"
broken=$(
  for root in "$WGS_WD" "$CHIP_WD"; do
    [ -d "$root" ] || continue
    for f in "$root"/*.classes; do
      [ -f "$f" ] || continue
      n=$(tr ' ' '\n' < "$f" | awk '$1 ~ /^[0-9]+$/ && $1+0>0' | sort -u | wc -l)
      if [ "$n" -lt 3 ]; then
        printf "  BROKEN\t%s\t%dcls\n" "$f" "$n"
      fi
    done
  done
)
if [ -n "$broken" ]; then
    echo "[finalize_pilot] ABORT: broken .classes still present:"
    echo "$broken"
    exit 2
fi
echo "  ok, no broken .classes files"

# ---------- 2. WGS accuracy pilot ----------
echo
echo ">>> [2/4] WGS pilot ($WGS_WD)"
CHRS="$CHRS" \
PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
WORKDIR="$WGS_WD" OUTDIR="${PROJECT_ROOT}/04_accuracy" \
REVIEW_DIR="${REPO_DIR}/results/accuracy_pilot" \
bash "$REPO_DIR/run_accuracy_pilot.sh"

# ---------- 3. Chip accuracy pilot ----------
echo
echo ">>> [3/4] Chip-GSA pilot ($CHIP_WD)"
CHIP_OUT="/tmp/accuracy_chip_gsa_final_$$"
CHIP_REVIEW="${REPO_DIR}/results/accuracy_pilot_chip_gsa"
rm -rf "$CHIP_REVIEW"
mkdir -p "$CHIP_REVIEW/plots"

CHRS="$CHRS" \
PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
WORKDIR="$CHIP_WD" OUTDIR="$CHIP_OUT" \
REVIEW_DIR="$CHIP_REVIEW" \
bash "$REPO_DIR/run_accuracy_pilot.sh"

# ---------- 3b. FLARE accuracy (WGS + chip), the same scorer ----------
TAGS="accuracy_pilot accuracy_pilot_chip_gsa"
if [ "$SCORE_FLARE" = "1" ]; then
    for density in wgs chip; do
        if [ "$density" = wgs ]; then FWD="$FLARE_WD"; FTAG=accuracy_pilot_flare; FOUT="${PROJECT_ROOT}/04_accuracy_flare"
        else FWD="$FLARE_CHIP_WD"; FTAG=accuracy_pilot_flare_chip_gsa; FOUT="/tmp/accuracy_flare_chip_final_$$"; fi
        echo
        echo ">>> [3b] FLARE ${density} ($FWD)"
        python3 -I "$REPO_DIR/scripts/flare-to-lat3.py" --flare-dir "$FWD" --out-dir "$FWD/lat3" \
            --gen "$GEN" --jobs "${SLURM_CPUS_PER_TASK:-4}" || die "FLARE ${density}: conversion failed (above)"
        rm -rf "${REPO_DIR}/results/${FTAG}"; mkdir -p "${REPO_DIR}/results/${FTAG}/plots"
        CHRS="$CHRS" \
        PROJECT_ROOT="$PROJECT_ROOT" ADMIX_POP="$ADMIX_POP" GEN="$GEN" \
        WORKDIR="$FWD/lat3" OUTDIR="$FOUT" \
        REVIEW_DIR="${REPO_DIR}/results/${FTAG}" \
        bash "$REPO_DIR/run_accuracy_pilot.sh"
        TAGS="$TAGS $FTAG"
    done
fi

# ---------- 3.5 Post-run assertion: every chr present, no silent skips ----------
echo
echo ">>> [3.5] Post-run assertion — per-chr coverage complete"
for tag in $TAGS; do
    f="${REPO_DIR}/results/${tag}/accuracy_per_chr.tsv"
    [ -s "$f" ] || die "post-run: expected TSV missing/empty: $f"
    # Extract chr column for NAT_HGDP/NAT/NAT (proxy for full coverage)
    chrs_present=$(awk -F'\t' 'NR>1 && $2=="NAT_HGDP" && $1=="NAT" && $3=="NAT" {print $4}' \
        "$f" | sort -nu | paste -sd,)
    n_chrs=$(echo "$chrs_present" | tr ',' '\n' | wc -l)
    if [ "$n_chrs" -lt 22 ]; then
        die "post-run: $tag has only $n_chrs / 22 chrs — got: $chrs_present"
    fi
    echo "  $tag: 22/22 chrs present ($chrs_present)"
done
# Also assert every panel in the summary has full-coverage n_sites
TAGS="$TAGS" Rscript -e '
    for (path in file.path("results", strsplit(Sys.getenv("TAGS"), " +")[[1]], "accuracy_summary.tsv")) {
        d <- read.delim(path, sep="\t", stringsAsFactors=FALSE)
        d <- d[d$ancestry=="NAT" & d$track=="NAT",]
        max_n <- max(d$total_sites)
        low   <- d[d$total_sites < 0.5 * max_n, c("panel","total_sites")]
        if (nrow(low)) {
            cat(sprintf("[post-run] %s: partial-coverage panels detected:\n", path))
            print(low); quit(status=3)
        }
    }
    cat("[post-run] all panels have >= 50% of max per-panel site coverage\n")
'

# ---------- 4. Commit + push ----------
echo
echo ">>> [4/4] git commit + push"
cd "$REPO_DIR"
# Sync with origin in case any hand-commits landed
git fetch origin claude/setup-shapeit5-conda-xlcdl 2>&1 | tail -3
git pull --rebase -X theirs origin claude/setup-shapeit5-conda-xlcdl 2>&1 | tail -5 || {
    # rebase conflicts: prefer the pilot's fresh TSVs
    git checkout --theirs results/accuracy_pilot*/*.tsv 2>/dev/null || true
    git add results/
    git rebase --continue || git rebase --abort
}

git add -A results/
if git diff --cached --quiet; then
    echo "  no result changes to commit — nothing to push"
else
    git commit -m "accuracy pilot: full 22-chr RFMix WGS + chip$([ "$SCORE_FLARE" = 1 ] && echo ' and FLARE WGS + chip') (finalize_pilot job $SLURM_JOB_ID)

Auto-committed by finalize_pilot.sbatch after the 3b arrays completed:
every (track, panel, chr) combo present, no 1- or 2-class .classes file,
22/22 chromosomes in every accuracy table ($TAGS)."
    for attempt in 1 2 3; do
        if git push origin claude/setup-shapeit5-conda-xlcdl; then
            echo "  push ok on attempt $attempt"
            break
        fi
        echo "  push failed on attempt $attempt, sleeping..."
        sleep $((attempt * 10))
    done
fi

echo
echo "==================================================================="
echo "[finalize_pilot] DONE at $(date -u +%FT%TZ)"
echo "==================================================================="
