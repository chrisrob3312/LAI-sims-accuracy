#!/usr/bin/env bash
#SBATCH --job-name=merge_sites
#SBATCH --partition=atkinson,mhgcp
# Same exclusion as every other stage: the mhgcp-t* CPUs lack instructions this
# conda env's bcftools is built with (SIGILL, exit 132, empty log).
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/merge_sites_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/merge_sites_%j.err
# =============================================================================
# merge-sites-collect.sh -- gather the merge's per-chromosome site tables
# (1c's provenance/, fill check, soft-union groups, ADMIXTURE site rule, and
# scripts/cohort-site-report.sh output when it exists) into a few TSVs under
# results/merge_sites/ in the repo, to commit and plot off the cluster.
# READ-ONLY on the panel. Columns: scripts/merge-sites-collect.py.
#
#   bash scripts/merge-sites-collect.sh                    # 1c's tables: seconds
#   SKIP_SITE_REPORT=1 bash scripts/merge-sites-collect.sh # same, never reads site reports
#   sbatch --dependency=afterok:<site_report job> scripts/merge-sites-collect.sh
#                                                          # after the per-site report:
#                                                          # ~27M records, minutes
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
OUT="${OUT:-${REPO_DIR}/results/merge_sites}"
EPOCH="${PROJECT_ROOT}/rebuild_epoch"
CHROMS="${CHROMS:-$(seq -s ' ' 1 22)}"
PY="${REPO_DIR}/scripts/merge-sites-collect.py"
[[ -f "$PY" ]] || { echo "ERROR: ${PY} not found (run from the repo, or set REPO_DIR)" >&2; exit 2; }

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    module load anaconda3/2024.06
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    set +u; conda activate shapeit5; set -u
    unset PYTHONHOME PYTHONPATH
fi
command -v bcftools > /dev/null || { echo "ERROR: bcftools not on PATH (conda activate shapeit5)" >&2; exit 2; }
[[ -d "${PANEL_DIR}/provenance" ]] || { echo "ERROR: no ${PANEL_DIR}/provenance (1c has not run)" >&2; exit 2; }
mkdir -p "$OUT" || exit 2
echo "[merge-sites] panel: ${PANEL_DIR}"
echo "[merge-sites] out:   ${OUT}"

# Per chromosome, from the published files themselves.
TOT=$(mktemp "${TMPDIR:-/tmp}/merge_sites_totals.XXXXXX") || exit 2
trap 'rm -f "$TOT"' EXIT
printf 'chrom\tsamples\tsoftunion_records\tadmix_positions\tsoftunion_new\trecipe\n' > "$TOT"
for CHR in $CHROMS; do
    CHR="${CHR#chr}"
    SU="${PANEL_DIR}/merged_chr${CHR}.softunion.bcf"
    AB="${PANEL_DIR}/merged_chr${CHR}.admix_sites.bed"
    if [[ -s "$SU" ]]; then
        ns=$(bcftools query -l "$SU" | wc -l)
        nr=$(bcftools index -s "$SU" | awk '{s += $3} END {print s + 0}')
    else
        ns=NA; nr=NA
    fi
    if [[ -s "$AB" ]]; then na=$(wc -l < "$AB"); else na=NA; fi
    new=0; [[ -s "$SU" && -f "$EPOCH" && "$SU" -nt "$EPOCH" ]] && new=1
    rc=NA; [[ -s "${SU}.recipe" ]] && rc=$(head -c 400 "${SU}.recipe" | tr '\t\n' '  ')
    printf 'chr%s\t%s\t%s\t%s\t%s\t%s\n' "$CHR" "$ns" "$nr" "$na" "$new" "$rc" >> "$TOT"
done
echo "[merge-sites] per-chromosome totals (raw, recipe cut):"
cut -f1-5 "$TOT" | sed 's/^/  /'

python3 -I "$PY" --panel-dir "$PANEL_DIR" --out-dir "$OUT" --totals "$TOT" --epoch "$EPOCH" \
    --chroms "$CHROMS" ${SKIP_SITE_REPORT:+--no-site-report} \
    || { echo "ERROR: merge-sites-collect.py failed (stderr above)" >&2; exit 1; }
echo "[merge-sites] done: $(ls "$OUT" | tr '\n' ' ')"
