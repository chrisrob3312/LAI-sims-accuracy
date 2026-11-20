#!/usr/bin/env bash
# ============================================================================
# 6c_extract_real_test_genomes_clm.sh  --  the real half of the benchmark cohort
#
#SBATCH --job-name=6c_real
#SBATCH --partition=atkinson,mhgcp
#SBATCH --array=1-22
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/6c_real_%a_%j.out
#SBATCH --error=logs/6c_real_%a_%j.err
# The t-nodes cannot execute these binaries: 270/270 tasks of an earlier array
# died there with SIGILL (exit 132) before writing a log. c02 is fine.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#
# Subsets the phased panel to the real individuals that
# assemble-benchmark-cohort.py put in the cohort, one BCF per chromosome.
# These are the genomes the harmonization pipeline is tested on alongside the
# simulated 2,160, so they need to exist as data, not just as a manifest row.
#
# WHY THE MANIFEST IS THE ONLY SOURCE OF THE SAMPLE LIST. The four exclusions
# that make a real genome test-eligible -- not a donor, not in either
# reference panel, not 1000 Genomes, not in sgdp_server_panel_ids.txt -- are
# applied in one place, by the assembler. Re-deriving the list here from the
# panel plus some filters would mean two implementations of the same rule, and
# the one that drifts is the one that silently puts a donor in the test set.
# So this script applies no policy at all: it reads column 1 of the manifest.
#
# Paths are all overridable (PROJECT_ROOT, PANEL_DIR, OUTDIR, MANIFEST).
# ============================================================================
set -euo pipefail

CHR_N="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job (1-22)}"
CHR="chr${CHR_N}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_%s.shapeit5_phased.softunion_maf005.bcf}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/07_benchmark_real}"
MANIFEST="${MANIFEST:-${REPO_DIR}/results/benchmark_cohort_manifest.tsv}"
CONDA_ENV="${CONDA_ENV:-shapeit5}"
MIN_PANEL="${MIN_PANEL:-3900}"

# Conda before any tool call. A check that ran ahead of activation once made
# 2 of 594 chunk tasks exit 127 on "bcftools: command not found".
if [[ -z "${CONDA_PREFIX:-}" || "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH in ${CONDA_ENV}" >&2; exit 127; }

PANEL=$(printf "$PANEL_TPL" "$CHR")
[[ -s "$PANEL" ]] || { echo "ERROR: missing panel $PANEL" >&2; exit 1; }
[[ -s "$MANIFEST" ]] || { echo "ERROR: missing ${MANIFEST}
       Run scripts/assemble-benchmark-cohort.py first. Without it there is no
       test-eligible sample list, and guessing one here would risk putting a
       donor in the test set." >&2; exit 1; }

N_PANEL=$(bcftools query -l "$PANEL" | wc -l)
if [[ "$N_PANEL" -lt "$MIN_PANEL" ]]; then
    echo "ERROR: panel has ${N_PANEL} samples, expected >= ${MIN_PANEL}.
       This is the pre-rerun build; SGDP, Almarri and the restored samples
       are missing, and most of the real test cohort lives in those." >&2
    exit 3
fi

mkdir -p "$OUTDIR"
KEEP="${OUTDIR}/.keep_${CHR}.txt"
awk -F'\t' 'NR>1 && NF && $1 != "" {print $1}' "$MANIFEST" | sort -u > "$KEEP"
N_WANT=$(wc -l < "$KEEP")
[[ "$N_WANT" -gt 0 ]] || { echo "ERROR: manifest lists no samples" >&2; exit 1; }

# Every manifest sample must be in the panel. A missing one means the manifest
# was built against a different panel than this, and a silent -S --force-samples
# would quietly ship a smaller cohort than the manifest claims.
bcftools query -l "$PANEL" | sort -u > "${OUTDIR}/.panel_${CHR}.txt"
MISSING=$(comm -23 "$KEEP" "${OUTDIR}/.panel_${CHR}.txt" | wc -l)
if [[ "$MISSING" -gt 0 ]]; then
    echo "ERROR: ${MISSING} of ${N_WANT} manifest samples are not in ${PANEL}:" >&2
    comm -23 "$KEEP" "${OUTDIR}/.panel_${CHR}.txt" | head -10 >&2
    echo "       The manifest and the panel disagree. Rebuild the manifest from
       this panel's Q table rather than shipping a partial cohort." >&2
    exit 4
fi

OUT="${OUTDIR}/benchmark_real.${CHR}.bcf"
if [[ -s "$OUT" && -s "${OUT}.csi" ]]; then
    N_HAVE=$(bcftools query -l "$OUT" | wc -l)
    if [[ "$N_HAVE" -eq "$N_WANT" ]]; then
        echo "[$(date +%T)] [${CHR}] already done (${N_HAVE} samples), skipping"
        rm -f "$KEEP" "${OUTDIR}/.panel_${CHR}.txt"
        exit 0
    fi
    echo "[$(date +%T)] [${CHR}] STALE: has ${N_HAVE} samples, manifest wants ${N_WANT} -- rebuilding"
    rm -f "$OUT" "${OUT}.csi"
fi

echo "[$(date +%T)] [${CHR}] extracting ${N_WANT} real test genomes from ${N_PANEL}-sample panel"
bcftools view -S "$KEEP" --force-samples=0 -Ob -o "$OUT" "$PANEL" 2>/dev/null \
  || bcftools view -S "$KEEP" -Ob -o "$OUT" "$PANEL"
bcftools index -f "$OUT"

N_OUT=$(bcftools query -l "$OUT" | wc -l)
if [[ "$N_OUT" -ne "$N_WANT" ]]; then
    echo "ERROR: wrote ${N_OUT} samples but the manifest lists ${N_WANT}" >&2
    exit 5
fi
rm -f "$KEEP" "${OUTDIR}/.panel_${CHR}.txt"
echo "[$(date +%T)] [${CHR}] DONE -> ${OUT} (${N_OUT} samples)"
