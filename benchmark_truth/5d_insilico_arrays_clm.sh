#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5d_bench_array
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5d_bench_array_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5d_bench_array_chr%a_%j.err

# ============================================================================
# DESCRIPTION
# ============================================================================
# Make in-silico genotyping-array files (SNP6 and GSA) from the truth VCFs:
# site restriction + strand reporting + palindromic scramble + genotype error
# + missingness + unphasing (see insilico_array.py). One file set per
# (cohort, platform) per chr, only for the platforms listed for that cohort in
# cohorts.tsv. Also emits PLINK bed/bim/fam per platform (what most labs hand
# the genotyping pipeline).
#
# Site lists (chr:pos, hg38), one file per chr:
#   GSA : reference_ids/chip_gsa_jessica/chr%s.snps   (already in repo)
#   SNP6: reference_ids/chip_snp6/chr%s.snps          (build once with
#         scripts/prepare-chip-snps.sh from the Affymetrix GenomeWideSNP_6
#         hg38 annotation CSV: --manifest <affy csv> --name snp6 ; the script
#         looks for Chr + MapInfo columns - for the Affy CSV rename
#         "Chromosome" -> Chr and "Physical Position" -> MapInfo first, or
#         pass a 2-column TSV you make yourself)
#
# Outputs under $BENCH_DIR/arrays_hg38/<PLATFORM>/:
#   <COHORT>.<PLATFORM>.chr${CHR}.vcf.gz{,.tbi}  + .changes.tsv (ground truth of every modification)
#   <COHORT>.<PLATFORM>.chr${CHR}.{bed,bim,fam}
# ============================================================================
CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
BENCH_DIR="${BENCH_DIR:-${PROJECT_ROOT}/06_benchmark_truth}"
COHORTS="${COHORTS:-${REPO_DIR}/benchmark_truth/cohorts.tsv}"
GSA_TPL="${GSA_TPL:-${REPO_DIR}/reference_ids/chip_gsa_jessica/chr%s.snps}"
SNP6_TPL="${SNP6_TPL:-${REPO_DIR}/reference_ids/chip_snp6/chr%s.snps}"
# platform error models (documented in the README; change here, not in the .py)
GSA_FLIP="${GSA_FLIP:-0.50}";  GSA_ERR="${GSA_ERR:-0.002}";  GSA_MISS="${GSA_MISS:-0.010}"
SNP6_FLIP="${SNP6_FLIP:-0.35}"; SNP6_ERR="${SNP6_ERR:-0.004}"; SNP6_MISS="${SNP6_MISS:-0.020}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?array job}
THREADS=${SLURM_CPUS_PER_TASK:-2}
TRUTH="${BENCH_DIR}/truth_hg38"; OUT="${BENCH_DIR}/arrays_hg38"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH; export PATH="${CONDA_PREFIX}/bin:${PATH}"

make_array () {   # cohort platform sites flip err miss
    local cohort="$1" plat="$2" sites="$3" flip="$4" err="$5" miss="$6"
    local d="${OUT}/${plat}"; mkdir -p "$d"
    local base="${d}/${cohort}.${plat}.chr${CHR}"
    [[ -s "${base}.bed" ]] && { echo "[$(date +%T)] [chr${CHR}] [$cohort/$plat] exists"; return 0; }
    [[ -s "$sites" ]] || { echo "ERROR: missing site list $sites for $plat"; exit 1; }
    python3 "${REPO_DIR}/benchmark_truth/insilico_array.py" \
        --vcf "${TRUTH}/${cohort}.chr${CHR}.truth.vcf.gz" --sites "$sites" --out "$base" \
        --flip-frac "$flip" --error "$err" --missing "$miss" --seed $((CHR * 1000 + ${#cohort} * 7 + ${#plat}))
    bcftools index -t -f "${base}.vcf.gz"
    plink2 --vcf "${base}.vcf.gz" --double-id --make-bed --threads "$THREADS" --out "$base" \
        > "${base}.plink.log" 2>&1
}

tail -n +2 "$COHORTS" | while IFS=$'\t' read -r cohort ctype n pops props gen build platforms notes; do
    [[ -z "$cohort" || "$cohort" == \#* ]] && continue
    IFS=';' read -r -a plats <<< "$platforms"
    for plat in "${plats[@]}"; do
        case "$plat" in
            gsa)  make_array "$cohort" gsa  "$(printf "$GSA_TPL" "$CHR")"  "$GSA_FLIP"  "$GSA_ERR"  "$GSA_MISS" ;;
            snp6) make_array "$cohort" snp6 "$(printf "$SNP6_TPL" "$CHR")" "$SNP6_FLIP" "$SNP6_ERR" "$SNP6_MISS" ;;
            *) echo "WARN: unknown platform '$plat' for $cohort - skipped" ;;
        esac
    done
done
echo "[$(date +%T)] [chr${CHR}] # Complete. Outputs in $OUT"
