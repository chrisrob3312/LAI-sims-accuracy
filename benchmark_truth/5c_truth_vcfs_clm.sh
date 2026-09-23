#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5c_bench_vcf
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5c_bench_vcf_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5c_bench_vcf_chr%a_%j.err

# ============================================================================
# DESCRIPTION
# ============================================================================
# Convert every simulated cohort's SHAPEIT .haps/.sample (from 5b) into
# phased, indexed GRCh38 "WGS truth" VCFs, one per cohort per chr, plus a
# merged all-cohort truth VCF per chr. Contigs are '1..22' (the panel's
# 'rechr' convention); 5e renames to 'chr1..' where a platform needs it.
#
# Outputs under $BENCH_DIR/truth_hg38/:
#   <COHORT>.chr${CHR}.truth.vcf.gz{,.tbi}
#   ALL.chr${CHR}.truth.vcf.gz{,.tbi}        (bcftools merge of all cohorts)
#   ALL.chr${CHR}.samples.tsv                 sample_id  cohort
# ============================================================================
CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
BENCH_DIR="${BENCH_DIR:-${PROJECT_ROOT}/06_benchmark_truth}"
COHORTS="${COHORTS:-${REPO_DIR}/benchmark_truth/cohorts.tsv}"
HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg38/hg38.fa}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?array job}
THREADS=${SLURM_CPUS_PER_TASK:-4}
SIMDIR="${BENCH_DIR}/sim"; OUT="${BENCH_DIR}/truth_hg38"; mkdir -p "$OUT"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH; export PATH="${CONDA_PREFIX}/bin:${PATH}"

: > "${OUT}/ALL.chr${CHR}.samples.tsv"
LIST="${OUT}/ALL.chr${CHR}.merge_list.txt"; : > "$LIST"
tail -n +2 "$COHORTS" | while IFS=$'\t' read -r cohort rest; do
    [[ -z "$cohort" || "$cohort" == \#* ]] && continue
    haps="${SIMDIR}/${cohort}/${cohort}.chr${CHR}.haps"; sample="${SIMDIR}/${cohort}/${cohort}.chr${CHR}.sample"
    [[ -s "$haps" && -s "$sample" ]] || { echo "ERROR: missing $haps (run 5b)"; exit 1; }
    vcf="${OUT}/${cohort}.chr${CHR}.truth.vcf.gz"
    if [[ ! -s "$vcf" ]]; then
        echo "[$(date +%T)] [chr${CHR}] [$cohort] haps -> VCF"
        # haps allele columns are (REF, ALT) as exported from the panel; keep that order.
        plink2 --haps "$haps" --sample "$sample" --export vcf bgz id-paste=iid \
               --threads "$THREADS" --out "${OUT}/${cohort}.chr${CHR}.tmp" \
               > "${OUT}/${cohort}.chr${CHR}.plink.log" 2>&1
        # normalise against the reference so REF/ALT are guaranteed correct for truth
        bcftools norm -f "$HG38_FA" --check-ref ws -Oz -o "$vcf" "${OUT}/${cohort}.chr${CHR}.tmp.vcf.gz" 2>> "${OUT}/${cohort}.chr${CHR}.norm.log" \
            || { echo "WARN: bcftools norm --check-ref failed for $cohort (contig names?) - keeping unnormalised"; mv "${OUT}/${cohort}.chr${CHR}.tmp.vcf.gz" "$vcf"; }
        rm -f "${OUT}/${cohort}.chr${CHR}.tmp.vcf.gz"
        bcftools index -t -f "$vcf"
    fi
    echo "$vcf" >> "$LIST"
    bcftools query -l "$vcf" | awk -v c="$cohort" '{print $1"\t"c}' >> "${OUT}/ALL.chr${CHR}.samples.tsv"
done

echo "[$(date +%T)] [chr${CHR}] merging cohorts"
bcftools merge --threads "$THREADS" -l "$LIST" -Oz -o "${OUT}/ALL.chr${CHR}.truth.vcf.gz"
bcftools index -t -f "${OUT}/ALL.chr${CHR}.truth.vcf.gz"
echo "[$(date +%T)] [chr${CHR}] truth: $(bcftools query -l "${OUT}/ALL.chr${CHR}.truth.vcf.gz" | wc -l) samples, $(bcftools view -H "${OUT}/ALL.chr${CHR}.truth.vcf.gz" | wc -l) sites"
