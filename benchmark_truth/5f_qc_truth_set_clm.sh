#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5f_bench_qc
#SBATCH --partition=atkinson,mhgcp
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5f_bench_qc_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5f_bench_qc_%j.err

# ============================================================================
# DESCRIPTION  (one-shot, after 5b-5e)
# ============================================================================
# Does the simulated truth set look like the design?
#   1. counts: samples per cohort, sites per platform per chr, flips/errors
#      injected (from *.changes.tsv), rejected liftover sites
#   2. structure: plink2 PCA on LD-pruned merged truth -> qc/pca.eigenvec;
#      qc_truth_set.R plots PC1-4 by cohort (homogeneous cohorts must cluster,
#      admixed cohorts must sit between their sources at the designed mix)
#   3. manifest: qc/MANIFEST.tsv listing every deliverable file with md5
#
# Outputs: $BENCH_DIR/qc/
# ============================================================================
CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
BENCH_DIR="${BENCH_DIR:-${PROJECT_ROOT}/06_benchmark_truth}"

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-8}
QC="${BENCH_DIR}/qc"; mkdir -p "$QC"
module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH; export PATH="${CONDA_PREFIX}/bin:${PATH}"

# 1. counts
echo -e "cohort\tn_samples" > "${QC}/samples_per_cohort.tsv"
cut -f2 "${BENCH_DIR}/truth_hg38/ALL.chr22.samples.tsv" | sort | uniq -c | awk '{print $2"\t"$1}' >> "${QC}/samples_per_cohort.tsv"
echo -e "platform\tcohort\tchr\tn_sites\tn_flip\tn_palflip\tn_palkept\tsum_err\tsum_mis" > "${QC}/array_sites.tsv"
for f in "${BENCH_DIR}"/arrays_hg38/*/*.changes.tsv; do
    plat=$(basename "$(dirname "$f")"); b=$(basename "$f" .changes.tsv); cohort=${b%%.*}; chr=${b##*.chr}
    vcf="${f%.changes.tsv}.vcf.gz"; n=$(bcftools view -H "$vcf" | wc -l)
    awk -F'\t' -v p="$plat" -v c="$cohort" -v ch="$chr" -v n="$n" 'NR>1{ if($2=="flip")f++; if($2=="palflip")pf++; if($2=="palkept")pk++;
         split($3,a,";"); split(a[1],e,"="); split(a[2],m,"="); err+=e[2]; mis+=m[2] }
         END{print p"\t"c"\t"ch"\t"n"\t"f+0"\t"pf+0"\t"pk+0"\t"err+0"\t"mis+0}' "$f" >> "${QC}/array_sites.tsv"
done
if compgen -G "${BENCH_DIR}/arrays_hg19/*/*.rejected.vcf.gz" > /dev/null; then
    echo -e "file\tn_rejected" > "${QC}/liftover_rejects.tsv"
    for f in "${BENCH_DIR}"/arrays_hg19/*/*.rejected.vcf.gz; do echo -e "$(basename "$f")\t$(zcat "$f" | grep -vc '^#' || echo 0)" >> "${QC}/liftover_rejects.tsv"; done
fi

# 2. PCA on merged truth (all chrs), LD-pruned
LIST="${QC}/truth_concat_list.txt"; : > "$LIST"
for c in {1..22}; do echo "${BENCH_DIR}/truth_hg38/ALL.chr${c}.truth.vcf.gz" >> "$LIST"; done
bcftools concat --threads "$THREADS" -f "$LIST" -Ob -o "${QC}/ALL.truth.bcf"; bcftools index -f "${QC}/ALL.truth.bcf"
plink2 --bcf "${QC}/ALL.truth.bcf" --set-missing-var-ids '@:#' --maf 0.01 --indep-pairwise 500 50 0.1 --threads "$THREADS" --out "${QC}/prune"
plink2 --bcf "${QC}/ALL.truth.bcf" --set-missing-var-ids '@:#' --extract "${QC}/prune.prune.in" --pca 10 approx --threads "$THREADS" --out "${QC}/pca"
Rscript "${REPO_DIR}/benchmark_truth/qc_truth_set.R" "${QC}/pca.eigenvec" "${BENCH_DIR}/truth_hg38/ALL.chr22.samples.tsv" "${QC}"

# 3. manifest
( cd "$BENCH_DIR" && find truth_hg38 arrays_hg38 arrays_hg19 -type f \( -name '*.vcf.gz' -o -name '*.bed' -o -name '*.bim' -o -name '*.fam' \) 2>/dev/null \
    | sort | xargs -r md5sum | awk '{print $2"\t"$1}' ) > "${QC}/MANIFEST.tsv"
echo "[$(date +%T)] QC done: $(wc -l < "${QC}/MANIFEST.tsv") deliverable files listed in ${QC}/MANIFEST.tsv"
