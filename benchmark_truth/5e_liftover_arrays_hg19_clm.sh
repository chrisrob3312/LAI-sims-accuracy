#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5e_bench_lift
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5e_bench_lift_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5e_bench_lift_chr%a_%j.err

# ============================================================================
# DESCRIPTION
# ============================================================================
# Deliver the GRCh37 arm of the benchmark: for every cohort whose `build` in
# cohorts.tsv is 37 (or both), lift its in-silico ARRAY files hg38 -> hg19 with
# Picard LiftoverVcf (same recipe as 5a_liftover_panel_hg38_to_hg19_clm.sh),
# and write hg19 PLINK sets. The TRUTH stays on GRCh38: the genotyping
# pipeline is expected to harmonize everything to GRCh38 before comparison,
# which is the point of having mixed-build inputs.
#
# Cohorts with build 38 are untouched; 'both' emits both.
# Outputs under $BENCH_DIR/arrays_hg19/<PLATFORM>/:
#   <COHORT>.<PLATFORM>.chr${CHR}.hg19.vcf.gz{,.tbi}  (+ .rejected.vcf.gz)
#   <COHORT>.<PLATFORM>.chr${CHR}.hg19.{bed,bim,fam}
# ============================================================================
CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
BENCH_DIR="${BENCH_DIR:-${PROJECT_ROOT}/06_benchmark_truth}"
COHORTS="${COHORTS:-${REPO_DIR}/benchmark_truth/cohorts.tsv}"
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg38ToHg19.over.chain.gz}"
HG19_FA="${HG19_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg19/hg19.fa}"
PICARD_XMX="${PICARD_XMX:-16g}"
RENAME_TO_ENSEMBL="${RENAME_TO_ENSEMBL:-1}"   # 'chr1' -> '1' after lift (arrays from labs are usually '1')

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?array job}
THREADS=${SLURM_CPUS_PER_TASK:-4}
IN="${BENCH_DIR}/arrays_hg38"; OUT="${BENCH_DIR}/arrays_hg19"; TMP="${OUT}/tmp/chr${CHR}"; mkdir -p "$TMP"
[[ -s "$CHAIN" ]] || { echo "ERROR: missing $CHAIN"; exit 1; }
[[ -s "$HG19_FA" && -s "${HG19_FA%.fa}.dict" ]] || { echo "ERROR: missing $HG19_FA or its .dict"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH; export PATH="${CONDA_PREFIX}/bin:${PATH}"

# Picard needs 'chr' contig names to match the UCSC chain; the panel uses '1..22'.
RENAME_UP="${TMP}/up.tsv"; RENAME_DOWN="${TMP}/down.tsv"
for c in {1..22} X Y; do echo -e "${c}\tchr${c}"; done > "$RENAME_UP"
for c in {1..22} X Y; do echo -e "chr${c}\t${c}"; done > "$RENAME_DOWN"

lift_one () {   # cohort platform
    local cohort="$1" plat="$2"
    local src="${IN}/${plat}/${cohort}.${plat}.chr${CHR}.vcf.gz"
    local d="${OUT}/${plat}"; mkdir -p "$d"
    local base="${d}/${cohort}.${plat}.chr${CHR}.hg19"
    [[ -s "${base}.bed" ]] && { echo "[$(date +%T)] [chr${CHR}] [$cohort/$plat] hg19 exists"; return 0; }
    [[ -s "$src" ]] || { echo "ERROR: missing $src (run 5d)"; exit 2; }
    local up="${TMP}/${cohort}.${plat}.up.vcf.gz" lifted="${TMP}/${cohort}.${plat}.lifted.vcf.gz"
    bcftools annotate --rename-chrs "$RENAME_UP" -Oz -o "$up" "$src"; bcftools index -t -f "$up"
    picard "-Xmx${PICARD_XMX}" LiftoverVcf I="$up" O="$lifted" CHAIN="$CHAIN" R="$HG19_FA" \
        REJECT="${base}.rejected.vcf.gz" RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true \
        CREATE_INDEX=false DISABLE_SORT=true MAX_RECORDS_IN_RAM=100000
    if [[ "$RENAME_TO_ENSEMBL" == "1" ]]; then
        bcftools annotate --rename-chrs "$RENAME_DOWN" -Oz -o "${lifted}.r" "$lifted" && mv "${lifted}.r" "$lifted"
    fi
    bcftools sort -m 4G -T "${TMP}/sort" -Oz -o "${base}.vcf.gz" "$lifted"; bcftools index -t -f "${base}.vcf.gz"
    plink2 --vcf "${base}.vcf.gz" --double-id --make-bed --threads "$THREADS" --out "$base" > "${base}.plink.log" 2>&1
    rm -f "$up" "${up}.tbi" "$lifted"
    echo "[$(date +%T)] [chr${CHR}] [$cohort/$plat] lifted $(bcftools view -H "${base}.vcf.gz" | wc -l) sites; rejected $(zcat "${base}.rejected.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)"
}

tail -n +2 "$COHORTS" | while IFS=$'\t' read -r cohort ctype n pops props gen build platforms notes; do
    [[ -z "$cohort" || "$cohort" == \#* ]] && continue
    [[ "$build" == "37" || "$build" == "both" ]] || continue
    IFS=';' read -r -a plats <<< "$platforms"
    for plat in "${plats[@]}"; do lift_one "$cohort" "$plat"; done
done
echo "[$(date +%T)] [chr${CHR}] # Complete. GRCh37 arrays in $OUT"
