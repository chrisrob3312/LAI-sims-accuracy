#!/usr/bin/env bash
# ============================================================================
# 6e_cohort_kinship_clm.sh  --  genome-wide kinship over the COMBINED cohort
#
#SBATCH --job-name=6e_kinship
#SBATCH --partition=atkinson,mhgcp
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/6e_kinship_%j.out
#SBATCH --error=logs/6e_kinship_%j.err
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#
# NOT AN ARRAY. Kinship is a genome-wide quantity: per-chromosome estimates
# cannot be averaged into one, and a per-chromosome run would under-call every
# relationship. One job over all 22.
#
# WHY THE COHORT NEEDS THIS AT ALL. ~500 real genomes is too few for a power
# calculation, so the simulated 2,160 are used for power and GWAS simulations
# too. But simulated individuals built from a shared donor pool are cryptic
# relatives by construction -- SAS at ~24 donors for 200 SASEUR individuals is
# the binding case -- so the roster count overstates the information in the
# cohort. cohort-power-summary.py turns this table into an effective sample
# size per group, which is the number a power calculation should use.
#
# This is a DIFFERENT run from 1e_king_dedup_clm.sh. That one de-duplicates
# PEOPLE in the reference panel before ADMIXTURE and drops only pairs above
# 0.354. This one measures relatedness across the finished cohort, simulated
# individuals included, and drops nothing: the relatedness is real, GENESIS
# models it with a GRM, and the point is to quantify what it costs rather than
# to remove anybody.
#
# It also catches a leak the other checks cannot: a REAL test genome related
# to a SIMULATED one means a donor, or a donor's close relative, is in the
# test set.
#
# Chip and build do not enter here. Kinship is a property of the PEOPLE, so
# it is computed once over the distinct individuals whatever design put them
# in which cells -- and under the default paired-builds design an individual
# appears in two cells, so filtering by cell would count some twice.
# ============================================================================
set -euo pipefail

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REAL_DIR="${REAL_DIR:-${PROJECT_ROOT}/07_benchmark_real}"
SIMSET_DIR="${SIMSET_DIR:-${PROJECT_ROOT}/06_benchmark_simset}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/09_cohort_kinship}"
BENCH_OUT="${BENCH_OUT:-${REPO_DIR}/results}"
SIM_MANIFEST="${SIM_MANIFEST:-${BENCH_OUT}/benchmark_cohort_simset_individuals.tsv}"
REAL_MANIFEST="${REAL_MANIFEST:-${BENCH_OUT}/benchmark_cohort_manifest.tsv}"
CONDA_ENV="${CONDA_ENV:-shapeit5}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MAF="${MAF:-0.05}"
PRUNE="${PRUNE:-200 50 0.2}"
KING_MIN="${KING_MIN:-0.0442}"   # third degree; below this is noise here

if [[ -z "${CONDA_PREFIX:-}" || "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
for t in bcftools plink2; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH in ${CONDA_ENV}" >&2; exit 127; }
done
[[ -s "$SIM_MANIFEST" ]] || { echo "ERROR: missing $SIM_MANIFEST (run the assembler)" >&2; exit 1; }

WORK="${OUTDIR}/.work"; rm -rf "$WORK"; mkdir -p "$WORK" "$OUTDIR"; cd "$WORK"

# --- per chromosome: real + every simulated model, merged ------------------
MERGED=()
for c in $(seq 1 22); do
    CHR="chr${c}"
    parts=()
    real="${REAL_DIR}/benchmark_real.${CHR}.bcf"
    if [[ -s "$real" ]]; then
        bcftools view -Oz -o "real.${CHR}.vcf.gz" "$real"
        bcftools index -f -t "real.${CHR}.vcf.gz"
        parts+=("real.${CHR}.vcf.gz")
    fi
    while read -r model; do
        [[ -n "$model" ]] || continue
        haps="${SIMSET_DIR}/${CHR}/${model}.${CHR}.haps"
        samp="${SIMSET_DIR}/${CHR}/${model}.${CHR}.sample"
        [[ -s "$haps" && -s "$samp" ]] || continue
        n_want=$(awk -F'\t' -v m="$model" 'NR>1 && $2==m' "$SIM_MANIFEST" | wc -l)
        n_have=$(awk 'NR==1{print int((NF-5)/2); exit}' "$haps")
        if [[ "$n_want" -ne "$n_have" ]]; then
            echo "ERROR: ${model} ${CHR}: .haps has ${n_have} individuals, roster says
       ${n_want}. reheader renames POSITIONALLY, so kinship would be reported
       under the wrong ids." >&2
            exit 4
        fi
        awk -F'\t' -v m="$model" 'NR>1 && $2==m {print $1}' "$SIM_MANIFEST" > "${model}.names"
        plink2 --haps "$haps" --sample "$samp" --export vcf bgz \
               --out "${model}.${CHR}.raw" --threads "$THREADS" --silent
        bcftools reheader -s "${model}.names" \
                 -o "${model}.${CHR}.vcf.gz" "${model}.${CHR}.raw.vcf.gz"
        bcftools index -f -t "${model}.${CHR}.vcf.gz"
        parts+=("${model}.${CHR}.vcf.gz")
    done < <(awk -F'\t' 'NR>1{print $2}' "$SIM_MANIFEST" | sort -u)

    [[ ${#parts[@]} -gt 0 ]] || { echo "[$(date +%T)] ${CHR}: no inputs, skipping"; continue; }
    if [[ ${#parts[@]} -eq 1 ]]; then
        cp "${parts[0]}" "cohort.${CHR}.vcf.gz"
    else
        # -0: a site missing from one part is reference there, not missing.
        bcftools merge -0 --threads "$THREADS" -Oz -o "cohort.${CHR}.vcf.gz" "${parts[@]}"
    fi
    bcftools index -f -t "cohort.${CHR}.vcf.gz"
    MERGED+=("cohort.${CHR}.vcf.gz")
    echo "[$(date +%T)] ${CHR}: merged $(( ${#parts[@]} )) part(s)"
done

[[ ${#MERGED[@]} -gt 0 ]] || { echo "ERROR: nothing merged; run 6b and 6c first" >&2; exit 1; }

echo "[$(date +%T)] concatenating ${#MERGED[@]} chromosomes"
bcftools concat --threads "$THREADS" -Oz -o cohort.all.vcf.gz "${MERGED[@]}"
bcftools index -f -t cohort.all.vcf.gz

N_SAMP=$(bcftools query -l cohort.all.vcf.gz | wc -l)
N_EXP=$(( $(awk 'NR>1' "$SIM_MANIFEST" | wc -l) \
          + $( [[ -s "$REAL_MANIFEST" ]] && awk 'NR>1' "$REAL_MANIFEST" | wc -l || echo 0 ) ))
echo "[$(date +%T)] cohort: ${N_SAMP} samples (rosters list ${N_EXP})"

# --- LD-prune, then KING --------------------------------------------------
# Kinship on LD-correlated markers is inflated; pruning is not optional.
echo "[$(date +%T)] MAF >= ${MAF}, LD prune ${PRUNE}"
plink2 --vcf cohort.all.vcf.gz --maf "$MAF" --snps-only --max-alleles 2 \
       --indep-pairwise $PRUNE --out prune --threads "$THREADS" --silent
echo "[$(date +%T)] $(wc -l < prune.prune.in) markers kept"

plink2 --vcf cohort.all.vcf.gz --extract prune.prune.in \
       --make-king-table --king-table-filter "$KING_MIN" \
       --out "${OUTDIR}/cohort_king" --threads "$THREADS" --silent

KT="${OUTDIR}/cohort_king.kin0"
[[ -s "$KT" ]] || { echo "ERROR: plink2 wrote no king table at ${KT}" >&2; exit 5; }
echo "[$(date +%T)] $(( $(wc -l < "$KT") - 1 )) pair(s) at kinship >= ${KING_MIN}"

SUM="${REPO_DIR}/scripts/cohort-power-summary.py"
if [[ -s "$SUM" ]]; then
    python3 -I "$SUM" --kinship "$KT" \
        --manifest "$REAL_MANIFEST" --individuals "$SIM_MANIFEST" \
        --out "${BENCH_OUT}/cohort_power_summary.txt" || rc=$?
    # exit 2 from the summary means a real-simulated leak, which is a finding
    # to surface loudly, not a reason to lose the kinship table we just spent
    # hours computing.
    if [[ "${rc:-0}" -eq 2 ]]; then
        echo "[$(date +%T)] *** real-simulated relatedness detected -- see"
        echo "[$(date +%T)] *** ${BENCH_OUT}/cohort_power_summary.txt"
    fi
fi

rm -rf "$WORK"
echo "[$(date +%T)] DONE -> ${KT}"
echo "[$(date +%T)]         ${BENCH_OUT}/cohort_power_summary.txt"
