#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=5b_bench_sim
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=00:30:00
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5b_bench_sim_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/5b_bench_sim_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Simulate the benchmarking truth-set cohorts (benchmark_truth/cohorts.tsv)
# from HGDP + MXB donor haplotypes with admix-simu, one chromosome per array
# task, every cohort in one job. Same machinery as 2b_simulation_clm.sh, but
# K-way (any number of source populations) and driven by a cohort table:
#
#   homog cohorts  : per constituent POPULATION of the superpop, run admix-simu
#                    with K=1 (within-population splicing), sizes proportional
#                    to donor counts, then pool to the superpop cohort.
#   admixed cohorts: superpop-level donor pools (all included donors of that
#                    superpop, pooled) and a K-way .dat with the given
#                    proportions and generations.
#
# Outputs per chr under $BENCH_DIR/sim/<COHORT>/:
#   <COHORT>.chr${CHR}.haps / .sample   simulated haplotypes (SHAPEIT format)
#   <COHORT>${CHR}.bp / .hanc2          local-ancestry TRUTH (admixed cohorts;
#                                       for homog cohorts it is trivially 1 pop)
#   <COHORT>.chr${CHR}.dat              the admix-simu config used
# Shared per chr under $BENCH_DIR/sim/_donors/:
#   <SUPERPOP>.<POP>.chr${CHR}.phgeno   donor haplotype matrices (cached)
#   bench.chr${CHR}.snp                 site list with cM (insert-map.pl)
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.rechr.bcf}"
BENCH_DIR="${BENCH_DIR:-${PROJECT_ROOT}/06_benchmark_truth}"
DONOR_DIR="${DONOR_DIR:-${BENCH_DIR}/donors}"
COHORTS="${COHORTS:-${REPO_DIR}/benchmark_truth/cohorts.tsv}"
ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/simu-jointcall/admix-simu-master}"
GMAP_TPL="${GMAP_TPL:-${REPO_DIR}/resources/gmap_hg38/chr%s.hg38.gmap.txt}"
SEED="${SEED:-20261005}"          # admix-simu takes no seed; recorded for the README
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-4}
# shellcheck disable=SC2059
PANEL=$(printf "$PANEL_TPL" "$CHR")
# shellcheck disable=SC2059
GMAP=$(printf "$GMAP_TPL" "$CHR")
[[ -s "$PANEL"         ]] || { echo "ERROR: missing $PANEL"; exit 1; }
[[ -s "$GMAP"          ]] || { echo "ERROR: missing $GMAP"; exit 1; }
[[ -d "$ADMIXSIMU_DIR" ]] || { echo "ERROR: missing $ADMIXSIMU_DIR"; exit 1; }
[[ -s "$COHORTS"       ]] || { echo "ERROR: missing $COHORTS"; exit 1; }
[[ -s "${DONOR_DIR}/donors.tsv" ]] || { echo "ERROR: missing ${DONOR_DIR}/donors.tsv (run 5a first)"; exit 1; }

SIMDIR="${BENCH_DIR}/sim"
DDIR="${SIMDIR}/_donors"
mkdir -p "$SIMDIR" "$DDIR" "$LOGDIR"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate "$CONDA_ENV"; set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ----------------------------------------------------------------------------
# 1. Site list with cM for this chr (once), from the whole panel.
# ----------------------------------------------------------------------------
SNP="${DDIR}/bench.chr${CHR}.snp"
if [[ ! -s "$SNP" ]]; then
    echo "[$(date +%T)] [chr${CHR}] building .snp"
    plink2 --bcf "$PANEL" --const-fid 0 --max-alleles 2 --snps-only just-acgt \
           --make-bed --threads "$THREADS" --out "${DDIR}/panel_${CHR}" \
           > "${DDIR}/panel_${CHR}.plink.log" 2>&1
    perl "${ADMIXSIMU_DIR}/insert-map.pl" "${DDIR}/panel_${CHR}.bim" "$GMAP" > "${DDIR}/chr${CHR}.pos"
    awk -F' ' '{ print $1":"$4"_"$5"_"$6, $1, $3, $4, $5, $6 }' "${DDIR}/chr${CHR}.pos" > "$SNP"
    # sites actually used (biallelic SNPs) - every phgeno below must use the same set
    awk '{print $2":"$4}' "$SNP" > "${DDIR}/bench.chr${CHR}.sites"
fi
NSITES=$(wc -l < "$SNP")

# ----------------------------------------------------------------------------
# 2. Donor haplotype matrices (.phgeno), cached per (superpop, population).
#    Restrict to the .snp site set so all phgeno rows align.
# ----------------------------------------------------------------------------
make_phgeno () {   # label keepfile -> ${DDIR}/${label}.chr${CHR}.phgeno
    local label="$1" keep="$2" out="${DDIR}/${1}.chr${CHR}.phgeno"
    [[ -s "$out" ]] && return 0
    plink2 --bcf "$PANEL" --keep <(awk '{print 0, $1}' "$keep") \
           --extract <(awk '{print $2":"$4}' "$SNP") --set-all-var-ids '@:#' \
           --export haps --threads "$THREADS" --out "${DDIR}/${label}.chr${CHR}" \
           > "${DDIR}/${label}.chr${CHR}.plink.log" 2>&1
    cut -d' ' -f 6- "${DDIR}/${label}.chr${CHR}.haps" | sed 's/\s//g' > "$out"
    [[ $(wc -l < "$out") -eq "$NSITES" ]] || { echo "ERROR: $label phgeno has $(wc -l < "$out") rows, expected $NSITES"; exit 1; }
    rm -f "${DDIR}/${label}.chr${CHR}.haps"
}

# superpop pooled donor list (all included donors of a superpop)
pooled_keep () {   # superpop -> file
    local sp="$1" f="${DDIR}/pool.${sp}.txt"
    [[ -s "$f" ]] || awk -F'\t' -v sp="$sp" 'NR>1 && $3==sp && $6==1 {print $1}' "${DONOR_DIR}/donors.tsv" > "$f"
    echo "$f"
}

# ----------------------------------------------------------------------------
# 3. admix-simu driver for one (cohort, tag, pops, props, gen, n)
# ----------------------------------------------------------------------------
run_sim () {
    local cohort="$1" tag="$2" n="$3" gen="$4"; shift 4
    local pops=("$@")                       # labels; each has ${DDIR}/${label}.chr${CHR}.phgeno
    local props="${PROPS_CSV}"              # set by caller
    local wd="${SIMDIR}/${cohort}"; mkdir -p "$wd"
    local prefix="${wd}/${tag}${CHR}"
    local dat="${wd}/${tag}.chr${CHR}.dat"
    IFS=';' read -r -a parr <<< "$props"
    python3 "${REPO_DIR}/benchmark_truth/make_dat.py" --n "$n" --pops "${pops[@]}" --props "${parr[@]}" --gen "$gen" --out "$dat"
    local args=()
    for p in "${pops[@]}"; do args+=("-${p}" "${DDIR}/${p}.chr${CHR}.phgeno"); done
    echo "[$(date +%T)] [chr${CHR}] [$cohort/$tag] simu-mix.pl n=$n gen=$gen pops=${pops[*]} props=$props"
    perl "${ADMIXSIMU_DIR}/simu-mix.pl" "$dat" "$SNP" "$prefix" "${args[@]}"
    perl "${ADMIXSIMU_DIR}/bp2anc.pl" "${prefix}.bp" > "${prefix}.hanc"
    sed 's/./& /g' "${prefix}.hanc" > "${prefix}.hanc1"
    python3 -c "import sys; rows=[l.split() for l in sys.stdin if l.strip()]; print('\n'.join(' '.join(c) for c in zip(*rows)))" \
        < "${prefix}.hanc1" > "${prefix}.hanc2"
    rm -f "${prefix}.hanc1"
    # simulated haps (SHAPEIT): 5 site cols + 2 cols per individual
    sed 's/./& /g' "${prefix}.phgeno" > "${prefix}.tmp_geno"
    awk -F' ' '{ print $2, $1, $4, $5, $6 }' "$SNP" > "${prefix}.tmp_snp"
    paste -d' ' "${prefix}.tmp_snp" "${prefix}.tmp_geno" > "${wd}/${tag}.chr${CHR}.haps"
    rm -f "${prefix}.tmp_geno" "${prefix}.tmp_snp"
    { echo "ID_1 ID_2 missing"; echo "0 0 0"; for i in $(seq 1 "$n"); do echo "SIM_${tag}_${i} SIM_${tag}_${i} 0"; done; } \
        > "${wd}/${tag}.chr${CHR}.sample"
}

# ----------------------------------------------------------------------------
# 4. Loop over cohorts.tsv
# ----------------------------------------------------------------------------
tail -n +2 "$COHORTS" | while IFS=$'\t' read -r cohort ctype n pops props gen build platforms notes; do
    [[ -z "$cohort" || "$cohort" == \#* ]] && continue
    wd="${SIMDIR}/${cohort}"
    if [[ -s "${wd}/${cohort}.chr${CHR}.haps" ]]; then
        echo "[$(date +%T)] [chr${CHR}] [$cohort] exists, skipping"; continue
    fi
    if [[ "$ctype" == "homog" ]]; then
        sp="$pops"
        # constituent populations with included donors, sizes proportional to donor counts
        mapfile -t popfiles < <(ls "${DONOR_DIR}"/donors."${sp}".*.txt 2>/dev/null)
        [[ ${#popfiles[@]} -gt 0 ]] || { echo "ERROR: no donor lists for superpop $sp"; exit 1; }
        total=0; declare -A cnt=()
        for f in "${popfiles[@]}"; do c=$(wc -l < "$f"); cnt["$f"]=$c; total=$((total + c)); done
        parts=()
        for f in "${popfiles[@]}"; do
            pop=$(basename "$f" .txt); pop=${pop#donors.}; label="${pop}"        # e.g. AFR.Yoruba
            ni=$(( (n * cnt["$f"] + total / 2) / total )); [[ $ni -lt 2 ]] && ni=2
            make_phgeno "$label" "$f"
            PROPS_CSV="1" run_sim "$cohort" "${cohort}__${label}" "$ni" "$gen" "$label"
            parts+=("${wd}/${cohort}__${label}.chr${CHR}")
        done
        # pool: haps columns side by side (same site rows), samples concatenated
        python3 - "${wd}/${cohort}.chr${CHR}" "${parts[@]}" <<'PYEOF'
import sys
out, parts = sys.argv[1], sys.argv[2:]
hs = [open(p + ".haps") for p in parts]
with open(out + ".haps", "w") as g:
    for lines in zip(*hs):
        first = lines[0].rstrip("\n").split(" ")
        rest = [l.rstrip("\n").split(" ")[5:] for l in lines[1:]]
        g.write(" ".join(first + sum(rest, [])) + "\n")
ids = []
for p in parts:
    ids += [l.split()[0] for l in open(p + ".sample").read().splitlines()[2:]]
with open(out + ".sample", "w") as g:
    g.write("ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{i} {i} 0\n" for i in ids))
PYEOF
        unset cnt
    else
        IFS=';' read -r -a parr <<< "$pops"
        for sp in "${parr[@]}"; do make_phgeno "$sp" "$(pooled_keep "$sp")"; done
        PROPS_CSV="$props" run_sim "$cohort" "$cohort" "$n" "$gen" "${parr[@]}"
    fi
done

echo "[$(date +%T)] [chr${CHR}] # Complete. Outputs in $SIMDIR  (seed note: $SEED)"
