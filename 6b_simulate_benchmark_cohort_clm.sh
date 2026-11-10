#!/usr/bin/env bash
# ============================================================================
# 6b_simulate_benchmark_cohort_clm.sh  --  build step 6: simulate the 2,160
#
#SBATCH --job-name=6b_simset
#SBATCH --partition=atkinson,mhgcp
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/6b_simset_%a_%j.out
#SBATCH --error=logs/6b_simset_%a_%j.err
# The t-nodes cannot execute these binaries: a prior array put 270/270 tasks
# on mhgcp-t01..t12 and every one died with SIGILL (exit 132) before writing
# a log. c02 is NOT excluded -- nothing ties it to that failure, and it has
# since run chunked phasing cleanly.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#
# ONE TASK PER CHROMOSOME, ALL 13 MODELS INSIDE. 2b extracted donor
# haplotypes once per simulated track, which meant re-running
# `plink2 --export haps` over a ~3,900-sample BCF for every model. The
# extraction dominates the runtime and the models share their component
# pools, so this extracts each component ONCE per chromosome and then loops
# the models. 22 tasks rather than 286, and the expensive step runs 6 times
# per chromosome instead of 30.
#
# WHAT IT READS
#   results/benchmark_cohort_simset_individuals.tsv   from
#       scripts/assemble-benchmark-cohort.py -- one row per simulated
#       individual, carrying its model and its chip x build cell. Step 10
#       masks against the SAME file, which is why the cell assignment lives
#       there and is not re-dealt here.
#   reference_ids/<pool>_simulation_<suffix>.txt      the benchmark donor half
#   01_merged_phased_panel/merged_chr<N>...bcf        the phased panel
#
# PATHS ARE ALL OVERRIDABLE. Nothing below assumes this repository's location,
# so the script runs unchanged from another working tree or cluster: set
# PROJECT_ROOT, PANEL_DIR, OUTDIR, REFS and ADMIXSIMU_DIR.
#
# REPRODUCIBILITY, and why it is preservation rather than seeding.
# admix-simu splits its randomness across two programs. The compiled `mixer`
# binary draws the recombination breakpoints into <prefix>.bp; it takes no
# seed and ships without source, so its RNG cannot be controlled. simu-mix.pl
# then draws which donor haplotype to copy, and scripts/simu-mix-seeded.pl
# makes that half seedable via SIMU_SEED. So every unit archives its .bp --
# the .bp file IS the recombination realisation, and simu-mix.pl's -bp mode
# replays from it exactly. Archiving .bp + SIMU_SEED makes a unit regenerable
# even though half its RNG is out of reach, which is stronger than seeding.
# ============================================================================
set -euo pipefail

# --dat-only inspects the roster and touches no chromosome, so it must not
# require a Slurm array context.
if [[ "${DAT_ONLY:-}" == "1" || "${1:-}" == "--dat-only" ]]; then
    CHR_N="${SLURM_ARRAY_TASK_ID:-1}"
else
    CHR_N="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job (1-22), or pass --dat-only}"
fi
CHR="chr${CHR_N}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_%s.shapeit5_phased.softunion_maf005.bcf}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/06_benchmark_simset}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"
BENCH_SUFFIX="${BENCH_SUFFIX:-bench}"
INDIV_TSV="${INDIV_TSV:-${REPO_DIR}/results/benchmark_cohort_simset_individuals.tsv}"
ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/simu-jointcall/admix-simu-master}"
SIMU_MIX="${SIMU_MIX:-${REPO_DIR}/scripts/simu-mix-seeded.pl}"
GMAP_DIR="${GMAP_DIR:-${REPO_DIR}/resources/gmap_hg38}"
CONDA_ENV="${CONDA_ENV:-shapeit5}"
FROM_BP="${FROM_BP:-}"
SEED_BASE="${SEED_BASE:-20261108}"

# NAT is simulated from the AMR donor pool; every other component's name is
# its own pool.
declare -A COMPONENT_POOL=(
  [NAT]=AMR [EUR]=EUR [AFR]=AFR [MEN]=MEN [EAS]=EAS [SAS]=SAS [OCE]=OCE
)

# --- the .dat file per model ----------------------------------------------
# admix-simu .dat format. Line 1 is the header that simu-mix.pl validates:
#
#     <n_samples> <ADMIXED_LABEL> <POP1> <POP2> ... <POPn>
#
# and each later line is a generation followed by the migrant fraction
# contributed by each population at that generation, in the header's order.
# A single founding pulse at generation 1 followed by G generations of drift
# is written as the pulse line plus a terminal all-zero line at G:
#
#     600 LAT3 NAT EUR AFR
#     1 0.5 0.45 0.05
#     16 0 0 0
#
# VERIFY THIS AGAINST YOUR EXISTING Brasa.dat BEFORE THE FULL RUN. The header
# is checked by simu-mix.pl itself (check_pops_match_dat), so a wrong pop list
# fails loudly; the generation lines are consumed by the `mixer` binary, which
# will accept a differently-shaped history WITHOUT complaining and silently
# simulate the wrong demography. That is the one assumption in this script
# that cannot fail loudly, so `--dat-only` prints the generated .dat for every
# model and exits -- runnable anywhere, needing no panel, no conda env and no
# cluster, so the two can be diffed before 2,160 genomes are built on it.
write_dat () {
    local model="$1" comps="$2" props="$3" gen="$4" n="$5" out="$6"
    local -a c p z
    IFS=',' read -r -a c <<< "$comps"
    IFS=',' read -r -a p <<< "$props"
    if [[ ${#c[@]} -ne ${#p[@]} ]]; then
        echo "ERROR: ${model} has ${#c[@]} components but ${#p[@]} proportions" >&2
        return 1
    fi
    printf '%s %s %s\n' "$n" "$model" "${c[*]}" > "$out"
    printf '1 %s\n' "${p[*]}" >> "$out"
    if [[ "$gen" -gt 1 ]]; then
        z=(); for _ in "${c[@]}"; do z+=("0"); done
        printf '%s %s\n' "$gen" "${z[*]}" >> "$out"
    fi
}

if [[ "${DAT_ONLY:-}" == "1" || "${1:-}" == "--dat-only" ]]; then
    [[ -s "$INDIV_TSV" ]] || { echo "ERROR: missing ${INDIV_TSV}" >&2; exit 1; }
    # write_dat truncates its output file, so it must not be handed
    # /dev/stdout here: inside a pipeline each call would rewind to offset 0
    # and overwrite the previous model, leaving only the last one. A real temp
    # file, then cat.
    _t=$(mktemp); trap 'rm -f "$_t"' EXIT
    awk -F'\t' 'NR>1{print $2"\t"$3"\t"$4"\t"$5}' "$INDIV_TSV" | sort -u |
    while IFS=$'\t' read -r model comps props gen; do
        n=$(awk -F'\t' -v m="$model" 'NR>1 && $2==m' "$INDIV_TSV" | wc -l)
        if [[ "$gen" -eq 0 ]]; then
            echo "--- ${model} (n=${n}, gen=0 CONTROL -- no .dat used) ---"
            echo "    single component ${comps}, ancestry constant; simulated by"
            echo "    sampling donor haplotype pairs, not by admix-simu."
            echo
            continue
        fi
        echo "--- ${model}.dat (n=${n}, gen=${gen}) ---"
        write_dat "$model" "$comps" "$props" "$gen" "$n" "$_t"
        cat "$_t"
        echo
    done
    exit 0
fi

# --- conda first, before any tool call -------------------------------------
# A staleness check that ran before `conda activate` once made 2 of 594 chunk
# tasks exit 127 on `bcftools: command not found`. Activation goes at the top,
# ahead of every tool invocation, not next to the first one that looks like it
# needs it.
if [[ -z "${CONDA_PREFIX:-}" || "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
for t in bcftools plink2 perl; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH in env ${CONDA_ENV}" >&2; exit 127; }
done

PANEL=$(printf "$PANEL_TPL" "$CHR")
[[ -s "$PANEL" ]] || { echo "ERROR: missing panel $PANEL" >&2; exit 1; }
[[ -s "$INDIV_TSV" ]] || { echo "ERROR: missing ${INDIV_TSV}
       Run scripts/assemble-benchmark-cohort.py first -- without it there is
       no roster, so there is nothing to simulate and no cell assignment." >&2; exit 1; }
[[ -s "$SIMU_MIX" ]] || SIMU_MIX="${ADMIXSIMU_DIR}/simu-mix.pl"
[[ -s "$SIMU_MIX" ]] || { echo "ERROR: no simu-mix.pl at ${SIMU_MIX}" >&2; exit 1; }
INSERT_MAP="${ADMIXSIMU_DIR}/insert-map.pl"
BP2ANC="${ADMIXSIMU_DIR}/bp2anc.pl"
for p in "$INSERT_MAP" "$BP2ANC"; do
    [[ -s "$p" ]] || { echo "ERROR: missing ${p} (set ADMIXSIMU_DIR)" >&2; exit 1; }
done

# The panel must be the 3,921-sample build, not the 3,448-sample July one.
# Simulating from the old panel would silently drop SGDP, Almarri and the 60
# restored samples from every donor pool.
N_PANEL=$(bcftools query -l "$PANEL" | wc -l)
echo "[$(date +%T)] [${CHR}] panel ${PANEL} -- ${N_PANEL} samples"
MIN_PANEL="${MIN_PANEL:-3900}"
if [[ "$N_PANEL" -lt "$MIN_PANEL" ]]; then
    echo "ERROR: panel has ${N_PANEL} samples, expected >= ${MIN_PANEL}.
       This looks like the pre-rerun panel. Refusing to simulate from it:
       the donor pools would be missing SGDP, Almarri and the restored
       samples, and nothing downstream would say so." >&2
    exit 3
fi

WORKDIR="${OUTDIR}/${CHR}"
PROV_DIR="${OUTDIR}/provenance/${CHR}"
mkdir -p "$WORKDIR" "$PROV_DIR"
cd "$WORKDIR"

# --- which components any model needs --------------------------------------
mapfile -t COMPONENTS < <(awk -F'\t' 'NR>1{n=split($3, a, ","); for(i=1;i<=n;i++) print a[i]}' \
                          "$INDIV_TSV" | sort -u)
[[ ${#COMPONENTS[@]} -gt 0 ]] || { echo "ERROR: no components in $INDIV_TSV" >&2; exit 1; }
echo "[$(date +%T)] [${CHR}] components needed: ${COMPONENTS[*]}"

# --- extract each component's donor haplotypes ONCE ------------------------
for comp in "${COMPONENTS[@]}"; do
    pool="${COMPONENT_POOL[$comp]:-$comp}"
    keep="${REFS}/${pool}_simulation_${BENCH_SUFFIX}.txt"
    if [[ ! -s "$keep" ]]; then
        echo "ERROR: no donor list ${keep} for component ${comp}.
       Run split-homog-donors-reference.py --mode benchmark --suffix ${BENCH_SUFFIX}" >&2
        exit 1
    fi
    n_don=$(awk 'NF && $1 !~ /^#/' "$keep" | wc -l)
    # A component with too few donors makes its simulated individuals cryptic
    # relatives of one another. That is a real finding to report, not a reason
    # to stop -- BENCHMARK_COHORT_PLAN.md section 5 leaves the SAS decision
    # open deliberately -- so warn loudly and continue.
    if [[ "$n_don" -lt 25 ]]; then
        echo "[$(date +%T)] [${CHR}] WARNING: component ${comp} (pool ${pool}) has only ${n_don} donors -> $((2*n_don)) haplotypes. Heavy reuse; check the kinship distribution."
    fi
    if [[ -s "${comp}.phgeno" ]]; then
        echo "[$(date +%T)] [${CHR}] ${comp}.phgeno present, reusing"
        continue
    fi
    echo "[$(date +%T)] [${CHR}] extracting ${comp} haps (${n_don} donors from ${pool})"
    plink2 --bcf "$PANEL" --keep <(awk 'NF && $1 !~ /^#/ {print $1"\t"$1}' "$keep") \
           --export haps --out "${comp}_${CHR}" --silent
    # .phgeno = one line per site, concatenated 0/1 with no separators.
    cut -d' ' -f 6- "${comp}_${CHR}.haps" | tr -d ' \t' > "${comp}.phgeno"
    [[ -s "${comp}.phgeno" ]] || { echo "ERROR: ${comp}.phgeno is empty" >&2; exit 1; }
done

# --- the .snp file with cM positions, built once ---------------------------
SNP="${CHR}.snp"
if [[ ! -s "$SNP" ]]; then
    first="${COMPONENTS[0]}"
    echo "[$(date +%T)] [${CHR}] building ${SNP}"
    plink2 --haps "${first}_${CHR}.haps" --sample "${first}_${CHR}.sample" \
           --make-just-bim --out "${CHR}_sites" --silent
    GMAP=$(ls "${GMAP_DIR}"/*"${CHR}"*.txt* 2>/dev/null | head -1 || true)
    [[ -n "$GMAP" ]] || { echo "ERROR: no genetic map for ${CHR} in ${GMAP_DIR}" >&2; exit 1; }
    perl "$INSERT_MAP" "${CHR}_sites.bim" "$GMAP" > "$SNP"
    [[ -s "$SNP" ]] || { echo "ERROR: ${SNP} is empty" >&2; exit 1; }
fi
N_SITES=$(wc -l < "$SNP")
echo "[$(date +%T)] [${CHR}] ${N_SITES} sites"

# ---------------------------------------------------------------------------
# DONOR-CONCURRENCY FEASIBILITY. The hard limit nobody documents.
#
# simu-mix.pl draws a donor haplotype uniformly at random but REJECTS any that
# is already being copied from by another simulated haplotype at that marker
# (scripts/simu-mix-seeded.pl ~line 310: the do/while skips
# $copying_sample[$pop][$id] when set). Every simulated haplotype holds its
# donor until its next breakpoint, so the draws are WITHOUT REPLACEMENT across
# the whole simulated set. When a population runs out it does not degrade --
# it dies:
#
#   die "Error: used all $num_haplotypes haplotypes in pop $pop_name ..."
#
# So a model needs, at every marker, at least as many donor haplotypes in
# component P as it has simulated haplotypes currently drawing from P --
# roughly 2 * n * proportion(P), and that is a random quantity, so headroom is
# needed, not equality.
#
# This is why 2b has never hit it and 6b would: 2b's model is small against
# its ~30-donor pools, while Table 1's sizes (LAT3 600, the rest 200) were set
# from the experimental design without reference to donor-pool size.
#
# Checked HERE, before any extraction, because the alternative is discovering
# it hours in, on the primary model, after the haplotypes are built.
FEAS_MARGIN="${FEAS_MARGIN:-0.80}"
feasible_batch () {   # <components csv> <proportions csv> -> max n per batch
    local comps="$1" props="$2" c p i=0 best=0 n_hap
    local -a ca pa
    IFS=',' read -r -a ca <<< "$comps"
    IFS=',' read -r -a pa <<< "$props"
    for i in "${!ca[@]}"; do
        c="${ca[$i]}"; p="${pa[$i]}"
        [[ "$(awk -v x="$p" 'BEGIN{print (x>0)?1:0}')" == "1" ]] || continue
        local pool="${COMPONENT_POOL[$c]:-$c}"
        local keep="${REFS}/${pool}_simulation_${BENCH_SUFFIX}.txt"
        local n_don; n_don=$(awk 'NF && $1 !~ /^#/' "$keep" 2>/dev/null | wc -l)
        n_hap=$(( n_don * 2 ))
        # n such that 2*n*p <= margin * n_hap
        local cap; cap=$(awk -v h="$n_hap" -v p="$p" -v m="$FEAS_MARGIN" \
                         'BEGIN{v=int(m*h/(2*p)); print (v<1)?1:v}')
        if [[ "$best" -eq 0 || "$cap" -lt "$best" ]]; then best="$cap"; fi
    done
    echo "$best"
}

# --- simulate each model ---------------------------------------------------
# One pass over the distinct models in the roster. Generations and proportions
# come from the roster, not from this script, so Table 1 lives in exactly one
# place and a change there does not need a change here.
awk -F'\t' 'NR>1{print $2"\t"$3"\t"$4"\t"$5}' "$INDIV_TSV" | sort -u > .models.tsv
N_MODELS=$(wc -l < .models.tsv)
echo "[$(date +%T)] [${CHR}] ${N_MODELS} models to simulate"

while IFS=$'\t' read -r model comps props gen; do
    n=$(awk -F'\t' -v m="$model" 'NR>1 && $2==m' "$INDIV_TSV" | wc -l)
    prefix="${model}.${CHR}"

    if [[ -s "${prefix}.hanc" && -s "${prefix}.haps" ]]; then
        echo "[$(date +%T)] [${CHR}] [${model}] already done (n=${n}), skipping"
        continue
    fi

    # Controls have gen 0 and a single component at proportion 1.0: there is
    # no admixture to simulate, so admix-simu is not the right tool. Sample
    # donor haplotype PAIRS directly instead; the ancestry is constant, which
    # is exactly what makes a control's error pure method error.
    if [[ "$gen" -eq 0 ]]; then
        echo "[$(date +%T)] [${CHR}] [${model}] gen 0 control (n=${n}) -- sampling donor haplotype pairs"
        comp="${comps%%,*}"
        SIMU_SEED="${SEED_BASE}${CHR_N}" \
        awk -v n="$n" -v seed="${SEED_BASE}${CHR_N}" -v model="$model" '
            BEGIN { srand(seed) }
            { nh = length($0); for (i = 1; i <= nh; i++) h[NR,i] = substr($0,i,1); ns = NR }
            END {
              for (j = 1; j <= n; j++) {
                a = int(rand()*nh)+1; b = int(rand()*nh)+1
                while (b == a && nh > 1) b = int(rand()*nh)+1
                pa[j] = a; pb[j] = b
              }
              for (s = 1; s <= ns; s++) {
                line = ""
                for (j = 1; j <= n; j++) line = line h[s,pa[j]] h[s,pb[j]]
                print line
              }
            }' "${comp}.phgeno" > "${prefix}.phgeno"
        # Constant ancestry: every site of every haplotype is this component.
        awk -v n="$n" -v ns="$N_SITES" 'BEGIN{
            line=""; for(j=1;j<=2*n;j++) line=line "0"
            for(s=1;s<=ns;s++) print line }' > "${prefix}.hanc"
        : > "${prefix}.haps"   # marker: this model produced no admix-simu haps
        echo "${model} gen0 control component=${comp} n=${n}" > "${prefix}.provenance"
        continue
    fi

    # -POP args in the .dat header's order.
    popargs=()
    IFS=',' read -r -a carr <<< "$comps"
    for comp in "${carr[@]}"; do
        [[ -s "${comp}.phgeno" ]] || { echo "ERROR: missing ${comp}.phgeno for ${model}" >&2; exit 1; }
        popargs+=("-${comp}" "${comp}.phgeno")
    done

    # --- BATCH SIZE, from the donor-concurrency limit --------------------
    CAP=$(feasible_batch "$comps" "$props")
    if [[ "$CAP" -ge "$n" ]]; then
        NBATCH=1; BSIZE="$n"
    else
        NBATCH=$(( (n + CAP - 1) / CAP )); BSIZE="$CAP"
        echo "[$(date +%T)] [${CHR}] [${model}] n=${n} exceeds the donor-concurrency"
        echo "            cap of ${CAP}; simulating in ${NBATCH} batches of <= ${BSIZE}."
        echo "            simu-mix holds one donor haplotype per simulated"
        echo "            haplotype and dies when a component runs out, so this"
        echo "            is a hard limit, not a quality knob. Individuals in"
        echo "            DIFFERENT batches may share donors -- that relatedness"
        echo "            is real and is what 6e measures."
    fi

    SEED_M="${SEED_BASE}$(printf '%02d%02d' "$CHR_N" "$(cksum <<< "$model" | cut -c1-2)")"
    rm -f "${prefix}.batchlist"
    left="$n"; b=0
    while [[ "$left" -gt 0 ]]; do
        b=$(( b + 1 ))
        take=$(( left < BSIZE ? left : BSIZE ))
        bp_prefix="${prefix}.b${b}"
        BDAT="${bp_prefix}.dat"
        write_dat "$model" "$comps" "$props" "$gen" "$take" "$BDAT"
        # A distinct seed per (model, chr, batch): independent realisations,
        # each individually regenerable.
        SEED="${SEED_M}${b}"
        archived="${FROM_BP}/${bp_prefix}.bp"
        if [[ -n "$FROM_BP" ]]; then
            if [[ ! -s "$archived" ]]; then
                echo "ERROR: FROM_BP set but ${archived} is missing. Refusing to draw
       fresh breakpoints under a name that implies replay -- that would be a
       different realisation presented as the same one." >&2
                exit 4
            fi
            echo "[$(date +%T)] [${CHR}] [${model}] batch ${b}/${NBATCH} REPLAY from ${archived}"
            cp "$archived" "${bp_prefix}.bp"
            SIMU_SEED="$SEED" perl "$SIMU_MIX" -bp "${bp_prefix}.bp" "$bp_prefix" "${popargs[@]}"
        else
            echo "[$(date +%T)] [${CHR}] [${model}] batch ${b}/${NBATCH} n=${take} gen=${gen} seed=${SEED}"
            SIMU_SEED="$SEED" perl "$SIMU_MIX" "$BDAT" "$SNP" "$bp_prefix" "${popargs[@]}"
        fi
        [[ -s "${bp_prefix}.bp" ]] || { echo "ERROR: ${bp_prefix}.bp was not written" >&2; exit 1; }
        perl "$BP2ANC" "${bp_prefix}.bp" > "${bp_prefix}.hanc"
        [[ -s "${bp_prefix}.hanc" ]] || { echo "ERROR: ${bp_prefix}.hanc is empty" >&2; exit 1; }
        cp "${bp_prefix}.bp"  "${PROV_DIR}/${bp_prefix}.bp"
        cp "$BDAT"            "${PROV_DIR}/${bp_prefix}.dat"
        echo "SIMU_SEED=${SEED}" > "${PROV_DIR}/${bp_prefix}.seed"
        echo "$bp_prefix" >> "${prefix}.batchlist"
        left=$(( left - take ))
    done

    # --- stitch the batches ----------------------------------------------
    # COLUMN-wise: .hanc and .phgeno carry one character per HAPLOTYPE per
    # site, so batches are joined side by side, not end to end. `paste -d ''`
    # on one line per site does exactly that. Row counts must match or the
    # sites are not the same sites.
    if [[ "$NBATCH" -eq 1 ]]; then
        mv "${prefix}.b1.hanc" "${prefix}.hanc"
        [[ -s "${prefix}.b1.haps" ]] && mv "${prefix}.b1.haps" "${prefix}.haps"
        [[ -s "${prefix}.b1.phgeno" ]] && mv "${prefix}.b1.phgeno" "${prefix}.phgeno"
        mv "${prefix}.b1.dat" "$DAT" 2>/dev/null || true
        cp "${PROV_DIR}/${prefix}.b1.dat" "${PROV_DIR}/${prefix}.dat" 2>/dev/null || true
    else
        rows=$(wc -l < "${prefix}.b1.hanc")
        while read -r bpx; do
            r=$(wc -l < "${bpx}.hanc")
            [[ "$r" -eq "$rows" ]] || { echo "ERROR: ${bpx}.hanc has ${r} sites, ${prefix}.b1.hanc has ${rows}. Batches must share a marker set." >&2; exit 6; }
        done < "${prefix}.batchlist"
        paste -d '' $(awk '{printf "%s.hanc ", $0}' "${prefix}.batchlist") > "${prefix}.hanc"
        if ls "${prefix}".b*.phgeno >/dev/null 2>&1; then
            paste -d '' $(awk '{printf "%s.phgeno ", $0}' "${prefix}.batchlist") > "${prefix}.phgeno"
        fi
        : > "${prefix}.haps"
        # The model's .dat records the FULL n; the per-batch ones are archived.
        write_dat "$model" "$comps" "$props" "$gen" "$n" "${prefix}.dat"
        cp "${prefix}.dat" "${PROV_DIR}/${prefix}.dat"
        echo "[$(date +%T)] [${CHR}] [${model}] stitched ${NBATCH} batches -> $(head -1 "${prefix}.hanc" | wc -c) haplotype columns"
    fi
    DAT="${prefix}.dat"
done < .models.tsv
rm -f .models.tsv

echo "[$(date +%T)] [${CHR}] DONE -> ${WORKDIR}"
echo "[$(date +%T)] [${CHR}] breakpoints + seeds archived in ${PROV_DIR}"

# The .dat generation lines are the one assumption here that cannot fail
# loudly: `mixer` simulates whatever history it parses and reports nothing if
# that is not the history that was meant. Tract length is the observable that
# betrays it, so the check runs as part of the job rather than waiting to be
# remembered. It does not fail the task -- the output is still a valid
# realisation of SOME demography, and one chromosome's deviation is noisy --
# but the verdict lands in the log next to the run that produced it.
#
# TWO CHECKS, DELIBERATELY DIFFERENT STATISTICS over the same truth, so a bug
# in one does not reproduce in the other:
#
#   validate-simset-tracts.py    counts ancestry SWITCHES per Morgan
#   cohort-ancestry-summary.py   measures ancestry FRACTIONS per individual
#
# The second exists because the ancestry MEANS stay correct under the failure
# -- the pulse proportions come off the .dat's first line and survive even
# when the terminal generation line is ignored -- so checking means alone is a
# false reassurance. It is the between-individual SPREAD that carries the
# signal: ~0.30 at gen 2 against ~0.077 at gen 16 per chromosome. A cohort
# that silently ran one generation shows right means with SDs several times
# too wide, uniformly.
VALIDATE="${VALIDATE:-${REPO_DIR}/scripts/validate-simset-tracts.py}"
ANCSUM="${ANCSUM:-${REPO_DIR}/scripts/cohort-ancestry-summary.py}"
if [[ -s "$VALIDATE" ]] && command -v python3 >/dev/null; then
    echo "[$(date +%T)] [${CHR}] validating tract lengths against intended generations"
    python3 -I "$VALIDATE" --simset "$OUTDIR" --individuals "$INDIV_TSV" \
        --chr "$CHR" || {
        echo "[$(date +%T)] [${CHR}] *** TRACT VALIDATION FAILED -- read the table above"
        echo "[$(date +%T)] [${CHR}] *** do not score this cohort until it passes"
    }
fi
if [[ -s "$ANCSUM" ]] && command -v python3 >/dev/null; then
    echo "[$(date +%T)] [${CHR}] ancestry proportions: mean AND between-individual spread"
    python3 -I "$ANCSUM" --simset "$OUTDIR" --individuals "$INDIV_TSV" \
        --chr "$CHR" --out "${PROV_DIR}/ancestry_summary.${CHR}.txt" || {
        echo "[$(date +%T)] [${CHR}] *** ANCESTRY SPREAD DISAGREES WITH THE INTENDED GENERATIONS"
        echo "[$(date +%T)] [${CHR}] *** note the MEANS can look right while this fails"
    }
fi
