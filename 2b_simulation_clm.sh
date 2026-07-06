#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=2b_simu
#SBATCH --partition=mhgcp
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=00:15:00
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/2b_simu_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/2b_simu_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Modernized replacement for simulation.sh. Builds admix-simu donor .phgeno
# files for both NAT donor configurations (HGDP-NAT-only and HGDP-NAT+MXB)
# directly from the new merged + jointly phased panel produced by
# 1b_phasing-jointcall_clm.sh, runs admix-simu, and emits simulated admixed
# haps + truth files per chrom.
#
# Differences from legacy simulation.sh:
#   - Source panel is merged_chr*.shapeit5_phased.softunion_maf005.rechr.bcf
#   - No SHAPEIT2 `shapeit -convert` -- uses plink2 --export haps directly.
#   - Two NAT donor tracks per run: NAT (HGDP-only) and NATMXB (HGDP + MXB).
#   - SLURM array 1-22 (parallel per-chrom). One (ADMIX_POP, GEN) per submit;
#     resubmit with different env vars to scan models.
#
# Outputs per chr (under $OUTDIR/$ADMIX_POP/gen$GEN/):
#   NAT${i}.phgeno           HGDP-NAT-simu donor haplotypes (track 1)
#   NATMXB${i}.phgeno        HGDP-NAT-simu + MXB-simu donor haplotypes (track 2)
#   EUR${i}.phgeno           IBS donor haplotypes
#   AFR${i}.phgeno           YRI donor haplotypes
#   <track>.${ADMIX_POP}${i}.bp / .hanc / .hanc2     truth ancestry per sim track
#   <track>.${ADMIX_POP}.chr${i}.haps / .sample      simulated admixed haps
# ============================================================================

# ----------------------------------------------------------------------------
# Variable configuration  --  paths and tunables (override via env vars at submit time)
# ----------------------------------------------------------------------------
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Project root -- all outputs anchored here so nothing collides with other
# lab members' work.
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"

# Phased panel from 1b_phasing-jointcall_clm.sh (chr-stripped, contigs 1..22)
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PANEL_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.rechr.bcf}"

# Sample-ID lists (in repo)
REFS="${REFS:-${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}/reference_ids}"
NAT_HGDP_SIMU="${NAT_HGDP_SIMU:-${REFS}/amr_simulation.txt}"
NAT_HGDPMXB_SIMU="${NAT_HGDPMXB_SIMU:-${REFS}/amr_hgdpmxb_simulation.txt}"
EUR_SIMU="${EUR_SIMU:-${REFS}/eur_simulation.txt}"
AFR_SIMU="${AFR_SIMU:-${REFS}/afr_simulation.txt}"

# admix-simu repo (contains insert-map.pl, simu-mix.pl, bp2anc.pl)
ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/simu-jointcall/admix-simu-master}"

# RFMix-format genetic map (3 cols: pos chr cM) used by insert-map.pl
GMAP_DIR="${GMAP_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/genetic_map/recomb-hg38}"
GMAP_TPL="${GMAP_TPL:-${GMAP_DIR}/genetic_map_chr%s_hg38.txt}"

# Simulation model
ADMIX_POP="${ADMIX_POP:-Brasa}"     # name of the simulated admixed pop (matches .dat file)
GEN="${GEN:-12}"                    # generations (matches .dat file)
DAT_FILE="${DAT_FILE:-${PROJECT_ROOT}/02_simulations/${ADMIX_POP}.dat}"
SAMPLE_TEMPLATE="${SAMPLE_TEMPLATE:-${PROJECT_ROOT}/02_simulations/${ADMIX_POP}.sample.txt}"

# Outputs
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/02_simulations}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-4}

# shellcheck disable=SC2059
PANEL=$(printf "$PANEL_TPL" "$CHR")
# shellcheck disable=SC2059
GMAP=$(printf "$GMAP_TPL" "$CHR")

[[ -s "$PANEL"           ]] || { echo "ERROR: missing $PANEL"; exit 1; }
[[ -s "$GMAP"            ]] || { echo "ERROR: missing $GMAP"; exit 1; }
[[ -d "$ADMIXSIMU_DIR"   ]] || { echo "ERROR: missing $ADMIXSIMU_DIR"; exit 1; }
[[ -s "$DAT_FILE"        ]] || { echo "ERROR: missing $DAT_FILE (admix-simu config)"; exit 1; }
[[ -s "$SAMPLE_TEMPLATE" ]] || { echo "ERROR: missing $SAMPLE_TEMPLATE"; exit 1; }

WORKDIR="${OUTDIR}/${ADMIX_POP}/gen${GEN}"
mkdir -p "$WORKDIR" "$LOGDIR"
cd "$WORKDIR"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u

# ----------------------------------------------------------------------------
# Helper: extract a per-pop .haps/.sample (SHAPEIT format) from the merged
# phased panel using plink2.
# ----------------------------------------------------------------------------
extract_haps () {
    local label="$1" keep_file="$2"
    plink2 --bcf "$PANEL" \
           --keep <(awk '{print 0, $1}' "$keep_file") \
           --export haps \
           --threads "$THREADS" \
           --out "${label}_chr${CHR}" \
           > "${label}_chr${CHR}.plink.log" 2>&1
}

echo "[$(date +%T)] [chr${CHR}] extracting per-pop haps from $PANEL"
extract_haps NAT_HGDP_SIMU       "$NAT_HGDP_SIMU"
extract_haps NAT_HGDPMXB_SIMU    "$NAT_HGDPMXB_SIMU"
extract_haps IBS_1KG_SIMU        "$EUR_SIMU"
extract_haps YRI_1KG_SIMU        "$AFR_SIMU"

# .phgeno = concatenated 0/1 haplotype matrix (one line per site, no separators)
echo "[$(date +%T)] [chr${CHR}] building .phgeno"
cut -d' ' -f 6- NAT_HGDP_SIMU_chr${CHR}.haps    | sed 's/\s//g' > NAT${CHR}.phgeno
cut -d' ' -f 6- NAT_HGDPMXB_SIMU_chr${CHR}.haps | sed 's/\s//g' > NATMXB${CHR}.phgeno
cut -d' ' -f 6- IBS_1KG_SIMU_chr${CHR}.haps     | sed 's/\s//g' > EUR${CHR}.phgeno
cut -d' ' -f 6- YRI_1KG_SIMU_chr${CHR}.haps     | sed 's/\s//g' > AFR${CHR}.phgeno

# .snp file with cM positions for admix-simu (one .snp per chr, shared across tracks).
# insert-map.pl reads a plink .bim; build it from a single donor .haps via plink2.
echo "[$(date +%T)] [chr${CHR}] building .snp (cM positions)"
plink2 --bcf "$PANEL" \
       --keep <(awk '{print 0, $1}' "$NAT_HGDP_SIMU") \
       --const-fid 0 --max-alleles 2 \
       --make-bed \
       --threads "$THREADS" \
       --out "1kg_hgdp_pBim2_${CHR}" \
       > "1kg_hgdp_pBim2_${CHR}.plink.log" 2>&1
perl "${ADMIXSIMU_DIR}/insert-map.pl" \
    "1kg_hgdp_pBim2_${CHR}.bim" "$GMAP" > "chr${CHR}.pos"
awk -F' ' '{ print $1":"$4"_"$5"_"$6, $1, $3, $4, $5, $6 }' "chr${CHR}.pos" \
    > "1kg_hgdp.chr${CHR}.snp"

# Per NAT donor track: simu-mix.pl -> .bp -> .hanc -> .hanc2 + simulated .haps/.sample
run_admix_simu () {
    local track="$1"        # "NAT" or "NATMXB"
    local prefix="${track}.${ADMIX_POP}${CHR}"

    echo "[$(date +%T)] [chr${CHR}] [$track] simu-mix.pl"
    perl "${ADMIXSIMU_DIR}/simu-mix.pl" \
        "$DAT_FILE" "1kg_hgdp.chr${CHR}.snp" "$prefix" \
        -NAT "${track}${CHR}.phgeno" \
        -EUR "EUR${CHR}.phgeno" \
        -AFR "AFR${CHR}.phgeno"

    echo "[$(date +%T)] [chr${CHR}] [$track] bp -> hanc -> hanc2"
    perl "${ADMIXSIMU_DIR}/bp2anc.pl" "${prefix}.bp" > "${prefix}.hanc"
    sed 's/./& /g' "${prefix}.hanc" > "${prefix}.hanc1"
    python3 -c "import sys; rows=[l.split() for l in sys.stdin if l.strip()]; print('\n'.join(' '.join(c) for c in zip(*rows)))" \
        < "${prefix}.hanc1" > "${prefix}.hanc2"
    rm "${prefix}.hanc1"

    # Build simulated-admixed-individual .haps/.sample
    sed 's/./& /g' "${prefix}.phgeno" > "${prefix}.tmp_geno"
    awk -F' ' '{ print $2, $1, $4, $5, $6 }' "1kg_hgdp.chr${CHR}.snp" > "${prefix}.tmp_snp"
    paste "${prefix}.tmp_snp" "${prefix}.tmp_geno" | sed 's/\t/ /g' \
        > "${track}.${ADMIX_POP}.chr${CHR}.haps"
    rm "${prefix}.tmp_geno" "${prefix}.tmp_snp"

    cp "$SAMPLE_TEMPLATE" "${track}.${ADMIX_POP}.chr${CHR}.sample"
}

run_admix_simu NAT
run_admix_simu NATMXB

echo "[$(date +%T)] [chr${CHR}] # Complete. Outputs in $WORKDIR"
