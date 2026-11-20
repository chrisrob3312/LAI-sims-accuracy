#!/usr/bin/env bash
#SBATCH --job-name=site_report
#SBATCH --partition=atkinson,mhgcp
# Same exclusion as every other stage: the mhgcp-t* CPUs lack instructions this
# conda env's bcftools is built with (SIGILL, exit 132, empty log).
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/site_report_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/site_report_chr%a_%j.err
# =============================================================================
# cohort-site-report.sh -- per site and per cohort, over 1c's published
# soft-union: where each cohort's OWN file stood at every record (present;
# absent, so 1c filled it 0/0; swapped / flipped / other allele), how much of
# the cohort is missing there before phasing (what SHAPEIT5 will impute), and
# the cohort's allele frequency. READ-ONLY on the panel: writes only its own
# report files.
#
# 1c already logs this as COUNTS, per cohort, at three stages (raw merge,
# after QC, soft-union): chrN.{merged,qc,softunion}.cohorts.tsv and
# .patterns.tsv in ${PANEL_DIR}/provenance/, plus fill_check.tsv,
# softunion_groups.tsv and merged_chrN.admix_sites.summary.tsv. This report
# adds the POSITIONS, and genotype missingness per cohort, which 1c does not
# write. It runs after 1c, so nothing in the queued chain changes.
#
#   bash scripts/cohort-site-report.sh 21 22            # interactive, shapeit5 env
#   sbatch --array=1-22 scripts/cohort-site-report.sh   # every chromosome
#
# Writes, per chromosome:
#   ${PANEL_DIR}/provenance/chrN.site_report.tsv.gz    one row per record:
#       pos ref alt admix, then per cohort src_<C> miss_<C> af_<C>
#   ${PANEL_DIR}/provenance/chrN.site_report.bins.tsv  per BIN_MB-Mb bin
# and a summary on stdout (codes in scripts/cohort-site-report.py).
# =============================================================================
set -uo pipefail   # no -e: report every chromosome it can

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
# The HGDP+1KG file 1c merged -- same default as 1c and the rebuild driver.
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chr__CHR__.filtered.SNV_INDEL.phased.shapeit5.bcf}"
BIN_MB="${BIN_MB:-1}"
PY="${REPO_DIR}/scripts/cohort-site-report.py"
[[ -f "$PY" ]] || { echo "ERROR: ${PY} not found (run from the repo, or set REPO_DIR)" >&2; exit 2; }

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    module load anaconda3/2024.06
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    set +u; conda activate shapeit5; set -u
    unset PYTHONHOME PYTHONPATH
fi
command -v bcftools > /dev/null || { echo "ERROR: bcftools not on PATH (conda activate shapeit5)" >&2; exit 2; }

if [[ $# -gt 0 ]]; then CHRS=("$@")
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then CHRS=("$SLURM_ARRAY_TASK_ID")
else echo "usage: bash $0 <chr> [chr...]   or   sbatch --array=1-22 $0" >&2; exit 2; fi

STATUS=0
for CHR in "${CHRS[@]}"; do
    CHR="${CHR#chr}"
    SU="${PANEL_DIR}/merged_chr${CHR}.softunion.bcf"
    echo "== chr${CHR}: ${SU}"
    if [[ ! -s "$SU" ]]; then
        echo "chr${CHR}: no soft-union yet (1c has not published it)"; STATUS=1; continue
    fi
    if [[ -f "${PROJECT_ROOT}/rebuild_epoch" && ! "$SU" -nt "${PROJECT_ROOT}/rebuild_epoch" ]]; then
        echo "chr${CHR}: WARNING: the soft-union predates the rebuild epoch -- this reports the previous build"
    fi
    [[ -s "${SU}.recipe" ]] && echo "chr${CHR}: 1c stamp (raw): $(cat "${SU}.recipe")"
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/site_report.chr${CHR}.XXXXXX") || { STATUS=1; continue; }

    # The cohort files 1c merged, in its merge order.
    LABELS=(); FILES=()
    add () {
        if [[ -s "$2" ]]; then LABELS+=("$1"); FILES+=("$2")
        else echo "chr${CHR}: ${1}: no file at ${2} -- left out of the report"; fi
    }
    add HGDP1KG "${HGDP1KG_PATTERN//__CHR__/$CHR}"
    add MXB     "${PANEL_DIR}/mxb_lifted.chr${CHR}.bcf"
    add SGDP    "${PANEL_DIR}/sgdp_prepped.chr${CHR}.bcf"
    add ALMARRI "${PANEL_DIR}/almarri_prepped.chr${CHR}.bcf"

    # Each soft-union sample's cohort: the first file (merge order) listing it.
    bcftools query -l "$SU" | sort -u > "${WORK}/su.samples"
    : > "${WORK}/groups.tsv"
    COHORT_ARGS=(); N_ARGS=(); USE=(); FMT='%POS\t%REF\t%ALT'
    for i in "${!LABELS[@]}"; do
        L="${LABELS[$i]}"
        bcftools query -l "${FILES[$i]}" | sort -u | comm -12 - "${WORK}/su.samples" \
            | comm -23 - <(cut -f1 "${WORK}/groups.tsv" | sort -u) \
            | awk -v l="$L" '{print $1 "\t" l}' >> "${WORK}/groups.tsv"
        n=$(awk -v l="$L" '$2 == l' "${WORK}/groups.tsv" | wc -l)
        if [[ "$n" -eq 0 ]]; then echo "chr${CHR}: ${L}: none of its samples are in the soft-union"; continue; fi
        bcftools query -f '%POS\t%REF\t%ALT\n' "${FILES[$i]}" > "${WORK}/sites.${L}.tsv"
        USE+=("$L"); N_ARGS+=(--n "${L}=${n}"); COHORT_ARGS+=(--cohort "${L}=${WORK}/sites.${L}.tsv")
        FMT+="\t%INFO/AN_${L}\t%INFO/AC_${L}"
    done
    n_su=$(wc -l < "${WORK}/su.samples"); n_assigned=$(wc -l < "${WORK}/groups.tsv")
    echo "chr${CHR}: ${n_su} soft-union samples, ${n_assigned} assigned to a cohort ($(cut -f2 "${WORK}/groups.tsv" | sort | uniq -c | awk '{printf "%s %s; ", $2, $1}'))"
    [[ "$n_assigned" -eq "$n_su" ]] || echo "chr${CHR}: WARNING: $((n_su - n_assigned)) soft-union sample(s) in no cohort file (not counted)"
    if [[ ${#USE[@]} -eq 0 ]]; then echo "chr${CHR}: no cohort to report"; STATUS=1; rm -rf "$WORK"; continue; fi

    # Per-cohort AN and AC at every soft-union record, from the genotypes.
    bcftools +fill-tags "$SU" -Ou -- -S "${WORK}/groups.tsv" -t AN,AC \
        | bcftools query -f "${FMT}\n" > "${WORK}/su.tsv" \
        || { echo "chr${CHR}: ERROR: fill-tags / query failed (stderr above)"; STATUS=1; rm -rf "$WORK"; continue; }

    AB="${PANEL_DIR}/merged_chr${CHR}.admix_sites.bed"; AB_ARG=()
    if [[ -s "$AB" ]]; then AB_ARG=(--admix-bed "$AB"); else echo "chr${CHR}: no ${AB}: admix column all 0"; fi
    mkdir -p "${PANEL_DIR}/provenance"
    python3 -I "$PY" --chrom "chr${CHR}" --su "${WORK}/su.tsv" --labels "${USE[*]}" \
        "${N_ARGS[@]}" "${COHORT_ARGS[@]}" ${AB_ARG[@]+"${AB_ARG[@]}"} --bin-mb "$BIN_MB" \
        --out-sites "${PANEL_DIR}/provenance/chr${CHR}.site_report.tsv.gz" \
        --out-bins "${PANEL_DIR}/provenance/chr${CHR}.site_report.bins.tsv" \
        || { echo "chr${CHR}: ERROR: cohort-site-report.py failed (stderr above)"; STATUS=1; }
    rm -rf "$WORK"
done
exit "$STATUS"
