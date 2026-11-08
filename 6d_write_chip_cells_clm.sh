#!/usr/bin/env bash
# ============================================================================
# 6d_write_chip_cells_clm.sh  --  build step 10: the chip x build x chr files
#
#SBATCH --job-name=6d_cells
#SBATCH --partition=atkinson,mhgcp
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/6d_cells_%a_%j.out
#SBATCH --error=logs/6d_cells_%a_%j.err
# The t-nodes cannot execute these binaries: 270/270 tasks of an earlier array
# died there with SIGILL (exit 132) before writing a log. c02 is fine.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#
# STATUS: the bcftools/plink2 pipeline below has NOT been run against real
# data. The allele-orientation logic it depends on -- scripts/lift-cell-to-hg19.py
# -- IS unit-tested (18 cases: same / flip / swap / flip+swap / palindromic /
# unresolvable, plus GT-field handling), because that is where a silent wrong
# answer would come from. Run one chromosome and read the reports before
# launching all 22.
#
# WHAT A CELL IS. Three chips x two builds = six cells. Every individual, real
# and simulated, is in exactly one; the assignment lives in the manifests and
# is NOT recomputed here, so step 6 and step 10 cannot disagree about who is
# where.
#
# THE hg19 CELLS USE THE CHIP'S NATIVE hg19 MANIFEST, never a down-lift of the
# hg38 sites. A down-lift contains only the probes that map back cleanly, by
# construction, so it would hide every probe the pipeline's own CrossMap
# drops, moves or strand-flips -- which is part of what this benchmark is for.
# Genotypes are read at the GRCh38 site and written at the hg19 site the same
# probe interrogates, with GT re-oriented where the builds disagree about the
# reference allele. The truth stays GRCh38 for everyone.
#
# ONE TASK PER CHROMOSOME, SIX CELLS INSIDE: the per-chromosome inputs (the
# real BCF and each model's simulated haplotypes) are read once and reused
# across cells, as in 6b.
# ============================================================================
set -euo pipefail

CHR_N="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job (1-22)}"
CHR="chr${CHR_N}"

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REAL_DIR="${REAL_DIR:-${PROJECT_ROOT}/07_benchmark_real}"
SIMSET_DIR="${SIMSET_DIR:-${PROJECT_ROOT}/06_benchmark_simset}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/08_chip_cells}"
MANIFEST_DIR="${MANIFEST_DIR:-${REPO_DIR}/chip_manifests/manifests}"
BENCH_OUT="${BENCH_OUT:-${REPO_DIR}/results}"
REAL_MANIFEST="${REAL_MANIFEST:-${BENCH_OUT}/benchmark_cohort_manifest.tsv}"
SIM_MANIFEST="${SIM_MANIFEST:-${BENCH_OUT}/benchmark_cohort_simset_individuals.tsv}"
LIFT="${LIFT:-${REPO_DIR}/scripts/lift-cell-to-hg19.py}"
CHIPS="${CHIPS:-snp6 affy500k gsa}"
BUILDS="${BUILDS:-hg38 hg19}"
CONDA_ENV="${CONDA_ENV:-shapeit5}"

# Conda before any tool call -- a check that ran ahead of activation once made
# 2 of 594 chunk tasks exit 127 on "bcftools: command not found".
if [[ -z "${CONDA_PREFIX:-}" || "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
for t in bcftools plink2 python3; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH in ${CONDA_ENV}" >&2; exit 127; }
done

for f in "$REAL_MANIFEST" "$SIM_MANIFEST" "$LIFT"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

WORK="${OUTDIR}/.work_${CHR}"
REP="${OUTDIR}/reports"
rm -rf "$WORK"; mkdir -p "$WORK" "$REP"
cd "$WORK"

col () {   # col <file> <header-name>  -> 1-based index, resolved BY NAME
    awk -F'\t' -v want="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==want){print i; exit}}' "$1"
}

# ---------------------------------------------------------------------------
# Simulated genotypes -> VCF, once per model per chromosome.
#
# SAMPLE NAMING IS THE JOINT. 6b writes one .haps per model holding that
# model's individuals in roster order, and the assembler names them
# <MODEL>_0001.. in that same order. bcftools reheader -s replaces names
# POSITIONALLY, so the two conventions line up by construction -- but a count
# mismatch would silently shift every id by some offset and put genotypes
# under the wrong individual, so it is checked, not assumed.
# ---------------------------------------------------------------------------
declare -A MODEL_VCF=()
prepare_model () {
    local model="$1"
    [[ -n "${MODEL_VCF[$model]:-}" ]] && return 0
    local haps="${SIMSET_DIR}/${CHR}/${model}.${CHR}.haps"
    local samp="${SIMSET_DIR}/${CHR}/${model}.${CHR}.sample"
    if [[ ! -s "$haps" ]]; then
        echo "  [${model}] no .haps for ${CHR} -- run 6b first; skipping model" >&2
        MODEL_VCF[$model]="MISSING"
        return 0
    fi
    local n_want n_have
    n_want=$(awk -F'\t' -v c="$C_SIM_MODEL" -v m="$model" 'NR>1 && $c==m' "$SIM_MANIFEST" | wc -l)
    # .haps columns: 5 leading fields then 2 per individual.
    n_have=$(awk 'NR==1{print int((NF-5)/2); exit}' "$haps")
    if [[ "$n_want" -ne "$n_have" ]]; then
        echo "ERROR: ${model} ${CHR}: .haps carries ${n_have} individuals but the
       roster lists ${n_want}. reheader -s renames POSITIONALLY, so a
       mismatch would put genotypes under the wrong sim_id. Re-run 6b for
       this model rather than shipping a shifted cohort." >&2
        exit 4
    fi
    awk -F'\t' -v c="$C_SIM_MODEL" -v i="$C_SIM_ID" -v m="$model" \
        'NR>1 && $c==m {print $i}' "$SIM_MANIFEST" > "${model}.names"
    plink2 --haps "$haps" --sample "$samp" --export vcf bgz \
           --out "${model}.raw" --silent
    bcftools reheader -s "${model}.names" -o "${model}.vcf.gz" "${model}.raw.vcf.gz"
    bcftools index -f -t "${model}.vcf.gz"
    MODEL_VCF[$model]="${WORK}/${model}.vcf.gz"
}

C_REAL_ID=$(col "$REAL_MANIFEST" sample)
C_REAL_CHIP=$(col "$REAL_MANIFEST" chip)
C_REAL_BUILD=$(col "$REAL_MANIFEST" build)
C_SIM_ID=$(col "$SIM_MANIFEST" sim_id)
C_SIM_MODEL=$(col "$SIM_MANIFEST" model)
C_SIM_CHIP=$(col "$SIM_MANIFEST" chip)
C_SIM_BUILD=$(col "$SIM_MANIFEST" build)
for v in C_REAL_ID C_REAL_CHIP C_REAL_BUILD C_SIM_ID C_SIM_MODEL C_SIM_CHIP C_SIM_BUILD; do
    [[ -n "${!v}" ]] || { echo "ERROR: could not resolve column for $v" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
for CHIP in $CHIPS; do
  for BUILD in $BUILDS; do
    CELL="${CHIP}.${BUILD}"
    SITES="${MANIFEST_DIR}/${CHIP}.${BUILD}.sites.tsv"
    if [[ ! -s "$SITES" ]]; then
        echo "[$(date +%T)] [${CHR}] [${CELL}] no site list ${SITES} -- skipping"
        echo "            (only GSA v1.0 is built; Affy 6.0 and 500K are pending)"
        continue
    fi
    OUT="${OUTDIR}/benchmark.${CHIP}.${BUILD}.${CHR}.bcf"
    if [[ -s "$OUT" && -s "${OUT}.csi" ]]; then
        echo "[$(date +%T)] [${CHR}] [${CELL}] already done, skipping"
        continue
    fi

    # Extraction positions are always GRCh38: that is what the panel and the
    # simulations are on. For an hg19 cell the hg38 coordinate comes from the
    # manifest's grch38_site, which is the vendor's own probe->hg38 mapping.
    POS="${CELL}.pos38.txt"
    if [[ "$BUILD" == "hg38" ]]; then
        awk -F'\t' 'NR>1 && NF>=2 {c=$1; sub(/^chr/,"",c); print c"\t"$2}' "$SITES" \
          | sort -k1,1 -k2,2n -u > "$POS"
    else
        g=$(col "$SITES" grch38_site)
        [[ -n "$g" ]] || { echo "ERROR: ${SITES} has no grch38_site column" >&2; exit 1; }
        awk -F'\t' -v g="$g" 'NR>1 && $g ~ /:/ {split($g,a,":"); c=a[1];
             sub(/^chr/,"",c); print c"\t"a[2]}' "$SITES" \
          | sort -k1,1 -k2,2n -u > "$POS"
    fi
    N_SITES=$(wc -l < "$POS")
    [[ "$N_SITES" -gt 0 ]] || { echo "ERROR: no ${CHR} sites for ${CELL}" >&2; exit 1; }
    awk -v c="$CHR_N" '$1==c' "$POS" > "${CELL}.pos38.chr.txt"
    N_CHR=$(wc -l < "${CELL}.pos38.chr.txt")
    echo "[$(date +%T)] [${CHR}] [${CELL}] ${N_CHR} chip sites on this chromosome"
    [[ "$N_CHR" -gt 0 ]] || { echo "  no sites on ${CHR}; skipping"; continue; }

    PARTS=()

    # --- real individuals in this cell ---
    awk -F'\t' -v i="$C_REAL_ID" -v c="$C_REAL_CHIP" -v b="$C_REAL_BUILD" \
        -v chip="$CHIP" -v build="$BUILD" \
        'NR>1 && $c==chip && $b==build {print $i}' "$REAL_MANIFEST" > "${CELL}.real.ids"
    N_REAL=$(wc -l < "${CELL}.real.ids")
    REAL_BCF="${REAL_DIR}/benchmark_real.${CHR}.bcf"
    if [[ "$N_REAL" -gt 0 && -s "$REAL_BCF" ]]; then
        bcftools view -S "${CELL}.real.ids" -T "${CELL}.pos38.chr.txt" \
                 -Oz -o "${CELL}.real.vcf.gz" "$REAL_BCF"
        bcftools index -f -t "${CELL}.real.vcf.gz"
        PARTS+=("${CELL}.real.vcf.gz")
    elif [[ "$N_REAL" -gt 0 ]]; then
        echo "  WARNING: ${N_REAL} real individuals assigned here but ${REAL_BCF} is missing (run 6c)"
    fi

    # --- simulated individuals in this cell, per model ---
    awk -F'\t' -v m="$C_SIM_MODEL" -v c="$C_SIM_CHIP" -v b="$C_SIM_BUILD" \
        -v chip="$CHIP" -v build="$BUILD" \
        'NR>1 && $c==chip && $b==build {print $m}' "$SIM_MANIFEST" | sort -u > "${CELL}.models"
    while read -r model; do
        [[ -n "$model" ]] || continue
        prepare_model "$model"
        [[ "${MODEL_VCF[$model]}" == "MISSING" ]] && continue
        awk -F'\t' -v i="$C_SIM_ID" -v mm="$C_SIM_MODEL" -v c="$C_SIM_CHIP" \
            -v b="$C_SIM_BUILD" -v chip="$CHIP" -v build="$BUILD" -v m="$model" \
            'NR>1 && $mm==m && $c==chip && $b==build {print $i}' \
            "$SIM_MANIFEST" > "${CELL}.${model}.ids"
        [[ -s "${CELL}.${model}.ids" ]] || continue
        bcftools view -S "${CELL}.${model}.ids" -T "${CELL}.pos38.chr.txt" \
                 -Oz -o "${CELL}.${model}.vcf.gz" "${MODEL_VCF[$model]}"
        bcftools index -f -t "${CELL}.${model}.vcf.gz"
        PARTS+=("${CELL}.${model}.vcf.gz")
    done < "${CELL}.models"

    if [[ ${#PARTS[@]} -eq 0 ]]; then
        echo "  nothing to write for ${CELL} on ${CHR} (no inputs present yet)"
        continue
    fi

    if [[ ${#PARTS[@]} -eq 1 ]]; then
        cp "${PARTS[0]}" "${CELL}.merged.vcf.gz"
    else
        # -0 : a site absent from one part is reference, not missing. These are
        # chip sites present in every input by construction, so a gap means the
        # part simply had no ALT carrier there.
        bcftools merge -0 -Oz -o "${CELL}.merged.vcf.gz" "${PARTS[@]}"
    fi
    bcftools index -f -t "${CELL}.merged.vcf.gz"

    if [[ "$BUILD" == "hg38" ]]; then
        bcftools view -Ob -o "$OUT" "${CELL}.merged.vcf.gz"
    else
        # Rewrite onto the chip's native hg19 coordinates, re-orienting GT
        # where the builds disagree about the reference allele. Read the
        # report: a large "unresolvable" count means the manifest and the VCF
        # disagree about which reference the hg38 side was called against.
        bcftools view -Ov "${CELL}.merged.vcf.gz" \
          | python3 -I "$LIFT" --manifest "$SITES" \
                --report "${REP}/lift.${CHIP}.${CHR}.txt" \
          | bcftools sort -Ob -o "$OUT" -
    fi
    bcftools index -f "$OUT"

    N_OUT_S=$(bcftools query -l "$OUT" | wc -l)
    N_OUT_V=$(bcftools index -n "$OUT")
    # Every id file for this cell, deduplicated: real.ids plus one per model.
    # Counting N_REAL separately and adding would double-count it.
    N_EXP=$(cat "${CELL}".real.ids "${CELL}".*.ids 2>/dev/null | sort -u | grep -c . || true)
    if [[ "$N_OUT_S" -ne "$N_EXP" ]]; then
        echo "ERROR: ${CELL} ${CHR} holds ${N_OUT_S} samples but ${N_EXP} were
       assigned to this cell. A cell that silently lost or gained individuals
       would unbalance the chip x build design, and a chip effect and an
       ancestry effect are not separable once that happens." >&2
        exit 5
    fi
    echo "[$(date +%T)] [${CHR}] [${CELL}] -> $(basename "$OUT"): ${N_OUT_S} samples, ${N_OUT_V} sites"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CHR" "$CHIP" "$BUILD" "$N_OUT_S" "$N_OUT_V" "$N_CHR" \
      >> "${REP}/cells.${CHR}.tsv"
  done
done

echo "[$(date +%T)] [${CHR}] DONE -> ${OUTDIR}"
echo "[$(date +%T)] [${CHR}] reports in ${REP}"
rm -rf "$WORK"
