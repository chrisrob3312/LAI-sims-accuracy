#!/usr/bin/env bash
# =============================================================================
# check-sgdp-raw.sh -- where the raw SGDP WGS data stand. READ-ONLY.
#
#   1. the download job (squeue) and the last lines of its newest log
#   2. the download: bytes on disk vs the server's size, the DONE stamp
#   3. the per-sample VCFs: how many, with indexes
#   4. which manifest samples map to a per-sample VCF (and which do not), and
#      which VCFs on disk no manifest sample uses
#   5. the prepped chromosomes: recipe stamp, samples, records
#
# Run from the repo checkout inside the shapeit5 env (bcftools, python3):
#   bash scripts/check-sgdp-raw.sh
# SGDP_METADATA=<path> points at a copy of the SGDP metadata TSV if the prep
# has not downloaded one yet (the Block A probe copy works).
# Nothing here writes outside a temp file, and no stderr is hidden.
# =============================================================================
set -uo pipefail   # no -e: one failed check must not hide the others

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
MANIFEST="${SGDP_MANIFEST:-${REPO_DIR}/reference_ids/sgdp_manifest.tsv}"
SGDP_METADATA="${SGDP_METADATA:-${SGDP_RAW_DIR}/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt}"

echo "== 1. download job =="
squeue -u "${USER:-$(id -un)}" -n sgdp_fetch -o '%i %T %M %R'
LOG=$(ls -t "$LOGDIR"/sgdp_fetch_*.out | awk 'NR == 1')
if [[ -n "$LOG" ]]; then
    echo "-- newest log: $LOG"; tail -8 "$LOG"
    ERR="${LOG%.out}.err"; [[ -s "$ERR" ]] && { echo "-- its .err:"; tail -8 "$ERR"; }
fi

echo "== 2. download =="
ls -la "${SGDP_RAW_DIR}/download"
if [[ -f "${SGDP_RAW_DIR}/vcfs/DONE" ]]; then
    echo "-- DONE stamp:"; cat "${SGDP_RAW_DIR}/vcfs/DONE"
else
    echo "-- no DONE stamp yet (${SGDP_RAW_DIR}/vcfs/DONE)"
fi

echo "== 3. per-sample VCFs =="
VCFDIR="${SGDP_RAW_DIR}/vcfs"
[[ -d "$VCFDIR" ]] || VCFDIR="${SGDP_RAW_DIR}/vcfs.partial"
if [[ -d "$VCFDIR" ]]; then
    echo "-- ${VCFDIR}: $(ls "$VCFDIR" | grep -c '\.vcf\.gz$') VCFs, $(ls "$VCFDIR" | grep -c '\.vcf\.gz\.tbi$') indexes"
    ls "$VCFDIR" | grep '\.vcf\.gz$' | head -3
else
    echo "-- not extracted yet"
fi

echo "== 4. manifest -> per-sample VCF =="
if [[ ! -s "$SGDP_METADATA" ]]; then
    echo "-- no SGDP metadata at ${SGDP_METADATA} (the prep downloads it); point SGDP_METADATA at a copy"
elif [[ ! -d "$VCFDIR" ]]; then
    echo "-- skipped: no VCF directory yet"
else
    MAP=$(mktemp)
    python3 -I "${REPO_DIR}/scripts/sgdp-raw-sample-map.py" --manifest "$MANIFEST" \
        --metadata "$SGDP_METADATA" --vcf-dir "$VCFDIR" --out "$MAP"
    echo "-- exit $? ; first rows (Illumina ID, panel name, file):"
    head -3 "$MAP"
    echo "-- super-pops of the mapped samples:"
    awk -F'\t' 'NR == FNR {m[$2] = 1; next} FNR > 1 && ($2 in m) {print $1}' "$MAP" "$MANIFEST" | sort | uniq -c
    rm -f "$MAP"
fi

echo "== 5. prepped chromosomes (${OUTDIR}/sgdp_prepped.chrN.bcf) =="
printf 'chr\tsamples\trecords\trecipe\n'
for c in $(seq 1 22); do
    f="${OUTDIR}/sgdp_prepped.chr${c}.bcf"
    if [[ -s "$f" ]]; then
        printf '%s\t%s\t%s\t%s\n' "$c" "$(bcftools query -l "$f" | wc -l)" "$(bcftools index -n "$f")" \
            "$(cat "${f}.recipe" 2>&1)"
    else
        printf '%s\t-\t-\t(not built)\n' "$c"
    fi
done
echo "-- imputed (phased_data2021) outputs moved aside by the raw prep: $(ls "$OUTDIR" | grep -c 'sgdp_prepped\.chr[0-9]*\.phased_data2021\.bcf$')"
