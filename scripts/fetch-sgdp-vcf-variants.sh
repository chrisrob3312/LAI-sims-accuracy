#!/usr/bin/env bash
#SBATCH --job-name=sgdp_fetch
#SBATCH --partition=atkinson,mhgcp
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/sgdp_fetch_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/sgdp_fetch_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu
# =============================================================================
# fetch-sgdp-vcf-variants.sh -- download the SGDP per-sample WGS calls
#
# Source: the Reich lab share, vcf_variants/ -- one GATK UnifiedGenotyper VCF
# per public sample (279), variant sites only, hs37d5, with the per-call
# filter level FL; shipped as ONE 57 GB tar with an .md5sum beside every file:
#   https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/
#     vcfs.variants.public_samples.279samples.tar        (the data)
#     vcfs.variants.public_samples.279samples.file.list  (its `tar -tv` listing)
# This replaces phased_data2021 (GLIMPSE-imputed against 1KG phase 3) as the
# SGDP source of the panel. This job only fetches and verifies; the prep
# (filter, merge with absent = 0/0, liftover, normalise) is a separate step.
#
# What it does, and refuses to do:
#   1. downloads the tar to ${SGDP_RAW_DIR}/download, RESUMABLE (curl -C -):
#      a stall (< 10 kB/s for 10 min) or a dropped connection resumes where
#      it stopped, and so does resubmitting the job after a timeout;
#   2. requires the exact byte count the server reports;
#   3. lists the members BEFORE extracting and refuses the archive if any is
#      not a regular file with a bare name (no '/', no '..', no link);
#   4. extracts into a fresh directory of its own, checks every file against
#      the .md5sum shipped with it, and only then moves it to
#      ${SGDP_RAW_DIR}/vcfs with a DONE stamp;
#   5. never deletes the tar -- remove it yourself once the prep has run.
# Only system tools (curl, tar, md5sum, awk): no conda env needed.
#
# TLS: verified by default. SGDP_CACERT=<bundle> uses another CA bundle (e.g.
# the conda env's ssl/cacert.pem); SGDP_INSECURE=1 is curl -k, which
# 1b_sgdp_liftover_clm.sh already uses for this host because the cluster's CA
# bundle lacks its chain. Integrity then rests on the exact byte count and the
# md5sums inside the tar. Pick whichever the Block A check showed works.
#
# Space: tar + extracted files ~ 2 x 57 GB at the peak. Checked up front
# against df (SGDP_SKIP_SPACE_CHECK=1 to override, e.g. on a quota'd FS where
# df does not show the quota).
#
#   sbatch scripts/fetch-sgdp-vcf-variants.sh                 # verified TLS
#   SGDP_INSECURE=1 sbatch scripts/fetch-sgdp-vcf-variants.sh  # curl -k
# =============================================================================
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
BASE_URL="${SGDP_VCF_URL:-https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants}"
TAR_NAME="vcfs.variants.public_samples.279samples.tar"
LIST_NAME="vcfs.variants.public_samples.279samples.file.list"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-60}"

TLS=()
if [[ -n "${SGDP_CACERT:-}" ]]; then
    TLS=(--cacert "$SGDP_CACERT")
elif [[ "${SGDP_INSECURE:-0}" == "1" ]]; then
    TLS=(-k)
fi

DL="${SGDP_RAW_DIR}/download"
OUT="${SGDP_RAW_DIR}/vcfs"
PART="${SGDP_RAW_DIR}/vcfs.partial"
TAR="${DL}/${TAR_NAME}"
log () { echo "[$(date '+%F %T')] [sgdp-fetch] $*"; }

if [[ -f "${OUT}/DONE" ]]; then
    log "already fetched and verified:"; cat "${OUT}/DONE"; exit 0
fi
if [[ -e "$OUT" ]]; then
    log "ERROR: ${OUT} exists without a DONE stamp -- not this job's to overwrite. Look at it, then move it away."
    exit 1
fi
mkdir -p "$DL"
log "TLS: ${TLS[*]:-verified with the default CA bundle}; target ${SGDP_RAW_DIR}"

# 1. The published listing, and the size the server reports for the tar.
curl -fsS ${TLS[@]+"${TLS[@]}"} --retry 5 --retry-delay 30 -o "${DL}/${LIST_NAME}" "${BASE_URL}/${LIST_NAME}"
HEAD=$(curl -fsSI ${TLS[@]+"${TLS[@]}"} --retry 5 --retry-delay 30 "${BASE_URL}/${TAR_NAME}")
echo "$HEAD"
EXPECTED=$(awk 'tolower($1) == "content-length:" {gsub("\r", ""); n = $2} END {print n + 0}' <<<"$HEAD")
LISTED=$(awk '$3 ~ /^[0-9]+$/ {s += $3; n++} END {printf "%d bytes in %d files", s, n}' "${DL}/${LIST_NAME}")
log "server: ${EXPECTED} bytes; listing: ${LISTED}"
[[ "$EXPECTED" -gt "${SGDP_MIN_BYTES:-1000000000}" ]] || { log "ERROR: implausible size ${EXPECTED} -- is the URL right?"; exit 1; }

# Space for what is still to come: the rest of the tar, plus the extraction.
HAVE=0
if [[ -f "$TAR" ]]; then HAVE=$(stat -c%s "$TAR"); fi
NEED_KB=$(( (EXPECTED - HAVE + EXPECTED) / 1024 + 1048576 ))
FREE_KB=$(df -Pk "$SGDP_RAW_DIR" | awk 'NR == 2 {print $4}')
log "space: need ~$(( NEED_KB / 1048576 )) GiB, df shows $(( FREE_KB / 1048576 )) GiB free"
if [[ "$FREE_KB" -lt "$NEED_KB" && "${SGDP_SKIP_SPACE_CHECK:-0}" != "1" ]]; then
    log "ERROR: not enough space (SGDP_SKIP_SPACE_CHECK=1 to go ahead anyway)"; exit 1
fi

# 2. Resumable download.
attempt=0
while :; do
    have=0
    if [[ -f "$TAR" ]]; then have=$(stat -c%s "$TAR"); fi
    [[ "$have" -lt "$EXPECTED" ]] || break
    attempt=$(( attempt + 1 ))
    if [[ "$attempt" -gt "$MAX_ATTEMPTS" ]]; then
        log "ERROR: stopped after ${MAX_ATTEMPTS} attempts at ${have} of ${EXPECTED} bytes; resubmit to resume"
        exit 1
    fi
    log "attempt ${attempt}: have ${have} of ${EXPECTED} bytes"
    t0=$SECONDS
    if curl -fsS ${TLS[@]+"${TLS[@]}"} -C - --speed-limit 10000 --speed-time 600 \
            -o "$TAR" "${BASE_URL}/${TAR_NAME}"; then
        :
    else
        log "curl exited $? -- resuming in 60 s"
        sleep 60
    fi
    now=$(stat -c%s "$TAR"); dt=$(( SECONDS - t0 ))
    log "  +$(( now - have )) bytes in ${dt} s ($(( (now - have) / (dt > 0 ? dt : 1) / 1024 )) KiB/s)"
done
SIZE=$(stat -c%s "$TAR")
[[ "$SIZE" -eq "$EXPECTED" ]] || { log "ERROR: ${TAR} is ${SIZE} bytes, server says ${EXPECTED}"; exit 1; }
log "download complete: ${SIZE} bytes"

# 3. Vet every member before anything is written outside download/.
tar -tf "$TAR" > "${DL}/members.txt"
tar -tvf "$TAR" > "${DL}/members.tv.txt"
BADTYPE=$(awk 'substr($1, 1, 1) != "-"' "${DL}/members.tv.txt" | wc -l)
BADNAME=$(awk '!/^[A-Za-z0-9][A-Za-z0-9._-]*$/ || /\.\./' "${DL}/members.txt" | wc -l)
if [[ "$BADTYPE" -gt 0 || "$BADNAME" -gt 0 ]]; then
    log "ERROR: refusing to extract: ${BADTYPE} member(s) not regular files, ${BADNAME} with unsafe names:"
    awk 'substr($1, 1, 1) != "-"' "${DL}/members.tv.txt" | head
    awk '!/^[A-Za-z0-9][A-Za-z0-9._-]*$/ || /\.\./' "${DL}/members.txt" | head
    exit 1
fi
if ! diff <(awk '{print $NF}' "${DL}/${LIST_NAME}" | sort) <(sort "${DL}/members.txt") > "${DL}/list_vs_tar.diff"; then
    log "WARNING: the tar's members differ from the published listing (${DL}/list_vs_tar.diff):"
    head -20 "${DL}/list_vs_tar.diff"
fi
log "$(wc -l < "${DL}/members.txt") members, all regular files with bare names"

# 4. Extract into a fresh directory, verify, then publish.
rm -rf "${PART:?}"
mkdir -p "$PART"
tar -xf "$TAR" -C "$PART" --no-same-owner --no-same-permissions
ok=0; bad=0; nosum=0
while read -r f; do
    case "$f" in *.md5sum) continue ;; esac
    if [[ -f "${PART}/${f}.md5sum" ]]; then
        want=$(awk 'NR == 1 {print $1}' "${PART}/${f}.md5sum")
        got=$(md5sum "${PART}/${f}" | awk '{print $1}')
        if [[ "$want" == "$got" ]]; then
            ok=$(( ok + 1 ))
        else
            bad=$(( bad + 1 )); log "MD5 MISMATCH ${f}: expected ${want}, got ${got}"
        fi
    else
        nosum=$(( nosum + 1 )); log "no .md5sum shipped for ${f}"
    fi
done < "${DL}/members.txt"
NVCF=$(awk '/\.vcf\.gz$/' "${DL}/members.txt" | wc -l)
log "md5: ${ok} match, ${bad} mismatch, ${nosum} without a checksum; ${NVCF} per-sample VCFs"
if [[ "$bad" -gt 0 ]]; then
    log "ERROR: checksum mismatch -- ${PART} and the tar are kept for inspection"; exit 1
fi
{
    echo "fetched   $(date '+%F %T') from ${BASE_URL}/${TAR_NAME}"
    echo "tar_bytes ${SIZE}"
    echo "members   $(wc -l < "${DL}/members.txt")"
    echo "vcfs      ${NVCF}"
    echo "md5_ok    ${ok} (mismatch ${bad}, no checksum ${nosum})"
} > "${PART}/DONE"
mv "$PART" "$OUT"
log "done -> ${OUT}"
log "the tar is kept at ${TAR}; delete it yourself once the SGDP prep has read ${OUT}"
