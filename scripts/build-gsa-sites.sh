#!/usr/bin/env bash
# ============================================================================
# build-gsa-sites.sh
#
# Thin cluster wrapper over build-chip-manifest-sites.py for the Illumina GSA.
# Run this on mhgcp, where the two reference FASTAs live; without them the
# python script writes only provisional lists and exits 3 on purpose.
#
# WHICH MANIFEST FILES THIS NEEDS, AND WHY THESE ONES.
# Document 20 step 9: the two builds must come from the vendor's own manifests
# for the SAME array version, joined on the probe ID. Of the files downloaded
# from Illumina, only GSA-24 v1.0 has both:
#
#   GSA-24v1-0_C2.hg38.annotated.txt   v1.0, GRCh38, carries Name/Chr/MapInfo/Alleles
#   GSA-24v1-0_C1.bed                  v1.0, GRCh37, carries Name + position
#   GSA-24v1-0_C1_b150_rsids.txt       v1.0, probe id -> dbSNP b150 rsID
#
# GSA-24 v2.0 (A1) was downloaded in GRCh37 only (GSA-24v2-0_A1.hg19.annotated
# .txt and GSA-24v2-0_A1.bed are both hg19), so v2.0 CANNOT be built as a pair
# until the vendor's GRCh38 A1 manifest is downloaded too. Do not substitute a
# liftover of the hg19 A1 file: that is exactly the down-lift Document 20
# forbids, because it would hand the pipeline only the sites that map back
# cleanly and hide the ones its own liftover drops or flips.
#
# The .bpm (GSA-24v1-0_C1.bpm, build 37) is Illumina's binary manifest. It is
# not read here: it carries RefStrand, but this script determines strand
# empirically from the reference sequence instead, which is both independent of
# the vendor's convention and the same method 1b/5a already use.
#
# Usage:
#   MANIFEST_DIR=~/manifests scripts/build-gsa-sites.sh
#   MANIFEST_DIR=~/manifests OUT=reference_ids/chip_gsa scripts/build-gsa-sites.sh
# ============================================================================
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST_DIR="${MANIFEST_DIR:-${HOME}/manifests}"
OUT="${OUT:-${REPO_DIR}/reference_ids/chip_gsa}"

HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
HG19_FA="${HG19_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg19/hg19.fa}"

# The manifest zips unpack to different directory names depending on how they
# were extracted, so locate each file by name rather than by path.
find_one () {   # pattern -> the single matching path, or abort
    local pat="$1" hits n
    hits=$(find "$MANIFEST_DIR" -type f -name "$pat" 2>/dev/null | sort)
    n=$(printf '%s\n' "$hits" | awk 'NF' | wc -l)
    if [[ "$n" -eq 0 ]]; then
        echo "ERROR: no file matching '$pat' under $MANIFEST_DIR" >&2
        echo "       unzip the Illumina GSA manifest archives there first." >&2
        exit 2
    fi
    if [[ "$n" -gt 1 ]]; then
        echo "ERROR: $n files match '$pat' under $MANIFEST_DIR:" >&2
        printf '%s\n' "$hits" | sed 's/^/         /' >&2
        echo "       keep one, or set MANIFEST_DIR to a directory with one copy." >&2
        exit 2
    fi
    printf '%s\n' "$hits"
}

HG38_ANNOT=$(find_one 'GSA-24v1-0_C2.hg38.annotated.txt')
HG19_BED=$(find_one 'GSA-24v1-0_C1.bed')
RSIDS=$(find_one 'GSA-24v1-0_C1_b150_rsids.txt')

for f in "$HG38_FA" "${HG38_FA}.fai" "$HG19_FA" "${HG19_FA}.fai"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

echo "[gsa] manifests : $MANIFEST_DIR"
echo "[gsa]   hg38     : $(basename "$HG38_ANNOT")"
echo "[gsa]   hg19     : $(basename "$HG19_BED")"
echo "[gsa]   rsids    : $(basename "$RSIDS")"
echo "[gsa] out       : $OUT"

# -I so an attacker-planted module beside the manifests cannot be imported.
python3 -I "${REPO_DIR}/scripts/build-chip-manifest-sites.py" \
    --chip gsa \
    --hg38-annot "$HG38_ANNOT" \
    --hg19-bed   "$HG19_BED" \
    --rsids      "$RSIDS" \
    --hg38-fasta "$HG38_FA" \
    --hg19-fasta "$HG19_FA" \
    --out        "$OUT"

cat <<EOF

Crop the simulated panel to GSA density with the existing chip machinery:

  # hg38 cell
  CHIP_NAME=gsa_hg38 CHIP_SNPS_TPL=${OUT}/hg38/chr%d.snps \\
    bash 3b_chip_wgs-rfmix-jointcall_clm.sh

  # hg19 cell -- positions are the chip's NATIVE hg19, never a down-lift, and
  # ${OUT}/manifests/gsa.hg19.sites.tsv carries grch38_site so the
  # intake can check each hg19 genotype against the GRCh38 truth.
  CHIP_NAME=gsa_hg19 CHIP_SNPS_TPL=${OUT}/hg19/chr%d.snps \\
    bash 3b_chip_wgs-rfmix-jointcall_clm.sh

Read ${OUT}/build_report.tsv before using either: it counts every
dropped probe by reason, the direct/revcomp split (which measures the vendor's
strand convention rather than assuming one) and the palindromic census.
EOF
