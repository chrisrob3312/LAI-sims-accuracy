#!/usr/bin/env bash
# ============================================================================
# fetch-ucsc-chip-sites.sh
#
# Build chip site lists from UCSC's genotyping-array tracks, for cropping
# simulated WGS down to array density in the genotyping_pipeline benchmark.
#
# COORDINATE TRAP -- the reason this script exists rather than a one-liner.
# UCSC database tables are 0-BASED HALF-OPEN. Our .snps lists, VCF/BCF and
# PLINK .bim are all 1-BASED. For a SNP, UCSC stores:
#     chromStart=564620  chromEnd=564621
# and the 1-based position is chromEnd, NOT chromStart. Using chromStart
# shifts every site by one base. Nothing downstream errors: the crop simply
# intersects almost nothing, or intersects neighbouring bases and produces a
# site set that looks superficially fine. This script always takes chromEnd.
#
# WHY NOT POSITION-ONLY. The existing reference_ids/chip_gsa_jessica lists are
# bare chr:pos. That is enough to crop, but it cannot support strand or
# palindromic-SNP handling -- which is precisely what a build-harmonisation
# benchmark is meant to stress. UCSC's tables carry `strand` and `observed`
# (e.g. "C/T"), so this emits BOTH:
#     chr<N>.snps        chr:pos          (drop-in for CHIP_SNPS_TPL)
#     chr<N>.sites.tsv   chrom pos name strand observed rsId   (allele-aware)
# The TSV is what you need to decide A/T and C/G sites, which cannot be
# strand-resolved by allele alone and are the classic silent-corruption source
# when merging data from two different arrays.
#
# TABLES (UCSC database names):
#   snpArrayAffy6        Affymetrix Genome-Wide Human SNP Array 6.0 (~909k)
#   snpArrayAffy5        Affymetrix SNP Array 5.0
#   snpArrayAffy250Nsp   Affymetrix Mapping 500K, Nsp half  (~262k)
#   snpArrayAffy250Sty   Affymetrix Mapping 500K, Sty half  (~238k)
#   (500K = Nsp + Sty together; pass BOTH and merge, see --merge-500k)
#
# BUILD. These tracks are curated on hg19. An hg38 copy may or may not exist
# for a given table; the script probes and tells you rather than assuming. If
# hg38 is absent, lift the hg19 list with the chain already on this cluster --
# but lift the ALLELE-BEARING TSV, never a bare position list, or you lose the
# ability to detect ref/alt swaps introduced by the lift itself.
#
# Usage:
#   scripts/fetch-ucsc-chip-sites.sh --list                      # what exists
#   scripts/fetch-ucsc-chip-sites.sh snpArrayAffy6               # hg19
#   DB=hg38 scripts/fetch-ucsc-chip-sites.sh snpArrayAffy6
#   scripts/fetch-ucsc-chip-sites.sh --merge-500k                # Nsp + Sty
# ============================================================================
set -euo pipefail

DB="${DB:-hg19}"
BASE="${BASE:-https://hgdownload.soe.ucsc.edu/goldenPath}"
REFS="${REFS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/reference_ids}"
CACHE="${CACHE:-${HOME}/.cache/ucsc_chip}"
mkdir -p "$CACHE"

probe_table () {   # table -> prints OK/MISSING without downloading the body
    local t="$1"
    local code
    code=$(curl -sS -o /dev/null -w '%{http_code}' --max-time 30 \
             "${BASE}/${DB}/database/${t}.txt.gz" 2>/dev/null || echo 000)
    printf '  %-22s %s (HTTP %s)\n' "$t" \
        "$([[ "$code" == 200 ]] && echo OK || echo MISSING)" "$code"
    [[ "$code" == 200 ]]
}

if [[ "${1:-}" == "--list" ]]; then
    echo "[ucsc] probing array tables in ${DB}:"
    for t in snpArrayAffy6 snpArrayAffy5 snpArrayAffy250Nsp snpArrayAffy250Sty \
             snpArrayIllumina650 snpArrayIllumina550 snpArrayIllumina300 \
             snpArrayIlluminaHuman1MDuo snpArrayIlluminaHumanOmni1Quad; do
        probe_table "$t" || true
    done
    echo
    echo "[ucsc] a MISSING table in ${DB} usually means the track is curated on"
    echo "[ucsc] hg19 only. Re-run with DB=hg19, then lift the .sites.tsv."
    exit 0
fi

# --merge-500k: the Affy 500K set is two physical arrays; the site list is
# their union, so neither half alone is "the 500K array".
MERGE500K=""
if [[ "${1:-}" == "--merge-500k" ]]; then
    MERGE500K=1; set -- snpArrayAffy250Nsp snpArrayAffy250Sty
fi
[[ $# -gt 0 ]] || { echo "usage: $0 [--list|--merge-500k|<table> ...]" >&2; exit 1; }

OUTNAME="${OUTNAME:-}"
if [[ -n "$MERGE500K" ]]; then OUTNAME="affy500k_${DB}"
elif [[ -z "$OUTNAME" ]]; then OUTNAME="$(echo "$1" | tr 'A-Z' 'a-z')_${DB}"; fi
OUT="${REFS}/chip_${OUTNAME}"
mkdir -p "$OUT"
echo "[ucsc] db=${DB}  tables=$*  ->  ${OUT}"

ALL="${CACHE}/${OUTNAME}.combined.tsv"
: > "$ALL"
for T in "$@"; do
    GZ="${CACHE}/${DB}.${T}.txt.gz"
    if [[ ! -s "$GZ" ]]; then
        echo "[ucsc] downloading ${T} (${DB})"
        curl -sSLf --retry 3 -o "$GZ" "${BASE}/${DB}/database/${T}.txt.gz" || {
            echo "ERROR: could not download ${T} from ${DB}." >&2
            echo "       Run '$0 --list' to see which tables exist in this build." >&2
            rm -f "$GZ"; exit 2; }
    else
        echo "[ucsc] using cached ${GZ}"
    fi
    # Columns: bin chrom chromStart chromEnd name score strand observed rsId
    # chromEnd is the 1-based position. Autosomes only; the panel is chr1-22.
    zcat "$GZ" | awk -F'\t' -v OFS='\t' '
        $2 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
            c=$2; sub(/^chr/,"",c)
            print c, $4, $5, $7, $8, $9        # chrom pos name strand observed rsId
        }' >> "$ALL"
    echo "[ucsc]   ${T}: $(zcat "$GZ" | wc -l) rows in table"
done

# Dedup on (chrom,pos): the 500K halves can share sites, and a duplicated
# position would be cropped twice.
sort -u -k1,1n -k2,2n "$ALL" -o "$ALL"
echo "[ucsc] combined autosomal sites: $(wc -l < "$ALL")"

TOTAL=0
for c in $(seq 1 22); do
    awk -F'\t' -v c="$c" '$1==c {print $1":"$2}' "$ALL" | sort -u > "${OUT}/chr${c}.snps"
    awk -F'\t' -v c="$c" -v OFS='\t' 'BEGIN{print "chrom","pos","name","strand","observed","rsId"}
        $1==c {print}' "$ALL" > "${OUT}/chr${c}.sites.tsv"
    n=$(wc -l < "${OUT}/chr${c}.snps"); TOTAL=$((TOTAL+n))
    printf '  chr%-2s %8d\n' "$c" "$n"
done
echo "[ucsc] wrote ${TOTAL} sites across 22 autosomes -> ${OUT}"

# Palindromic census. A/T and C/G sites cannot be strand-resolved from alleles
# alone, so they are where a two-array merge silently corrupts genotypes.
PAL=$(awk -F'\t' 'NR>1 && ($5=="A/T" || $5=="T/A" || $5=="C/G" || $5=="G/C")' \
        "${OUT}"/chr*.sites.tsv | wc -l)
echo "[ucsc] palindromic (A/T, C/G) sites: ${PAL} ($(awk -v p="$PAL" -v t="$TOTAL" \
        'BEGIN{printf "%.1f", (t>0?100*p/t:0)}')%)"
echo "[ucsc]   ^ these need an explicit policy (drop, or resolve by allele"
echo "[ucsc]     frequency) before merging two arrays. They are the classic"
echo "[ucsc]     silent-corruption source in cross-platform merges."

cat <<EOF

Use with the existing chip machinery:
  CHIP_NAME=${OUTNAME} \\
  CHIP_SNPS_TPL=${OUT}/chr%d.snps \\
  bash 3b_chip_wgs-rfmix-jointcall_clm.sh

NOTE ON BUILD: these coordinates are ${DB}. The simulated panel is hg38. If
${DB} is not hg38, lift chr*.sites.tsv (not chr*.snps) with
  /storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg19ToHg38.over.chain.gz
keeping the allele columns, so a ref/alt swap introduced by the lift is
detectable rather than silent.
EOF
