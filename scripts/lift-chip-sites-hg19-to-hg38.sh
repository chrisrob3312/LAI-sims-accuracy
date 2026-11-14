#!/usr/bin/env bash
# ============================================================================
# lift-chip-sites-hg19-to-hg38.sh
#
# Lift a UCSC-derived chip site list from hg19 to hg38, keeping alleles.
#
# WHY NOT `liftOver` ON A BED. A BED carries positions only, so a bare
# position lift cannot tell you when the lift landed on a base whose
# REF/ALT relationship changed -- and that is common enough to matter.
# Instead this builds a minimal VCF (REF read from the hg19 FASTA, ALT taken
# from UCSC's `observed`), lifts with Picard RECOVER_SWAPPED_REF_ALT, and
# reports what was swapped or rejected. Same tool and chain the rest of the
# pipeline already uses (1a MXB, 1b SGDP, 5a panel), so behaviour matches.
#
# STRAND. UCSC gives `strand` and `observed` per the source array manifest,
# so `observed` is NOT guaranteed to be on the genome + strand. This script
# resolves each site empirically: read the hg19 reference base, and if
# neither observed allele matches it, try the reverse complement. Every site
# is classified and counted, so the strand convention is measured rather than
# assumed:
#     direct      one observed allele equals the hg19 reference base
#     revcomp     matched only after reverse-complementing
#     unresolved  neither matched -- excluded, and reported
#
# PALINDROMES. A/T and C/G sites cannot be strand-resolved by allele identity
# at all: both the allele pair and its reverse complement match. They are
# carried through but FLAGGED in the output (palindromic=1) so a downstream
# merge can apply an explicit policy. ~15-16% of Affy sites are palindromic,
# so silently dropping them would discard a sixth of the array.
#
# COORDINATES. The input .sites.tsv is already 1-based (fetch-ucsc-chip-sites
# converted from UCSC's 0-based chromEnd). VCF is 1-based, so no further
# conversion happens here. Do not re-apply an off-by-one "fix".
#
# Usage:
#   scripts/lift-chip-sites-hg19-to-hg38.sh reference_ids/chip_snparrayaffy6_hg19
#   scripts/lift-chip-sites-hg19-to-hg38.sh reference_ids/chip_affy500k_hg19
# ============================================================================
set -euo pipefail

SRC="${1:-}"
[[ -n "$SRC" && -d "$SRC" ]] || { echo "usage: $0 <chip_<name>_hg19 dir>" >&2; exit 1; }

HG19_FA="${HG19_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/hg19/hg19.fa}"
HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg19ToHg38.over.chain.gz}"
PICARD_XMX="${PICARD_XMX:-8g}"

DST="${DST:-${SRC%_hg19}_hg38}"
WORK="${WORK:-${TMPDIR:-/tmp}/liftchip.$$}"
mkdir -p "$DST" "$WORK"
trap 'rm -rf "$WORK"' EXIT

for f in "$HG19_FA" "${HG19_FA}.fai" "$HG38_FA" "${HG38_FA}.fai" "$CHAIN"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
HG38_DICT="${HG38_FA%.*}.dict"
[[ -s "$HG38_DICT" ]] || { echo "ERROR: missing $HG38_DICT" >&2; exit 1; }

echo "[lift] $SRC  ->  $DST"
echo "[lift] chain: $CHAIN"

# --- 1. one minimal VCF for all autosomes --------------------------------
# UCSC hg19 contigs are chr-prefixed, matching hg19.fa, so no renaming.
VCF="${WORK}/sites.hg19.vcf"
{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=PROBE,Number=1,Type=String,Description="Array probeset name">'
  echo '##INFO=<ID=RSID,Number=1,Type=String,Description="dbSNP rsID">'
  echo '##INFO=<ID=PAL,Number=1,Type=Integer,Description="1 if palindromic (A/T or C/G)">'
  echo '##INFO=<ID=STR,Number=1,Type=String,Description="strand resolution: direct or revcomp">'
  awk '{print "##contig=<ID=chr"$1">"}' <(seq 1 22)
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$VCF"

TOTAL=0; DIRECT=0; REVCOMP=0; UNRES=0; PAL=0
for c in $(seq 1 22); do
    T="${SRC}/chr${c}.sites.tsv"
    [[ -s "$T" ]] || { echo "[lift] WARNING: no $T"; continue; }
    # Pull every hg19 reference base for this chromosome in one samtools call.
    awk -F'\t' -v c="$c" 'NR>1 {print "chr"c":"$2"-"$2}' "$T" > "${WORK}/regions.chr${c}.txt"
    samtools faidx "$HG19_FA" -r "${WORK}/regions.chr${c}.txt" 2>/dev/null \
      | awk '/^>/ {next} {print toupper($0)}' > "${WORK}/refbase.chr${c}.txt"

    paste <(awk -F'\t' 'NR>1' "$T") "${WORK}/refbase.chr${c}.txt" \
      | awk -F'\t' -v OFS='\t' '
        function rc(b) { return b=="A"?"T": b=="T"?"A": b=="C"?"G": b=="G"?"C":"N" }
        {
          chrom=$1; pos=$2; name=$3; obs=$5; rsid=$6; ref=$7
          split(obs, a, "/")
          if (length(a) != 2 || ref=="" || ref=="N") { unres++; next }
          A1=a[1]; A2=a[2]
          # Palindromic: the pair equals its own reverse complement.
          pal = ((A1=="A"&&A2=="T")||(A1=="T"&&A2=="A")||(A1=="C"&&A2=="G")||(A1=="G"&&A2=="C")) ? 1 : 0
          strres=""
          if (A1==ref)      { alt=A2; strres="direct" }
          else if (A2==ref) { alt=A1; strres="direct" }
          else {
            R1=rc(A1); R2=rc(A2)
            if (R1==ref)      { alt=R2; strres="revcomp" }
            else if (R2==ref) { alt=R1; strres="revcomp" }
            else { unres++; next }
          }
          if (alt==ref) { unres++; next }
          print "chr"chrom, pos, (rsid=="" ? name : rsid), ref, alt, ".", ".", \
                "PROBE="name";RSID="rsid";PAL="pal";STR="strres
          tot++; if (strres=="direct") d++; else rv++
          if (pal) p++
        }
        END { print "#STATS", tot+0, d+0, rv+0, unres+0, p+0 > "/dev/stderr" }
      ' 2>"${WORK}/stats.chr${c}" >> "$VCF"

    read -r _ t d rv u p < "${WORK}/stats.chr${c}"
    TOTAL=$((TOTAL+t)); DIRECT=$((DIRECT+d)); REVCOMP=$((REVCOMP+rv))
    UNRES=$((UNRES+u)); PAL=$((PAL+p))
done

echo "[lift] built hg19 VCF: ${TOTAL} sites"
echo "[lift]   strand direct      : ${DIRECT}"
echo "[lift]   strand revcomp     : ${REVCOMP}"
echo "[lift]   unresolved (dropped): ${UNRES}"
echo "[lift]   palindromic        : ${PAL}"
if [[ "$REVCOMP" -gt "$DIRECT" ]]; then
    echo "[lift]   NOTE: most sites needed reverse-complementing, i.e. this"
    echo "[lift]   array's 'observed' is reported on the minus strand."
fi

bcftools sort -Oz -o "${WORK}/sites.hg19.vcf.gz" "$VCF"
bcftools index -t "${WORK}/sites.hg19.vcf.gz"

# --- 2. lift --------------------------------------------------------------
echo "[lift] Picard LiftoverVcf hg19 -> hg38"
picard "-Xmx${PICARD_XMX}" LiftoverVcf \
    I="${WORK}/sites.hg19.vcf.gz" \
    O="${WORK}/sites.hg38.vcf.gz" \
    CHAIN="$CHAIN" R="$HG38_FA" \
    REJECT="${DST}/rejected.hg19.vcf.gz" \
    RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true CREATE_INDEX=false

N_OUT=$(bcftools view -H "${WORK}/sites.hg38.vcf.gz" | wc -l)
N_REJ=$(zcat "${DST}/rejected.hg19.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
PCT=$(awk -v r="$N_REJ" -v t="$TOTAL" 'BEGIN{printf "%.2f", (t>0?100*r/t:0)}')
echo "[lift] lifted ${N_OUT}, rejected ${N_REJ} (${PCT}%)"
[[ $(awk -v p="$PCT" 'BEGIN{print (p>5)?1:0}') == 1 ]] && \
  echo "[lift]   WARNING: >5% rejected -- inspect ${DST}/rejected.hg19.vcf.gz"

# --- 3. emit per-chr lists in the same shape as the hg19 input -----------
TOT38=0
for c in $(seq 1 22); do
    bcftools view -H -r "chr${c}" "${WORK}/sites.hg38.vcf.gz" 2>/dev/null \
      | awk -F'\t' -v OFS='\t' -v c="$c" '
        BEGIN{ print "chrom","pos","name","ref","alt","rsId","palindromic","strand_res" > ("'"${DST}"'/chr" c ".sites.tsv") }
        { probe=rsid=pal=str=""
          n=split($8, kv, ";")
          for(i=1;i<=n;i++){ split(kv[i], x, "=")
            if(x[1]=="PROBE") probe=x[2]; else if(x[1]=="RSID") rsid=x[2]
            else if(x[1]=="PAL") pal=x[2]; else if(x[1]=="STR") str=x[2] }
          print c, $2 >> ("'"${DST}"'/chr" c ".snps.tmp")
          print c, $2, probe, $4, $5, rsid, pal, str >> ("'"${DST}"'/chr" c ".sites.tsv")
        }'
    if [[ -s "${DST}/chr${c}.snps.tmp" ]]; then
        awk '{print $1":"$2}' "${DST}/chr${c}.snps.tmp" | sort -u > "${DST}/chr${c}.snps"
        rm -f "${DST}/chr${c}.snps.tmp"
    else
        : > "${DST}/chr${c}.snps"
    fi
    n=$(wc -l < "${DST}/chr${c}.snps"); TOT38=$((TOT38+n))
    printf '  chr%-2s %8d\n' "$c" "$n"
done

PAL38=$(awk -F'\t' 'FNR>1 && $7==1' "${DST}"/chr*.sites.tsv | wc -l)
echo "[lift] hg38 sites: ${TOT38}  (palindromic ${PAL38})"
cat <<EOF

Use with the chip machinery:
  CHIP_NAME=$(basename "$DST") \\
  CHIP_SNPS_TPL=${DST}/chr%d.snps \\
  bash 3b_chip_wgs-rfmix-jointcall_clm.sh

The .sites.tsv now carries ref/alt verified against each build plus a
palindromic flag, so a cross-array merge can apply an explicit policy instead
of discovering flipped genotypes downstream.
EOF
