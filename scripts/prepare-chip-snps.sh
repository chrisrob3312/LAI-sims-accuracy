#!/usr/bin/env bash
# Convert an Illumina manifest CSV (or any TSV with chr + hg38 pos) into
# per-chr SNP-position lists usable by --extract in 3b_chip.
#
# Illumina publishes manifest CSVs at their product pages; download the
# hg38 build version for GSA / Omni / MEGA etc. Example fetch:
#   wget "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/infinium-global-screening-array-24-v3-0-a1-manifest-file-csv.zip"
#
# Manifest CSVs are ~30-60MB with a header block, then per-marker rows.
# Fields we care about: Chr, MapInfo (== hg38 position for the b38 file).
#
# Output: reference_ids/chip_${CHIP_NAME}/chr${1..22}.snps  (chr:pos per line)
#
# Usage:
#   bash scripts/prepare-chip-snps.sh --manifest /path/to/GSA_manifest.csv --name gsa_v3
#   bash scripts/prepare-chip-snps.sh --manifest /path/to/Omni_manifest.csv --name omni25

set -euo pipefail

MANIFEST=""
CHIP_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest) MANIFEST="$2"; shift 2 ;;
        --name)     CHIP_NAME="$2"; shift 2 ;;
        *) echo "unknown arg: $1"; exit 1 ;;
    esac
done

[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || { echo "ERROR: --manifest FILE required"; exit 1; }
[[ -n "$CHIP_NAME" ]] || { echo "ERROR: --name STRING required"; exit 1; }

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${REPO_DIR}/reference_ids/chip_${CHIP_NAME}"
mkdir -p "$OUTDIR"

echo "[prepare-chip] manifest=$MANIFEST  name=$CHIP_NAME  out=$OUTDIR"

# Detect where the manifest data starts (Illumina headers vary in length)
# The [Assay] section is followed by column headers including "Chr" and "MapInfo".
DATA_START=$(awk -F',' 'toupper($1)=="ILMNID" || toupper($1)=="ILLMNID" {print NR; exit}' "$MANIFEST")
if [[ -z "$DATA_START" ]]; then
    DATA_START=$(grep -in "^Chr,\|,Chr," "$MANIFEST" | head -1 | cut -d: -f1)
fi
[[ -n "$DATA_START" ]] || { echo "ERROR: could not find data-header row in manifest"; exit 1; }
echo "[prepare-chip] data header on line $DATA_START"

# Locate Chr and MapInfo columns
HDR=$(sed -n "${DATA_START}p" "$MANIFEST")
CHR_COL=$(echo "$HDR" | awk -F',' '{for(i=1;i<=NF;i++) if(toupper($i)=="CHR") print i; exit}')
POS_COL=$(echo "$HDR" | awk -F',' '{for(i=1;i<=NF;i++) if(toupper($i)=="MAPINFO") print i; exit}')
[[ -n "$CHR_COL" && -n "$POS_COL" ]] || { echo "ERROR: no Chr/MapInfo columns"; exit 1; }
echo "[prepare-chip] Chr col=$CHR_COL  MapInfo col=$POS_COL"

# Emit per-chr SNP lists (chr:pos, matching the snp_id convention in 3b's
# panel BCFs). Only autosomes 1..22.
awk -F',' -v cc="$CHR_COL" -v pc="$POS_COL" -v out="$OUTDIR" -v start="$DATA_START" '
    NR <= start { next }
    {
        chr = $cc
        pos = $pc
        # strip whitespace + quotes
        gsub(/["\r ]/, "", chr)
        gsub(/["\r ]/, "", pos)
        if (chr ~ /^[0-9]+$/ && pos ~ /^[0-9]+$/ && chr+0 >= 1 && chr+0 <= 22) {
            print chr ":" pos > (out "/chr" chr+0 ".snps")
        }
    }
' "$MANIFEST"

# Report per-chr counts
echo "[prepare-chip] per-chr SNP counts:"
for c in {1..22}; do
    f="$OUTDIR/chr${c}.snps"
    [[ -f "$f" ]] && printf "  chr%-2d  %d\n" $c $(wc -l < "$f") || printf "  chr%-2d  (none)\n" $c
done

total=$(wc -l $OUTDIR/chr*.snps 2>/dev/null | tail -1 | awk '{print $1}')
echo "[prepare-chip] total SNPs across autosomes: $total"
echo "[prepare-chip] use with 3b_chip:  CHIP_SNPS_TPL=$OUTDIR/chr%s.snps  CHIP_NAME=$CHIP_NAME  sbatch 3b_chip_..."
