#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=1b_sgdp_liftover
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=04:00:00
#SBATCH --time-min=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_sgdp_liftover_chr%a_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1b_sgdp_liftover_chr%a_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- Path C stage 1: SGDP acquisition + liftover + normalization.
# ============================================================================
# Downloads SGDP phased per-chr VCFs from the Reich lab public share
# (phased_data2021, which is hg19 / GRCh37), filters to the candidate samples
# produced by scripts/fetch-sgdp-anchors.py (union of all 7 super-pops),
# LIFTS hg19 -> hg38, splits multi-allelics, normalizes, and emits
# sgdp_prepped.chr${CHR}.bcf ready to feed into 1c_phasing-jointcall_clm.sh's
# bcftools merge alongside HGDP+1KG and MXB.
#
# IMPORTANT -- SGDP phased_data2021 is hg19, our panel is hg38:
#   The Reich lab "phased_data2021" release is aligned to hs37d5 (hg19/GRCh37).
#   Our HGDP+1KG+MXB panel is hg38. We therefore lift SGDP hg19 -> hg38 using
#   the SAME Picard recipe + reference/chain files already used for the MXB
#   liftover in 1a_prep-mxb-liftover_clm.sh. Filtering to the ~277 candidates
#   happens BEFORE liftover so Picard only processes a tiny panel (fast, low
#   heap). The order is: download -> sample-filter -> (rename contigs) ->
#   Picard LiftoverVcf hg19->hg38 -> norm -m- -> +fixref -m flip -d -> sort.
#
# NOTE -- "candidates", not the final homogeneous set:
#   SGDP_CANDIDATES is a population-majority candidate pool from SGDP metadata.
#   Which of these actually clear Q>=0.95 is decided later by
#   4a_build_homogeneity_panel_clm.sh over the merged panel.
#
# Why re-phase from scratch (downstream, in 1c): RFMix assumes reference
# haplotypes are jointly phased. Appending an already-phased SGDP release to
# an already-phased HGDP+1KG+MXB panel creates cohort-specific switch-error
# patterns that RFMix would misinterpret as ancestry switches. Re-phasing
# after merge is the only way to keep LAI accuracy honest.
#
# URLs -- Reich lab public share (hg19 phased_data2021):
#     https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/
#   Set SGDP_URL_TPL to the per-chr pattern using __CHR__ for 1..22. CONFIRM
#   the exact filename stem from the directory listing before first run:
#     curl -skL ".../phased_data2021/" | grep -oiE 'href="[^"]+"'
#     curl -skIL "${SGDP_URL_TPL//__CHR__/22}"      # expect HTTP 200
#   If the release ships as ONE genome-wide file instead of per-chr, set
#   SGDP_COMBINED_URL (no __CHR__); the script downloads it once and subsets
#   per chr with `bcftools view -r`.
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"

# Where SGDP prepped BCFs land (analogous to mxb_lifted.chr${CHR}.bcf).
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
SGDP_STAGE="${SGDP_STAGE:-${PROJECT_ROOT}/00_sgdp_raw}"   # cache for downloads

# Reich lab public share URL template (hg19 phased_data2021). __CHR__ -> 1..22.
# CONFIRM the filename stem from the directory listing before running:
#   curl -skL "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/" \
#     | grep -oiE 'href="[^"]+"'
# CONFIRMED 2026-10-02 from the directory listing: phased_data2021 ships as
# per-chr BCFs named chr.sgdp.pub.{1..22}.bcf (+ .csi + .stats), hg19/hs37d5.
SGDP_URL_TPL="${SGDP_URL_TPL:-https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/chr.sgdp.pub.__CHR__.bcf}"

# If the release ships as ONE genome-wide file (no per-chr split), set this to
# its URL (leave __CHR__ out). The script then downloads it once (shared across
# array tasks) and subsets per chr with `bcftools view -r`. Leave empty to use
# the per-chr SGDP_URL_TPL above.
SGDP_COMBINED_URL="${SGDP_COMBINED_URL:-}"

# Combined candidate list (union of all 7 super-pops from fetch-sgdp-anchors.py).
# One sample ID per line, in BCF-native IDs. This gates who gets pulled from
# SGDP into the panel. NOTE: candidates, NOT anchors and NOT the final
# homogeneous set -- 4a_k7 discovers anchors and decides homogeneity downstream.
# (SGDP_ANCHORS kept as a back-compat alias for the old env var name.)
SGDP_CANDIDATES="${SGDP_CANDIDATES:-${SGDP_ANCHORS:-${REPO_DIR}/reference_ids/sgdp_combined_candidates.txt}}"

# --- hg19 -> hg38 liftover (identical recipe + files to 1a MXB liftover) ----
# hg38 target reference for Picard + bcftools norm.
HG38_FA="${HG38_FA:-/storage/atkinson/shared_resources/reference/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
# UCSC hg19->hg38 chain (chr-prefixed contigs expected on both sides).
CHAIN="${CHAIN:-/storage/atkinson/shared_resources/reference/genetic_maps/liftover/hg19ToHg38.over.chain.gz}"
# Picard heap. SGDP candidate panel is tiny (~277 samples) so 1a's 40g is plenty.
PICARD_XMX="${PICARD_XMX:-24g}"

set -euo pipefail
CHR=${SLURM_ARRAY_TASK_ID:?must run as SLURM array job (sbatch --array=1-22 ...)}
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUTDIR" "$SGDP_STAGE" "${OUTDIR}/tmp/chr${CHR}"
TMPDIR="${OUTDIR}/tmp/chr${CHR}"

[[ -s "$SGDP_CANDIDATES" ]] || {
    echo "ERROR: missing $SGDP_CANDIDATES"
    echo "       Generate it with:"
    echo "         python3 scripts/fetch-sgdp-anchors.py --outdir reference_ids"
    echo "         cat reference_ids/sgdp_{sas,oceania,men,amr}_candidates.txt \\"
    echo "             | sort -u > reference_ids/sgdp_combined_candidates.txt"
    exit 1
}
[[ -s "$HG38_FA" ]] || { echo "ERROR: missing $HG38_FA"; exit 1; }
[[ -s "$CHAIN"   ]] || { echo "ERROR: missing hg19->hg38 chain $CHAIN"; exit 1; }
# Picard LiftoverVcf needs a .dict sidecar next to the target FASTA. Without it
# Picard fails only AFTER the slow BCF->VCF conversion, so check up front.
HG38_DICT="${HG38_FA%.*}.dict"
[[ -s "$HG38_DICT" ]] || { echo "ERROR: missing Picard sequence dictionary $HG38_DICT (make with: picard CreateSequenceDictionary R=$HG38_FA)"; exit 1; }

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
# Prevent system anaconda3 from leaking Python into subprocesses on some nodes
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ---------------------------------------------------------------------------
# 1. Download the per-chr SGDP phased BCF (hg19). Idempotent -- skip if cached.
#    phased_data2021 ships .bcf + .csi; grab both so we don't re-index.
# ---------------------------------------------------------------------------
PREPPED="${OUTDIR}/sgdp_prepped.chr${CHR}.bcf"
if [[ -s "$PREPPED" && -s "${PREPPED}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] $PREPPED exists, skipping"
    exit 0
fi

if [[ -n "$SGDP_COMBINED_URL" ]]; then
    # Single genome-wide file mode: download once (array tasks share it), then
    # subset to this chr. SGDP hg19 contigs are bare-numeric, so region = CHR.
    SGDP_URL="$SGDP_COMBINED_URL"
    RAW="${SGDP_STAGE}/sgdp.phased.hg19.allchr.bcf"
else
    SGDP_URL="${SGDP_URL_TPL//__CHR__/${CHR}}"
    RAW="${SGDP_STAGE}/chr.sgdp.pub.${CHR}.bcf"
fi

if [[ -s "$RAW" ]]; then
    # CACHE HIT. Say so loudly and print size+mtime. A silent skip here is
    # exactly how the chr13 PEL cache-poisoning bug hid: a truncated file from
    # a killed job was reused forever and every downstream artefact was wrong,
    # with nothing in the log to show the file had not been re-fetched.
    echo "[$(date +%T)] [chr${CHR}] CACHE HIT: reusing $RAW"
    echo "[$(date +%T)] [chr${CHR}]   size=$(stat -c%s "$RAW") bytes" \
         "($(du -h "$RAW" | cut -f1))  mtime=$(stat -c%y "$RAW" | cut -d. -f1)"
    echo "[$(date +%T)] [chr${CHR}]   delete it to force a fresh download"
else
    echo "[$(date +%T)] [chr${CHR}] downloading SGDP (hg19): $SGDP_URL"
    DL_T0=$SECONDS
    # -k: the cluster CA bundle doesn't cover the Harvard chain. Payload is
    # public academic data; integrity relies on size sanity + bcftools header
    # validation below, not TLS.
    curl -kSLf --retry 3 --retry-delay 30 -o "$RAW" "$SGDP_URL" || {
        echo "ERROR: SGDP download failed. Confirm URL:"
        echo "       $SGDP_URL"
        echo "       curl -skIL '$SGDP_URL'   # expect HTTP 200"
        rm -f "$RAW"
        exit 1
    }
    DL_SECS=$(( SECONDS - DL_T0 ))
    DL_SIZE=$(stat -c%s "$RAW")
    echo "[$(date +%T)] [chr${CHR}] downloaded ${DL_SIZE} bytes" \
         "($(du -h "$RAW" | cut -f1)) in ${DL_SECS}s"
    if [[ "$DL_SIZE" -lt 1000000 ]]; then
        echo "ERROR: downloaded file suspiciously small (${DL_SIZE} bytes)."
        echo "       First 200 bytes (usually an HTML error page):"
        head -c 200 "$RAW"
        echo
        rm -f "$RAW"
        exit 1
    fi
fi
# Validate the BCF is readable BEFORE any downstream step trusts it. A
# truncated-but-large file passes the size check; only bcftools can tell.
if ! bcftools view -h "$RAW" >/dev/null 2>&1; then
    echo "ERROR: $RAW exists but bcftools cannot read its header -- the file is"
    echo "       truncated or corrupt (a killed download, or a poisoned cache)."
    echo "       Delete it and re-run:  rm -f '$RAW' '${RAW}.csi'"
    exit 1
fi
echo "[$(date +%T)] [chr${CHR}] BCF header OK, $(bcftools query -l "$RAW" | wc -l) samples in source"
# Index (fetch remote .csi for per-chr; build locally for combined/subsetted).
if [[ -s "${RAW}.csi" ]]; then
    echo "[$(date +%T)] [chr${CHR}] index: reusing cached ${RAW}.csi"
elif [[ -n "$SGDP_COMBINED_URL" ]] \
     || ! curl -kSLf --retry 3 -o "${RAW}.csi" "${SGDP_URL}.csi" 2>/dev/null; then
    echo "[$(date +%T)] [chr${CHR}] index: no remote .csi, building locally"
    rm -f "${RAW}.csi"
    bcftools index --threads "$THREADS" "$RAW"
else
    echo "[$(date +%T)] [chr${CHR}] index: fetched remote ${SGDP_URL}.csi"
fi

# ---------------------------------------------------------------------------
# 2. Filter to candidate samples (hg19), then lift hg19 -> hg38.
#    Filtering FIRST keeps Picard fast (only ~277 samples to lift).
# ---------------------------------------------------------------------------
KEEP_IDS="${TMPDIR}/sgdp_keep_ids.chr${CHR}.txt"
awk 'NF' "$SGDP_CANDIDATES" > "$KEEP_IDS"
N_ANCHORS=$(wc -l < "$KEEP_IDS")
echo "[$(date +%T)] [chr${CHR}] candidate list: ${N_ANCHORS} SGDP samples requested"

# Guard against the silent-empty-filter trap: if NONE of the candidate IDs are
# present in the BCF header, the sample-ID scheme doesn't match (e.g. anchors
# are S_Papuan-1 but the BCF uses LP6005... Illumina IDs). Fail loud with a
# hint rather than emitting an empty panel that 1c would silently merge.
PRESENT=$(bcftools query -l "$RAW" | grep -Fxf "$KEEP_IDS" | wc -l || true)
echo "[$(date +%T)] [chr${CHR}] ${PRESENT}/${N_ANCHORS} candidate IDs present in BCF header"
if [[ "$PRESENT" -eq 0 ]]; then
    echo "ERROR: none of the candidate IDs match the SGDP BCF sample names."
    echo "       The BCF uses a different ID scheme than reference_ids/sgdp_combined_candidates.txt."
    echo "       Inspect with:  bcftools query -l '$RAW' | head"
    echo "       SGDP carries 3 parallel IDs (ena SAMEA..., Illumina LP6005...,"
    echo "       reich S_Population-N). Remap the anchor list to the header scheme."
    exit 3
fi

# 2a. Region-subset (combined mode only) + sample filter, as a bgzipped VCF.
#     Picard's htsjdk rejects bcftools-emitted BCFs ('BCF magic header info
#     not found'), so we feed it VCF.gz (same workaround as 5a).
FILT_VCF="${TMPDIR}/sgdp.chr${CHR}.hg19.filt.vcf.gz"
REGION_ARGS=()
[[ -n "$SGDP_COMBINED_URL" ]] && REGION_ARGS=(-r "${CHR}")
bcftools view "${REGION_ARGS[@]}" -S "$KEEP_IDS" --force-samples \
    --threads "$THREADS" -Ou "$RAW" \
  | bcftools norm --threads "$THREADS" -m -any -Oz -o "$FILT_VCF"
bcftools index -t "$FILT_VCF"

# 2b. hg19 contigs in SGDP/hs37d5 are bare-numeric (1,2,..); the UCSC
#     hg19ToHg38 chain + Picard expect chr-prefixed. Rename if needed.
first_contig=$(bcftools view -h "$FILT_VCF" \
    | awk '/^##contig=<ID=/ {sub(/^##contig=<ID=/,""); sub(/[,>].*$/,""); print; exit}')
echo "[$(date +%T)] [chr${CHR}] SGDP first contig: ${first_contig:-<none>}"
LIFT_IN="$FILT_VCF"
if [[ "$first_contig" != chr* ]]; then
    echo "[$(date +%T)] [chr${CHR}] renaming numeric contigs -> chr-prefixed for liftover"
    RENAME="${TMPDIR}/rename_to_chr.chr${CHR}.txt"
    : > "$RENAME"
    for c in {1..22} X Y MT; do echo "$c chr$c" >> "$RENAME"; done
    LIFT_IN="${TMPDIR}/sgdp.chr${CHR}.hg19.chrnamed.vcf.gz"
    bcftools annotate --rename-chrs "$RENAME" --threads "$THREADS" \
        -Oz -o "$LIFT_IN" "$FILT_VCF"
    bcftools index -t "$LIFT_IN"
fi

# 2c. Picard LiftoverVcf hg19 -> hg38 (same recipe as 1a MXB liftover).
echo "[$(date +%T)] [chr${CHR}] Picard LiftoverVcf hg19 -> hg38"
LIFTED_VCF="${TMPDIR}/sgdp.chr${CHR}.hg38.lifted.vcf.gz"
picard "-Xmx${PICARD_XMX}" LiftoverVcf \
    I="$LIFT_IN" \
    O="$LIFTED_VCF" \
    CHAIN="$CHAIN" \
    R="$HG38_FA" \
    REJECT="${OUTDIR}/sgdp_prepped.chr${CHR}.rejected.vcf.gz" \
    RECOVER_SWAPPED_REF_ALT=true \
    WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=false \
    MAX_RECORDS_IN_RAM=200000

# 2d. norm (split multiallelics, left-align) + fixref (flip strand, drop
#     unfixable) against hg38, then sort -> PREPPED.bcf (hg38, ready for 1c).
echo "[$(date +%T)] [chr${CHR}] norm -m- + fixref + sort"
bcftools norm -f "$HG38_FA" -m- --threads "$THREADS" \
    -Ob -o "${TMPDIR}/sgdp.chr${CHR}.hg38.norm.bcf" "$LIFTED_VCF"
bcftools index --threads "$THREADS" "${TMPDIR}/sgdp.chr${CHR}.hg38.norm.bcf"
bcftools +fixref "${TMPDIR}/sgdp.chr${CHR}.hg38.norm.bcf" \
    --threads "$THREADS" -Ob -o "${TMPDIR}/sgdp.chr${CHR}.hg38.fixref.bcf" \
    -- -f "$HG38_FA" -m flip -d
# Keep only this chr (liftover can scatter a handful of variants to other
# contigs); then sort.
bcftools view "${TMPDIR}/sgdp.chr${CHR}.hg38.fixref.bcf" -r "chr${CHR}" -Ou \
  | bcftools sort -m 8G -T "${TMPDIR}/sort" -Ob -o "$PREPPED"
bcftools index "$PREPPED"

# Clean large intermediates (keep the reject log for provenance).
rm -f "$FILT_VCF" "${FILT_VCF}.tbi" "$LIFTED_VCF" \
      "${TMPDIR}/sgdp.chr${CHR}.hg38.norm.bcf"* \
      "${TMPDIR}/sgdp.chr${CHR}.hg38.fixref.bcf"*
[[ "$LIFT_IN" != "$FILT_VCF" ]] && rm -f "$LIFT_IN" "${LIFT_IN}.tbi"

# ---------------------------------------------------------------------------
# 3. Report contents.
# ---------------------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "$PREPPED" | wc -l)
N_VARIANTS=$(bcftools index -n "$PREPPED")
N_REJECT=$(zcat "${OUTDIR}/sgdp_prepped.chr${CHR}.rejected.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
# Report the reject RATE, not just the count: a liftover that silently drops
# 30% of sites still "succeeds" and the raw count alone does not show it.
REJ_PCT=$(awk -v r="$N_REJECT" -v k="$N_VARIANTS" \
    'BEGIN{ t=r+k; printf "%.2f", (t>0 ? 100*r/t : 0) }')
echo "[$(date +%T)] [chr${CHR}] DONE. ${N_SAMPLES} samples, ${N_VARIANTS} variants lifted, ${N_REJECT} rejected (${REJ_PCT}% of attempted) -> $PREPPED"
if awk -v p="$REJ_PCT" 'BEGIN{exit !(p>5)}'; then
    echo "[$(date +%T)] [chr${CHR}] WARNING: reject rate ${REJ_PCT}% exceeds 5%." \
         "Check the chain/reference pair and the rejected VCF:" \
         "${OUTDIR}/sgdp_prepped.chr${CHR}.rejected.vcf.gz"
fi

# Emit a per-chr sample list for provenance / debugging.
bcftools query -l "$PREPPED" > "${PREPPED%.bcf}.samples.txt"
