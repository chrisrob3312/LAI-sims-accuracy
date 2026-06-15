#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=4a_homog
#SBATCH --partition=mhgcp
#SBATCH --nodelist=mhgcp-a01,mhgcp-m00,mhgcp-c00,mhgcp-c01,mhgcp-c02
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION
# ============================================================================
# Build the homogeneity-filtered 5-superpop reference panel from the phased
# HGDP+1KG+MXB BCFs. Workflow:
#   1. Concat all 22 per-chr phased BCFs (chr-prefixed, hg38).
#   2. plink2 import + LD-prune (--indep-pairwise 50 5 0.2) -> ~500-700k SNPs.
#   3. Build supervised .pop file using high-confidence anchor sub-pops:
#        EUR = CEU, GBR, IBS, TSI, FIN
#        AFR = YRI, ESN, MSL, GWD, LWK             (exclude ASW/ACB - admixed)
#        EAS = CHB, JPT, CHS, CDX, KHV
#        SAS = GIH, PJL, BEB, STU, ITU             (gnomAD lumps these into CSA)
#        AMR = Karitiana, Surui, Pima, Maya         (HGDP-NAT only, NOT 1KG-AMR)
#      Other samples (incl. ASW/ACB, all 1KG-AMR, MXB, OCE, MEN, CSA non-SAS):
#      labelled "-" (unsupervised). Anchors define K=5 components; everyone
#      else gets inferred ancestry.
#   4. admixture --supervised -j16 K=5
#   5. Parse .Q -> per-sample Q_<comp> proportions
#   6. Write lai_ref_panel_samples.tsv with cohort/pop/super_pop/Q/included
#      Included flag:
#        - MXB samples: 1 (forced in -- they are the study cohort)
#        - HGDP+1KG: 1 if max(Q)>=0.95 AND assigned component matches
#          self-reported super_pop label (EUR/AFR/AMR/EAS/SAS), else 0
#      OCE/MEN/CSA-non-SAS: excluded from panel (super_pop = NA, included=0)
#   7. Subset each per-chr rechr.bcf to the included samples -> the homog
#      5-superpop panel.
#
# Inputs:
#   $PROJECT_ROOT/01_merged_phased_panel/merged_chr*.shapeit5_phased.softunion_maf005.bcf
#   $GNOMAD_META (sample-level metadata from gnomAD HGDP+1KG release)
#   $MXB_SAMPLES (50 MXB sample IDs, one per line)
#
# Outputs in $HOMOG_DIR:
#   admixture/lai_ref.{bed,bim,fam,Q,P,pop}
#   lai_ref_panel_samples.tsv  (sample, cohort, pop, super_pop, Q_*, included)
#   homog_5pop_samples.txt     (just the included sample IDs)
#   homog_5pop/merged_chr${CHR}.homog_5pop.bcf{,.csi}  per chr
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
GNOMAD_META="${GNOMAD_META:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
MXB_SAMPLES="${MXB_SAMPLES:-${PROJECT_ROOT}/reference_ids/mxb_all_50.txt}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.95}"

# Anchor sub-pops per super-pop for supervised K=5.
EUR_ANCHORS="${EUR_ANCHORS:-CEU,GBR,IBS,TSI,FIN}"
AFR_ANCHORS="${AFR_ANCHORS:-YRI,ESN,MSL,GWD,LWK}"
EAS_ANCHORS="${EAS_ANCHORS:-CHB,JPT,CHS,CDX,KHV}"
SAS_ANCHORS="${SAS_ANCHORS:-GIH,PJL,BEB,STU,ITU}"
AMR_ANCHORS="${AMR_ANCHORS:-Karitiana,Surui,Pima,Maya}"

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$HOMOG_DIR" "$HOMOG_DIR/admixture" "$HOMOG_DIR/homog_5pop"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Verify all 22 phased BCFs exist.
echo "[$(date +%T)] verifying phased inputs"
for CHR in {1..22}; do
    F="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    [[ -s "$F" && -s "${F}.csi" ]] || { echo "ERROR: missing $F (phasing not complete)"; exit 1; }
done

# ---------------------------------------------------------------------------
# 1. Concat all 22 chrs to a single BCF.
# ---------------------------------------------------------------------------
CONCAT="${HOMOG_DIR}/admixture/all_chr.phased.bcf"
if [[ ! -s "$CONCAT" ]]; then
    echo "[$(date +%T)] concat 22 chrs"
    : > "${HOMOG_DIR}/admixture/concat_list.txt"
    for CHR in {1..22}; do
        echo "${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf" \
            >> "${HOMOG_DIR}/admixture/concat_list.txt"
    done
    bcftools concat --threads "$THREADS" -Ob -o "$CONCAT" \
        -f "${HOMOG_DIR}/admixture/concat_list.txt"
    bcftools index "$CONCAT"
fi

# ---------------------------------------------------------------------------
# 2. plink2 import + LD-prune.
# ---------------------------------------------------------------------------
PRUNED="${HOMOG_DIR}/admixture/lai_ref"
if [[ ! -s "${PRUNED}.bed" ]]; then
    echo "[$(date +%T)] plink2 LD-prune (50 5 0.2)"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' \
           --indep-pairwise 50 5 0.2 \
           --out "${HOMOG_DIR}/admixture/lai_ref.prune"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' \
           --extract "${HOMOG_DIR}/admixture/lai_ref.prune.prune.in" \
           --max-alleles 2 --snps-only just-acgt \
           --make-bed --out "$PRUNED"
    N_PRUNED=$(wc -l < "${PRUNED}.bim")
    N_SAMPLES=$(wc -l < "${PRUNED}.fam")
    echo "[$(date +%T)] pruned: ${N_PRUNED} SNPs, ${N_SAMPLES} samples"
fi

# ---------------------------------------------------------------------------
# 3. Build supervised .pop file from gnomAD meta + anchor sub-pops + MXB list.
#    Format: one label per sample (matches .fam row order), "-" for unsup.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] build supervised .pop file"
POP_FILE="${PRUNED}.pop"
META_LOOKUP="${HOMOG_DIR}/admixture/sample_meta.tsv"

[[ -s "$GNOMAD_META" ]] || { echo "ERROR: missing $GNOMAD_META"; exit 1; }
[[ -s "$MXB_SAMPLES" ]] || { echo "ERROR: missing $MXB_SAMPLES"; exit 1; }

# Header lookup
SAMPLE_COL=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="project_meta.sample_id") {print i; exit}}' "$GNOMAD_META")
POP_COL=$(   awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="hgdp_tgp_meta.Population") {print i; exit}}' "$GNOMAD_META")
REGION_COL=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="hgdp_tgp_meta.Genetic.region") {print i; exit}}' "$GNOMAD_META")
[[ -n "$SAMPLE_COL" && -n "$POP_COL" && -n "$REGION_COL" ]] || \
    { echo "ERROR: gnomAD meta missing required columns"; exit 1; }

# Build sample_id -> (cohort, pop, super_pop) lookup. cohort heuristic:
#  HGDP* / LP6005* / SS6004471 -> HGDP
#  HG* / NA*                   -> 1KG
#  in $MXB_SAMPLES             -> MXB
awk -F'\t' -v sc="$SAMPLE_COL" -v pc="$POP_COL" -v rc="$REGION_COL" \
    -v mxb_file="$MXB_SAMPLES" '
    BEGIN {
        while ((getline line < mxb_file) > 0) mxb[line] = 1
        close(mxb_file)
    }
    NR > 1 {
        s = $sc; p = $pc; r = $rc
        if (s in mxb) { cohort = "MXB"; super_pop = "AMR" }
        else if (s ~ /^(HGDP|LP6005|SS6004471|SS)/) cohort = "HGDP"
        else if (s ~ /^(HG|NA)/) cohort = "1KG"
        else cohort = "OTHER"
        # super_pop mapping (only for non-MXB; MXB handled above)
        if (cohort != "MXB") {
            if (r == "EUR") super_pop = "EUR"
            else if (r == "AFR") super_pop = "AFR"
            else if (r == "EAS") super_pop = "EAS"
            else if (r == "AMR") super_pop = "AMR"
            else if (r == "CSA") {
                # gnomAD lumps 1KG-SAS into CSA; restrict SAS to true 1KG-SAS subpops
                if (p == "GIH" || p == "PJL" || p == "BEB" || p == "STU" || p == "ITU")
                    super_pop = "SAS"
                else super_pop = "NA"
            }
            else super_pop = "NA"   # OCE, MEN, etc.
        }
        # MXB samples have no gnomAD meta row; handled below.
        print s "\t" cohort "\t" p "\t" super_pop
    }
    ' "$GNOMAD_META" > "$META_LOOKUP"

# Append MXB samples (not in gnomAD meta).
awk '{print $1 "\tMXB\tMXB\tAMR"}' "$MXB_SAMPLES" >> "$META_LOOKUP"

# Build the .pop file in .fam row order.
awk 'BEGIN {FS=OFS="\t"} NR==FNR {m[$1]=$2"\t"$3"\t"$4; next}
     {
        s=$1
        if (s in m) {
            split(m[s], a, "\t")
            cohort=a[1]; pop=a[2]; super_pop=a[3]
            # supervised label only if in anchor sub-pop list
            # anchors set via env vars (csv): EUR_ANCHORS, ..., AMR_ANCHORS
            label = "-"
            if (super_pop != "NA") {
                cmd = "echo \"" ENVIRON["EUR_ANCHORS"] "," ENVIRON["AFR_ANCHORS"] "," ENVIRON["EAS_ANCHORS"] "," ENVIRON["SAS_ANCHORS"] "," ENVIRON["AMR_ANCHORS"] "\" | tr , \"\\n\" | grep -xF \"" pop "\""
                # cheaper: build maps once below
            }
        }
        print s, cohort, pop, super_pop  # we will post-process to .pop in next step
    }' "$META_LOOKUP" "${PRUNED}.fam" \
    > "${HOMOG_DIR}/admixture/fam_with_meta.tsv"

# Build anchor set as a single awk-friendly check.
python3 - "$EUR_ANCHORS" "$AFR_ANCHORS" "$EAS_ANCHORS" "$SAS_ANCHORS" "$AMR_ANCHORS" \
    "${HOMOG_DIR}/admixture/fam_with_meta.tsv" "$POP_FILE" <<'PYEOF'
import sys
eur, afr, eas, sas, amr = [set(a.split(",")) for a in sys.argv[1:6]]
inp, outp = sys.argv[6], sys.argv[7]
anchor_map = {}
for p in eur: anchor_map[p] = "EUR"
for p in afr: anchor_map[p] = "AFR"
for p in eas: anchor_map[p] = "EAS"
for p in sas: anchor_map[p] = "SAS"
for p in amr: anchor_map[p] = "AMR"
with open(inp) as f, open(outp, "w") as g:
    for line in f:
        sample, cohort, pop, super_pop = line.rstrip("\n").split("\t")
        # supervised label: only if sub-pop is in anchor set AND aligns with super_pop
        label = anchor_map.get(pop, "-")
        if label != "-" and label != super_pop:
            # Pop name collision (shouldn't happen with default anchors); use "-"
            label = "-"
        g.write(label + "\n")
PYEOF
echo "[$(date +%T)] anchor counts:"
sort "$POP_FILE" | uniq -c | sort -rn | sed 's/^/    /'

# ---------------------------------------------------------------------------
# 4. Run supervised ADMIXTURE K=5.
# ---------------------------------------------------------------------------
ADM_DIR="${HOMOG_DIR}/admixture"
if [[ ! -s "${ADM_DIR}/lai_ref.5.Q" ]]; then
    echo "[$(date +%T)] admixture --supervised K=5 (this is the slow step)"
    cd "$ADM_DIR"
    admixture --supervised -j${THREADS} lai_ref.bed 5
    cd - >/dev/null
fi
[[ -s "${ADM_DIR}/lai_ref.5.Q" ]] || { echo "ERROR: ADMIXTURE failed"; exit 1; }

# ---------------------------------------------------------------------------
# 5. Parse .Q -> per-sample Q + decide included flag.
#    Component order in .Q rows depends on supervised label order. With our
#    .pop file, ADMIXTURE writes columns in alphabetical order of the
#    distinct supervised labels seen: AFR, AMR, EAS, EUR, SAS.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] write lai_ref_panel_samples.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.tsv"
python3 - "${ADM_DIR}/lai_ref.fam" "${ADM_DIR}/lai_ref.5.Q" "$META_LOOKUP" \
    "$MXB_SAMPLES" "$HOMOG_THRESHOLD" "$SAMPLES_TSV" <<'PYEOF'
import sys
fam, qfile, meta, mxb_file, thresh, out = sys.argv[1:7]
thresh = float(thresh)
mxb = set(open(mxb_file).read().split())
meta_map = {}
with open(meta) as f:
    for line in f:
        s, cohort, pop, super_pop = line.rstrip("\n").split("\t")
        meta_map[s] = (cohort, pop, super_pop)
COMPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]  # alpha order
samples = []
with open(fam) as f:
    for line in f:
        samples.append(line.split()[1])
Q = [list(map(float, l.split())) for l in open(qfile)]
assert len(samples) == len(Q), f"{len(samples)} samples vs {len(Q)} Q rows"
with open(out, "w") as g:
    g.write("sample\tcohort\tpop\tsuper_pop\t" + "\t".join("Q_"+c for c in COMPS) +
            "\tassigned\tmax_Q\tincluded\n")
    for s, q in zip(samples, Q):
        cohort, pop, super_pop = meta_map.get(s, ("OTHER", "NA", "NA"))
        max_q = max(q)
        assigned = COMPS[q.index(max_q)]
        if s in mxb:
            included = 1
        elif super_pop not in COMPS:
            included = 0
        else:
            included = 1 if (max_q >= thresh and assigned == super_pop) else 0
        g.write(f"{s}\t{cohort}\t{pop}\t{super_pop}\t" +
                "\t".join(f"{x:.4f}" for x in q) +
                f"\t{assigned}\t{max_q:.4f}\t{included}\n")
PYEOF
N_INCL=$(awk -F'\t' 'NR>1 && $NF==1' "$SAMPLES_TSV" | wc -l)
N_TOTAL=$(awk 'NR>1' "$SAMPLES_TSV" | wc -l)
echo "[$(date +%T)] included: ${N_INCL}/${N_TOTAL} samples in homog 5-pop panel"
awk -F'\t' 'NR>1 && $NF==1 {print $4}' "$SAMPLES_TSV" | sort | uniq -c | sed 's/^/    /'

# Just-IDs file for downstream filter steps.
awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_5pop_samples.txt"

# ---------------------------------------------------------------------------
# 6. Subset each per-chr rechr.bcf to the homogeneous samples.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] subset per-chr BCFs to homog samples"
HOMOG_KEEP="${HOMOG_DIR}/homog_5pop_samples.txt"
for CHR in {1..22}; do
    SRC="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    DST="${HOMOG_DIR}/homog_5pop/merged_chr${CHR}.homog_5pop.bcf"
    if [[ -s "$DST" && -s "${DST}.csi" ]]; then continue; fi
    bcftools view "$SRC" -S "$HOMOG_KEEP" --force-samples --threads "$THREADS" \
        -Ob -o "$DST"
    bcftools index "$DST"
done

echo "[$(date +%T)] DONE."
echo "  homog samples TSV: $SAMPLES_TSV"
echo "  homog 5-pop BCFs : $HOMOG_DIR/homog_5pop/merged_chr*.homog_5pop.bcf"
