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
# HGDP+1KG+MXB BCFs. Uses Jessica's existing per-super-pop anchor lists as
# the supervised-ADMIXTURE prior rather than re-deriving from sub-pop names.
#
# Anchor sources (reference_ids/):
#   EUR = eur_rfmix.txt   (77 samples, 1KG IBS)
#   AFR = afr_rfmix.txt   (77 samples, 1KG YRI)
#   AMR = amr_rfmix.txt   (31 samples, HGDP-NAT) + all 50 MXB samples
#         The 50 MXB samples are ~98% Amerindigenous by curation
#         (selected by MXB scientists from various Mexican geographic
#         locations). Anchoring them alongside HGDP-NAT roughly triples
#         the AMR exemplar count and stabilizes the AMR component
#         centroid for the rest of the inference.
#   EAS = eas_rfmix.txt   (build with build-eas-sas-anchors.sh)
#   SAS = sas_rfmix.txt   (build with build-eas-sas-anchors.sh)
#   Every sample NOT in one of these sets (1KG-AMR PEL/MXL/CLM/PUR,
#   ASW/ACB, OCE, MEN, etc.) gets supervised label "-" so the anchors
#   define the K=5 components without admixed samples dragging centroids.
#
# Super-pop assignment for the samples TSV's included-flag check:
#   - sample_groups.tsv (from make_sample_groups.sh) gives super_pop. CSA gets
#     refined: 1KG-SAS sub-pops (in sas_rfmix.txt) -> SAS; other CSA -> drop.
#   - MXB samples -> AMR (forced included=1 regardless of Q)
#   - OCE / MEN / non-SAS CSA -> excluded (included=0, super_pop=NA)
#
# Pipeline:
#   1. Concat all 22 phased chrs -> one BCF
#   2. plink2 LD-prune (--indep-pairwise 50 5 0.2) -> ~500-700k SNPs
#   3. Build supervised .pop file from anchor lists (in .fam row order)
#   4. admixture --supervised -j16 K=5
#   5. Parse .Q -> lai_ref_panel_samples.tsv (cohort, pop, super_pop, Q_*, included)
#   6. Subset each per-chr phased BCF to included samples -> homog_5pop/
#
# Outputs in $HOMOG_DIR:
#   admixture/lai_ref.{bed,bim,fam,5.Q,5.P,pop}
#   lai_ref_panel_samples.tsv
#   homog_5pop_samples.txt
#   homog_5pop/merged_chr${CHR}.homog_5pop.bcf{,.csi}
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$PROJECT_ROOT}"   # adjust if your repo lives elsewhere
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
GNOMAD_META="${GNOMAD_META:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
MXB_POPINFO="${MXB_POPINFO:-${REPO_DIR}/reference_ids/MXB50genomes_popinfo.tsv}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.95}"

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$HOMOG_DIR" "$HOMOG_DIR/admixture" "$HOMOG_DIR/homog_5pop"

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Verify required inputs.
echo "[$(date +%T)] verifying inputs"
for CHR in {1..22}; do
    F="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    [[ -s "$F" && -s "${F}.csi" ]] || { echo "ERROR: missing $F (phasing not done)"; exit 1; }
done
for F in eur_rfmix.txt afr_rfmix.txt amr_rfmix.txt eas_rfmix.txt sas_rfmix.txt; do
    [[ -s "${REFS}/${F}" ]] || { echo "ERROR: missing ${REFS}/${F} (run build-eas-sas-anchors.sh first)"; exit 1; }
done
[[ -s "$SAMPLE_GROUPS" ]] || { echo "ERROR: missing $SAMPLE_GROUPS"; exit 1; }
[[ -s "$GNOMAD_META"   ]] || { echo "ERROR: missing $GNOMAD_META"; exit 1; }
[[ -s "$MXB_POPINFO"   ]] || { echo "ERROR: missing $MXB_POPINFO"; exit 1; }

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
# 3. Build supervised .pop file in .fam row order from anchor lists.
# ---------------------------------------------------------------------------
POP_FILE="${PRUNED}.pop"
echo "[$(date +%T)] build supervised .pop file"
python3 - "${PRUNED}.fam" "$POP_FILE" \
    "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" "${REFS}/amr_rfmix.txt" \
    "${REFS}/eas_rfmix.txt" "${REFS}/sas_rfmix.txt" "$MXB_POPINFO" <<'PYEOF'
import sys
fam, outp, eur_f, afr_f, amr_f, eas_f, sas_f, mxb_pop_f = sys.argv[1:9]
def load(p): return set(open(p).read().split())
# 50 MXB samples (~98% Amerindigenous, curated) join the AMR anchor set.
mxb = set()
with open(mxb_pop_f) as f:
    next(f)  # header
    for line in f:
        mxb.add(line.split("\t")[0])
anchors = {"EUR": load(eur_f), "AFR": load(afr_f),
           "AMR": load(amr_f) | mxb,
           "EAS": load(eas_f), "SAS": load(sas_f)}
with open(fam) as f, open(outp, "w") as g:
    counts = {k: 0 for k in anchors}
    counts["-"] = 0
    for line in f:
        sample = line.split()[1]  # IID
        label = "-"
        for sp, ids in anchors.items():
            if sample in ids:
                label = sp
                break
        g.write(label + "\n")
        counts[label] = counts.get(label, 0) + 1
    for k in ("EUR", "AFR", "AMR", "EAS", "SAS", "-"):
        print(f"  {k}: {counts.get(k, 0)}", file=sys.stderr)
PYEOF

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
# 5. Parse .Q -> lai_ref_panel_samples.tsv.
#    Component order in .Q follows alpha order of supervised labels: AFR AMR EAS EUR SAS.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] write lai_ref_panel_samples.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.tsv"
python3 - "${ADM_DIR}/lai_ref.fam" "${ADM_DIR}/lai_ref.5.Q" \
    "$SAMPLE_GROUPS" "$GNOMAD_META" "$MXB_POPINFO" \
    "${REFS}/sas_rfmix.txt" "$HOMOG_THRESHOLD" "$SAMPLES_TSV" <<'PYEOF'
import sys
fam, qfile, groups_f, gnomad_f, mxb_pop_f, sas_f, thresh, out = sys.argv[1:9]
thresh = float(thresh)

# sample -> super_pop from make_sample_groups.sh (AFR/AMR/EUR/EAS/CSA/OCE/MEN).
sg = {}
with open(groups_f) as f:
    for line in f:
        s, sp = line.rstrip("\n").split("\t")[:2]
        sg[s] = sp

# sas_rfmix: 1KG-SAS sub-pops (refine CSA -> SAS for these samples).
sas_ids = set(open(sas_f).read().split())

# gnomAD meta -> sample -> sub-pop (hgdp_tgp_meta.Population).
gnomad_pop = {}
with open(gnomad_f) as f:
    header = next(f).rstrip("\n").split("\t")
    sc = header.index("project_meta.sample_id")
    pc = header.index("hgdp_tgp_meta.Population")
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) > max(sc, pc):
            gnomad_pop[fields[sc]] = fields[pc]

# MXB sub-pop from MXB50genomes_popinfo.tsv (Sample_ID -> inferred_genetic_cluster).
mxb_pop = {}
with open(mxb_pop_f) as f:
    h = next(f).rstrip("\n").split("\t")
    sc = h.index("Sample_ID")
    pc = h.index("inferred_genetic_cluster")
    for line in f:
        ff = line.rstrip("\n").split("\t")
        mxb_pop[ff[sc]] = ff[pc]

COMPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
samples = [l.split()[1] for l in open(fam)]
Q = [list(map(float, l.split())) for l in open(qfile)]
assert len(samples) == len(Q), f"sample/Q row mismatch: {len(samples)} vs {len(Q)}"

def cohort_of(s):
    if s.startswith("MXB"):                          return "MXB"
    if s.startswith(("HGDP", "LP6005", "SS")):       return "HGDP"
    if s.startswith(("HG", "NA")):                   return "1KG"
    return "OTHER"

def super_pop_of(s):
    if s in mxb_pop:                  return "AMR"
    raw = sg.get(s, "NA")
    if raw == "CSA":                  return "SAS" if s in sas_ids else "NA"
    if raw in ("AFR","AMR","EUR","EAS"): return raw
    return "NA"

def pop_of(s):
    if s in mxb_pop:           return mxb_pop[s]
    return gnomad_pop.get(s, "NA")

with open(out, "w") as g:
    hdr = ["sample", "cohort", "pop", "super_pop"] + [f"Q_{c}" for c in COMPS] + \
          ["assigned", "max_Q", "included"]
    g.write("\t".join(hdr) + "\n")
    n_incl = 0
    for s, q in zip(samples, Q):
        co = cohort_of(s); sp = super_pop_of(s); p = pop_of(s)
        max_q = max(q); assigned = COMPS[q.index(max_q)]
        if co == "MXB":
            included = 1                # always include MXB
        elif sp == "NA":
            included = 0                # OCE / MEN / non-SAS CSA
        else:
            included = 1 if (max_q >= thresh and assigned == sp) else 0
        n_incl += included
        g.write("\t".join([s, co, p, sp] + [f"{x:.4f}" for x in q] +
                          [assigned, f"{max_q:.4f}", str(included)]) + "\n")
    print(f"included {n_incl}/{len(samples)} samples", file=sys.stderr)
PYEOF

N_INCL=$(awk -F'\t' 'NR>1 && $NF==1' "$SAMPLES_TSV" | wc -l)
N_TOT=$(awk 'NR>1' "$SAMPLES_TSV" | wc -l)
echo "[$(date +%T)] included: ${N_INCL}/${N_TOT} samples"
echo "  per super_pop (included only):"
awk -F'\t' 'NR>1 && $NF==1 {print "    " $4}' "$SAMPLES_TSV" | sort | uniq -c

# Just-IDs file for the per-chr subset step + downstream tools.
awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_5pop_samples.txt"

# ---------------------------------------------------------------------------
# 6. Subset each per-chr phased BCF to the homogeneous samples.
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
echo "  samples TSV : $SAMPLES_TSV"
echo "  homog BCFs  : $HOMOG_DIR/homog_5pop/merged_chr*.homog_5pop.bcf"
