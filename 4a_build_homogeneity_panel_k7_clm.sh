#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=4a_homog_k7
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=01:00:00
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_k7_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_k7_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- K=7 homogeneity panel with SGDP OCE/MEN anchor DISCOVERY.
# ============================================================================
# K=7 companion to 4a_build_homogeneity_panel_clm.sh. The K=5 script uses
# PRE-BUILT anchor lists (eur/afr/amr/eas/sas_rfmix.txt) and cannot resolve
# Oceanian or Middle Eastern, because those super-pops have no 1KG sub-pop
# label and no gnomAD inference prob to pick anchors from.
#
# This script implements the two-stage anchor-discovery method used to build
# the original anchor sets, applied to the SGDP additions:
#
#   STAGE A (unsupervised, anchor DISCOVERY):
#     Run UNSUPERVISED ADMIXTURE at K=7 on the merged panel (which now
#     includes the SGDP candidates from 1c/1b). Map each of the 7 clusters to
#     a super-pop by exemplar majority vote (the 5 existing anchor sets + the
#     SGDP OCE/MEN candidate pools). For OCE and MEN, take the SGDP candidates
#     whose Q in their own cluster is >= HOMOG_THRESHOLD, SORTED by that Q
#     descending -> reference_ids/oce_rfmix.txt, men_rfmix.txt.
#     => "order by the unsupervised sorted order and identify the anchors
#         from the unsupervised" (Christina's spec).
#
#   STAGE B (supervised, anchor-INFORMED):
#     Rerun ADMIXTURE --supervised at K=7 with the full 7-way .pop prior
#     (EUR/AFR/AMR/EAS/SAS existing + OCE/MEN just discovered). Supervised Q
#     is the truth for the homogeneity filter.
#
#   STAGE C (homogeneity filter):
#     included = 1 iff supervised max_Q >= HOMOG_THRESHOLD AND the supervised
#     argmax super-pop matches the sample's assigned super-pop. MXB forced in.
#
# The SGDP fetcher's 277 outputs remain strictly a CANDIDATE POOL; this script
# is what decides which clear Q>=0.95 and enter the published panel.
#
# Outputs in $HOMOG_DIR (k7 subdirs, so the K=5 outputs are never clobbered):
#   admixture_k7/lai_ref.{bed,bim,fam,7.Q,7.P,pop}
#   admixture_k7/unsup/lai_ref.7.Q
#   reference_ids/oce_rfmix.txt, men_rfmix.txt         (DISCOVERED anchors)
#   lai_ref_panel_samples.k7.tsv
#   homog_7pop_samples.txt
#   homog_7pop/merged_chr${CHR}.homog_7pop.bcf{,.csi}
# ============================================================================

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
GNOMAD_META="${GNOMAD_META:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
MXB_POPINFO="${MXB_POPINFO:-${REPO_DIR}/reference_ids/MXB50genomes_popinfo.tsv}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"

# K=7 and strict 0.95 homogeneity (publication). Override via env if needed.
K="${K:-7}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.95}"

# SGDP candidate pools (from scripts/fetch-sgdp-anchors.py). Used as exemplars
# for cluster labeling + as the pool OCE/MEN anchors are discovered from.
SGDP_MANIFEST="${SGDP_MANIFEST:-${REFS}/sgdp_manifest.tsv}"
SGDP_OCE_CAND="${SGDP_OCE_CAND:-${REFS}/sgdp_oceania_anchors.txt}"
SGDP_MEN_CAND="${SGDP_MEN_CAND:-${REFS}/sgdp_men_anchors.txt}"

# Discovered anchor outputs (written by STAGE A).
OCE_ANCHORS="${OCE_ANCHORS:-${REFS}/oce_rfmix.txt}"
MEN_ANCHORS="${MEN_ANCHORS:-${REFS}/men_rfmix.txt}"

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-16}
ADM_DIR="${HOMOG_DIR}/admixture_k7"
UNSUP_DIR="${ADM_DIR}/unsup"
HOMOG_SUB="${HOMOG_DIR}/homog_7pop"
mkdir -p "$HOMOG_DIR" "$ADM_DIR" "$UNSUP_DIR" "$HOMOG_SUB"

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
# Verify required inputs.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] K=${K}, HOMOG_THRESHOLD=${HOMOG_THRESHOLD}"
echo "[$(date +%T)] verifying inputs"
for CHR in {1..22}; do
    F="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    [[ -s "$F" && -s "${F}.csi" ]] || { echo "ERROR: missing $F (phasing not done)"; exit 1; }
done
for F in eur_rfmix.txt afr_rfmix.txt amr_rfmix.txt eas_rfmix.txt sas_rfmix.txt; do
    [[ -s "${REFS}/${F}" ]] || { echo "ERROR: missing ${REFS}/${F}"; exit 1; }
done
[[ -s "$SAMPLE_GROUPS" ]] || { echo "ERROR: missing $SAMPLE_GROUPS"; exit 1; }
[[ -s "$GNOMAD_META"   ]] || { echo "ERROR: missing $GNOMAD_META"; exit 1; }
[[ -s "$MXB_POPINFO"   ]] || { echo "ERROR: missing $MXB_POPINFO"; exit 1; }
[[ -s "$SGDP_MANIFEST" ]] || { echo "ERROR: missing $SGDP_MANIFEST (run fetch-sgdp-anchors.py)"; exit 1; }
[[ -s "$SGDP_OCE_CAND" ]] || { echo "ERROR: missing $SGDP_OCE_CAND"; exit 1; }
[[ -s "$SGDP_MEN_CAND" ]] || { echo "ERROR: missing $SGDP_MEN_CAND"; exit 1; }

# ---------------------------------------------------------------------------
# 1. Concat all 22 chrs to a single BCF (k7 workspace).
# ---------------------------------------------------------------------------
CONCAT="${ADM_DIR}/all_chr.phased.bcf"
if [[ ! -s "$CONCAT" ]]; then
    echo "[$(date +%T)] concat 22 chrs"
    : > "${ADM_DIR}/concat_list.txt"
    for CHR in {1..22}; do
        echo "${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf" \
            >> "${ADM_DIR}/concat_list.txt"
    done
    bcftools concat --threads "$THREADS" -Ob -o "$CONCAT" -f "${ADM_DIR}/concat_list.txt"
    bcftools index "$CONCAT"
fi

# ---------------------------------------------------------------------------
# 2. plink2 import + LD-prune.
# ---------------------------------------------------------------------------
PRUNED="${ADM_DIR}/lai_ref"
if [[ ! -s "${PRUNED}.bed" ]]; then
    echo "[$(date +%T)] plink2 LD-prune (500 50 0.1 maf 0.01)"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' --maf 0.01 \
           --indep-pairwise 500 50 0.1 --out "${ADM_DIR}/lai_ref.prune"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' \
           --extract "${ADM_DIR}/lai_ref.prune.prune.in" \
           --max-alleles 2 --snps-only just-acgt \
           --make-bed --out "$PRUNED"
    echo "[$(date +%T)] pruned: $(wc -l < "${PRUNED}.bim") SNPs, $(wc -l < "${PRUNED}.fam") samples"
fi

# ---------------------------------------------------------------------------
# 3. STAGE A -- UNSUPERVISED ADMIXTURE K (runs FIRST; discovery depends on it).
# ---------------------------------------------------------------------------
if [[ ! -s "${UNSUP_DIR}/lai_ref.${K}.Q" ]]; then
    echo "[$(date +%T)] admixture (unsupervised) K=${K} (slow)"
    for ext in bed bim fam; do ln -sf "${PRUNED}.${ext}" "${UNSUP_DIR}/lai_ref.${ext}"; done
    cd "$UNSUP_DIR"
    admixture -j${THREADS} lai_ref.bed "${K}"
    cd - >/dev/null
fi
[[ -s "${UNSUP_DIR}/lai_ref.${K}.Q" ]] || { echo "ERROR: unsupervised ADMIXTURE failed"; exit 1; }

# ---------------------------------------------------------------------------
# 4. STAGE A -- discover OCE/MEN anchors from the unsupervised Q.
#    Map 7 clusters -> super-pop by exemplar majority vote (5 existing anchor
#    sets + SGDP OCE/MEN candidate pools), then take the SGDP OCE/MEN
#    candidates with cluster-Q >= threshold, sorted descending.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] discover OCE/MEN anchors from unsupervised K=${K}"
python3 - "${PRUNED}.fam" "${UNSUP_DIR}/lai_ref.${K}.Q" "$K" "$HOMOG_THRESHOLD" \
    "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" "${REFS}/amr_rfmix.txt" \
    "${REFS}/eas_rfmix.txt" "${REFS}/sas_rfmix.txt" \
    "$SGDP_OCE_CAND" "$SGDP_MEN_CAND" "$MXB_POPINFO" \
    "$OCE_ANCHORS" "$MEN_ANCHORS" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam_f, qu_f, K, thresh, eur_f, afr_f, amr_f, eas_f, sas_f,
 oce_cand_f, men_cand_f, mxb_f, oce_out, men_out) = sys.argv[1:15]
K = int(K); thresh = float(thresh)

def load(p): return set(open(p).read().split())
mxb = set()
with open(mxb_f) as f:
    next(f)
    for line in f:
        if line.strip(): mxb.add(line.split("\t")[0])

# Exemplar -> super-pop. Existing 5 come from the anchor lists; OCE/MEN
# exemplars are the SGDP candidate pools (used ONLY to label the clusters).
exemplars = {}
for sp, ids in [("EUR", load(eur_f)), ("AFR", load(afr_f)),
                ("AMR", load(amr_f) | mxb), ("EAS", load(eas_f)),
                ("SAS", load(sas_f)), ("OCE", load(oce_cand_f)),
                ("MEN", load(men_cand_f))]:
    for s in ids: exemplars[s] = sp

samples = [l.split()[1] for l in open(fam_f)]
Qu = [list(map(float, l.split())) for l in open(qu_f)]
assert len(samples) == len(Qu), f"{len(samples)} fam != {len(Qu)} Q rows"

# Majority vote: each cluster -> the super-pop whose exemplars peak there.
# Greedy by descending vote mass so no two super-pops claim one cluster.
votes = defaultdict(Counter)  # super-pop -> Counter(cluster)
for s, q in zip(samples, Qu):
    sp = exemplars.get(s)
    if sp is None: continue
    votes[sp][q.index(max(q))] += 1
cluster_of_sp = {}; used = set()
for sp, ctr in sorted(votes.items(), key=lambda kv: -sum(kv[1].values())):
    for cl, _ in ctr.most_common():
        if cl not in used:
            cluster_of_sp[sp] = cl; used.add(cl); break
print(f"[discover] unsup cluster map: {cluster_of_sp}", file=sys.stderr)

def discover(cand_f, sp, out):
    if sp not in cluster_of_sp:
        print(f"[discover] WARNING: {sp} did not claim a cluster at K={K}; "
              f"writing empty {out}", file=sys.stderr)
        open(out, "w").close(); return 0
    cl = cluster_of_sp[sp]
    cand = load(cand_f)
    picked = []
    for s, q in zip(samples, Qu):
        if s in cand and q[cl] >= thresh:
            picked.append((s, q[cl]))
    picked.sort(key=lambda x: -x[1])  # unsupervised sorted order, descending
    with open(out, "w") as g:
        for s, _ in picked: g.write(s + "\n")
    n_cand_in_panel = sum(1 for s in samples if s in cand)
    print(f"[discover] {sp}: cluster {cl}, {len(picked)}/{n_cand_in_panel} "
          f"candidates-in-panel pass Q>={thresh} -> {out}", file=sys.stderr)
    if picked:
        print(f"           top: {picked[0][0]} Q={picked[0][1]:.4f} | "
              f"min kept Q={picked[-1][1]:.4f}", file=sys.stderr)
    return len(picked)

n_oce = discover(oce_cand_f, "OCE", oce_out)
n_men = discover(men_cand_f, "MEN", men_out)
if n_oce == 0 or n_men == 0:
    print(f"[discover] NOTE: OCE={n_oce} MEN={n_men} anchors discovered. "
          f"Zero means that super-pop didn't form a clean cluster at K={K} "
          f"or no candidate cleared Q>={thresh}. Supervised stage will still "
          f"run but that label will be weak/absent.", file=sys.stderr)
PYEOF

N_OCE=$(awk 'NF' "$OCE_ANCHORS" 2>/dev/null | wc -l)
N_MEN=$(awk 'NF' "$MEN_ANCHORS" 2>/dev/null | wc -l)
echo "[$(date +%T)] discovered anchors: OCE=${N_OCE}, MEN=${N_MEN}"

# ---------------------------------------------------------------------------
# 5. STAGE B -- build supervised .pop file (7 labels) in .fam row order.
# ---------------------------------------------------------------------------
POP_FILE="${PRUNED}.pop"
echo "[$(date +%T)] build supervised .pop file (K=${K})"
python3 - "${PRUNED}.fam" "$POP_FILE" \
    "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" "${REFS}/amr_rfmix.txt" \
    "${REFS}/eas_rfmix.txt" "${REFS}/sas_rfmix.txt" \
    "$OCE_ANCHORS" "$MEN_ANCHORS" "$MXB_POPINFO" <<'PYEOF'
import sys
(fam, outp, eur_f, afr_f, amr_f, eas_f, sas_f, oce_f, men_f, mxb_pop_f) = sys.argv[1:11]
def load(p):
    try: return set(open(p).read().split())
    except FileNotFoundError: return set()
mxb = set()
with open(mxb_pop_f) as f:
    next(f)
    for line in f:
        if line.strip(): mxb.add(line.split("\t")[0])
# Order matters: ADMIXTURE Q columns follow first-appearance of labels in .pop.
anchors = [("EUR", load(eur_f)), ("AFR", load(afr_f)), ("AMR", load(amr_f) | mxb),
           ("EAS", load(eas_f)), ("SAS", load(sas_f)),
           ("OCE", load(oce_f)), ("MEN", load(men_f))]
counts = {sp: 0 for sp, _ in anchors}; counts["-"] = 0
with open(fam) as f, open(outp, "w") as g:
    for line in f:
        sample = line.split()[1]
        label = "-"
        for sp, ids in anchors:
            if sample in ids: label = sp; break
        g.write(label + "\n"); counts[label] = counts.get(label, 0) + 1
for sp, _ in anchors: print(f"  {sp}: {counts.get(sp,0)}", file=sys.stderr)
print(f"  -: {counts.get('-',0)}", file=sys.stderr)
PYEOF

# ---------------------------------------------------------------------------
# 6. STAGE B -- supervised ADMIXTURE K.
# ---------------------------------------------------------------------------
if [[ ! -s "${ADM_DIR}/lai_ref.${K}.Q" ]]; then
    echo "[$(date +%T)] admixture --supervised K=${K} (slow)"
    cd "$ADM_DIR"
    admixture --supervised -j${THREADS} lai_ref.bed "${K}"
    cd - >/dev/null
fi
[[ -s "${ADM_DIR}/lai_ref.${K}.Q" ]] || { echo "ERROR: supervised ADMIXTURE failed"; exit 1; }

# ---------------------------------------------------------------------------
# 7. STAGE C -- parse both .Q -> lai_ref_panel_samples.k7.tsv + homogeneity.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] write lai_ref_panel_samples.k7.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.k7.tsv"
python3 - "${ADM_DIR}/lai_ref.fam" "${ADM_DIR}/lai_ref.${K}.Q" \
    "${UNSUP_DIR}/lai_ref.${K}.Q" "${ADM_DIR}/lai_ref.pop" \
    "$SAMPLE_GROUPS" "$GNOMAD_META" "$MXB_POPINFO" "$SGDP_MANIFEST" \
    "${REFS}/sas_rfmix.txt" "$HOMOG_THRESHOLD" "$K" "$SAMPLES_TSV" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam, qs_f, qu_f, pop_f, groups_f, gnomad_f, mxb_pop_f, sgdp_manifest_f,
 sas_f, thresh, K, out) = sys.argv[1:13]
thresh = float(thresh); K = int(K)

sg = {}
with open(groups_f) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2: sg[parts[0]] = parts[1]
sas_ids = set(open(sas_f).read().split())

gnomad_pop = {}
with open(gnomad_f) as f:
    h = next(f).rstrip("\n").split("\t")
    sc = h.index("project_meta.sample_id"); pc = h.index("hgdp_tgp_meta.Population")
    for line in f:
        ff = line.rstrip("\n").split("\t")
        if len(ff) > max(sc, pc): gnomad_pop[ff[sc]] = ff[pc]

mxb_pop = {}
with open(mxb_pop_f) as f:
    h = next(f).rstrip("\n").split("\t")
    sc = h.index("Sample_ID"); pc = h.index("inferred_genetic_cluster")
    for line in f:
        ff = line.rstrip("\n").split("\t")
        if len(ff) > max(sc, pc): mxb_pop[ff[sc]] = ff[pc]

# SGDP: manifest gives super_pop + population for each SGDP candidate.
sgdp_sp = {}; sgdp_pop = {}
with open(sgdp_manifest_f) as f:
    h = next(f).rstrip("\n").split("\t")
    si = h.index("sample"); spi = h.index("super_pop"); pi = h.index("population")
    for line in f:
        ff = line.rstrip("\n").split("\t")
        if len(ff) > max(si, spi, pi):
            sgdp_sp[ff[si]] = ff[spi]; sgdp_pop[ff[si]] = ff[pi]

samples = [l.split()[1] for l in open(fam)]
Qs = [list(map(float, l.split())) for l in open(qs_f)]
Qu = [list(map(float, l.split())) for l in open(qu_f)]
anchors = [l.strip() for l in open(pop_f)]
assert len(samples) == len(Qs) == len(Qu) == len(anchors)

def majority_vote_map(anchors, Q, K):
    votes = defaultdict(Counter)
    for a, q in zip(anchors, Q):
        if a == "-": continue
        votes[a][q.index(max(q))] += 1
    m = {}; used = set()
    for sp, ctr in sorted(votes.items(), key=lambda kv: -sum(kv[1].values())):
        for cl, _ in ctr.most_common():
            if cl not in used: m[cl] = sp; used.add(cl); break
    for cl in range(K): m.setdefault(cl, f"UNKc{cl+1}")
    return m

sup_cluster_to_sp   = majority_vote_map(anchors, Qs, K)
unsup_cluster_to_sp = majority_vote_map(anchors, Qu, K)
sup_sp_to_cluster   = {v: k for k, v in sup_cluster_to_sp.items()}
unsup_sp_to_cluster = {v: k for k, v in unsup_cluster_to_sp.items()}
print(f"sup_cluster_pop_map   = {sup_cluster_to_sp}",   file=sys.stderr)
print(f"unsup_cluster_pop_map = {unsup_cluster_to_sp}", file=sys.stderr)

POP_ORDER = ["AFR", "AMR", "EAS", "EUR", "MEN", "OCE", "SAS"]

def cohort_of(s):
    if s in sgdp_sp: return "SGDP"
    if s.startswith("MXB"): return "MXB"
    if s.startswith(("HGDP", "LP6005", "SS")): return "HGDP"
    if s.startswith(("HG", "NA")): return "1KG"
    return "OTHER"

def super_pop_of(s):
    if s in mxb_pop: return "AMR"
    if s in sgdp_sp: return sgdp_sp[s]          # OCE/MEN/AFR/EUR/EAS/SAS/AMR
    raw = sg.get(s, "NA")
    if raw == "CSA": return "SAS" if s in sas_ids else "NA"
    if raw in ("AFR", "AMR", "EUR", "EAS"): return raw
    return "NA"

def pop_of(s):
    if s in mxb_pop: return mxb_pop[s]
    if s in sgdp_pop: return sgdp_pop[s]
    return gnomad_pop.get(s, "NA")

with open(out, "w") as g:
    hdr = ["sample", "cohort", "pop", "super_pop", "assigned"] \
        + [f"Q_sup_{c}" for c in POP_ORDER] + ["max_Q_sup", "argmax_pop_sup"] \
        + [f"Q_unsup_{c}" for c in POP_ORDER] + ["max_Q_unsup", "argmax_pop_unsup",
                                                 "sup_unsup_agree", "included"]
    g.write("\t".join(hdr) + "\n")
    n_incl = 0
    for s, qs, qu, a in zip(samples, Qs, Qu, anchors):
        co = cohort_of(s); sp = super_pop_of(s); p = pop_of(s)
        max_qs = max(qs); ap_sup   = sup_cluster_to_sp[qs.index(max_qs)]
        max_qu = max(qu); ap_unsup = unsup_cluster_to_sp[qu.index(max_qu)]
        agree = 1 if ap_sup == ap_unsup else 0
        qs_ordered = [qs[sup_sp_to_cluster[c]]   if c in sup_sp_to_cluster   else 0.0 for c in POP_ORDER]
        qu_ordered = [qu[unsup_sp_to_cluster[c]] if c in unsup_sp_to_cluster else 0.0 for c in POP_ORDER]
        if co == "MXB":
            included = 1
        elif sp == "NA":
            included = 0
        else:
            included = 1 if (max_qs >= thresh and ap_sup == sp) else 0
        n_incl += included
        g.write("\t".join(
            [s, co, p, sp, a]
            + [f"{x:.4f}" for x in qs_ordered] + [f"{max_qs:.4f}", ap_sup]
            + [f"{x:.4f}" for x in qu_ordered] + [f"{max_qu:.4f}", ap_unsup,
                                                  str(agree), str(included)]
        ) + "\n")
    print(f"included {n_incl}/{len(samples)}", file=sys.stderr)
PYEOF

N_INCL=$(awk -F'\t' 'NR>1 && $NF==1' "$SAMPLES_TSV" | wc -l)
N_TOT=$(awk 'NR>1' "$SAMPLES_TSV" | wc -l)
echo "[$(date +%T)] included: ${N_INCL}/${N_TOT} samples"
echo "  per super_pop (included only):"
awk -F'\t' 'NR>1 && $NF==1 {print "    " $4}' "$SAMPLES_TSV" | sort | uniq -c

awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_7pop_samples.txt"

# ---------------------------------------------------------------------------
# 8. Subset each per-chr phased BCF to the homogeneous samples (K=7 set).
# ---------------------------------------------------------------------------
echo "[$(date +%T)] subset per-chr BCFs to homog samples"
HOMOG_KEEP="${HOMOG_DIR}/homog_7pop_samples.txt"
for CHR in {1..22}; do
    SRC="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    DST="${HOMOG_SUB}/merged_chr${CHR}.homog_7pop.bcf"
    if [[ -s "$DST" && -s "${DST}.csi" ]]; then continue; fi
    bcftools view "$SRC" -S "$HOMOG_KEEP" --force-samples --threads "$THREADS" -Ob -o "$DST"
    bcftools index "$DST"
done

echo "[$(date +%T)] DONE."
echo "  samples TSV      : $SAMPLES_TSV"
echo "  discovered OCE   : $OCE_ANCHORS ($N_OCE)"
echo "  discovered MEN   : $MEN_ANCHORS ($N_MEN)"
echo "  homog BCFs       : ${HOMOG_SUB}/merged_chr*.homog_7pop.bcf"
