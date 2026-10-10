#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=4a_homog_kN
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time-min=01:00:00
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_kN_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_kN_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- K=8 homogeneity panel: 7 named super-pops + empirical OTH.
# ============================================================================
# K-parametric companion to 4a_build_homogeneity_panel_clm.sh (K=5). Default
# K=8: the 7 named super-pops (AFR AMR EAS EUR MEN OCE SAS) plus OTH, the
# leftover cluster that captures samples not fitting an HGDP super-pop.
#
# The K=5 script uses PRE-BUILT anchor lists and cannot resolve Oceanian or
# Middle Eastern, because those have no 1KG sub-pop label and no gnomAD
# inference prob to pick anchors from. This script discovers them instead:
#
#   STAGE A (unsupervised, anchor DISCOVERY):
#     Run UNSUPERVISED ADMIXTURE at K on the merged panel (which now includes
#     the SGDP candidates from 1b/1c). Match clusters to super-pop NAMES by
#     MEAN Q-loading of labelled exemplars, one-to-one. The 7 named groups
#     claim 7 clusters; the LEFTOVER cluster becomes OTH -- defined by the
#     data, not by a hand-written population list.
#     Anchors = candidates sorted by their cluster-Q descending, kept if
#     >= ANCHOR_MIN_Q, capped at ANCHOR_CAP -> reference_ids/{oce,men,oth}_rfmix.txt.
#     OTH has no metadata candidate pool, so its anchors are drawn from the
#     whole panel; this also makes supervised K=8 well-posed, since ADMIXTURE
#     --supervised requires exactly K distinct labels in the .pop file.
#
#   STAGE B (supervised, anchor-INFORMED):
#     Rerun ADMIXTURE --supervised at K with the full K-way .pop prior
#     (5 existing *_rfmix.txt + OCE/MEN/OTH discovered in stage A).
#
#   STAGE C (homogeneity filter + refinement evidence):
#     included = 1 iff supervised max_Q >= HOMOG_THRESHOLD AND the supervised
#     argmax super-pop matches the sample's assigned super-pop. MXB forced in.
#     Emits population_cluster_meanQ.k${K}.tsv -- mean unsupervised Q per
#     POPULATION per cluster. That table is the evidence for the refinement
#     pass: where Mozabite actually loads (AFR vs MEN), whether the Beringian
#     set (Tlingit/Aleut/Eskimo_*) sits with EAS or AMR or OTH, whether the
#     Central Asians (Tajik/Kyrgyz/Uygur/Hazara) concentrate or split, and
#     where the SGDP populations with no HGDP counterpart (Australian, Maori,
#     Hawaiian, Dusun, Igorot) fall. Populations whose argmax disagrees with
#     their assigned label, or that fail to concentrate (best mean-Q < 0.60),
#     are printed as [refine] lines.
#
# The SGDP fetcher's 277 outputs remain strictly a CANDIDATE POOL; this script
# is what decides which clear the homogeneity threshold and enter the panel.
#
# Outputs in $HOMOG_DIR (all paths carry ${K}, so K=5/7/8 never collide):
#   admixture_k${K}/lai_ref.{bed,bim,fam,${K}.Q,${K}.P,pop}
#   admixture_k${K}/unsup/lai_ref.${K}.Q
#   reference_ids/{oce,men,oth}_rfmix.txt              (DISCOVERED anchors)
#   lai_ref_panel_samples.k${K}.tsv
#   population_cluster_meanQ.k${K}.tsv                 (refinement evidence)
#   homog_${K}pop_samples.txt
#   homog_${K}pop/merged_chr${CHR}.homog_${K}pop.bcf{,.csi}
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
#     includes the SGDP candidates from 1b/1c). Map each of the 7 clusters to
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
#   lai_ref_panel_samples.k${K}.tsv
#   homog_${K}pop_samples.txt
#   homog_7pop/merged_chr${CHR}.homog_${K}pop.bcf{,.csi}
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

# K=7. Two DISTINCT thresholds (do not conflate):
#   ANCHOR_MIN_Q   = floor to QUALIFY as a discovered anchor in stage A
#                    (unsupervised cluster-Q). Looser so small/clinal OCE/MEN
#                    yield enough training exemplars.
#   ANCHOR_CAP     = max anchors kept per super-pop (top-N by Q). Matches
#                    build-anchors-from-koenig.py (cap 50).
#   HOMOG_THRESHOLD= floor for FINAL inclusion in stage C (0.93, matching the
#                    committed K=5 panel -- kept deliberately below 0.95 so
#                    under-represented super-pops retain enough members)
#                    (supervised Q). This is the published homogeneity cut.
K="${K:-8}"
ANCHOR_MIN_Q="${ANCHOR_MIN_Q:-0.93}"
ANCHOR_CAP="${ANCHOR_CAP:-50}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.93}"

# SGDP candidate pools (from scripts/fetch-sgdp-anchors.py). Used as exemplars
# for cluster labeling + as the pool OCE/MEN anchors are discovered from.
# These are CANDIDATES, not anchors -- anchors are the top-N discovered here.
SGDP_MANIFEST="${SGDP_MANIFEST:-${REFS}/sgdp_manifest.tsv}"
SGDP_OCE_CAND="${SGDP_OCE_CAND:-${REFS}/sgdp_oceania_candidates.txt}"
SGDP_MEN_CAND="${SGDP_MEN_CAND:-${REFS}/sgdp_men_candidates.txt}"

# Discovered anchor outputs (written by STAGE A).
OCE_ANCHORS="${OCE_ANCHORS:-${REFS}/oce_rfmix.txt}"
MEN_ANCHORS="${MEN_ANCHORS:-${REFS}/men_rfmix.txt}"
# OTH = the empirical leftover cluster at K=8 (no metadata candidate pool).
OTH_ANCHORS="${OTH_ANCHORS:-${REFS}/oth_rfmix.txt}"

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-16}
ADM_DIR="${HOMOG_DIR}/admixture_k${K}"
UNSUP_DIR="${ADM_DIR}/unsup"
HOMOG_SUB="${HOMOG_DIR}/homog_${K}pop"
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
echo "[$(date +%T)] K=${K}  ANCHOR_MIN_Q=${ANCHOR_MIN_Q}  ANCHOR_CAP=${ANCHOR_CAP}  HOMOG_THRESHOLD=${HOMOG_THRESHOLD}"
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
#    Map 7 clusters -> super-pop by mean-Q (known labels); then for OCE/MEN
#    take their candidates with cluster-Q >= ANCHOR_MIN_Q, SORTED descending,
#    CAPPED at ANCHOR_CAP. This is the "top 25-50 >= 0.93-0.95" recipe used to
#    build the original anchor sets, applied to the SGDP additions.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] discover OCE/MEN anchors from unsupervised K=${K} (top ${ANCHOR_CAP} >= ${ANCHOR_MIN_Q})"
python3 - "${PRUNED}.fam" "${UNSUP_DIR}/lai_ref.${K}.Q" "$K" "$ANCHOR_MIN_Q" "$ANCHOR_CAP" \
    "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" "${REFS}/amr_rfmix.txt" \
    "${REFS}/eas_rfmix.txt" "${REFS}/sas_rfmix.txt" \
    "$SGDP_OCE_CAND" "$SGDP_MEN_CAND" "$MXB_POPINFO" \
    "$OCE_ANCHORS" "$MEN_ANCHORS" "$OTH_ANCHORS" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam_f, qu_f, K, thresh, cap, eur_f, afr_f, amr_f, eas_f, sas_f,
 oce_cand_f, men_cand_f, mxb_f, oce_out, men_out, oth_out) = sys.argv[1:17]
K = int(K); thresh = float(thresh); cap = int(cap)

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

# Match clusters to super-pop NAMES using the KNOWN labels (exemplars), by
# MEAN Q-loading with a one-to-one (greedy-by-max-cell) assignment. This
# replaces the old argmax-majority-vote, which processed super-pops by
# descending exemplar count and let the biggest (EUR/AMR) grab clusters first.
# That misassigns the poorly-separated OCE/MEN: a MEN sample can argmax onto
# the EUR cluster even when MEN forms its own weaker cluster, so a count-based
# plurality lands MEN on the wrong (taken) cluster. Mean-Q is robust to that
# argmax-splitting: MEN claims the cluster where MEN exemplars load highest on
# AVERAGE, even if individual MEN argmaxes scatter into adjacent EUR/SAS.
# The 7 NAMED super-pops have metadata-defined exemplar pools. At K=8 they
# claim 7 clusters and ONE cluster is left over -- that leftover is OTH: the
# empirical "does not fit an HGDP super-pop" group. It is defined by the data,
# not by a hand-written population list, which is the whole point of K=8.
SPS = ["EUR", "AFR", "AMR", "EAS", "SAS", "OCE", "MEN"]
sums = {sp: [0.0] * K for sp in SPS}; cnts = {sp: 0 for sp in SPS}
for s, q in zip(samples, Qu):
    sp = exemplars.get(s)
    if sp is None: continue
    cnts[sp] += 1
    for c in range(K): sums[sp][c] += q[c]
meanQ = {sp: [(sums[sp][c] / cnts[sp] if cnts[sp] else 0.0) for c in range(K)] for sp in SPS}
cells = sorted(((meanQ[sp][c], sp, c) for sp in SPS for c in range(K)), reverse=True)
cluster_of_sp = {}; used_cl = set(); used_sp = set()
for val, sp, c in cells:
    if sp in used_sp or c in used_cl: continue
    cluster_of_sp[sp] = c; used_sp.add(sp); used_cl.add(c)
print(f"[discover] cluster map (by mean-Q): {cluster_of_sp}", file=sys.stderr)
for sp in SPS:
    if sp in cluster_of_sp:
        c = cluster_of_sp[sp]
        runner = sorted((meanQ[sp][cc], cc) for cc in range(K))[-2]
        print(f"[discover]   {sp:3s} -> cluster {c} "
              f"(meanQ={meanQ[sp][c]:.3f}, n_exemplar={cnts[sp]}, "
              f"2nd-best cluster {runner[1]} @ {runner[0]:.3f})", file=sys.stderr)

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
    n_pass = len(picked)
    picked = picked[:cap]             # cap at top-N per super-pop
    with open(out, "w") as g:
        for s, _ in picked: g.write(s + "\n")
    n_cand_in_panel = sum(1 for s in samples if s in cand)
    print(f"[discover] {sp}: cluster {cl}, {n_pass}/{n_cand_in_panel} "
          f"candidates-in-panel pass Q>={thresh}; kept top {len(picked)} "
          f"(cap {cap}) -> {out}", file=sys.stderr)
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

# ---------------------------------------------------------------------------
# OTH -- the leftover cluster(s). With K=8 and 7 named super-pops exactly one
# cluster is unclaimed; it is the empirical "does not fit an HGDP super-pop"
# group. Unlike the other 7 it has NO metadata candidate pool, so its anchors
# are simply the top-N samples loading on that cluster, drawn from the whole
# panel. This is what makes supervised K=8 well-posed (ADMIXTURE --supervised
# needs exactly K distinct labels in the .pop file).
# ---------------------------------------------------------------------------
leftover = [c for c in range(K) if c not in used_cl]
n_oth = 0
if not leftover:
    print(f"[discover] no leftover cluster at K={K} (all {K} claimed by named "
          f"super-pops); writing empty {oth_out}", file=sys.stderr)
    open(oth_out, "w").close()
else:
    if len(leftover) > 1:
        print(f"[discover] NOTE: {len(leftover)} unclaimed clusters {leftover}; "
              f"using the one with the most high-Q members as OTH.", file=sys.stderr)
    # Pick the leftover cluster with the most members above threshold.
    best = max(leftover, key=lambda c: sum(1 for q in Qu if q[c] >= thresh))
    picked = sorted(((s, q[best]) for s, q in zip(samples, Qu) if q[best] >= thresh),
                    key=lambda x: -x[1])
    n_pass = len(picked); picked = picked[:cap]
    with open(oth_out, "w") as g:
        for s, _ in picked: g.write(s + "\n")
    n_oth = len(picked)
    cluster_of_sp["OTH"] = best
    print(f"[discover] OTH: leftover cluster {best}, {n_pass} samples pass "
          f"Q>={thresh}; kept top {n_oth} (cap {cap}) -> {oth_out}", file=sys.stderr)
    if picked:
        print(f"           top: {picked[0][0]} Q={picked[0][1]:.4f} | "
              f"min kept Q={picked[-1][1]:.4f}", file=sys.stderr)
    if n_oth == 0:
        print(f"[discover] NOTE: OTH cluster {best} exists but no sample cleared "
              f"Q>={thresh}. Supervised K={K} will have only "
              f"{len(cluster_of_sp)-1} labels and ADMIXTURE will refuse to run; "
              f"either lower ANCHOR_MIN_Q or drop K to {K-1}.", file=sys.stderr)
PYEOF

N_OCE=$(awk 'NF' "$OCE_ANCHORS" 2>/dev/null | wc -l)
N_MEN=$(awk 'NF' "$MEN_ANCHORS" 2>/dev/null | wc -l)
N_OTH=$(awk 'NF' "$OTH_ANCHORS" 2>/dev/null | wc -l)
echo "[$(date +%T)] discovered anchors: OCE=${N_OCE}, MEN=${N_MEN}, OTH=${N_OTH}"

# ---------------------------------------------------------------------------
# 5. STAGE B -- build supervised .pop file (7 labels) in .fam row order.
# ---------------------------------------------------------------------------
POP_FILE="${PRUNED}.pop"
echo "[$(date +%T)] build supervised .pop file (K=${K})"
python3 - "${PRUNED}.fam" "$POP_FILE" \
    "${REFS}/eur_rfmix.txt" "${REFS}/afr_rfmix.txt" "${REFS}/amr_rfmix.txt" \
    "${REFS}/eas_rfmix.txt" "${REFS}/sas_rfmix.txt" \
    "$OCE_ANCHORS" "$MEN_ANCHORS" "$OTH_ANCHORS" "$MXB_POPINFO" <<'PYEOF'
import sys
(fam, outp, eur_f, afr_f, amr_f, eas_f, sas_f, oce_f, men_f, oth_f, mxb_pop_f) = sys.argv[1:12]
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
           ("OCE", load(oce_f)), ("MEN", load(men_f)), ("OTH", load(oth_f))]
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
# 7. STAGE C -- parse both .Q -> lai_ref_panel_samples.k${K}.tsv + homogeneity.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] write lai_ref_panel_samples.k${K}.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv"
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

def assign_map(anchor_labels, Q, K):
    # Match clusters -> super-pop names by MEAN Q-loading of labeled anchors,
    # one-to-one greedy-by-max-cell. Robust to argmax-splitting for the
    # poorly-separated OCE/MEN clusters (see discovery-stage note). Same method
    # for supervised (anchors fixed ~1.0 on own component -> trivially correct)
    # and unsupervised Q.
    sps = sorted({a for a in anchor_labels if a != "-"})
    sums = {sp: [0.0] * K for sp in sps}; cnts = {sp: 0 for sp in sps}
    for a, q in zip(anchor_labels, Q):
        if a == "-": continue
        cnts[a] += 1
        for c in range(K): sums[a][c] += q[c]
    meanQ = {sp: [(sums[sp][c] / cnts[sp] if cnts[sp] else 0.0) for c in range(K)] for sp in sps}
    cells = sorted(((meanQ[sp][c], sp, c) for sp in sps for c in range(K)), reverse=True)
    m = {}; used_cl = set(); used_sp = set()
    for val, sp, c in cells:
        if sp in used_sp or c in used_cl: continue
        m[c] = sp; used_sp.add(sp); used_cl.add(c)
    for cl in range(K): m.setdefault(cl, f"UNKc{cl+1}")
    return m

sup_cluster_to_sp   = assign_map(anchors, Qs, K)
unsup_cluster_to_sp = assign_map(anchors, Qu, K)
sup_sp_to_cluster   = {v: k for k, v in sup_cluster_to_sp.items()}
unsup_sp_to_cluster = {v: k for k, v in unsup_cluster_to_sp.items()}
print(f"sup_cluster_pop_map   = {sup_cluster_to_sp}",   file=sys.stderr)
print(f"unsup_cluster_pop_map = {unsup_cluster_to_sp}", file=sys.stderr)

POP_ORDER = ["AFR", "AMR", "EAS", "EUR", "MEN", "OCE", "OTH", "SAS"]
VALID_SUPERPOPS = set(POP_ORDER)

def cohort_of(s):
    if s in sgdp_sp: return "SGDP"
    if s.startswith("MXB"): return "MXB"
    if s.startswith(("HGDP", "LP6005", "SS")): return "HGDP"
    if s.startswith(("HG", "NA")): return "1KG"
    return "OTHER"

def super_pop_of(s):
    # BUGFIX: the K=5 allowlist was ("AFR","AMR","EUR","EAS") -- only 4 values --
    # while make_sample_groups.sh emits 8 (AFR AMR EUR EAS SAS CSA OCE MID->MEN).
    # OCE and MEN therefore fell through to "NA" and were hard-excluded by
    # `elif sp == "NA": included = 0`. At K=7 that silently discarded EVERY
    # HGDP+1KG Oceanian and Middle Eastern sample AFTER discovering anchors for
    # them -- the "argmax == assigned super_pop matched nothing" failure.
    # The allowlist must cover the full K-way vocabulary.
    if s in mxb_pop: return "AMR"
    if s in sgdp_sp: return sgdp_sp[s]          # OCE/MEN/AFR/EUR/EAS/SAS/AMR
    raw = (sg.get(s, "NA") or "NA").strip().upper()
    # gnomAD folds 1KG-SAS into CSA. Map CSA -> SAS for ALL CSA samples, not
    # just those in sas_rfmix.txt: gating membership on the ANCHOR file
    # conflated "is a training exemplar" with "belongs to this super-pop" and
    # dropped all 504 1KG SAS samples (BEB/GIH/ITU/PJL/STU) to NA.
    if raw == "CSA": return "SAS"
    if raw == "MID": return "MEN"               # defensive: pre-rename value
    if raw in VALID_SUPERPOPS: return raw
    return "NA"

def pop_of(s):
    if s in mxb_pop: return mxb_pop[s]
    if s in sgdp_pop: return sgdp_pop[s]
    return gnomad_pop.get(s, "NA")

# --------------------------------------------------------------------------
# LABEL-VOCABULARY CROSS-CHECK (the diagnostic that was missing).
# included=1 requires `argmax_super_pop == assigned_super_pop`, a STRING
# comparison between two independently-built vocabularies. If they diverge the
# count silently collapses toward zero and the old per-super-pop summary could
# not show it, because it filtered to included==1 -- the NA bucket was
# invisible by construction. Print both sides, the NA census, and abort if the
# overlap is empty or if a whole super-pop would be zeroed out.
# --------------------------------------------------------------------------
assigned_counts = Counter(super_pop_of(s) for s in samples)
argmax_vocab = sorted(set(sup_cluster_to_sp.values()))
assigned_vocab = sorted(k for k in assigned_counts if k != "NA")
print("[labels] argmax vocabulary   : " + ", ".join(argmax_vocab), file=sys.stderr)
print("[labels] assigned vocabulary : " + ", ".join(assigned_vocab), file=sys.stderr)
print("[labels] assigned census     : "
      + ", ".join(f"{k}={assigned_counts[k]}" for k in sorted(assigned_counts)),
      file=sys.stderr)
overlap = set(argmax_vocab) & set(assigned_vocab)
only_argmax   = sorted(set(argmax_vocab)   - set(assigned_vocab))
only_assigned = sorted(set(assigned_vocab) - set(argmax_vocab))
if only_argmax:
    print(f"[labels] WARNING: argmax-only labels (no sample assigned to them): "
          f"{only_argmax}", file=sys.stderr)
if only_assigned:
    print(f"[labels] WARNING: assigned-only labels (no cluster carries them, so "
          f"every such sample is excluded): {only_assigned}", file=sys.stderr)
if not overlap:
    print("[labels] FATAL: argmax and assigned vocabularies are DISJOINT. Every "
          "sample would be excluded. Refusing to write a degenerate panel.",
          file=sys.stderr)
    sys.exit(4)
n_na = assigned_counts.get("NA", 0)
if n_na > 0.5 * len(samples):
    print(f"[labels] FATAL: {n_na}/{len(samples)} samples have assigned "
          f"super_pop=NA (>50%). super_pop_of() is dropping a vocabulary the "
          f"upstream metadata emits. Refusing to continue.", file=sys.stderr)
    sys.exit(5)
if any(c in (f"UNKc{i+1}" for i in range(K)) for c in argmax_vocab):
    print(f"[labels] WARNING: unclaimed clusters present {argmax_vocab} -- a "
          f"super-pop failed to claim a cluster at K={K}.", file=sys.stderr)

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

# --------------------------------------------------------------------------
# PER-POPULATION MEAN-Q TABLE (the refinement artifact).
# Rows = population, columns = mean UNSUPERVISED Q on each named cluster.
# This is the evidence for deciding, from data rather than from a hand-written
# list, where each population actually sits: whether Mozabite loads AFR or MEN,
# whether the Beringian set (Tlingit/Aleut/Eskimo_*) loads EAS or AMR or OTH,
# whether Central Asians (Tajik/Kyrgyz/Uygur/Hazara) split rather than cluster,
# and where the SGDP populations with no HGDP counterpart (Australian, Maori,
# Hawaiian, Dusun, Igorot) fall. Review this, then re-run with overrides or
# exclusions as needed.
# --------------------------------------------------------------------------
import os
pop_rows = defaultdict(list)
for s, qu in zip(samples, Qu):
    pop_rows[(cohort_of(s), pop_of(s), super_pop_of(s))].append(qu)
pop_tbl = os.path.join(os.path.dirname(out),
                       f"population_cluster_meanQ.k{K}.tsv")
with open(pop_tbl, "w") as g:
    g.write("cohort\tpopulation\tassigned_super_pop\tn\t"
            + "\t".join(f"meanQ_{c}" for c in POP_ORDER)
            + "\targmax_cluster\tmax_meanQ\tsecond_cluster\tsecond_meanQ\n")
    for (co, p, sp), rows in sorted(pop_rows.items()):
        n = len(rows)
        mean_by_cluster = [sum(r[c] for r in rows) / n for c in range(K)]
        # Re-express in POP_ORDER using the unsupervised cluster->name map.
        named = {unsup_cluster_to_sp[c]: mean_by_cluster[c] for c in range(K)}
        ordered = [named.get(c, 0.0) for c in POP_ORDER]
        rank = sorted(named.items(), key=lambda kv: -kv[1])
        top, top_v = rank[0]
        sec, sec_v = (rank[1] if len(rank) > 1 else ("NA", 0.0))
        g.write(f"{co}\t{p}\t{sp}\t{n}\t"
                + "\t".join(f"{x:.4f}" for x in ordered)
                + f"\t{top}\t{top_v:.4f}\t{sec}\t{sec_v:.4f}\n")
print(f"wrote per-population mean-Q table -> {pop_tbl}", file=sys.stderr)

# Flag populations whose empirical cluster disagrees with their metadata label,
# and populations that fail to concentrate anywhere (admixed / clinal).
print("[refine] populations whose argmax cluster != assigned super_pop:",
      file=sys.stderr)
n_flag = 0
for (co, p, sp), rows in sorted(pop_rows.items()):
    n = len(rows)
    mean_by_cluster = [sum(r[c] for r in rows) / n for c in range(K)]
    named = {unsup_cluster_to_sp[c]: mean_by_cluster[c] for c in range(K)}
    top, top_v = max(named.items(), key=lambda kv: kv[1])
    if sp != "NA" and top != sp:
        print(f"[refine]   {co:5s} {p:20s} assigned={sp:4s} -> loads {top} "
              f"({top_v:.3f})  n={n}", file=sys.stderr); n_flag += 1
    elif top_v < 0.60:
        print(f"[refine]   {co:5s} {p:20s} assigned={sp:4s} DIFFUSE: best "
              f"{top} only {top_v:.3f} -- admixed/clinal, poor anchor material",
              file=sys.stderr); n_flag += 1
if n_flag == 0:
    print("[refine]   (none -- every population loads on its assigned cluster)",
          file=sys.stderr)
PYEOF

N_INCL=$(awk -F'\t' 'NR>1 && $NF==1' "$SAMPLES_TSV" | wc -l)
N_TOT=$(awk 'NR>1' "$SAMPLES_TSV" | wc -l)
echo "[$(date +%T)] included: ${N_INCL}/${N_TOT} samples"
echo "  per super_pop (INCLUDED):"
awk -F'\t' 'NR>1 && $NF==1 {print "    " $4}' "$SAMPLES_TSV" | sort | uniq -c
# The old script printed only the included rows, which made a total label
# mismatch invisible (NA never showed up). Always print the excluded census too.
echo "  per super_pop (EXCLUDED):"
awk -F'\t' 'NR>1 && $NF==0 {print "    " $4}' "$SAMPLES_TSV" | sort | uniq -c
echo "  excluded reason breakdown (assigned vs argmax, top 15):"
# Resolve columns by HEADER NAME, not by offset -- POP_ORDER length changes the
# offsets, and an off-by-one here silently reports the wrong column.
awk -F'\t' '
  NR==1 { for (i=1; i<=NF; i++) h[$i]=i; next }
  $h["included"]==0 {
      # group on the label PAIR only -- adding max_Q would make every row
      # unique and defeat the uniq -c aggregation.
      print "    assigned=" $h["super_pop"] "  argmax=" $h["argmax_pop_sup"]
  }' "$SAMPLES_TSV" | sort | uniq -c | sort -rn | head -15

awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_${K}pop_samples.txt"

# ---------------------------------------------------------------------------
# 8. Subset each per-chr phased BCF to the homogeneous samples (K=7 set).
# ---------------------------------------------------------------------------
echo "[$(date +%T)] subset per-chr BCFs to homog samples"
HOMOG_KEEP="${HOMOG_DIR}/homog_${K}pop_samples.txt"
for CHR in {1..22}; do
    SRC="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    DST="${HOMOG_SUB}/merged_chr${CHR}.homog_${K}pop.bcf"
    if [[ -s "$DST" && -s "${DST}.csi" ]]; then continue; fi
    bcftools view "$SRC" -S "$HOMOG_KEEP" --force-samples --threads "$THREADS" -Ob -o "$DST"
    bcftools index "$DST"
done

echo "[$(date +%T)] DONE."
echo "  samples TSV      : $SAMPLES_TSV"
echo "  discovered OCE   : $OCE_ANCHORS ($N_OCE)"
echo "  discovered MEN   : $MEN_ANCHORS ($N_MEN)"
echo "  discovered OTH   : $OTH_ANCHORS ($N_OTH)"
echo "  homog BCFs       : ${HOMOG_SUB}/merged_chr*.homog_${K}pop.bcf"
