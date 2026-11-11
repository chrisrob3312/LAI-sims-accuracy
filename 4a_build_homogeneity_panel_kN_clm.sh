#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# SLURM directives
# ----------------------------------------------------------------------------
#SBATCH --job-name=4a_homog_kN
#SBATCH --partition=atkinson,mhgcp
# DO NOT REMOVE mhgcp-t* FROM THIS EXCLUSION TO GAIN CONCURRENCY.
#
# Those nodes' CPUs lack instructions the SHAPEIT5 / bcftools builds in this
# conda env are compiled with. A task dispatched to one dies instantly with
# SIGILL -- exit 132, ~3 seconds wall clock, ~2 MB RSS, and an EMPTY log,
# because the binary never reaches its first write.
#
# Measured 2026-11-07 on a 594-task chunk array run with the exclusion lifted:
#   270 / 270 tasks that landed on t02,t05,t06,t07,t08,t10,t11  FAILED
#     0 / 57  tasks that landed on a01, m00                     failed
#
# They look attractive because sinfo shows them idle with free CPUs. They are
# idle because nothing that needs these instruction sets can run there. An
# earlier 1d_almarri failure on t02 -- empty logs, 4 MB RSS, never explained
# -- was the same thing.
#
# mhgcp-c02 WAS on this list and has been removed deliberately. It has the
# same specs as c00 and c01 (56 CPUs, 232 GB), both of which run this software
# without trouble and were never excluded, and no evidence ties c02 to the
# instruction-set problem above -- no task ever landed on it. The exclusion
# looks like it dates from c02 simply being busy. Re-add it only if a task
# actually fails there, and write down what happened.
#SBATCH --exclude=mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
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
# K=8: the 7 named super-pops (AFR AMR EAS EUR MEN OCE SAS) plus OTH, a cluster
# the data forms that is not one of them.
#
#   0. MARKERS. Concat the 22 phased chromosomes; restrict to the sites 1c
#      found adequately called in EVERY cohort x super-pop group and overall
#      (ADMIX_SITE_FILTER); MAF >= 0.01; LD-prune 500 50 0.1 over all samples.
#
#   STAGE A (unsupervised fit, then DATA-FIRST anchor discovery):
#     Unsupervised ADMIXTURE at K. For each cluster, rank every eligible
#     sample by its loading, keep the top ANCHOR_CAP at Q >= ANCHOR_MIN_Q, and
#     NAME the cluster from who those anchors are (plurality metadata
#     super-pop, purity >= ANCHOR_NAME_PURITY; the stronger claim wins a
#     contested super-pop; one unnamed cluster becomes OTH)
#     -> reference_ids/<sp>_anchors.k${K}.txt for the 7, and OTH.
#     Every cluster's composition, and why it was or was not named, goes to
#     discovery_clusters.k${K}.tsv. A cluster carried by one minority cohort
#     across super-pops is flagged as a likely processing batch.
#
#   STAGE B (supervised, anchor-INFORMED):
#     ADMIXTURE --supervised at K with the .pop built from those anchors.
#     Rerun only when the anchors (the .pop) or the bed changed.
#
#   STAGE C (homogeneity filter + refinement evidence):
#     included = 1 iff supervised max_Q >= the tiered threshold AND the
#     supervised argmax super-pop matches the sample's assigned super-pop.
#     Emits population_cluster_meanQ.k${K}.tsv (mean unsupervised Q per
#     POPULATION per cluster) and [refine] lines for populations whose argmax
#     disagrees with their label or that fail to concentrate.
#
# Every stage reuses its output only if it is newer than all of its inputs, so
# a re-phased panel can never be met by a previous run's concat, bed or Q.
#
# Outputs in $HOMOG_DIR (all paths carry ${K}, so K=5/7/8 never collide):
#   admixture_k${K}/lai_ref.{bed,bim,fam,${K}.Q,${K}.P,pop}
#   admixture_k${K}/admix_sites.all.bed                (ADMIXTURE-eligible sites)
#   admixture_k${K}/unsup/lai_ref.${K}.Q
#   discovery_clusters.k${K}.tsv                       (every cluster: named or why not)
#   reference_ids/<sp>_anchors.k${K}.txt, oth_anchors.k${K}.txt
#   anchor_summary.k${K}.tsv
#   lai_ref_panel_samples.k${K}.tsv
#   population_cluster_meanQ.k${K}.tsv                 (refinement evidence)
#   homog_${K}pop_samples.txt
#   homog_${K}pop/merged_chr${CHR}.homog_${K}pop.bcf{,.csi}
#
# NOT read or written by discovery: reference_ids/*_rfmix.txt, the curated
# RFMix reference panels read by the 3b_*.sh scripts. Their overlap with the
# discovered anchors is reported.
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
# A COUNT floor to go with the Q floor. ANCHOR_MIN_Q stops an ADMIXED sample
# from becoming an anchor; this stops a super-pop's whole prior from resting on
# a handful of samples. Both matter and they fail differently:
#
#   no Q floor      -> an admixed sample DEFINES that ancestry in the prior,
#                      and every Q estimate for the component inherits it
#   no count floor  -> the label is populated, so ADMIXTURE runs happily, but
#                      the component is estimated from 2-3 exemplars
#
# The second is the quiet one. 4a only ever failed on a label with ZERO
# anchors, because that is what makes ADMIXTURE --supervised refuse to run; a
# super-pop with two passed silently and then fed the homogeneity call, the
# donor/reference split, the test-genome sample_type and the GWAS ancestry
# covariate. OCE, MEN and SAS are where this is likely to bite.
#
# A warning rather than an error by default: a thin super-pop may be a fact to
# report in the methods rather than a reason to stop. ANCHOR_FLOOR_FATAL=1
# makes it stop.
ANCHOR_FLOOR="${ANCHOR_FLOOR:-10}"
ANCHOR_FLOOR_FATAL="${ANCHOR_FLOOR_FATAL:-0}"
# Optional tiering of the Q floor, mirroring HOMOG_THRESHOLD_STRICT. OFF by
# default: 0.93 uniform is what the NAT/Mexican accuracy pipeline used, and
# changing it silently would move every anchor set. Set ANCHOR_MIN_Q_STRICT to
# apply a stricter floor to the deeply sampled super-pops, which can meet it.
ANCHOR_MIN_Q_STRICT="${ANCHOR_MIN_Q_STRICT:-}"
ANCHOR_STRICT_SUPERPOPS="${ANCHOR_STRICT_SUPERPOPS:-EUR,EAS,AFR}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.93}"

# TIERED homogeneity cut. HOMOG_THRESHOLD above is the FLOOR and the default
# for any super-pop not named here; HOMOG_THRESHOLD_STRICT applies to the
# groups that have enough samples to afford it.
#
# WHY TWO NUMBERS. The cut trades purity against panel size, and the two sides
# of that trade are not symmetric across super-pops. EUR, EAS and AFR are
# richly sampled in HGDP+1KG, so 0.95 costs them little and buys cleaner
# reference haplotypes. AMR, OCE, MEN and SAS are not: at 0.95 they lose
# members the panel cannot spare, and a thin reference panel for exactly the
# ancestries whose local-ancestry calls are hardest is the failure this whole
# project is measuring. So the under-represented groups sit at 0.93.
#
# This is a deliberate asymmetry and it belongs in the methods: the panel is
# not "samples above one Q", it is "samples above a per-super-pop Q", and the
# per-super-pop counts and cut are printed in the stage C log and carried in
# homog_composition.k${K}.tsv.
HOMOG_THRESHOLD_STRICT="${HOMOG_THRESHOLD_STRICT:-0.95}"
HOMOG_STRICT_SUPERPOPS="${HOMOG_STRICT_SUPERPOPS:-EUR,EAS,AFR}"

# SGDP candidate pools (from scripts/fetch-sgdp-anchors.py). Used as exemplars
# for cluster labeling + as the pool OCE/MEN anchors are discovered from.
# These are CANDIDATES, not anchors -- anchors are the top-N discovered here.
SGDP_MANIFEST="${SGDP_MANIFEST:-${REFS}/sgdp_manifest.tsv}"

# Almarri 2021 manifest (sample<TAB>pop). Defined HERE, not next to the labeling
# code, because the MEN draw pool below reads it and bash needs it set first.
# Keyed on the PANEL sample name, which is appg.meta's long_ID
# (longranger222_wgs_<n>_APPG<id>_GRCh38_gatk) -- not the pop_ID (UAE19) that
# the distributed .fam uses. Build it on the cluster:
#   bash scripts/build-almarri-manifest.sh > reference_ids/almarri_manifest.tsv
ALMARRI_META="${ALMARRI_META:-${REFS}/almarri_manifest.tsv}"
ALMARRI_COUNTRIES="${ALMARRI_COUNTRIES:-IRAQ JORDAN OMAN SAUDI SYRIA UAE YEMEN}"

# DATA-FIRST ANCHOR DISCOVERY (stage A).
#
# Anchors come from the UNSUPERVISED fit, and nothing else decides them:
#   1. for each of the K clusters, rank EVERY eligible panel sample by its
#      loading on that cluster;
#   2. keep those at Q >= ANCHOR_MIN_Q, the top ANCHOR_CAP of them;
#   3. name the cluster from WHO those anchors are: the plurality metadata
#      super-pop among them, which must hold at least ANCHOR_NAME_PURITY of
#      the set. When two clusters claim one super-pop, the one with more
#      supporting samples (passing members with that label) keeps it. One
#      cluster left unnamed is OTH.
# Eligible = every panel sample except named exclusions (panel_exclude_ids)
# and the copy 1e dropped from each duplicate pair.
#
# This replaces the pool design, in which curated label lists named the
# clusters by mean Q and per-super-pop draw pools limited who could be drawn.
# That named clusters a priori; and when seed_for() expanded `$sp` inside the
# same `local` that assigned it (job 3747721), the curated lists came back
# empty and the pools fell back to SGDP alone -- 276 of 3,921 samples were
# the only candidates for every super-pop. The *_rfmix.txt lists are no longer
# read by discovery; their overlap with the discovered anchors is reported.
#
# The metadata (sample_groups.tsv with CSA->SAS and MID->MEN, MXB->AMR, the
# SGDP manifest, the Almarri manifest) is used ONLY to name a cluster the data
# already formed, and to report on it.
ANCHOR_NAME_PURITY="${ANCHOR_NAME_PURITY:-0.5}"
# A cluster whose anchors come mostly (>= this share) from one cohort that is a
# minority of the panel is flagged as a likely processing batch, not ancestry.
ANCHOR_COHORT_WARN="${ANCHOR_COHORT_WARN:-0.8}"
PANEL_EXCLUDE="${PANEL_EXCLUDE:-${REFS}/panel_exclude_ids.txt}"
DUP_DROP="${DUP_DROP:-${REFS}/panel_duplicate_drop_ids.txt}"

# ADMIXTURE MARKER FILTER -- per cohort x super-pop missingness, then overall.
#
# 1c writes, per chromosome, the soft-union sites at which every cohort x
# super-pop group (HGDP1KG_EUR, SGDP_AFR, MXB_AMR, ALMARRI_MEN, ...) is at most
# ADMIX_GROUP_MAXMISS missing and the whole panel at most ADMIX_TOTAL_MAXMISS,
# both measured on the PRE-phasing genotypes. 4a restricts the ADMIXTURE marker
# set to those sites before the MAF filter and LD pruning.
#
# Why here and not with --geno: SHAPEIT5 fills every missing genotype, so the
# phased panel 4a reads has none left to filter. And why per group: a cohort
# that is a small share of the panel can be missing entirely and still pass a
# panel-wide cut (Almarri is 3.5% of samples). On chr22 in November every SGDP
# sample had no call at 38% of the ADMIXTURE markers and every MXB sample at
# 62%: those genotypes were filled in by imputation, not called.
#
# Only the ADMIXTURE marker set is restricted. The phased panel used for LAI,
# RFMix and simulation donors keeps the full soft-union (and with it the
# population-specific variants the soft-union exists to keep).
ADMIX_SITE_FILTER="${ADMIX_SITE_FILTER:-1}"
ADMIX_SITES_TPL="${ADMIX_SITES_TPL:-${PHASED_DIR}/merged_chr%s.admix_sites.bed}"
SUPERPOPS="EUR AFR AMR EAS SAS OCE MEN"

# DISCOVERED anchor outputs, written by STAGE A, consumed by STAGE B only.
#
# They are deliberately NOT named *_rfmix.txt. Two reasons, both load-bearing:
#   1. *_rfmix.txt are the curated RFMix reference panels. Writing discovery
#      output over them would silently change which haplotypes RFMix is given.
#   2. A discovered list is also its own super-pop's seed on the next run. If
#      discovery overwrote the seed, every rerun would re-draw from its own
#      capped output and the pool would erode by the cap on each pass.
# OCE/MEN previously wrote to oce_rfmix.txt/men_rfmix.txt, which nothing else
# reads; they move here too so a *_rfmix.txt name always means "curated".
anchors_for () { echo "${REFS}/$(echo "$1" | tr 'A-Z' 'a-z')_anchors.k${K}.txt"; }

# OTH = the empirical leftover cluster at K=8 (no metadata candidate pool).
OTH_ANCHORS="${OTH_ANCHORS:-${REFS}/oth_anchors.k${K}.txt}"

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
echo "[$(date +%T)] K=${K}  ANCHOR_MIN_Q=${ANCHOR_MIN_Q}${ANCHOR_MIN_Q_STRICT:+ (strict ${ANCHOR_MIN_Q_STRICT} for ${ANCHOR_STRICT_SUPERPOPS})}  ANCHOR_CAP=${ANCHOR_CAP}  ANCHOR_FLOOR=${ANCHOR_FLOOR}  HOMOG_THRESHOLD=${HOMOG_THRESHOLD}"
echo "[$(date +%T)] verifying inputs"
for CHR in {1..22}; do
    F="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    [[ -s "$F" && -s "${F}.csi" ]] || { echo "ERROR: missing $F (phasing not done)"; exit 1; }
done
if [[ "$ADMIX_SITE_FILTER" == "1" ]]; then
    for CHR in {1..22}; do
        F="$(printf "$ADMIX_SITES_TPL" "$CHR")"
        [[ -s "$F" ]] || {
            echo "ERROR: missing ADMIXTURE site list $F" >&2
            echo "       1c writes it (per-cohort x super-pop missingness on the" >&2
            echo "       pre-phasing soft-union). Re-run the current 1c, or set" >&2
            echo "       ADMIX_SITE_FILTER=0 to use every soft-union site." >&2
            exit 1; }
    done
fi
[[ -s "$SAMPLE_GROUPS" ]] || { echo "ERROR: missing $SAMPLE_GROUPS"; exit 1; }
[[ -s "$GNOMAD_META"   ]] || { echo "ERROR: missing $GNOMAD_META"; exit 1; }
[[ -s "$MXB_POPINFO"   ]] || { echo "ERROR: missing $MXB_POPINFO"; exit 1; }
[[ -s "$SGDP_MANIFEST" ]] || { echo "ERROR: missing $SGDP_MANIFEST (run fetch-sgdp-anchors.py)"; exit 1; }

# REBUILD WHEN AN INPUT IS NEWER. Every stage below used to skip whenever its
# output merely existed, so a re-phased panel would have been met by the old
# concat, the old LD-pruned bed and the old unsupervised Q, silently. A stage
# now reuses its output only if that output is newer than every input.
# newer_than <output> <input>... : true (= rebuild) if the output is missing or
# any input is newer than it.
newer_than () {
    local out="$1"; shift
    [[ -s "$out" ]] || return 0
    local f
    for f in "$@"; do [[ "$f" -nt "$out" ]] && return 0; done
    return 1
}

# ---------------------------------------------------------------------------
# 1. Concat all 22 chrs to a single BCF. Written beside the target and moved
#    into place, so a killed task cannot leave a truncated concat that is
#    newer than its inputs.
# ---------------------------------------------------------------------------
CONCAT="${ADM_DIR}/all_chr.phased.bcf"
PHASED_LIST=()
for CHR in {1..22}; do
    PHASED_LIST+=("${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf")
done
if newer_than "$CONCAT" "${PHASED_LIST[@]}"; then
    echo "[$(date +%T)] concat 22 chrs"
    printf '%s\n' "${PHASED_LIST[@]}" > "${ADM_DIR}/concat_list.txt"
    rm -f "${CONCAT}.partial" "${CONCAT}.partial.csi"
    bcftools concat --threads "$THREADS" -Ob -o "${CONCAT}.partial" -f "${ADM_DIR}/concat_list.txt"
    bcftools index "${CONCAT}.partial"
    mv "${CONCAT}.partial.csi" "${CONCAT}.csi"
    mv "${CONCAT}.partial" "$CONCAT"
else
    echo "[$(date +%T)] concat newer than every phased chromosome -- reusing ${CONCAT}"
fi

# ---------------------------------------------------------------------------
# 1b. ADMIXTURE marker filter (see ADMIX_SITE_FILTER at the top).
# ---------------------------------------------------------------------------
EXTRACT_SITES=()
PRUNE_INPUTS=("$CONCAT")
ADMIX_SITES_ALL="${ADM_DIR}/admix_sites.all.bed"
if [[ "$ADMIX_SITE_FILTER" == "1" ]]; then
    SITE_LIST=()
    for CHR in {1..22}; do SITE_LIST+=("$(printf "$ADMIX_SITES_TPL" "$CHR")"); done
    if newer_than "$ADMIX_SITES_ALL" "${SITE_LIST[@]}"; then
        cat "${SITE_LIST[@]}" > "${ADMIX_SITES_ALL}.partial"
        mv "${ADMIX_SITES_ALL}.partial" "$ADMIX_SITES_ALL"
    fi
    echo "[$(date +%T)] ADMIXTURE marker filter ON: $(wc -l < "$ADMIX_SITES_ALL") soft-union" \
         "sites pass the per-cohort x super-pop and overall missingness rule" \
         "(thresholds and per-group counts: $(printf "$ADMIX_SITES_TPL" N | sed 's/\.bed$//').summary.tsv)"
    EXTRACT_SITES=(--extract bed1 "$ADMIX_SITES_ALL")
    PRUNE_INPUTS+=("$ADMIX_SITES_ALL")
else
    echo "[$(date +%T)] ADMIXTURE marker filter OFF (ADMIX_SITE_FILTER=0): all soft-union sites eligible"
fi

# ---------------------------------------------------------------------------
# 2. plink2 import + LD-prune, over all samples jointly (as in July).
# ---------------------------------------------------------------------------
PRUNED="${ADM_DIR}/lai_ref"
if newer_than "${PRUNED}.bed" "${PRUNE_INPUTS[@]}"; then
    echo "[$(date +%T)] plink2 LD-prune (500 50 0.1 maf 0.01)"
    rm -f "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam" "${ADM_DIR}/lai_ref.prune."*
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' ${EXTRACT_SITES[@]+"${EXTRACT_SITES[@]}"} \
           --maf 0.01 --indep-pairwise 500 50 0.1 --out "${ADM_DIR}/lai_ref.prune"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' \
           --extract "${ADM_DIR}/lai_ref.prune.prune.in" \
           --max-alleles 2 --snps-only just-acgt \
           --make-bed --out "$PRUNED"
    echo "[$(date +%T)] pruned: $(wc -l < "${PRUNED}.bim") SNPs, $(wc -l < "${PRUNED}.fam") samples"
else
    echo "[$(date +%T)] LD-pruned bed newer than its inputs -- reusing ${PRUNED}.bed" \
         "($(wc -l < "${PRUNED}.bim") SNPs)"
fi

# ---------------------------------------------------------------------------
# 3. STAGE A -- UNSUPERVISED ADMIXTURE K (runs FIRST; discovery depends on it).
# ---------------------------------------------------------------------------
if newer_than "${UNSUP_DIR}/lai_ref.${K}.Q" "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam"; then
    echo "[$(date +%T)] admixture (unsupervised) K=${K} (slow)"
    rm -f "${UNSUP_DIR}/lai_ref.${K}.Q" "${UNSUP_DIR}/lai_ref.${K}.P"
    for ext in bed bim fam; do ln -sf "${PRUNED}.${ext}" "${UNSUP_DIR}/lai_ref.${ext}"; done
    cd "$UNSUP_DIR"
    admixture -j${THREADS} lai_ref.bed "${K}"
    cd - >/dev/null
else
    echo "[$(date +%T)] unsupervised Q newer than the bed -- reusing ${UNSUP_DIR}/lai_ref.${K}.Q"
fi
[[ -s "${UNSUP_DIR}/lai_ref.${K}.Q" ]] || { echo "ERROR: unsupervised ADMIXTURE failed"; exit 1; }
# ---------------------------------------------------------------------------
# 4. STAGE A -- data-first anchor discovery from the unsupervised Q.
#    Per cluster: rank every eligible sample, keep the top ANCHOR_CAP at
#    Q >= ANCHOR_MIN_Q, name the cluster from who those anchors are.
#    See DATA-FIRST ANCHOR DISCOVERY at the top for the rule and the history.
# ---------------------------------------------------------------------------
DISCOVERY_TABLE="${HOMOG_DIR}/discovery_clusters.k${K}.tsv"
echo "[$(date +%T)] data-first discovery from unsupervised K=${K}: per cluster, top ${ANCHOR_CAP} at Q >= ${ANCHOR_MIN_Q}, named by anchor metadata (purity >= ${ANCHOR_NAME_PURITY})"
python3 - "${PRUNED}.fam" "${UNSUP_DIR}/lai_ref.${K}.Q" "$K" "$ANCHOR_MIN_Q" "$ANCHOR_CAP" \
    "${ANCHOR_MIN_Q_STRICT:-}" "$ANCHOR_STRICT_SUPERPOPS" "$ANCHOR_NAME_PURITY" "$ANCHOR_COHORT_WARN" \
    "$SAMPLE_GROUPS" "$MXB_POPINFO" "$SGDP_MANIFEST" "$ALMARRI_META" "$PANEL_EXCLUDE" "$DUP_DROP" \
    "$REFS" "$SUPERPOPS" "$OTH_ANCHORS" "$DISCOVERY_TABLE" <<'PYEOF'
import sys, re, csv
from collections import Counter

(fam_f, qu_f, K, thresh, cap, strict_q, strict_sps, purity_min, cohort_warn,
 groups_f, mxb_f, sgdp_f, almarri_f, excl_f, dup_f, refs, superpops,
 oth_out, table_out) = sys.argv[1:20]
K = int(K); thresh = float(thresh); cap = int(cap)
strict_q = float(strict_q) if strict_q else None
strict_sps = {x for x in strict_sps.replace(",", " ").split() if x}
purity_min = float(purity_min); cohort_warn = float(cohort_warn)
SUPERPOPS = superpops.split()

def ids(path):
    try:
        with open(path) as fh:
            return {l.split()[0] for l in fh if l.strip() and not l.startswith("#")}
    except FileNotFoundError:
        return set()

# ---- metadata: used ONLY to name a cluster and to report on it ------------
groups = {}
with open(groups_f) as fh:
    for l in fh:
        f = l.split()
        if len(f) >= 2:
            g = f[1].upper()
            groups[f[0]] = {"CSA": "SAS", "MID": "MEN"}.get(g, g)
mxb = set()
with open(mxb_f) as fh:
    next(fh, None)
    for l in fh:
        if l.strip(): mxb.add(l.split("\t")[0].strip())
sgdp = {}
with open(sgdp_f) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sgdp[r["sample"]] = r["super_pop"]
almarri = set()
try:
    with open(almarri_f) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if f and f[0] and f[0] != "sample": almarri.add(f[0])
except FileNotFoundError:
    pass
ALM_RE = re.compile(r"(^longranger\d*_wgs_|(^|_)APPG\d+(_|$))")

def superpop(s):     # same resolution order as super_pop_of() in stage C
    if s in mxb: return "AMR"
    if s in sgdp: return sgdp[s]
    if s in almarri or ALM_RE.search(s): return "MEN"
    return groups.get(s, "NA")

def cohort(s):
    if s in sgdp: return "SGDP"
    if s in almarri or ALM_RE.search(s): return "ALMARRI"
    if s in mxb or s.startswith("MXB"): return "MXB"
    if s.startswith(("HGDP", "LP6005", "SS600")): return "HGDP"
    if re.match(r"^(HG|NA)\d", s): return "1KG"
    if s.startswith("CHMI"): return "CHMI"
    return "OTHER"

samples = [l.split()[1] for l in open(fam_f) if l.strip()]
Q = [list(map(float, l.split())) for l in open(qu_f) if l.strip()]
assert len(samples) == len(Q), f"{len(samples)} fam rows != {len(Q)} Q rows"
assert all(len(q) == K for q in Q), f"Q rows do not all have K={K} columns"

excluded = ids(excl_f) | ids(dup_f)
elig = [i for i, s in enumerate(samples) if s not in excluded]
SP = [superpop(s) for s in samples]
CO = [cohort(s) for s in samples]
print(f"[discover] {len(samples)} samples; {len(samples) - len(elig)} not eligible "
      f"(named exclusions + 1e duplicate drops); {len(elig)} ranked in every cluster",
      file=sys.stderr)
n_na = sum(1 for i in elig if SP[i] == "NA")
if n_na:
    print(f"[discover] NOTE: {n_na} eligible sample(s) have no metadata super-pop; they can "
          f"be anchors but do not count toward naming a cluster", file=sys.stderr)
panel_co = Counter(CO[i] for i in elig)

rows = []
for c in range(K):
    ranked = sorted(elig, key=lambda i: -Q[i][c])
    passing = [i for i in ranked if Q[i][c] >= thresh]
    top = passing[:cap]
    sp_top = Counter(SP[i] for i in top)
    co_top = Counter(CO[i] for i in top)
    labelled = sorted(((n, sp) for sp, n in sp_top.items() if sp != "NA"), key=lambda x: (-x[0], x[1]))
    plur, n_plur = (labelled[0][1], labelled[0][0]) if labelled else ("NA", 0)
    purity = n_plur / len(top) if top else 0.0
    support = sum(1 for i in passing if SP[i] == plur)
    mean_top = sum(Q[i][c] for i in top) / len(top) if top else 0.0
    # Batch signature: one minority cohort supplies the anchors AND they span
    # several super-pops (purity < 0.8). One cohort with one super-pop is not
    # flagged -- a population sampled by one cohort only looks like that -- but
    # the same super-pop's samples from OTHER cohorts are reported, because if
    # they barely load on the cluster it is the cohort, not the ancestry.
    flag, note = "", ""
    if top:
        co, n_co = co_top.most_common(1)[0]
        if n_co / len(top) >= cohort_warn and panel_co[co] / len(elig) < 0.5:
            if purity < 0.8:
                flag = (f"COHORT-DOMINATED({co} {n_co}/{len(top)} across "
                        f"{sum(1 for x in sp_top if x != 'NA')} super-pops; "
                        f"{100 * panel_co[co] / len(elig):.0f}% of panel)")
            else:
                peers = [i for i in elig if SP[i] == plur and CO[i] != co]
                if peers:
                    note = (f"single-cohort cluster ({co} {n_co}/{len(top)}): {plur} samples from "
                            f"other cohorts (n={len(peers)}) load {sum(Q[i][c] for i in peers) / len(peers):.2f} "
                            f"on it on average, vs {mean_top:.2f} for its anchors")
    rows.append(dict(c=c, passing=passing, top=top, sp_top=sp_top, co_top=co_top, plur=plur,
                     purity=purity, support=support, mean_top=mean_top, flag=flag, note=note))

# ---- name clusters: strongest claim on a super-pop wins it -----------------
name, taken = {}, {}
for r in sorted((r for r in rows if r["top"] and r["plur"] != "NA" and r["purity"] >= purity_min),
                key=lambda r: (-r["support"], -r["mean_top"], r["c"])):
    if r["plur"] not in taken:
        taken[r["plur"]] = r["c"]; name[r["c"]] = r["plur"]
unnamed = [r for r in rows if r["c"] not in name]
oth = None
if unnamed:
    best = max(unnamed, key=lambda r: (len(r["passing"]), r["mean_top"]))
    if best["top"]:
        oth = best["c"]; name[oth] = "OTH"

def why(r):
    if r["c"] in name and name[r["c"]] != "OTH": return "named"
    if not r["top"]: return f"no sample >= {thresh}"
    if r["plur"] == "NA": return "anchors carry no metadata super-pop"
    if r["purity"] < purity_min: return f"mixed ({r['plur']} {r['purity']:.2f} < {purity_min})"
    return f"second cluster for {r['plur']} (cluster {taken.get(r['plur'])} has more support)"

def fmt(cnt, k=4):
    return ",".join(f"{a}:{b}" for a, b in cnt.most_common(k)) or "-"

print(f"[discover] {'cl':>3s} {'name':5s} {'n>=thr':>6s} {'top':>4s} {'meanQtop':>8s}  "
      f"{'plurality':10s} {'purity':>6s} {'support':>7s}  top-anchor super-pops | cohorts", file=sys.stderr)
with open(table_out, "w") as g:
    g.write("cluster\tname\treason\tn_pass\tn_top\tmean_q_top\tplurality\tpurity\tsupport\t"
            "top_superpops\ttop_cohorts\tflag\tnote\n")
    for r in rows:
        nm = name.get(r["c"], "-")
        print(f"[discover] {r['c']:3d} {nm:5s} {len(r['passing']):6d} {len(r['top']):4d} {r['mean_top']:8.4f}  "
              f"{r['plur']:10s} {r['purity']:6.2f} {r['support']:7d}  {fmt(r['sp_top'])} | {fmt(r['co_top'])}"
              f"{'  ' + r['flag'] if r['flag'] else ''}", file=sys.stderr)
        if nm in ("-", "OTH"):
            print(f"[discover]       {'OTH: ' if nm == 'OTH' else 'unnamed: '}{why(r)}", file=sys.stderr)
        g.write(f"{r['c']}\t{nm}\t{why(r)}\t{len(r['passing'])}\t{len(r['top'])}\t{r['mean_top']:.4f}\t"
                f"{r['plur']}\t{r['purity']:.3f}\t{r['support']}\t{fmt(r['sp_top'], 9)}\t"
                f"{fmt(r['co_top'], 9)}\t{r['flag']}\t{r['note']}\n")
for r in rows:
    if r["flag"]:
        print(f"[discover] WARNING: cluster {r['c']} is {r['flag']}. A cluster carried by one cohort "
              f"across super-pops is the signature of a processing batch, not an ancestry.", file=sys.stderr)
    if r["note"]:
        print(f"[discover] NOTE: cluster {r['c']} is a {r['note']}.", file=sys.stderr)

# ---- write anchors ---------------------------------------------------------
def write(path, idx):
    with open(path, "w") as g:
        for i in idx: g.write(samples[i] + "\n")

for sp in SUPERPOPS:
    out = f"{refs}/{sp.lower()}_anchors.k{K}.txt"
    c = taken.get(sp)
    if c is None:
        write(out, [])
        print(f"[discover] {sp}: no cluster at K={K} -> empty {out}", file=sys.stderr)
        continue
    top = rows[c]["top"]
    if strict_q is not None and sp in strict_sps:
        top = [i for i in top if Q[i][c] >= strict_q]
    write(out, top)
    off = sum(1 for i in top if SP[i] != sp)
    cur = ids(f"{refs}/{sp.lower()}_rfmix.txt")
    extra = f"; {len({samples[i] for i in top} & cur)} also in {sp.lower()}_rfmix.txt" if cur else ""
    print(f"[discover] {sp}: cluster {c}, {len(top)} anchors "
          f"(Q {Q[top[0]][c]:.4f}..{Q[top[-1]][c]:.4f}){'' if not off else f', {off} with a different metadata label'}"
          f"{extra} -> {out}" if top else f"[discover] {sp}: cluster {c}, 0 anchors after strict floor -> {out}",
          file=sys.stderr)
write(oth_out, rows[oth]["top"] if oth is not None else [])
print(f"[discover] OTH: " + (f"cluster {oth}, {len(rows[oth]['top'])} anchors -> {oth_out}"
                              if oth is not None else f"no unnamed cluster with anchors -> empty {oth_out}"),
      file=sys.stderr)
left = [r["c"] for r in unnamed if r["c"] != oth]
if left:
    print(f"[discover] NOTE: {len(left)} more unnamed cluster(s) {left}. Supervised K={K} needs "
          f"K distinct labels, so a super-pop will be missing and 4a stops below. The table "
          f"above says why each cluster went unnamed.", file=sys.stderr)
print(f"[discover] cluster table -> {table_out}", file=sys.stderr)
PYEOF

# Every one of the K labels must be populated or ADMIXTURE --supervised will
# refuse to run. Fail here, with the counts, rather than after the LD-pruned
# bed has been handed to a multi-hour supervised run.
echo "[$(date +%T)] discovered anchors:"
EMPTY_LABELS=""
for SP in $SUPERPOPS; do
    F="$(anchors_for "$SP")"
    N=$(awk 'NF' "$F" 2>/dev/null | wc -l)
    printf '[%s]   %-4s %4d  %s\n' "$(date +%T)" "$SP" "$N" "$(basename "$F")"
    [[ "$N" -eq 0 ]] && EMPTY_LABELS="${EMPTY_LABELS} ${SP}"
done
N_OTH=$(awk 'NF' "$OTH_ANCHORS" 2>/dev/null | wc -l)
printf '[%s]   %-4s %4d  %s\n' "$(date +%T)" "OTH" "$N_OTH" "$(basename "$OTH_ANCHORS")"
[[ "$N_OTH" -eq 0 ]] && EMPTY_LABELS="${EMPTY_LABELS} OTH"
if [[ -n "$EMPTY_LABELS" ]]; then
    echo "ERROR: no anchors discovered for:${EMPTY_LABELS}" >&2
    echo "       ADMIXTURE --supervised needs all ${K} labels populated." >&2
    echo "       At K=${K} those labels got no cluster of their own, or none of the" >&2
    echo "       cluster cleared Q>=${ANCHOR_MIN_Q}. Why each cluster was named or not:" >&2
    echo "         ${DISCOVERY_TABLE}" >&2
    echo "       Do not lower ANCHOR_MIN_Q to force it. The unsupervised Q is kept, so" >&2
    echo "       a rerun with a different K or rule repeats only what changed." >&2
    exit 6
fi

# ---------------------------------------------------------------------------
# ANCHOR SUFFICIENCY. Written to a file as well as the log, because the anchor
# set is the thing every downstream Q depends on and "how many exemplars
# defined this ancestry" has to be answerable later without re-reading logs.
# ---------------------------------------------------------------------------
ANCHOR_SUMMARY="${ANCHOR_SUMMARY:-${HOMOG_DIR}/anchor_summary.k${K}.tsv}"
mkdir -p "$(dirname "$ANCHOR_SUMMARY")"
{
    printf 'super_pop	n_anchors	min_Q_floor	cap	below_count_floor
'
    THIN=""
    for SP in $SUPERPOPS OTH; do
        if [[ "$SP" == "OTH" ]]; then F="$OTH_ANCHORS"; else F="$(anchors_for "$SP")"; fi
        N=$(awk 'NF' "$F" 2>/dev/null | wc -l)
        # The strict floor applies only if ANCHOR_MIN_Q_STRICT is set.
        FLOOR="$ANCHOR_MIN_Q"
        if [[ -n "$ANCHOR_MIN_Q_STRICT" ]]            && [[ ",${ANCHOR_STRICT_SUPERPOPS}," == *",${SP},"* ]]; then
            FLOOR="$ANCHOR_MIN_Q_STRICT"
        fi
        BELOW=0
        if [[ "$N" -lt "$ANCHOR_FLOOR" ]]; then BELOW=1; THIN="${THIN} ${SP}(${N})"; fi
        printf '%s	%d	%s	%s	%d
' "$SP" "$N" "$FLOOR" "$ANCHOR_CAP" "$BELOW"
    done
} > "$ANCHOR_SUMMARY"
echo "[$(date +%T)] anchor summary -> ${ANCHOR_SUMMARY}"
THIN=$(awk -F'	' -v fl="$ANCHOR_FLOOR" 'NR>1 && $5==1 {printf " %s(%d)", $1, $2}' "$ANCHOR_SUMMARY")
if [[ -n "$THIN" ]]; then
    echo "[$(date +%T)] ============================================================"
    echo "[$(date +%T)] THIN ANCHOR SETS (below ANCHOR_FLOOR=${ANCHOR_FLOOR}):${THIN}"
    echo "[$(date +%T)]"
    echo "[$(date +%T)] Each of these super-pops' ancestry is defined in the"
    echo "[$(date +%T)] supervised prior by that many exemplars. ADMIXTURE will"
    echo "[$(date +%T)] run -- the labels are populated -- but every Q estimate"
    echo "[$(date +%T)] for those components rests on them, and so do the"
    echo "[$(date +%T)] homogeneity calls, the donor/reference split, the test"
    echo "[$(date +%T)] genomes' sample_type and the GWAS ancestry covariate."
    echo "[$(date +%T)]"
    echo "[$(date +%T)] Report it in the methods, or widen the candidate pool."
    echo "[$(date +%T)] Do NOT lower ANCHOR_MIN_Q to fill the set: that admits"
    echo "[$(date +%T)] admixed samples as the definition of a pure ancestry,"
    echo "[$(date +%T)] which is worse than a small set of clean ones."
    echo "[$(date +%T)] ============================================================"
    if [[ "$ANCHOR_FLOOR_FATAL" == "1" ]]; then
        echo "ERROR: ANCHOR_FLOOR_FATAL=1 and thin anchor sets present:${THIN}" >&2
        exit 7
    fi
fi

# ---------------------------------------------------------------------------
# 5. STAGE B -- build supervised .pop file (K labels) in .fam row order.
#    Every label now comes from a DISCOVERED anchor file. Previously five of
#    them were read straight from the curated *_rfmix.txt, which is what made
#    the method non-uniform.
# ---------------------------------------------------------------------------
POP_FILE="${PRUNED}.pop"
echo "[$(date +%T)] build supervised .pop file (K=${K})"
# Label order fixes the ADMIXTURE Q column order, so it must be deterministic:
# the SUPERPOPS order, then OTH.
POP_SPEC="${ADM_DIR}/pop_spec.tsv"
: > "$POP_SPEC"
for SP in $SUPERPOPS; do printf '%s\t%s\n' "$SP" "$(anchors_for "$SP")" >> "$POP_SPEC"; done
printf 'OTH\t%s\n' "$OTH_ANCHORS" >> "$POP_SPEC"
python3 - "${PRUNED}.fam" "${POP_FILE}.new" "$POP_SPEC" <<'PYEOF'
import sys
from collections import defaultdict
(fam, outp, spec_f) = sys.argv[1:4]
def load(p):
    try: return set(open(p).read().split())
    except FileNotFoundError: return set()
# Order matters: ADMIXTURE Q columns follow first-appearance of labels in .pop.
anchors = []
with open(spec_f) as f:
    for line in f:
        if not line.strip(): continue
        sp, path = line.rstrip("\n").split("\t")
        anchors.append((sp, load(path)))
# A sample carrying two labels would be assigned by list order alone, silently.
# Discovery draws each super-pop from its own pool, but the pools can overlap
# (a sample may sit in two candidate lists), so check rather than assume.
seen = defaultdict(list)
for sp, ids in anchors:
    for s in ids: seen[s].append(sp)
conflicts = {s: v for s, v in seen.items() if len(v) > 1}
if conflicts:
    print(f"ERROR: {len(conflicts)} sample(s) were discovered as anchors for "
          f"more than one super-pop; ADMIXTURE would be trained on a "
          f"contradiction.", file=sys.stderr)
    for s, v in list(conflicts.items())[:10]:
        print(f"         {s}: {', '.join(v)}", file=sys.stderr)
    sys.exit(7)
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
# Replace the .pop only if the anchors changed: its timestamp decides whether
# the multi-hour supervised run below has to be redone.
if [[ -s "$POP_FILE" ]] && cmp -s "${POP_FILE}.new" "$POP_FILE"; then
    rm -f "${POP_FILE}.new"
    echo "[$(date +%T)] .pop unchanged (same anchors as the existing supervised run)"
else
    mv "${POP_FILE}.new" "$POP_FILE"
fi

# ---------------------------------------------------------------------------
# 6. STAGE B -- supervised ADMIXTURE K.
# ---------------------------------------------------------------------------
if newer_than "${ADM_DIR}/lai_ref.${K}.Q" "$POP_FILE" "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam"; then
    rm -f "${ADM_DIR}/lai_ref.${K}.Q" "${ADM_DIR}/lai_ref.${K}.P"
    echo "[$(date +%T)] admixture --supervised K=${K} (slow)"
    cd "$ADM_DIR"
    admixture --supervised -j${THREADS} lai_ref.bed "${K}"
    cd - >/dev/null
fi
[[ -s "${ADM_DIR}/lai_ref.${K}.Q" ]] || { echo "ERROR: supervised ADMIXTURE failed"; exit 1; }

# ---------------------------------------------------------------------------
# 7. STAGE C -- parse both .Q -> lai_ref_panel_samples.k${K}.tsv + homogeneity.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Simulation-donor exclusion. Any sample used as an admix-simu DONOR must not
# also appear in the reference panel RFMix is given, or the benchmark scores
# the model on haplotypes it was handed as reference -- inflating accuracy.
# Previously this was enforced only by hand-curation of reference_ids/*, and
# only really for NAT. Enforce it here for EVERY ancestry by unioning all
# *_simulation.txt donor lists and excluding them from the homogeneous panel.
# ---------------------------------------------------------------------------
# Duplicate individuals found by 1e (the same person under two IDs across
# HGDP, 1KG, MXB, SGDP and Almarri). They are folded in with the simulation
# donors because both are "must not be in the homogeneous panel", and because
# a duplicate kept here would be counted twice by every per-super-pop statistic
# downstream as well as handing one copy to a method that was given the other.
# Absent file = 1e has not run; that is allowed, and said so, rather than
# silently treated as "no duplicates exist".
# ALMARRI 2021 METADATA. Without this every Almarri sample resolves to
# super_pop=NA and is hard-excluded by `elif sp == "NA": included = 0` -- all
# 137 of them, which is the single largest MEN contribution, prepped by 1d and
# merged by 1c only to be dropped here in silence. make_sample_groups.sh reads
# the gnomAD HGDP+1KG metadata and the MXB pop file and knows nothing about
# Almarri; neither does the SGDP manifest. So it needs its own route.
#
# Two columns, sample and pop. Build it on the cluster from the ENA metadata:
#   scripts/build-almarri-manifest.sh > reference_ids/almarri_manifest.tsv
# If the file is absent the ID PATTERN is used instead, and that is reported.
# The panel IDs are appg.meta's long_ID, longranger222_wgs_<n>_APPG<id>_GRCh38_
# gatk; the <COUNTRY><n> form (UAE19) is the pop_ID of the distributed .fam and
# matches no panel sample. Both patterns are tried. ALMARRI_META and
# ALMARRI_COUNTRIES are set near the top of the script (the MEN draw pool reads
# them before this point).

DUP_DROP="${DUP_DROP:-${REFS}/panel_duplicate_drop_ids.txt}"
# The named exclusion list, read here as well as in 1c.
#
# Not redundancy for its own sake: 1c resolved this path from BASH_SOURCE,
# which under sbatch points at Slurm's spool directory, so the list was never
# found and CHMI_CHMI3_WGS2 went into the phased panel. 1c is fixed, but a
# name on this list must not depend on one script having resolved one path
# correctly -- the cost of it slipping through is a cell line in an ADMIXTURE
# fit and possibly in a donor pool, which no downstream check would catch.
#
# Reading it here also means a panel that already contains such a sample can
# be corrected WITHOUT re-phasing: the sample stays in the phased BCF, and is
# simply never admitted to the homogeneous panel, the reference lists or the
# donor pools.
PANEL_EXCLUDE="${PANEL_EXCLUDE:-${REFS}/panel_exclude_ids.txt}"
SIM_DONORS="${SIM_DONORS:-${ADM_DIR}/simulation_donors.all.txt}"
: > "$SIM_DONORS"
if [[ -s "$DUP_DROP" ]]; then
    N_DUP=$(awk 'NF && $1 !~ /^#/' "$DUP_DROP" | wc -l)
    echo "[$(date +%T)]   duplicate drop list: $(basename "$DUP_DROP") (${N_DUP} ids)"
    awk 'NF && $1 !~ /^#/ {print $1}' "$DUP_DROP" >> "$SIM_DONORS"
elif [[ -s "$PANEL_EXCLUDE" ]]; then
    :   # handled below; the note only matters when neither list exists
else
    echo "[$(date +%T)]   NOTE: no duplicate drop list at ${DUP_DROP}"
    echo "[$(date +%T)]   1e_king_dedup has not run. The panel may contain the"
    echo "[$(date +%T)]   same person twice under two IDs, which ADMIXTURE counts"
    echo "[$(date +%T)]   twice and which can put one copy in a donor set and the"
    echo "[$(date +%T)]   other in a reference panel."
fi
for f in "${REFS}"/*_simulation.txt; do
    [[ -s "$f" ]] || continue
    echo "[$(date +%T)]   donor list: $(basename "$f") ($(awk 'NF' "$f" | wc -l) ids)"
    awk 'NF' "$f" >> "$SIM_DONORS"
done
# Named exclusions join the same drop set. They are reported separately
# because "excluded by name" and "a simulation donor" are different facts
# about a sample and the methods state them differently.
if [[ -s "$PANEL_EXCLUDE" ]]; then
    N_EXCL=$(awk 'NF && $1 !~ /^#/' "$PANEL_EXCLUDE" | wc -l)
    echo "[$(date +%T)]   named exclusions: $(basename "$PANEL_EXCLUDE") (${N_EXCL} ids)"
    while read -r _id; do
        if grep -Fxq "$_id" "${PRUNED}.fam" 2>/dev/null \
           || awk '{print $2}' "${PRUNED}.fam" 2>/dev/null | grep -Fxq "$_id"; then
            echo "[$(date +%T)]     ${_id}: PRESENT in the panel -- excluding"
        else
            echo "[$(date +%T)]     ${_id}: not in the panel (nothing to exclude)"
        fi
    done < <(awk 'NF && $1 !~ /^#/ {print $1}' "$PANEL_EXCLUDE")
    awk 'NF && $1 !~ /^#/ {print $1}' "$PANEL_EXCLUDE" >> "$SIM_DONORS"
else
    echo "[$(date +%T)]   WARNING: no named-exclusion list at ${PANEL_EXCLUDE}"
fi
sort -u -o "$SIM_DONORS" "$SIM_DONORS"
echo "[$(date +%T)] simulation donors to exclude from the homog panel: $(wc -l < "$SIM_DONORS")"

echo "[$(date +%T)] write lai_ref_panel_samples.k${K}.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv"
python3 - "${ADM_DIR}/lai_ref.fam" "${ADM_DIR}/lai_ref.${K}.Q" \
    "${UNSUP_DIR}/lai_ref.${K}.Q" "${ADM_DIR}/lai_ref.pop" \
    "$SAMPLE_GROUPS" "$GNOMAD_META" "$MXB_POPINFO" "$SGDP_MANIFEST" \
    "${REFS}/sas_rfmix.txt" "$HOMOG_THRESHOLD" "$K" "$SIM_DONORS" "$SAMPLES_TSV" \
    "$HOMOG_THRESHOLD_STRICT" "$HOMOG_STRICT_SUPERPOPS" \
    "$ALMARRI_META" "$ALMARRI_COUNTRIES" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam, qs_f, qu_f, pop_f, groups_f, gnomad_f, mxb_pop_f, sgdp_manifest_f,
 sas_f, thresh, K, donors_f, out, thresh_strict, strict_sps,
 almarri_meta_f, almarri_countries) = sys.argv[1:18]
thresh = float(thresh); K = int(K); thresh_strict = float(thresh_strict)

# Tiered homogeneity cut: the strict value for the richly-sampled super-pops,
# the floor for everyone else. See HOMOG_THRESHOLD_STRICT in the config block
# for why the asymmetry is deliberate.
STRICT = {x.strip().upper() for x in strict_sps.replace(",", " ").split() if x.strip()}
if thresh_strict < thresh:
    sys.exit(f"ERROR: HOMOG_THRESHOLD_STRICT ({thresh_strict}) is below "
             f"HOMOG_THRESHOLD ({thresh}). The strict tier must be at least as "
             f"strict as the floor, or the two tiers are the wrong way round.")
def thresh_for(sp):
    return thresh_strict if sp in STRICT else thresh
print(f"[homog] tiered cut: {sorted(STRICT)} at Q>={thresh_strict}, "
      f"all others at Q>={thresh}", file=sys.stderr)
n_by_sp = {}; n_seen_by_sp = {}

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

sim_donors = set()
try:
    sim_donors = {x for x in open(donors_f).read().split() if x}
except FileNotFoundError:
    pass
print(f"simulation donors loaded: {len(sim_donors)}", file=sys.stderr)

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
# OTH is a MODEL component, not a reference ancestry. It is the leftover
# cluster -- by construction a mixture of whatever did not fit the 7 named
# HGDP super-pops -- so it must NOT enter the HOMOGENEOUS reference panel even
# when its members load tightly on their own cluster. OTH samples still get an
# OTH anchor set (so supervised K=8 is well-posed), still appear in the TSV
# with their Q values, and still ship in the POOLED all-samples panel; they are
# simply excluded from homog_${K}pop. The homogeneous panel therefore carries
# exactly the 7 named super-pops.
PANEL_SUPERPOPS = [c for c in POP_ORDER if c != "OTH"]

def cohort_of(s):
    # Order matters: SGDP and Almarri are recognised by membership/ID form
    # BEFORE the prefix rules, because an Almarri long_ID would otherwise fall
    # through to OTHER and an LP6005/SS600 ID is an HGDP individual (sequenced
    # by SGDP, called by gnomAD) that lives in the HGDP+1KG source.
    if s in sgdp_sp: return "SGDP"
    if almarri_of(s): return "ALMARRI"
    if s.startswith("MXB"): return "MXB"
    if s.startswith(("HGDP", "LP6005", "SS600")): return "HGDP"
    if _re.match(r"^(HG|NA)\d", s): return "1KG"
    if s.startswith("CHMI"): return "CHMI"
    return "OTHER"

# ALMARRI 2021. Read the manifest when it exists; otherwise fall back to the
# ID pattern, because the alternative is silently losing all 137 samples.
import os as _os, re as _re
almarri_pop = {}
_countries = [c for c in almarri_countries.replace(",", " ").split() if c]
_almarri_re = _re.compile(r"^(" + "|".join(_countries) + r")(\d+)$", _re.I)
# The form the panel actually carries (appg.meta long_ID). Population is not in
# the ID, so a pattern hit without the manifest is labelled ALMARRI_UNKNOWN_POP.
_almarri_long_re = _re.compile(r"(^longranger\d*_wgs_|(^|_)APPG\d+(_|$))")
if almarri_meta_f and _os.path.exists(almarri_meta_f):
    with open(almarri_meta_f) as f:
        for line in f:
            ff = line.rstrip("\n").split("\t")
            if len(ff) >= 2 and ff[0] and ff[0] != "sample":
                almarri_pop[ff[0].strip()] = ff[1].strip()
    print(f"[labels] Almarri manifest: {len(almarri_pop)} samples",
          file=sys.stderr)
else:
    print(f"[labels] NOTE: no Almarri manifest at {almarri_meta_f}; resolving "
          f"Almarri samples from the ID pattern <COUNTRY><n> over "
          f"{_countries}. Build the manifest to get exact population labels: "
          f"scripts/build-almarri-manifest.sh > {almarri_meta_f}",
          file=sys.stderr)

def almarri_of(s):
    """-> population label, or None when the sample is not Almarri."""
    if s in almarri_pop:
        return almarri_pop[s]
    m = _almarri_re.match(s.strip())
    if m:
        return m.group(1).upper()
    if _almarri_long_re.search(s.strip()):
        return "ALMARRI_UNKNOWN_POP"
    return None

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
    if almarri_of(s): return "MEN"              # Almarri 2021, 8 Arab populations
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
    _a = almarri_of(s)
    if _a: return _a
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
        if s in sim_donors:
            included = 0   # donor exclusion outranks the MXB force-include
        elif co == "MXB":
            included = 1
        elif sp == "NA":
            included = 0
        elif s in sim_donors:
            # Used as an admix-simu donor -> must never be a reference
            # haplotype, for ANY ancestry, or the benchmark is circular.
            included = 0
        elif sp == "OTH" or ap_sup == "OTH":
            # Leftover cluster -> pooled panel only, never the homogeneous one.
            included = 0
        else:
            included = 1 if (max_qs >= thresh_for(sp) and ap_sup == sp) else 0
        n_incl += included
        n_by_sp[sp] = n_by_sp.get(sp, 0) + included
        n_seen_by_sp[sp] = n_seen_by_sp.get(sp, 0) + 1
        g.write("\t".join(
            [s, co, p, sp, a]
            + [f"{x:.4f}" for x in qs_ordered] + [f"{max_qs:.4f}", ap_sup]
            + [f"{x:.4f}" for x in qu_ordered] + [f"{max_qu:.4f}", ap_unsup,
                                                  str(agree), str(included)]
        ) + "\n")
    print(f"included {n_incl}/{len(samples)}", file=sys.stderr)
    print("  per super-pop (threshold applied, kept/seen):", file=sys.stderr)
    for _sp in sorted(n_seen_by_sp):
        print(f"    {_sp:4s} Q>={thresh_for(_sp):.2f}  "
              f"{n_by_sp.get(_sp,0):5d}/{n_seen_by_sp[_sp]:-5d}", file=sys.stderr)

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
  }' "$SAMPLES_TSV" | sort | uniq -c | sort -rn | awk 'NR<=15'
# `awk 'NR<=15'` rather than `head -15`: head closes the pipe at line 15, which
# SIGPIPEs the upstream `sort -rn` and, under `set -o pipefail`, fails the
# whole pipeline and kills the task. awk reads its input to the end, so the
# producer always finishes. This is a diagnostic print -- it must not be able
# to abort the stage it is reporting on.

awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_${K}pop_samples.txt"

# ---------------------------------------------------------------------------
# Composition breakdown: how many homogeneous samples each super-pop has, and
# which sub-populations (regions) supply them. This is the table to read when
# deciding which under-represented groups need more data collected -- it shows
# both the absolute count and the pass rate per population, so a super-pop that
# is small because it has few candidates is distinguishable from one that is
# small because its candidates keep failing the homogeneity cut.
# ---------------------------------------------------------------------------
COMPOSITION="${HOMOG_DIR}/homog_composition.k${K}.tsv"
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  {
    sp=$h["super_pop"]; pop=$h["pop"]; co=$h["cohort"]; inc=$h["included"]
    key=sp "\t" pop "\t" co
    tot[key]++; if (inc==1) ok[key]++
    sp_tot[sp]++; if (inc==1) sp_ok[sp]++
  }
  END {
    print "super_pop\tpopulation\tcohort\tn_included\tn_total\tpct_included"
    for (k in tot)
      printf "%s\t%d\t%d\t%.1f\n", k, ok[k]+0, tot[k], 100*(ok[k]+0)/tot[k]
  }' "$SAMPLES_TSV" | { read -r hdr; echo "$hdr"; sort -k1,1 -k4,4nr; } > "$COMPOSITION"

echo "  homogeneous composition per super_pop (sub-populations listed in $COMPOSITION):"
awk -F'\t' 'NR>1 {ok[$1]+=$4; tot[$1]+=$5; npop[$1]++}
  END { printf "    %-5s %8s %8s %7s %6s\n","SP","incl","total","pct","npops"
        for (s in ok) printf "    %-5s %8d %8d %6.1f%% %6d\n", s, ok[s], tot[s], 100*ok[s]/tot[s], npop[s] }' \
  "$COMPOSITION"
# Per-super-pop sample-ID lists, so an under-represented group can be traced
# straight back to the individuals and sub-populations supplying it.
mkdir -p "${HOMOG_DIR}/homog_${K}pop_by_superpop"
awk -F'\t' -v d="${HOMOG_DIR}/homog_${K}pop_by_superpop" '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  $h["included"]==1 { print $h["sample"] > (d "/" $h["super_pop"] "_homog_samples.txt") }
' "$SAMPLES_TSV"
echo "  per-super_pop ID lists -> ${HOMOG_DIR}/homog_${K}pop_by_superpop/"

# ---------------------------------------------------------------------------
# Publish the homogeneous lists where the 3b HOMOG forks actually look.
#
# THE BUG THIS FIXES. 3b_homog, 3b_chip_homog and the FLARE homog arm all read
# reference_ids/amr_homog_rfmix.txt. Nothing wrote it: this script put its
# output in homog_${K}pop_by_superpop/AMR_homog_samples.txt and the two names
# were never connected. So the HOMOG fork silently kept using whatever
# amr_homog_rfmix.txt was left from an earlier run -- a K=5 list of 45 ids --
# while every other arm used the new K=8 panel. The homog comparison would
# have been measuring the wrong panel, and nothing in its output would say so.
#
# Versioned copies are written first and never overwritten, so each K keeps
# its own list. The unversioned name the 3b scripts read is then updated, with
# the previous contents preserved beside it.
for _SP in $SUPERPOPS; do
    _src="${HOMOG_DIR}/homog_${K}pop_by_superpop/${_SP}_homog_samples.txt"
    [[ -s "$_src" ]] || continue
    _lc=$(echo "$_SP" | tr 'A-Z' 'a-z')
    _ver="${REFS}/${_lc}_homog_rfmix.k${K}.txt"
    cp "$_src" "$_ver"
    _n=$(awk 'NF' "$_ver" | wc -l)
    echo "  ${_SP} homogeneous -> $(basename "$_ver") (${_n})"

    _cur="${REFS}/${_lc}_homog_rfmix.txt"
    if [[ "${UPDATE_HOMOG_RFMIX:-1}" == "1" ]]; then
        if [[ -s "$_cur" ]] && ! cmp -s "$_cur" "$_ver"; then
            _bak="${_cur%.txt}.superseded.txt"
            cp "$_cur" "$_bak"
            echo "    previous $(basename "$_cur") ($(awk 'NF' "$_cur" | wc -l) ids)"
            echo "    preserved as $(basename "$_bak")"
        fi
        cp "$_ver" "$_cur"
        echo "    $(basename "$_cur") updated to the K=${K} set (${_n})"
    fi
done
echo "  NOTE: the 3b HOMOG forks read reference_ids/<sp>_homog_rfmix.txt at RUN"
echo "  time, so an already-queued stage 5 picks up these lists."

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
echo "  discovered anchors (stage A, all ${K} labels):"
for SP in $SUPERPOPS; do
    F="$(anchors_for "$SP")"
    printf '    %-4s %4d  %s\n' "$SP" "$(awk 'NF' "$F" 2>/dev/null | wc -l)" "$F"
done
printf '    %-4s %4d  %s\n' "OTH" "$N_OTH" "$OTH_ANCHORS"
echo "  curated *_rfmix.txt were NOT modified (they are the RFMix reference panels)"
echo "  homog BCFs       : ${HOMOG_SUB}/merged_chr*.homog_${K}pop.bcf"
