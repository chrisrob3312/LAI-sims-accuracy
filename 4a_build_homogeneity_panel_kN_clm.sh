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
# inference prob to pick anchors from. This script discovers anchors instead,
# and discovers them for ALL SEVEN super-pops by the same rule:
#
#   STAGE A (unsupervised, anchor DISCOVERY):
#     Run UNSUPERVISED ADMIXTURE at K on the merged panel (which now includes
#     the SGDP candidates from 1b/1c). Match clusters to super-pop NAMES by
#     MEAN Q-loading of labelled exemplars, one-to-one. The 7 named groups
#     claim 7 clusters; the LEFTOVER cluster becomes OTH -- defined by the
#     data, not by a hand-written population list.
#     Anchors = that super-pop's DRAW POOL sorted by its cluster-Q descending,
#     kept if >= ANCHOR_MIN_Q, capped at ANCHOR_CAP
#     -> reference_ids/<sp>_anchors.k${K}.txt, for all 7 plus OTH.
#     OTH has no metadata candidate pool, so its anchors are drawn from the
#     whole panel; this also makes supervised K=8 well-posed, since ADMIXTURE
#     --supervised requires exactly K distinct labels in the .pop file.
#
#     WHY UNIFORM. Previously EUR/AFR/AMR/EAS/SAS were handed the curated
#     *_rfmix.txt lists directly as supervised anchors while only OCE/MEN/OTH
#     were discovered. Two methods produced the one .pop prior, and the SGDP
#     genomes downloaded for the other five super-pops could never be drawn
#     for them -- so widening the panel did nothing for five of seven groups.
#
#   STAGE B (supervised, anchor-INFORMED):
#     Rerun ADMIXTURE --supervised at K with the full K-way .pop prior, built
#     entirely from stage A's discovered anchors.
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
# the original anchor sets, applied UNIFORMLY to all 7 super-pop candidate
# pools rather than only to the SGDP additions:
#
#   STAGE A (unsupervised, anchor DISCOVERY):
#     Run UNSUPERVISED ADMIXTURE at K on the merged panel (which now includes
#     the SGDP candidates from 1b/1c). Map each cluster to a super-pop by
#     MEAN-Q with a one-to-one assignment, using each super-pop's LABEL pool.
#     Then for EVERY one of the 7, take its DRAW pool members whose Q in that
#     super-pop's own cluster is >= ANCHOR_MIN_Q, SORTED by that Q descending,
#     capped at ANCHOR_CAP -> reference_ids/<sp>_anchors.k${K}.txt.
#     => "order by the unsupervised sorted order and identify the anchors
#         from the unsupervised" (Christina's spec), now the same way for all
#         seven pools.
#
#   STAGE B (supervised, anchor-INFORMED):
#     Rerun ADMIXTURE --supervised at K with the full K-way .pop prior, built
#     entirely from the stage-A discovered anchors. Supervised Q is the truth
#     for the homogeneity filter.
#
#   STAGE C (homogeneity filter):
#     included = 1 iff supervised max_Q >= HOMOG_THRESHOLD AND the supervised
#     argmax super-pop matches the sample's assigned super-pop. MXB forced in.
#
# Every candidate pool remains strictly a CANDIDATE POOL; this script is what
# decides which members clear Q >= HOMOG_THRESHOLD and enter the panel.
#
# Outputs in $HOMOG_DIR (k${K} subdirs, so the K=5 outputs are never clobbered):
#   admixture_k${K}/lai_ref.{bed,bim,fam,${K}.Q,${K}.P,pop}
#   admixture_k${K}/unsup/lai_ref.${K}.Q
#   admixture_k${K}/discover_spec.tsv                  (pools used, per super-pop)
#   reference_ids/<sp>_anchors.k${K}.txt               (DISCOVERED anchors, 8)
#   lai_ref_panel_samples.k${K}.tsv
#   homog_${K}pop_samples.txt
#   homog_${K}pop/merged_chr${CHR}.homog_${K}pop.bcf{,.csi}
#
# NOT written by this script: reference_ids/*_rfmix.txt. Those are the curated
# RFMix reference panels (read by all four 3b_*.sh) and this script only ever
# reads them, as seeds and as label pools.
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

# UNIFORM DISCOVERY ACROSS ALL 7 CANDIDATE POOLS.
#
# Every named super-pop now goes through the SAME path: assemble a candidate
# pool, rank its members by their loading on that super-pop's unsupervised
# cluster, keep the top ANCHOR_CAP at or above ANCHOR_MIN_Q, and use THOSE as
# the supervised anchors. Previously only OCE, MEN and OTH were discovered;
# EUR/AFR/AMR/EAS/SAS were handed the curated *_rfmix.txt lists directly, so
# two different methods produced the 8 supervised labels and the SGDP samples
# downloaded for the other five super-pops could never become anchors for them.
#
# A candidate pool does TWO jobs, and they are deliberately kept separate:
#
#   DRAW pool  -- who may become an anchor. Union of every *_<token>_candidates
#                 .txt for that super-pop PLUS the curated seed list, so SGDP
#                 (and later Almarri) samples can be drawn for ALL seven.
#   LABEL pool -- who names the unsupervised clusters by mean-Q. The curated
#                 seed list where one exists, because those are the vetted
#                 exemplars and cluster naming is the step everything else
#                 depends on; the draw pool for OCE/MEN, which have no curated
#                 list to use instead.
#
# Membership in either confers nothing. Which samples become ANCHORS is decided
# by the top-N >= ANCHOR_MIN_Q cut in stage A, and which enter the HOMOGENEOUS
# panel by the supervised Q >= HOMOG_THRESHOLD cut in stage C. Adding a source
# here is never the same as anchoring on it.
#
# Drop in e.g. almarri_men_candidates.txt and it is picked up with no code
# change. The token differs from the super-pop code only for OCE (the SGDP
# fetcher writes *_oceania_candidates.txt).
SUPERPOPS="EUR AFR AMR EAS SAS OCE MEN"

cand_token () { [[ "$1" == "OCE" ]] && echo oceania || echo "$1" | tr 'A-Z' 'a-z'; }

# Per-super-pop glob override: EUR_CAND_GLOB, AFR_CAND_GLOB, ... all honoured.
cand_glob_for () {
    local sp="$1" var="${1}_CAND_GLOB"
    if [[ -n "${!var:-}" ]]; then echo "${!var}"
    else echo "${REFS}/*_$(cand_token "$sp")_candidates.txt"; fi
}

# Curated seed/label lists. These are the RFMix REFERENCE PANELS, read by all
# four 3b_*.sh scripts, and are NEVER written by this script -- see the naming
# note on the discovery outputs below. OCE and MEN have none.
seed_for () {
    local sp="$1" f="${REFS}/$(echo "$sp" | tr 'A-Z' 'a-z')_rfmix.txt"
    [[ -s "$f" ]] && echo "$f" || echo ""
}

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

# Build the per-super-pop candidate unions from the globs, reporting every
# source that contributes so the provenance is in the log rather than implied.
build_cand_union () {
    local sp="$1" glob="$2" seed="$3" out="$4" n=0 tot=0
    : > "$out"
    if [[ -n "$seed" && -s "$seed" ]]; then
        local c; c=$(awk 'NF' "$seed" | wc -l)
        echo "[$(date +%T)]   ${sp} draw pool += $(basename "$seed") (${c}, curated seed)"
        awk 'NF {print $1}' "$seed" >> "$out"
        n=$((n+1)); tot=$((tot+c))
    fi
    for f in $glob; do
        [[ -s "$f" ]] || continue
        local c; c=$(awk 'NF' "$f" | wc -l)
        echo "[$(date +%T)]   ${sp} draw pool += $(basename "$f") (${c})"
        awk 'NF {print $1}' "$f" >> "$out"
        n=$((n+1)); tot=$((tot+c))
    done
    # AMR carries the Mexican Biobank genomes, which are part of the AMR
    # exemplar set everywhere else in this pipeline and must be drawable here.
    if [[ "$sp" == "AMR" && -s "$MXB_POPINFO" ]]; then
        local c; c=$(awk 'NR>1 && NF' "$MXB_POPINFO" | wc -l)
        echo "[$(date +%T)]   ${sp} draw pool += $(basename "$MXB_POPINFO") (${c}, MXB)"
        awk -F'\t' 'NR>1 && NF {print $1}' "$MXB_POPINFO" >> "$out"
        n=$((n+1)); tot=$((tot+c))
    fi
    sort -u -o "$out" "$out"
    local u; u=$(wc -l < "$out")
    if [[ "$n" -eq 0 ]]; then
        echo "ERROR: no candidate files matched for ${sp}: $glob" >&2
        echo "       expected at least one, e.g. ${REFS}/sgdp_$(cand_token "$sp")_candidates.txt" >&2
        echo "       or the curated seed ${REFS}/$(echo "$sp" | tr 'A-Z' 'a-z')_rfmix.txt" >&2
        return 1
    fi
    echo "[$(date +%T)]   ${sp}: ${n} source(s), ${tot} ids -> ${u} unique"
    return 0
}

# The label pool names the clusters; the draw pool is who may be drawn. Where a
# curated seed exists the label pool is just that seed, so cluster naming keeps
# using only vetted exemplars even though the draw pool is now wider.
echo "[$(date +%T)] assembling candidate pools (draw + label) for: ${SUPERPOPS}"
DISCOVER_SPEC="${ADM_DIR}/discover_spec.tsv"
: > "$DISCOVER_SPEC"
for SP in $SUPERPOPS; do
    DRAW="${ADM_DIR}/$(echo "$SP" | tr 'A-Z' 'a-z')_candidates.union.txt"
    build_cand_union "$SP" "$(cand_glob_for "$SP")" "$(seed_for "$SP")" "$DRAW" || exit 1
    SEED="$(seed_for "$SP")"
    if [[ -n "$SEED" ]]; then
        LABEL="${ADM_DIR}/$(echo "$SP" | tr 'A-Z' 'a-z')_label.txt"
        awk 'NF {print $1}' "$SEED" > "$LABEL"
        # MXB belongs to the AMR exemplar set for cluster naming too.
        [[ "$SP" == "AMR" && -s "$MXB_POPINFO" ]] && \
            awk -F'\t' 'NR>1 && NF {print $1}' "$MXB_POPINFO" >> "$LABEL"
        sort -u -o "$LABEL" "$LABEL"
        echo "[$(date +%T)]   ${SP} label pool: $(wc -l < "$LABEL") (curated)"
    else
        LABEL="$DRAW"
        echo "[$(date +%T)]   ${SP} label pool: $(wc -l < "$LABEL") (= draw pool; no curated list)"
    fi
    printf '%s\t%s\t%s\t%s\n' "$SP" "$DRAW" "$LABEL" "$(anchors_for "$SP")" >> "$DISCOVER_SPEC"
done

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
# 4. STAGE A -- discover anchors for ALL 7 super-pops from the unsupervised Q.
#    Map K clusters -> super-pop by mean-Q (label pools); then for EACH of the
#    7 take its DRAW pool members with cluster-Q >= ANCHOR_MIN_Q, SORTED
#    descending, CAPPED at ANCHOR_CAP. This is the "top 25-50 >= 0.93-0.95"
#    recipe that built the original anchor sets, now applied uniformly instead
#    of to OCE/MEN alone.
# ---------------------------------------------------------------------------
echo "[$(date +%T)] discover anchors for all 7 super-pops from unsupervised K=${K} (top ${ANCHOR_CAP} >= ${ANCHOR_MIN_Q})"
python3 - "${PRUNED}.fam" "${UNSUP_DIR}/lai_ref.${K}.Q" "$K" "$ANCHOR_MIN_Q" "$ANCHOR_CAP" \
    "$DISCOVER_SPEC" "$OTH_ANCHORS" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam_f, qu_f, K, thresh, cap, spec_f, oth_out) = sys.argv[1:8]
K = int(K); thresh = float(thresh); cap = int(cap)

def load(p):
    try: return set(open(p).read().split())
    except FileNotFoundError: return set()

# spec: super_pop <TAB> draw_pool <TAB> label_pool <TAB> anchors_out
SPEC = []
with open(spec_f) as f:
    for line in f:
        if not line.strip(): continue
        sp, draw, label, out = line.rstrip("\n").split("\t")
        SPEC.append((sp, draw, label, out))
SPS = [s[0] for s in SPEC]

# Exemplar -> super-pop, from the LABEL pools only. A sample that appears in
# two super-pops' label pools is ambiguous and is dropped from labelling
# entirely rather than silently going to whichever was read last.
owners = defaultdict(set)
for sp, _draw, label_f, _out in SPEC:
    for s in load(label_f): owners[s].add(sp)
exemplars = {}
dupes = 0
for s, sps in owners.items():
    if len(sps) == 1: exemplars[s] = next(iter(sps))
    else: dupes += 1
if dupes:
    print(f"[discover] NOTE: {dupes} sample(s) appear in more than one label "
          f"pool and were excluded from cluster naming.", file=sys.stderr)

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

def discover(cand_f, label_f, sp, out):
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
          f"draw-pool-in-panel pass Q>={thresh}; kept top {len(picked)} "
          f"(cap {cap}) -> {out}", file=sys.stderr)
    if picked:
        print(f"           top: {picked[0][0]} Q={picked[0][1]:.4f} | "
              f"min kept Q={picked[-1][1]:.4f}", file=sys.stderr)
    # How the discovered set compares with the curated seed it replaces. This
    # is the number to read before trusting a rerun: if discovery keeps almost
    # none of the curated anchors for a super-pop that used to work, the
    # cluster map or the threshold is wrong, not the curated list.
    curated = load(label_f)
    if curated:
        kept_ids = {s for s, _ in picked}
        overlap = len(kept_ids & curated)
        print(f"           vs curated label pool ({len(curated)}): "
              f"{overlap} retained, {len(kept_ids) - overlap} newly drawn",
              file=sys.stderr)
    return len(picked)

counts = {}
for sp, draw_f, label_f, out_f in SPEC:
    counts[sp] = discover(draw_f, label_f, sp, out_f)
empty = [sp for sp, n in counts.items() if n == 0]
if empty:
    print(f"[discover] NOTE: no anchors for {', '.join(empty)}. Zero means that "
          f"super-pop didn't form a clean cluster at K={K} or no draw-pool "
          f"member cleared Q>={thresh}. ADMIXTURE --supervised needs every "
          f"label populated, so this will fail downstream.", file=sys.stderr)

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
    echo "       Either lower ANCHOR_MIN_Q (currently ${ANCHOR_MIN_Q}), widen" >&2
    echo "       that super-pop's candidate pool, or drop K." >&2
    echo "       The cluster map and per-super-pop counts are above." >&2
    exit 6
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
python3 - "${PRUNED}.fam" "$POP_FILE" "$POP_SPEC" <<'PYEOF'
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
DUP_DROP="${DUP_DROP:-${REFS}/panel_duplicate_drop_ids.txt}"
SIM_DONORS="${SIM_DONORS:-${ADM_DIR}/simulation_donors.all.txt}"
: > "$SIM_DONORS"
if [[ -s "$DUP_DROP" ]]; then
    N_DUP=$(awk 'NF && $1 !~ /^#/' "$DUP_DROP" | wc -l)
    echo "[$(date +%T)]   duplicate drop list: $(basename "$DUP_DROP") (${N_DUP} ids)"
    awk 'NF && $1 !~ /^#/ {print $1}' "$DUP_DROP" >> "$SIM_DONORS"
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
sort -u -o "$SIM_DONORS" "$SIM_DONORS"
echo "[$(date +%T)] simulation donors to exclude from the homog panel: $(wc -l < "$SIM_DONORS")"

echo "[$(date +%T)] write lai_ref_panel_samples.k${K}.tsv"
SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv"
python3 - "${ADM_DIR}/lai_ref.fam" "${ADM_DIR}/lai_ref.${K}.Q" \
    "${UNSUP_DIR}/lai_ref.${K}.Q" "${ADM_DIR}/lai_ref.pop" \
    "$SAMPLE_GROUPS" "$GNOMAD_META" "$MXB_POPINFO" "$SGDP_MANIFEST" \
    "${REFS}/sas_rfmix.txt" "$HOMOG_THRESHOLD" "$K" "$SIM_DONORS" "$SAMPLES_TSV" \
    "$HOMOG_THRESHOLD_STRICT" "$HOMOG_STRICT_SUPERPOPS" <<'PYEOF'
import sys
from collections import defaultdict, Counter
(fam, qs_f, qu_f, pop_f, groups_f, gnomad_f, mxb_pop_f, sgdp_manifest_f,
 sas_f, thresh, K, donors_f, out, thresh_strict, strict_sps) = sys.argv[1:16]
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
