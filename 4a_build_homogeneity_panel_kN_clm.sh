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
# 72 h, and no --time-min: with --time-min Slurm may START the job with any
# limit down to that floor to backfill it, and an ADMIXTURE fit cut off part
# way loses all of its progress. Two fits took 8 h 23 min on 352,550 markers
# (16 Nov); without the global MAF floor there are more markers -- run
# PRUNE_ONLY=1 first, which prints the count and the expected hours.
#SBATCH --time=72:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_kN_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/4a_homog_kN_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu

# ============================================================================
# DESCRIPTION -- K=7 homogeneity panel: the 7 named super-pops.
# ============================================================================
# K-parametric companion to 4a_build_homogeneity_panel_clm.sh (K=5). Default
# K=7 (from 17 Nov): AFR AMR EAS EUR MEN OCE SAS. At K=8 the eighth cluster
# (OTH) had 2 anchors and became a second African component in the supervised
# fit (San 1.00, Mbuti 0.4-0.8, LWK 0.29), so it is dropped; K above the
# number of named super-pops brings OTH back as the leftover cluster.
#
#   0. MARKERS. Concat the 22 phased chromosomes. Site list (ADMIX_SITE_SOURCE
#      =4a): 4a measures, on 1c's PRE-phasing soft-union, which sites have
#      every super-pop <= 7% missing and every cohort overall <= 15%
#      (scripts/admix-site-list.sh), and restricts to them. LD-prune 500 50 0.1
#      on the UNRELATED fit set only (relatives and duplicate copies inflate
#      r2). No global MAF floor (ADMIX_PRUNE_MAF): every site already has
#      MAF >= 0.005 in at least one super-pop (1c's soft-union).
#      PRUNE_ONLY=1 stops here and reports the marker count and expected
#      hours; a count that would not fit in ADMIX_MAX_HOURS stops the job
#      before the fits rather than at the time limit.
#
#   FIT SET (FIT_UNRELATED=1): the panel keeps relatives, ADMIXTURE assumes
#      none. Both runs below are FITTED on an unrelated subset (1e's duplicate
#      copies out, then greedily one of each pair at KING >= FIT_KINSHIP
#      0.0884) and every relative is PROJECTED onto the fitted allele
#      frequencies with admixture -P. The merged Q has every sample in .fam
#      order, so what follows reads it exactly as before; relatives cannot be
#      anchors. Per-sample record: admixture_fit_set.k${K}.tsv.
#
#   STAGE A (unsupervised fit, then DATA-FIRST anchor discovery):
#     Unsupervised ADMIXTURE at K. For each cluster, rank every eligible
#     sample by its loading and NAME the cluster from its top ANCHOR_CAP at
#     Q >= ANCHOR_MIN_Q (plurality metadata super-pop, purity >=
#     ANCHOR_NAME_PURITY; the stronger claim wins a contested super-pop; when
#     K exceeds the named super-pops, one unnamed cluster becomes OTH). Then
#     CHOOSE the anchors (ANCHOR_SELECT):
#     central (default) = of the samples at Q >= ANCHOR_MIN_Q that carry the
#     cluster's super-pop label, the ANCHOR_CAP with the best combined rank of
#     Q (highest) and distance to the group's median in the top ANCHOR_PCS
#     principal components (nearest); top_q = the top ANCHOR_CAP by Q alone,
#     which picks each ancestry's most drifted end (16 Nov run: EUR anchors
#     all northern, AFR anchors half rainforest hunter-gatherers).
#     -> reference_ids/<sp>_anchors.k${K}.txt for the 7, and OTH;
#        every candidate's ranks in anchor_selection.k${K}.tsv.
#     Every cluster's composition, and why it was or was not named, goes to
#     discovery_clusters.k${K}.tsv. A cluster carried by one minority cohort
#     across super-pops is flagged as a likely processing batch.
#
#   STAGE B (supervised, anchor-INFORMED):
#     ADMIXTURE --supervised at K with the .pop built from those anchors.
#     Rerun only when the anchors (the .pop) or the bed changed.
#
#   STAGE C (homogeneity filter + refinement evidence):
#     included = 1 iff supervised max_Q >= HOMOG_THRESHOLD (0.94 for every
#     super-pop from 17 Nov) AND the supervised argmax super-pop matches the
#     sample's assigned super-pop. Publishes <sp>_homog_rfmix.k${K}.txt, and
#     amr_homog_rfmix.txt for the 3-way LAI HOMOG arms = July's K=5 AMR panel
#     plus the AMR samples new since July that pass here, minus donors and
#     their relatives (scripts/build-3way-amr-homog.py).
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
#   admixture_k${K}/lai_ref.recipe, lai_ref.mac_bins.tsv (pruning rule; markers by allele count)
#   admixture_k${K}/fit_pca.{eigenvec,eigenval}        (centrality for anchor choice)
#   admixture_k${K}/unsup/lai_ref.${K}.Q
#   discovery_clusters.k${K}.tsv                       (every cluster: named or why not)
#   anchor_selection.k${K}.tsv                         (every anchor candidate: Q, distance, ranks)
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

# Two DISTINCT thresholds (do not conflate), both 0.94 for EVERY super-pop
# from 17 Nov -- one number for all, between the 0.93 the under-represented
# groups had and the 0.95 the rest had:
#   ANCHOR_MIN_Q   = floor to be an anchor CANDIDATE in stage A (unsupervised
#                    cluster-Q, and the candidate must carry the cluster's
#                    super-pop label).
#   ANCHOR_CAP     = max anchors kept per super-pop. Matches
#                    build-anchors-from-koenig.py (cap 50).
#   HOMOG_THRESHOLD= floor for FINAL inclusion in stage C (supervised Q).
#                    This is the published homogeneity cut.
# K = 7 named super-pops (see DESCRIPTION for why not 8).
K="${K:-7}"
ANCHOR_MIN_Q="${ANCHOR_MIN_Q:-0.94}"
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
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.94}"

# WHICH ANCHORS: highest Q AND central to the super-pop (decided 16 Nov).
#
# The top ANCHOR_CAP by unsupervised Q are each ancestry's most DRIFTED
# members, not its typical ones, because a sample loads 0.999 on a cluster
# exactly when it is furthest from every other cluster. In the 16 Nov run that
# made the EUR anchors all northern (GBR, Orcadian, CEU, FIN, Estonian,
# Icelandic, Norwegian), the AFR anchors 24 of 50 rainforest hunter-gatherers
# (Biaka, Mbuti), the EAS anchors 45 CDX, the SAS anchors 49 ITU/STU -- and the
# supervised components, which the anchors define, centred on those ends.
#
#   ANCHOR_SELECT=central (default). Candidates = samples at Q >= ANCHOR_MIN_Q
#     on the cluster that carry the cluster's metadata super-pop. Rank them by
#     Q (highest first) and by distance to the candidates' median position in
#     the top ANCHOR_PCS principal components of the unrelated fit set (PCs
#     scaled by the square root of their eigenvalues; nearest first); keep the
#     ANCHOR_CAP with the smallest sum of the two ranks.
#     The centre is the median of the CANDIDATES, not of every sample with the
#     label. For AMR the median of every label-carrier sits among admixed
#     Latin Americans, and the "nearest" unadmixed samples would then be the
#     ones with the most European ancestry -- the opposite of an anchor.
#     ANCHOR_CENTER_MIN_Q lowers the Q a sample needs to count toward the
#     centre (default ANCHOR_MIN_Q) if a broader centre is wanted.
#   ANCHOR_SELECT=top_q: the top ANCHOR_CAP by Q (the rule before 16 Nov).
#
# PCs come from the markers' common subset (MAF >= ANCHOR_PCA_MAF over the fit
# set): PCA weights a variant by 1/sqrt(p(1-p)), so rare variants would make
# the axes track families and private haplotypes rather than populations.
# This MAF floor shapes the PCA only; it removes no ADMIXTURE marker.
ANCHOR_SELECT="${ANCHOR_SELECT:-central}"
ANCHOR_PCS="${ANCHOR_PCS:-10}"
ANCHOR_PCA_MAF="${ANCHOR_PCA_MAF:-0.05}"
ANCHOR_CENTER_MIN_Q="${ANCHOR_CENTER_MIN_Q:-${ANCHOR_MIN_Q}}"
case "$ANCHOR_SELECT" in central|top_q) ;; *)
    echo "ERROR: ANCHOR_SELECT must be central or top_q, got '${ANCHOR_SELECT}'" >&2; exit 1 ;;
esac

# Homogeneity cut: ONE number for every super-pop from 17 Nov (0.94).
#
# Until then it was tiered -- 0.95 for the richly sampled EUR/EAS/AFR (and SAS
# from 16 Nov), 0.93 for AMR/OCE/MEN so the thin ancestries kept members. 0.94
# for all is the compromise: the under-represented groups lose the few
# samples between 0.93 and 0.94, and every super-pop meets the same standard.
# The tier mechanism stays: HOMOG_STRICT_SUPERPOPS=EUR,EAS,AFR
# HOMOG_THRESHOLD_STRICT=0.95 restores a stricter cut for the named ones.
# The cut acts in stage C only, so changing it later re-runs minutes of
# work, not the ADMIXTURE fits.
HOMOG_THRESHOLD_STRICT="${HOMOG_THRESHOLD_STRICT:-${HOMOG_THRESHOLD}}"
HOMOG_STRICT_SUPERPOPS="${HOMOG_STRICT_SUPERPOPS:-}"

# SGDP candidate pools (from scripts/fetch-sgdp-anchors.py). Used as exemplars
# for cluster labeling + as the pool OCE/MEN anchors are discovered from.
# These are CANDIDATES, not anchors -- anchors are the top-N discovered here.
SGDP_MANIFEST="${SGDP_MANIFEST:-${REFS}/sgdp_manifest.tsv}"

# Almarri 2021 manifest (sample<TAB>pop). Defined HERE, not next to the labeling
# code, because the MEN draw pool below reads it and bash needs it set first.
# Keyed on the PANEL sample name, which is appg.meta's long_ID
# (longranger222_wgs_<n>_APPG<id>_GRCh38_gatk) -- not the pop_ID (UAE19) that
# the distributed .fam uses. Build it on the cluster:
#   bash scripts/build-almarri-manifest.sh <almarri_prepped.chrN.bcf> > reference_ids/almarri_manifest.tsv
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

# ADMIXTURE MARKER FILTER -- every super-pop <= 7% missing, every cohort <= 15%.
#
# THE RULE (17 Nov): at each soft-union site, every super-pop -- its cohorts
# pooled -- is at most ADMIX_SUPERPOP_MAXMISS (0.07) missing, and every cohort
# (HGDP1KG, MXB, SGDP, ALMARRI) -- its super-pops pooled -- at most
# ADMIX_COHORT_MAXMISS (0.15), measured on the PRE-phasing genotypes of 1c's
# published soft-union. ADMIX_SITE_SOURCE=4a (default): 4a computes these lists
# itself, per chromosome, with scripts/admix-site-list.sh, into
# admixture_k${K}/admix_sites/, stamped with the rule, and restricts the
# ADMIXTURE marker set to them before LD pruning.
#
# ADMIX_SITE_SOURCE=1c reads 1c's own lists instead. Those apply the 12 Nov
# rule -- each cohort x super-pop group <= 15%, panel <= 7% -- which read the
# two numbers the wrong way round. They stay on disk (1c is not re-run: the
# rule only decides which markers ADMIXTURE reads, and changing 1c re-phases).
# On chr22: 1c's rule kept 279,155 of 367,930 sites; "every super-pop <= 7%"
# alone kept 286,901 (scripts/admix-rule-compare.sh).
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
ADMIX_SITE_SOURCE="${ADMIX_SITE_SOURCE:-4a}"
ADMIX_SUPERPOP_MAXMISS="${ADMIX_SUPERPOP_MAXMISS:-0.07}"
ADMIX_COHORT_MAXMISS="${ADMIX_COHORT_MAXMISS:-0.15}"
ADMIX_MISS_MIN_N="${ADMIX_MISS_MIN_N:-5}"
ADMIX_SITE_JOBS="${ADMIX_SITE_JOBS:-8}"        # chromosomes measured at once
case "$ADMIX_SITE_SOURCE" in 4a|1c) ;; *)
    echo "ERROR: ADMIX_SITE_SOURCE must be 4a or 1c, got '${ADMIX_SITE_SOURCE}'" >&2; exit 1 ;;
esac

# PRUNE_ON_FIT=1 (default): LD is measured on the unrelated ADMIXTURE fit set
# only. Until 17 Nov the prune ran over all samples, so 638 relatives and 153
# duplicate copies of the same people entered every r2 -- a family or a person
# counted twice makes markers look more correlated than they are in the
# population. The bed is still made for every sample (relatives are projected).
PRUNE_ON_FIT="${PRUNE_ON_FIT:-1}"

# ADMIX_MAX_HOURS -- the job stops BEFORE the fits, with the count and the
# resubmit line, when the pruned marker count implies more than this many hours
# for the two ADMIXTURE fits (the job limit is 72 h). Set ADMIX_MAX_SNPS to
# thin, or raise this to accept the time.
ADMIX_MAX_HOURS="${ADMIX_MAX_HOURS:-60}"

# ADMIX_PRUNE_MAF -- an extra GLOBAL minor-allele-frequency floor before LD
# pruning. Empty (the default) = none. Every candidate site already passed
# 1c's soft-union rule -- MAF >= LAI_MAF (0.005) in at least one super-pop of
# >= MIN_SUBPOP_N samples -- and the missingness rule above. A global floor on
# top removes exactly the variants that are common in one super-pop and rare
# in the panel as a whole, i.e. the ones that tell the under-represented
# ancestries apart: on 16 Nov `--maf 0.01` removed 11,285,750 of the
# 18,880,026 sites the missingness rule kept.
# Set ADMIX_PRUNE_MAF=0.01 to reproduce the 16 Nov marker set.
#
# ADMIX_MAX_SNPS -- optional random thin (fixed seed) of the pruned markers to
# at most this many, if their number makes ADMIXTURE too slow. Empty = none.
#
# PRUNE_ONLY=1 -- build the pruned marker set, report its size, its allele
# count spectrum and the expected ADMIXTURE hours, and stop before the fits.
#
# The rule is stamped in lai_ref.recipe and a pruned set made under another
# rule is rebuilt, never reused: the timestamps alone cannot see a rule change.
ADMIX_PRUNE_MAF="${ADMIX_PRUNE_MAF:-}"
ADMIX_MAX_SNPS="${ADMIX_MAX_SNPS:-}"
PRUNE_ONLY="${PRUNE_ONLY:-0}"

# ADMIX_INPUT -- which genotypes ADMIXTURE reads.
#   phased     the phased panel. SHAPEIT5 has filled every gap, so a site a
#              cohort's file lacked carries imputed genotypes for that cohort.
#   softunion  the pre-phasing soft-union 1c publishes: real calls only, and a
#              site a cohort's file lacked stays MISSING for that cohort.
#              ADMIXTURE handles missing genotypes natively (they do not enter
#              its likelihood), so no site is dropped and none is invented --
#              the MXB-specific variants the soft-union keeps stay usable.
#              Combine with ADMIX_SITE_FILTER=0 to keep every site.
ADMIX_INPUT="${ADMIX_INPUT:-phased}"
case "$ADMIX_INPUT" in
    phased)    ADMIX_TPL="${PHASED_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.bcf" ;;
    softunion) ADMIX_TPL="${PHASED_DIR}/merged_chr%s.softunion.bcf" ;;
    *)         echo "ERROR: ADMIX_INPUT must be phased or softunion, got '${ADMIX_INPUT}'" >&2; exit 1 ;;
esac
# Where the per-chromosome site lists are: 4a's own (ADMIX_SITE_SOURCE=4a) or
# 1c's. Any other folder of lists can be named directly, e.g. one written by
# scripts/admix-rule-compare.sh: ADMIX_SITES_TPL=<out>/<rule>/merged_chr%s.admix_sites.bed
ADMIX_SITES_4A_DIR="${HOMOG_DIR}/admix_sites_4a"     # the rule does not depend on K
if [[ "$ADMIX_SITE_SOURCE" == "4a" ]]; then
    ADMIX_SITES_TPL="${ADMIX_SITES_TPL:-${ADMIX_SITES_4A_DIR}/merged_chr%s.admix_sites.bed}"
else
    ADMIX_SITES_TPL="${ADMIX_SITES_TPL:-${PHASED_DIR}/merged_chr%s.admix_sites.bed}"
fi
SUPERPOPS="EUR AFR AMR EAS SAS OCE MEN"
# OTH -- the leftover cluster -- exists only when K exceeds the named super-pops.
N_NAMED=$(wc -w <<< "$SUPERPOPS")
WANT_OTH=0; [[ "$K" -gt "$N_NAMED" ]] && WANT_OTH=1
if [[ "$K" -lt "$N_NAMED" ]]; then
    echo "ERROR: K=${K} is below the ${N_NAMED} named super-pops (${SUPERPOPS}); supervised" >&2
    echo "       ADMIXTURE needs one label per component. Use K >= ${N_NAMED}." >&2
    exit 1
fi

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

# OTH = the empirical leftover cluster when K exceeds the named super-pops (no
# metadata candidate pool). Not written, required or labelled at K=7.
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
echo "[$(date +%T)] K=${K}  ANCHOR_MIN_Q=${ANCHOR_MIN_Q}${ANCHOR_MIN_Q_STRICT:+ (strict ${ANCHOR_MIN_Q_STRICT} for ${ANCHOR_STRICT_SUPERPOPS})}  ANCHOR_CAP=${ANCHOR_CAP}  ANCHOR_FLOOR=${ANCHOR_FLOOR}  HOMOG_THRESHOLD=${HOMOG_THRESHOLD}${HOMOG_STRICT_SUPERPOPS:+ (${HOMOG_THRESHOLD_STRICT} for ${HOMOG_STRICT_SUPERPOPS})}"
echo "[$(date +%T)] ANCHOR_SELECT=${ANCHOR_SELECT}$([[ "$ANCHOR_SELECT" == central ]] && echo " (PCs ${ANCHOR_PCS}, PCA MAF ${ANCHOR_PCA_MAF}, centre Q >= ${ANCHOR_CENTER_MIN_Q})")  ADMIX_PRUNE_MAF=${ADMIX_PRUNE_MAF:-none}  ADMIX_MAX_SNPS=${ADMIX_MAX_SNPS:-none}  PRUNE_ONLY=${PRUNE_ONLY}  PRUNE_ON_FIT=${PRUNE_ON_FIT}  ADMIX_MAX_HOURS=${ADMIX_MAX_HOURS}"
echo "[$(date +%T)] ADMIXTURE sites: filter=${ADMIX_SITE_FILTER} source=${ADMIX_SITE_SOURCE}$([[ "$ADMIX_SITE_SOURCE" == 4a ]] && echo " (every super-pop <= ${ADMIX_SUPERPOP_MAXMISS} missing, every cohort <= ${ADMIX_COHORT_MAXMISS}, margins of n >= ${ADMIX_MISS_MIN_N})")  lists: ${ADMIX_SITES_TPL}"
echo "[$(date +%T)] OTH cluster: $([[ "$WANT_OTH" == 1 ]] && echo "yes (K=${K} > ${N_NAMED} named)" || echo "none (K=${K} = the ${N_NAMED} named super-pops)")"
echo "[$(date +%T)] verifying inputs"
for CHR in {1..22}; do
    F="${PHASED_DIR}/merged_chr${CHR}.shapeit5_phased.softunion_maf005.bcf"
    [[ -s "$F" && -s "${F}.csi" ]] || { echo "ERROR: missing $F (phasing not done)"; exit 1; }
done
for CHR in {1..22}; do
    F="$(printf "$ADMIX_TPL" "$CHR")"
    [[ -s "$F" ]] || { echo "ERROR: missing ADMIXTURE input $F (ADMIX_INPUT=${ADMIX_INPUT})" >&2; exit 1; }
done
if [[ "$ADMIX_SITE_FILTER" == "1" && "$ADMIX_SITE_SOURCE" == "4a" ]]; then
    for CHR in {1..22}; do
        for F in "${PHASED_DIR}/merged_chr${CHR}.softunion.bcf" "${PHASED_DIR}/merged_chr${CHR}.softunion.bcf.recipe" \
                 "${PHASED_DIR}/mxb_lifted.chr${CHR}.bcf"; do
            [[ -s "$F" ]] || {
                echo "ERROR: missing $F" >&2
                echo "       4a measures the ADMIXTURE site rule on 1c's published pre-phasing" >&2
                echo "       soft-union (PHASE_MODE=chunked publishes it). Set ADMIX_SITE_SOURCE=1c" >&2
                echo "       to use 1c's own lists (the 12 Nov rule) instead." >&2
                exit 1; }
        done
    done
    for F in "${REPO_DIR}/scripts/admix-site-list.sh" "${REPO_DIR}/scripts/admix-site-margins.py" \
             "${REPO_DIR}/scripts/build-missingness-groups.py"; do
        [[ -s "$F" ]] || { echo "ERROR: missing $F (REPO_DIR=${REPO_DIR} is not the repo?)" >&2; exit 1; }
    done
elif [[ "$ADMIX_SITE_FILTER" == "1" ]]; then
    for CHR in {1..22}; do
        F="$(printf "$ADMIX_SITES_TPL" "$CHR")"
        [[ -s "$F" ]] || {
            echo "ERROR: missing ADMIXTURE site list $F" >&2
            echo "       (ADMIX_SITE_SOURCE=${ADMIX_SITE_SOURCE}). Set ADMIX_SITE_SOURCE=4a to measure" >&2
            echo "       the lists here, or ADMIX_SITE_FILTER=0 to use every soft-union site." >&2
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
CONCAT="${ADM_DIR}/all_chr.${ADMIX_INPUT}.bcf"
PHASED_LIST=()
for CHR in {1..22}; do
    PHASED_LIST+=("$(printf "$ADMIX_TPL" "$CHR")")
done
echo "[$(date +%T)] ADMIXTURE input: ${ADMIX_INPUT} (${ADMIX_TPL})"
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
if [[ "$ADMIX_SITE_FILTER" == "1" && "$ADMIX_SITE_SOURCE" == "4a" ]]; then
    # Measure the rule per chromosome on 1c's pre-phasing soft-union. A list is
    # reused only if it is newer than its inputs AND was made under this rule
    # by this code (its .recipe stamp).
    mkdir -p "$ADMIX_SITES_4A_DIR"
    SITE_RECIPE="superpop<=${ADMIX_SUPERPOP_MAXMISS} cohort<=${ADMIX_COHORT_MAXMISS} min_n=${ADMIX_MISS_MIN_N} code=$(cat \
        "${REPO_DIR}/scripts/admix-site-list.sh" "${REPO_DIR}/scripts/admix-site-margins.py" \
        "${REPO_DIR}/scripts/build-missingness-groups.py" | md5sum | cut -c1-12)"
    SITE_TODO=()
    for CHR in {1..22}; do
        B="$(printf "$ADMIX_SITES_TPL" "$CHR")"
        if newer_than "$B" "${PHASED_DIR}/merged_chr${CHR}.softunion.bcf" "$SAMPLE_GROUPS" "$SGDP_MANIFEST" \
               ${PANEL_EXCLUDE:+"$PANEL_EXCLUDE"} \
           || [[ "$(cat "${B}.recipe" 2>/dev/null)" != "$SITE_RECIPE" ]]; then
            SITE_TODO+=("$CHR")
        fi
    done
    if [[ ${#SITE_TODO[@]} -gt 0 ]]; then
        echo "[$(date +%T)] ADMIXTURE site lists: measuring chr ${SITE_TODO[*]} (${ADMIX_SITE_JOBS} at a time): ${SITE_RECIPE}"
        site_list_one () {
            local chr="$1" b
            b="$(printf "$ADMIX_SITES_TPL" "$chr")"
            rm -f "$b" "${b}.recipe"
            bash "${REPO_DIR}/scripts/admix-site-list.sh" "$chr" "$b" > "${b}.log" 2>&1 || return 1
            echo "$SITE_RECIPE" > "${b}.recipe"
        }
        export -f site_list_one
        export ADMIX_SITES_TPL REPO_DIR SITE_RECIPE PROJECT_ROOT PHASED_DIR SAMPLE_GROUPS SGDP_MANIFEST \
               ADMIX_SUPERPOP_MAXMISS ADMIX_COHORT_MAXMISS ADMIX_MISS_MIN_N
        export EXCLUDE_IDS="$PANEL_EXCLUDE"
        if ! printf '%s\n' "${SITE_TODO[@]}" | xargs -P "$ADMIX_SITE_JOBS" -I{} bash -c 'site_list_one "$1"' _ {}; then
            echo "ERROR: the ADMIXTURE site list failed for:" >&2
            for CHR in "${SITE_TODO[@]}"; do
                B="$(printf "$ADMIX_SITES_TPL" "$CHR")"
                [[ -s "${B}.recipe" ]] && continue
                echo "------ chr${CHR}: ${B}.log" >&2; cat "${B}.log" >&2
            done
            exit 1
        fi
    else
        echo "[$(date +%T)] ADMIXTURE site lists: all 22 current under this rule -- reusing ${ADMIX_SITES_4A_DIR}"
    fi
    echo "[$(date +%T)] per chromosome: soft-union positions -> kept (every super-pop <= ${ADMIX_SUPERPOP_MAXMISS}, every cohort <= ${ADMIX_COHORT_MAXMISS})"
    for CHR in {1..22}; do
        B="$(printf "$ADMIX_SITES_TPL" "$CHR")"
        awk -v c="$CHR" '/^# records=/ { for (i = 2; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
                          printf "    chr%-3s %10d -> %10d (%5.1f%%)\n", c, v["positions"], v["positions_kept"],
                                 v["positions"] ? 100 * v["positions_kept"] / v["positions"] : 0 }' "${B}.summary.tsv"
    done
fi
if [[ "$ADMIX_SITE_FILTER" == "1" ]]; then
    SITE_LIST=()
    for CHR in {1..22}; do SITE_LIST+=("$(printf "$ADMIX_SITES_TPL" "$CHR")"); done
    # Rebuilt when a list is newer OR when ADMIX_SITES_TPL points somewhere else
    # (e.g. scripts/admix-rule-compare.sh's lists for another missingness rule,
    # which may be older than the combined file).
    if newer_than "$ADMIX_SITES_ALL" "${SITE_LIST[@]}" \
       || [[ "$(cat "${ADMIX_SITES_ALL}.recipe" 2>/dev/null)" != "tpl=${ADMIX_SITES_TPL}" ]]; then
        cat "${SITE_LIST[@]}" > "${ADMIX_SITES_ALL}.partial"
        mv "${ADMIX_SITES_ALL}.partial" "$ADMIX_SITES_ALL"
        echo "tpl=${ADMIX_SITES_TPL}" > "${ADMIX_SITES_ALL}.recipe"
    fi
    echo "[$(date +%T)] ADMIXTURE marker filter ON: $(wc -l < "$ADMIX_SITES_ALL") soft-union" \
         "sites pass the missingness rule (source ${ADMIX_SITE_SOURCE}; per-chromosome lists and" \
         "summaries: $(printf "$ADMIX_SITES_TPL" N)*)"
    EXTRACT_SITES=(--extract bed1 "$ADMIX_SITES_ALL")
    PRUNE_INPUTS+=("$ADMIX_SITES_ALL")
else
    echo "[$(date +%T)] ADMIXTURE marker filter OFF (ADMIX_SITE_FILTER=0): all soft-union sites eligible"
fi

# ---------------------------------------------------------------------------
# ADMIXTURE FIT SET. The panel keeps relatives (only duplicates are removed,
# by 1e), but ADMIXTURE assumes unrelated individuals: a family enters the
# likelihood as several draws of what is largely one genome and pulls the
# allele frequencies toward it -- and the rebuild puts ~600 HGDP+1KG relatives
# (mostly 1KG parent-child families) back in the panel. So both ADMIXTURE runs
# are FITTED on an unrelated subset and every relative is PROJECTED onto the
# fitted allele frequencies (admixture -P): relatives still get a Q, can still
# be homogeneous and still enter the panel, but do not shape the clusters.
# (Standard for ADMIXTURE and PCA: fit on unrelated, project the rest.)
# Fit set: 1e's duplicate copies out, then greedily the sample with the most
# remaining relatives at KING >= FIT_KINSHIP (0.0884, second degree -- gnomAD's
# "related") until no such pair is left among the fitted.
# FIT_UNRELATED=0 fits on every sample (the July/November behaviour).
# The same set is used for LD pruning (PRUNE_ON_FIT, step 2).
# ---------------------------------------------------------------------------
FIT_UNRELATED="${FIT_UNRELATED:-1}"
FIT_KINSHIP="${FIT_KINSHIP:-0.0884}"
KING_KIN0="${KING_KIN0:-${HOMOG_DIR}/king/king.kin0}"
if [[ "$FIT_UNRELATED" == "1" ]] && newer_than "$KING_KIN0" "${PHASED_LIST[@]}"; then
    echo "ERROR: ${KING_KIN0} is missing or older than the phased panel." >&2
    echo "       4a fits ADMIXTURE on unrelated samples chosen from 1e's KING table:" >&2
    echo "       run 1e first, or set FIT_UNRELATED=0 to fit on every sample." >&2
    exit 1
fi
# fit_set <.fam-like sample file> <fit out> <projected out> <table out>:
# writes "FID IID" rows (as given) of the fitted and the projected samples.
fit_set () {
python3 -I - "$1" "$KING_KIN0" "$FIT_KINSHIP" "$FIT_UNRELATED" "$DUP_DROP" "$2" "$3" "$4" <<'PYEOF'
import sys
from collections import defaultdict
fam_f, kin_f, thr, fit_on, dup_f, fit_out, proj_out, table_out = sys.argv[1:9]
thr = float(thr); fit_on = fit_on == "1"
fam = [l.rstrip("\n") for l in open(fam_f) if l.strip()]
iid = [l.split()[1] for l in fam]
inpanel = set(iid)
def ids(p):
    try:
        return {l.split()[0] for l in open(p) if l.strip() and not l.startswith("#")}
    except FileNotFoundError:
        return set()
removed = {}                      # sample -> (reason, related_to, kinship)
n_pairs = 0
if fit_on:
    for s in sorted(ids(dup_f) & inpanel):
        removed[s] = ("duplicate copy on 1e's drop list", "-", "-")
    adj = defaultdict(dict)
    with open(kin_f) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(head)}
        miss = [c for c in ("IID1", "IID2", "KINSHIP") if c not in ix]
        if miss:
            sys.exit(f"ERROR: {kin_f} header {head} lacks {miss}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b, k = f[ix["IID1"]], f[ix["IID2"]], float(f[ix["KINSHIP"]])
            if k >= thr and a != b and a in inpanel and b in inpanel:
                adj[a][b] = adj[b][a] = max(k, adj[a].get(b, k))
                n_pairs += 1
    def live(s):
        return [t for t in adj[s] if t not in removed]
    # Greedy: drop the sample with the most relatives still in the fit (ties:
    # the larger ID, so the result is the same on every run) until none left.
    while True:
        cand = [(len(live(s)), s) for s in adj if s not in removed and live(s)]
        if not cand:
            break
        deg, s = max(cand)
        t = max(live(s), key=lambda u: (adj[s][u], u))
        removed[s] = (f"relative of {deg} sample(s) left in the fit", t, f"{adj[s][t]:.4f}")
with open(fit_out, "w") as g, open(proj_out, "w") as p, open(table_out, "w") as tb:
    tb.write("sample\trole\treason\trelated_to\tkinship\n")
    for line, s in zip(fam, iid):
        w = line.split()
        if s in removed:
            p.write(f"{w[0]}\t{w[1]}\n")
            tb.write(f"{s}\tprojected\t" + "\t".join(removed[s]) + "\n")
        else:
            g.write(f"{w[0]}\t{w[1]}\n")
            tb.write(f"{s}\tfit\t-\t-\t-\n")
n_rel = sum(1 for r in removed.values() if r[0].startswith("relative"))
n_dup = len(removed) - n_rel
print(f"  pairs at KING >= {thr} within the panel: {n_pairs}; fitted {len(iid) - len(removed)}, "
      f"projected {len(removed)} ({n_rel} relatives, {n_dup} duplicate copies)", file=sys.stderr)
PYEOF
}

# ---------------------------------------------------------------------------
# 2. plink2 import + LD-prune on the unrelated fit set (PRUNE_ON_FIT=1), then
#    the bed for EVERY sample on the kept markers.
# ---------------------------------------------------------------------------
PRUNED="${ADM_DIR}/lai_ref"
PRUNE_KEEP_ARG=(); _keep_md5=all
if [[ "$PRUNE_ON_FIT" == "1" && "$FIT_UNRELATED" == "1" ]]; then
    PRUNE_KEEP="${ADM_DIR}/ld_prune_unrelated.iid"
    bcftools query -l "$CONCAT" | awk '{print $1, $1, 0, 0, 0, -9}' > "${ADM_DIR}/concat_samples.fam"
    echo "[$(date +%T)] LD-prune sample set: the unrelated ADMIXTURE fit set (KING >= ${FIT_KINSHIP} from ${KING_KIN0})"
    fit_set "${ADM_DIR}/concat_samples.fam" "${PRUNE_KEEP}.fit" "${PRUNE_KEEP}.proj" "${PRUNE_KEEP}.table"
    # "0 IID": the form the 3b arms --keep a --bcf import with (no FID on import).
    awk '{print 0, $2}' "${PRUNE_KEEP}.fit" > "$PRUNE_KEEP"
    PRUNE_KEEP_ARG=(--keep "$PRUNE_KEEP")
    _keep_md5=$(md5sum < "$PRUNE_KEEP" | cut -c1-12)
    echo "[$(date +%T)]   $(wc -l < "$PRUNE_KEEP") of $(wc -l < "${ADM_DIR}/concat_samples.fam") samples measure LD"
fi
PRUNE_RECIPE="ld=500,50,0.1 maf=${ADMIX_PRUNE_MAF:-none} thin=${ADMIX_MAX_SNPS:-none} sitefilter=${ADMIX_SITE_FILTER} sites=${ADMIX_SITES_TPL} input=${ADMIX_INPUT} ldsamples=${_keep_md5}"
MAF_ARG=()
[[ -n "$ADMIX_PRUNE_MAF" ]] && MAF_ARG=(--maf "$ADMIX_PRUNE_MAF")
if newer_than "${PRUNED}.bed" "${PRUNE_INPUTS[@]}" \
   || [[ "$(cat "${PRUNED}.recipe" 2>/dev/null)" != "$PRUNE_RECIPE" ]]; then
    echo "[$(date +%T)] plink2 LD-prune: ${PRUNE_RECIPE}"
    [[ -s "${PRUNED}.recipe" ]] && echo "[$(date +%T)]   (previous pruned set: $(cat "${PRUNED}.recipe"))"
    rm -f "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam" "${PRUNED}.recipe" "${ADM_DIR}/lai_ref.prune."* \
          "${PRUNED}.unthinned".{bed,bim,fam,log}
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' ${EXTRACT_SITES[@]+"${EXTRACT_SITES[@]}"} \
           ${PRUNE_KEEP_ARG[@]+"${PRUNE_KEEP_ARG[@]}"} \
           ${MAF_ARG[@]+"${MAF_ARG[@]}"} --indep-pairwise 500 50 0.1 --out "${ADM_DIR}/lai_ref.prune"
    _MADE="$PRUNED"; [[ -n "$ADMIX_MAX_SNPS" ]] && _MADE="${PRUNED}.unthinned"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' \
           --extract "${ADM_DIR}/lai_ref.prune.prune.in" \
           --max-alleles 2 --snps-only just-acgt \
           --make-bed --out "$_MADE"
    if [[ -n "$ADMIX_MAX_SNPS" ]]; then
        if [[ $(wc -l < "${_MADE}.bim") -gt "$ADMIX_MAX_SNPS" ]]; then
            echo "[$(date +%T)] thinning $(wc -l < "${_MADE}.bim") pruned SNPs to ${ADMIX_MAX_SNPS} at random (seed 20261116)"
            plink2 --bfile "$_MADE" --threads "$THREADS" --thin-count "$ADMIX_MAX_SNPS" --seed 20261116 \
                   --make-bed --out "$PRUNED"
        else
            for ext in bed bim fam; do mv "${_MADE}.${ext}" "${PRUNED}.${ext}"; done
        fi
        rm -f "${PRUNED}.unthinned".{bed,bim,fam,log}
    fi
    echo "$PRUNE_RECIPE" > "${PRUNED}.recipe"
    echo "[$(date +%T)] pruned: $(wc -l < "${PRUNED}.bim") SNPs, $(wc -l < "${PRUNED}.fam") samples"
else
    echo "[$(date +%T)] LD-pruned bed newer than its inputs and made under this rule -- reusing" \
         "${PRUNED}.bed ($(wc -l < "${PRUNED}.bim") SNPs; ${PRUNE_RECIPE})"
fi

# The pruned markers by minor allele count over every panel sample, and the
# ADMIXTURE time that number implies. With no global MAF floor the markers
# include variants common in one super-pop and rare overall -- and, in a
# super-pop as small as OCE (55 samples), single copies: MAF >= 0.005 of 110
# alleles is one allele. This table shows how many of each there are.
MAC_BINS="${PRUNED}.mac_bins.tsv"
if newer_than "$MAC_BINS" "${PRUNED}.bed" "${PRUNED}.bim"; then
    plink2 --bfile "$PRUNED" --threads "$THREADS" --freq counts --out "${PRUNED}.mac"
    awk 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
         { a = $h["ALT_CTS"]; n = $h["OBS_CT"]; m = (a < n - a) ? a : n - a; f = m / n
           b = (m <= 1) ? "1_MAC_1" : (m <= 5) ? "2_MAC_2-5" : (f < 0.005) ? "3_MAC_6+_MAF<0.005" \
             : (f < 0.01) ? "4_MAF_0.005-0.01" : (f < 0.05) ? "5_MAF_0.01-0.05" : "6_MAF>=0.05"
           c[b]++; t++ }
         END { print "bin\tSNPs\tpct"; for (k in c) printf "%s\t%d\t%.1f\n", k, c[k], 100 * c[k] / t }' \
        "${PRUNED}.mac.acount" | { read -r _h; echo "$_h"; sort; } > "${MAC_BINS}.partial"
    mv "${MAC_BINS}.partial" "$MAC_BINS"
fi
echo "[$(date +%T)] pruned markers by minor allele count over all $(wc -l < "${PRUNED}.fam") samples (${MAC_BINS}):"
sed 's/^/    /' "$MAC_BINS"
_NSNP=$(wc -l < "${PRUNED}.bim")
_HOURS=$(( (_NSNP * 503 * K / 8 / 352550 + 30) / 60 ))
echo "[$(date +%T)] expected ADMIXTURE time, both fits at K=${K}: about ${_HOURS} h" \
     "(scaled from 8 h 23 min for 352,550 markers, 3,766 fitted samples, K=8, 16 threads on 16 Nov)"
if [[ "$PRUNE_ONLY" == "1" ]]; then
    echo "[$(date +%T)] PRUNE_ONLY=1: stopping before the ADMIXTURE fits. The pruned set is kept"
    echo "[$(date +%T)] and the next run reuses it while the rule (${PRUNE_RECIPE}) is unchanged."
    exit 0
fi
if [[ "$_HOURS" -gt "$ADMIX_MAX_HOURS" ]]; then
    _FIT=$(( ADMIX_MAX_HOURS * 60 * 352550 * 8 / (503 * K) ))
    echo "ERROR: ${_NSNP} pruned markers imply about ${_HOURS} h of ADMIXTURE, over ADMIX_MAX_HOURS=${ADMIX_MAX_HOURS}" >&2
    echo "       (the job limit is 72 h). Stopping BEFORE the fits rather than at the limit." >&2
    echo "       The pruned set is kept. Resubmit with one of:" >&2
    echo "         ADMIX_MAX_SNPS=${_FIT}   random thin (seed 20261116) to what fits in ${ADMIX_MAX_HOURS} h" >&2
    echo "         ADMIX_MAX_HOURS=${_HOURS}  accept the time (a fit cut off by the limit loses all its progress;" >&2
    echo "                             a finished unsupervised fit is reused by the resubmission)" >&2
    exit 8
fi

# ---------------------------------------------------------------------------
# 2b. The fit set on the pruned bed (the same samples LD was measured on; see
#     ADMIXTURE FIT SET above). Per-sample record: admixture_fit_set.k${K}.tsv.
# ---------------------------------------------------------------------------
FIT_LIST="${ADM_DIR}/admixture_fit.keep"          # FID IID rows of lai_ref.fam
PROJ_LIST="${ADM_DIR}/admixture_projected.keep"
FIT_TABLE="${ADM_DIR}/admixture_fit_set.k${K}.tsv"
FIT_BED="${ADM_DIR}/lai_ref_fit"; PROJ_BED="${ADM_DIR}/lai_ref_proj"
echo "[$(date +%T)] ADMIXTURE fit set (FIT_UNRELATED=${FIT_UNRELATED}, KING >= ${FIT_KINSHIP} from ${KING_KIN0})"
fit_set "${PRUNED}.fam" "${FIT_LIST}.new" "${PROJ_LIST}.new" "${FIT_TABLE}.new"
# Replace only on change: these files' timestamps decide what is re-run below.
for f in "$FIT_LIST" "$PROJ_LIST" "$FIT_TABLE"; do
    if [[ -f "$f" ]] && cmp -s "${f}.new" "$f"; then rm -f "${f}.new"; else mv "${f}.new" "$f"; fi
done
N_FIT=$(awk 'NF' "$FIT_LIST" | wc -l); N_PROJ=$(awk 'NF' "$PROJ_LIST" | wc -l)
echo "[$(date +%T)] ADMIXTURE: ${N_FIT} samples fitted, ${N_PROJ} projected -> ${FIT_TABLE}"

if newer_than "${FIT_BED}.bed" "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam" "$FIT_LIST"; then
    rm -f "${FIT_BED}.bed" "${FIT_BED}.bim" "${FIT_BED}.fam"
    plink2 --bfile "$PRUNED" --keep "$FIT_LIST" --threads "$THREADS" --make-bed --out "$FIT_BED"
fi
[[ $(wc -l < "${FIT_BED}.fam") -eq "$N_FIT" ]] || { echo "ERROR: ${FIT_BED}.fam is not the ${N_FIT}-sample fit set" >&2; exit 1; }
if [[ "$N_PROJ" -gt 0 ]]; then
    if newer_than "${PROJ_BED}.bed" "${PRUNED}.bed" "${PRUNED}.bim" "${PRUNED}.fam" "$PROJ_LIST"; then
        rm -f "${PROJ_BED}.bed" "${PROJ_BED}.bim" "${PROJ_BED}.fam"
        plink2 --bfile "$PRUNED" --keep "$PROJ_LIST" --threads "$THREADS" --make-bed --out "$PROJ_BED"
    fi
    [[ $(wc -l < "${PROJ_BED}.fam") -eq "$N_PROJ" ]] || { echo "ERROR: ${PROJ_BED}.fam is not the ${N_PROJ} projected samples" >&2; exit 1; }
    # Projection needs the SAME markers in the same order as the fit.
    cmp -s "${FIT_BED}.bim" "${PROJ_BED}.bim" || { echo "ERROR: fit and projection .bim differ" >&2; exit 1; }
fi

# project <fitted .P> <dir>: the relatives' Q on fixed allele frequencies.
project () {
    local P="$1" wd="$2"
    mkdir -p "$wd"
    if [[ "$N_PROJ" -eq 0 ]]; then : > "${wd}/lai_ref_proj.${K}.Q"; return 0; fi
    if newer_than "${wd}/lai_ref_proj.${K}.Q" "$P" "${PROJ_BED}.bed" "${PROJ_BED}.bim" "${PROJ_BED}.fam"; then
        echo "[$(date +%T)] admixture -P: projecting ${N_PROJ} samples onto ${P}"
        rm -f "${wd}/lai_ref_proj.${K}.Q" "${wd}/lai_ref_proj.${K}.P"
        for ext in bed bim fam; do ln -sf "${PROJ_BED}.${ext}" "${wd}/lai_ref_proj.${ext}"; done
        cp "$P" "${wd}/lai_ref_proj.${K}.P.in"
        ( cd "$wd" && admixture -j${THREADS} -P lai_ref_proj.bed "${K}" )
    fi
    [[ -s "${wd}/lai_ref_proj.${K}.Q" ]] || { echo "ERROR: ADMIXTURE projection failed in ${wd}" >&2; exit 1; }
}
# merge_q <fit Q> <projected Q> <out Q>: one row per lai_ref.fam sample, in its
# order -- what discovery and stage C read, exactly as before.
merge_q () {
    python3 -I - "${PRUNED}.fam" "${FIT_BED}.fam" "$1" "${PROJ_BED}.fam" "$2" "$K" "${3}.new" <<'PYEOF'
import os, sys
fam, ffam, fq, pfam, pq, K, out = sys.argv[1:8]
K = int(K)
def load(fam_f, q_f):
    if not os.path.exists(q_f) or os.path.getsize(q_f) == 0:
        return {}
    ids = [l.split()[1] for l in open(fam_f) if l.strip()]
    Q = [l.rstrip("\n") for l in open(q_f) if l.strip()]
    if len(ids) != len(Q):
        sys.exit(f"ERROR: {fam_f} has {len(ids)} samples but {q_f} has {len(Q)} Q rows")
    bad = [q for q in Q if len(q.split()) != K]
    if bad:
        sys.exit(f"ERROR: {q_f} has rows without K={K} columns")
    return dict(zip(ids, Q))
qf, qp = load(ffam, fq), load(pfam, pq)
if set(qf) & set(qp):
    sys.exit(f"ERROR: {len(set(qf) & set(qp))} sample(s) both fitted and projected")
order = [l.split()[1] for l in open(fam) if l.strip()]
miss = [s for s in order if s not in qf and s not in qp]
if miss:
    sys.exit(f"ERROR: {len(miss)} sample(s) have no Q (first: {miss[:3]})")
with open(out, "w") as g:
    for s in order:
        g.write((qf[s] if s in qf else qp[s]) + "\n")
print(f"  merged Q: {len(qf)} fitted + {len(qp)} projected = {len(order)} rows", file=sys.stderr)
PYEOF
    if [[ -f "$3" ]] && cmp -s "${3}.new" "$3"; then rm -f "${3}.new"; else mv "${3}.new" "$3"; fi
}

# ---------------------------------------------------------------------------
# 3. STAGE A -- UNSUPERVISED ADMIXTURE K (runs FIRST; discovery depends on it):
#    fitted on the unrelated set, relatives projected, merged in .fam order.
# ---------------------------------------------------------------------------
UFIT="${UNSUP_DIR}/fit"; mkdir -p "$UFIT"
if newer_than "${UFIT}/lai_ref_fit.${K}.Q" "${FIT_BED}.bed" "${FIT_BED}.bim" "${FIT_BED}.fam"; then
    echo "[$(date +%T)] admixture (unsupervised) K=${K} on the ${N_FIT} fitted samples (slow)"
    rm -f "${UFIT}/lai_ref_fit.${K}.Q" "${UFIT}/lai_ref_fit.${K}.P"
    for ext in bed bim fam; do ln -sf "${FIT_BED}.${ext}" "${UFIT}/lai_ref_fit.${ext}"; done
    ( cd "$UFIT" && admixture -j${THREADS} lai_ref_fit.bed "${K}" )
else
    echo "[$(date +%T)] unsupervised fit newer than its bed -- reusing ${UFIT}/lai_ref_fit.${K}.Q"
fi
[[ -s "${UFIT}/lai_ref_fit.${K}.Q" ]] || { echo "ERROR: unsupervised ADMIXTURE failed"; exit 1; }
project "${UFIT}/lai_ref_fit.${K}.P" "${UNSUP_DIR}/proj"
merge_q "${UFIT}/lai_ref_fit.${K}.Q" "${UNSUP_DIR}/proj/lai_ref_proj.${K}.Q" "${UNSUP_DIR}/lai_ref.${K}.Q"
[[ -s "${UNSUP_DIR}/lai_ref.${K}.Q" ]] || { echo "ERROR: unsupervised ADMIXTURE failed"; exit 1; }

# A projected sample cannot be an anchor: anchors label the supervised FIT.
ANCHOR_INELIG="${ADM_DIR}/anchor_ineligible.txt"
{ [[ -s "$DUP_DROP" ]] && awk 'NF && $1 !~ /^#/ {print $1}' "$DUP_DROP"
  awk 'NF {print $2}' "$PROJ_LIST"; } | sort -u > "$ANCHOR_INELIG"

# ---------------------------------------------------------------------------
# 3b. Principal components of the fitted (unrelated) samples, for
#     ANCHOR_SELECT=central. Anchors are always fitted samples, so every
#     candidate has PCs. See WHICH ANCHORS at the top.
# ---------------------------------------------------------------------------
PCA_PREFIX="${ADM_DIR}/fit_pca"
PCA_RECIPE="pcs=${ANCHOR_PCS} maf=${ANCHOR_PCA_MAF}"
if [[ "$ANCHOR_SELECT" == "central" ]]; then
    if newer_than "${PCA_PREFIX}.eigenvec" "${FIT_BED}.bed" "${FIT_BED}.bim" "${FIT_BED}.fam" \
       || [[ "$(cat "${PCA_PREFIX}.recipe" 2>/dev/null)" != "$PCA_RECIPE" ]]; then
        echo "[$(date +%T)] PCA of the ${N_FIT} fitted samples (${PCA_RECIPE}) for the central anchor choice"
        rm -f "${PCA_PREFIX}.eigenvec" "${PCA_PREFIX}.eigenval" "${PCA_PREFIX}.recipe"
        plink2 --bfile "$FIT_BED" --threads "$THREADS" --maf "$ANCHOR_PCA_MAF" \
               --pca "$ANCHOR_PCS" --out "$PCA_PREFIX"
        echo "$PCA_RECIPE" > "${PCA_PREFIX}.recipe"
    else
        echo "[$(date +%T)] PCA newer than the fit set and made under this rule -- reusing ${PCA_PREFIX}.eigenvec"
    fi
    [[ -s "${PCA_PREFIX}.eigenvec" && -s "${PCA_PREFIX}.eigenval" ]] \
        || { echo "ERROR: PCA of the fit set failed (${PCA_PREFIX}.eigenvec/.eigenval missing)" >&2; exit 1; }
fi

# ---------------------------------------------------------------------------
# 4. STAGE A -- data-first anchor discovery from the unsupervised Q.
#    Per cluster: rank every eligible sample, name the cluster from its top
#    ANCHOR_CAP at Q >= ANCHOR_MIN_Q, then choose the anchors (ANCHOR_SELECT).
#    See DATA-FIRST ANCHOR DISCOVERY and WHICH ANCHORS at the top.
# ---------------------------------------------------------------------------
DISCOVERY_TABLE="${HOMOG_DIR}/discovery_clusters.k${K}.tsv"
ANCHOR_TABLE="${HOMOG_DIR}/anchor_selection.k${K}.tsv"
echo "[$(date +%T)] data-first discovery from unsupervised K=${K}: per cluster, named from its top ${ANCHOR_CAP} at Q >= ${ANCHOR_MIN_Q} (purity >= ${ANCHOR_NAME_PURITY}); anchors chosen ${ANCHOR_SELECT}"
python3 -I - "${PRUNED}.fam" "${UNSUP_DIR}/lai_ref.${K}.Q" "$K" "$ANCHOR_MIN_Q" "$ANCHOR_CAP" \
    "${ANCHOR_MIN_Q_STRICT:-}" "$ANCHOR_STRICT_SUPERPOPS" "$ANCHOR_NAME_PURITY" "$ANCHOR_COHORT_WARN" \
    "$SAMPLE_GROUPS" "$MXB_POPINFO" "$SGDP_MANIFEST" "$ALMARRI_META" "$PANEL_EXCLUDE" "$ANCHOR_INELIG" \
    "$REFS" "$SUPERPOPS" "$OTH_ANCHORS" "$DISCOVERY_TABLE" \
    "$ANCHOR_SELECT" "${PCA_PREFIX}.eigenvec" "${PCA_PREFIX}.eigenval" "$ANCHOR_CENTER_MIN_Q" \
    "$GNOMAD_META" "$ANCHOR_TABLE" "$WANT_OTH" <<'PYEOF'
import sys, re, csv, math
from collections import Counter

(fam_f, qu_f, K, thresh, cap, strict_q, strict_sps, purity_min, cohort_warn,
 groups_f, mxb_f, sgdp_f, almarri_f, excl_f, dup_f, refs, superpops,
 oth_out, table_out, select, evec_f, eval_f, center_min_q, gnomad_f,
 anchor_table_out, want_oth) = sys.argv[1:27]
want_oth = want_oth == "1"
K = int(K); thresh = float(thresh); cap = int(cap)
strict_q = float(strict_q) if strict_q else None
strict_sps = {x for x in strict_sps.replace(",", " ").split() if x}
purity_min = float(purity_min); cohort_warn = float(cohort_warn)
center_min_q = float(center_min_q)
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

# Population labels, ONLY to report who the anchors are. Same precedence as
# pop_of() in stage C: MXB, then SGDP, then Almarri, then gnomAD.
def tsv_map(path, key, val):
    out = {}
    try:
        with open(path) as fh:
            h = next(fh).rstrip("\n").split("\t")
            if key in h and val in h:
                ki, vi = h.index(key), h.index(val)
                for l in fh:
                    f = l.rstrip("\n").split("\t")
                    if len(f) > max(ki, vi): out[f[ki]] = f[vi]
    except (FileNotFoundError, StopIteration):
        pass
    return out
popl = tsv_map(gnomad_f, "project_meta.sample_id", "hgdp_tgp_meta.Population")
try:
    with open(almarri_f) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[0] != "sample": popl[f[0]] = f[1]
except FileNotFoundError:
    pass
popl.update(tsv_map(sgdp_f, "sample", "population"))
popl.update(tsv_map(mxb_f, "Sample_ID", "inferred_genetic_cluster"))

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
POPL = [popl.get(s, "NA") for s in samples]

# ---- principal components for ANCHOR_SELECT=central -------------------------
# Each PC scaled by the square root of its eigenvalue, so Euclidean distance
# weighs the axes by the variation they carry.
PC = {}
if select == "central":
    ev = [float(l.split()[0]) for l in open(eval_f) if l.strip()]
    with open(evec_f) as fh:
        h = fh.readline().lstrip("#").split()
        if "IID" not in h:
            sys.exit(f"ERROR: {evec_f} header {h} has no IID column")
        ii = h.index("IID")
        pcc = [j for j, col in enumerate(h) if re.fullmatch(r"PC\d+", col)]
        if not pcc or len(pcc) > len(ev):
            sys.exit(f"ERROR: {evec_f} has {len(pcc)} PC columns but {eval_f} has {len(ev)} eigenvalues")
        wt = [math.sqrt(max(ev[k], 0.0)) for k in range(len(pcc))]
        for l in fh:
            f = l.split()
            if f:
                PC[f[ii]] = [float(f[j]) * wt[k] for k, j in enumerate(pcc)]
    print(f"[discover] central anchor choice: {len(PC)} fitted samples in {len(pcc)} PCs "
          f"(eigenvalues {', '.join(f'{x:.1f}' for x in ev[:len(pcc)])}); centre = median of the "
          f"candidates at Q >= {center_min_q}", file=sys.stderr)

def avg_ranks(vals):
    """1 = smallest value; tied values share their average rank."""
    order = sorted(range(len(vals)), key=lambda j: vals[j])
    r = [0.0] * len(vals)
    j = 0
    while j < len(order):
        k = j
        while k + 1 < len(order) and vals[order[k + 1]] == vals[order[j]]:
            k += 1
        for m in range(j, k + 1):
            r[order[m]] = (j + k) / 2 + 1
        j = k + 1
    return r

def median(xs):
    xs = sorted(xs); n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

SEL = []   # every candidate considered: the anchor_selection table

def choose(c, name, label, pool):
    """The anchors for cluster c from pool (indices at the Q floor, Q-descending).
    label: the metadata super-pop a candidate must carry, or None (OTH)."""
    if select == "top_q":
        chosen = set(pool[:cap])
        rq = avg_ranks([-Q[i][c] for i in pool])
        for j, i in enumerate(pool):
            SEL.append([c, name, i, Q[i][c], None, rq[j], None, rq[j], int(i in chosen)])
        return pool[:cap], ""
    cand = [i for i in pool if label is None or SP[i] == label]
    nopc = [samples[i] for i in cand if samples[i] not in PC]
    if nopc:
        sys.exit(f"ERROR: {len(nopc)} {name} anchor candidate(s) have no PCs (first: {nopc[:3]}); "
                 f"candidates are fitted samples and every fitted sample is in {evec_f}")
    centre = [i for i in elig if Q[i][c] >= center_min_q and (label is None or SP[i] == label)
              and samples[i] in PC]
    if not cand or not centre:
        return [], "no candidate carries the label"
    D = len(PC[samples[cand[0]]])
    mid = [median([PC[samples[i]][d] for i in centre]) for d in range(D)]
    dist = {i: math.sqrt(sum((PC[samples[i]][d] - mid[d]) ** 2 for d in range(D))) for i in cand}
    rq = avg_ranks([-Q[i][c] for i in cand])
    rd = avg_ranks([dist[i] for i in cand])
    score = {i: rq[j] + rd[j] for j, i in enumerate(cand)}
    chosen = sorted(cand, key=lambda i: (score[i], -Q[i][c], dist[i], samples[i]))[:cap]
    ch = set(chosen)
    for j, i in enumerate(cand):
        SEL.append([c, name, i, Q[i][c], dist[i], rq[j], rd[j], score[i], int(i in ch)])
    dropped = len(pool) - len(cand)
    note = (f"central: {len(cand)} candidates at the Q floor carrying {label or 'any'} label"
            f"{f' ({dropped} with another label left out)' if dropped else ''}; centre = median of "
            f"{len(centre)} at Q >= {center_min_q}; chosen distance "
            f"{min(dist[i] for i in chosen):.4f}-{max(dist[i] for i in chosen):.4f} "
            f"(candidates {min(dist.values()):.4f}-{max(dist.values()):.4f})")
    return chosen, note
print(f"[discover] {len(samples)} samples; {len(samples) - len(elig)} not eligible "
      f"(named exclusions, 1e duplicate drops, relatives projected rather than fitted); "
      f"{len(elig)} ranked in every cluster",
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
if unnamed and want_oth:       # only when K exceeds the named super-pops
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

def report(name, c, top, note, out, extra=""):
    if not top:
        print(f"[discover] {name}: cluster {c}, 0 anchors ({note or 'none at the Q floor'}) -> {out}",
              file=sys.stderr)
        return
    qs = [Q[i][c] for i in top]
    off = sum(1 for i in top if name != "OTH" and SP[i] != name)
    print(f"[discover] {name}: cluster {c}, {len(top)} anchors (Q {min(qs):.4f}-{max(qs):.4f})"
          f"{'' if not off else f', {off} with a different metadata label'}{extra} -> {out}",
          file=sys.stderr)
    if note:
        print(f"[discover]     {note}", file=sys.stderr)
    print(f"[discover]     by population: "
          + ", ".join(f"{p} {n}" for p, n in Counter(POPL[i] for i in top).most_common(10)),
          file=sys.stderr)

for sp in SUPERPOPS:
    out = f"{refs}/{sp.lower()}_anchors.k{K}.txt"
    c = taken.get(sp)
    if c is None:
        write(out, [])
        print(f"[discover] {sp}: no cluster at K={K} -> empty {out}", file=sys.stderr)
        continue
    if select == "top_q":
        # The rule before 16 Nov, unchanged: the top ANCHOR_CAP, then the strict floor.
        top, note = choose(c, sp, sp, rows[c]["passing"])
        if strict_q is not None and sp in strict_sps:
            top = [i for i in top if Q[i][c] >= strict_q]
            kept = set(top)
            for row in SEL:
                if row[0] == c:
                    row[8] = int(row[2] in kept)
    else:
        pool = rows[c]["passing"]
        if strict_q is not None and sp in strict_sps:
            pool = [i for i in pool if Q[i][c] >= strict_q]
        top, note = choose(c, sp, sp, pool)
    write(out, top)
    cur = ids(f"{refs}/{sp.lower()}_rfmix.txt")
    report(sp, c, top, note, out,
           f"; {len({samples[i] for i in top} & cur)} also in {sp.lower()}_rfmix.txt" if cur else "")
if oth is not None:
    oth_top, note = choose(oth, "OTH", None, rows[oth]["passing"])
    write(oth_out, oth_top)
    report("OTH", oth, oth_top, note, oth_out)
elif want_oth:
    write(oth_out, [])
    print(f"[discover] OTH: no unnamed cluster with anchors -> empty {oth_out}", file=sys.stderr)
else:
    print(f"[discover] OTH: none at K={K} (one cluster per named super-pop)", file=sys.stderr)
left = [r["c"] for r in unnamed if r["c"] != oth]
if left:
    print(f"[discover] NOTE: {len(left)} unnamed cluster(s) {left}. Supervised K={K} needs "
          f"K distinct labels, so a super-pop will be missing and 4a stops below. The table "
          f"above says why each cluster went unnamed (e.g. two clusters claiming one super-pop).",
          file=sys.stderr)
with open(anchor_table_out, "w") as g:
    g.write("cluster\tname\tsample\tcohort\tpop\tmeta_super_pop\tQ\tdist_to_centre\t"
            "rank_Q\trank_dist\trank_sum\tselected\n")
    for c, nm, i, q, d, rq, rd, sc, s in SEL:
        g.write(f"{c}\t{nm}\t{samples[i]}\t{CO[i]}\t{POPL[i]}\t{SP[i]}\t{q:.6f}\t"
                f"{'NA' if d is None else f'{d:.5f}'}\t{rq:g}\t{'NA' if rd is None else f'{rd:g}'}\t"
                f"{sc:g}\t{s}\n")
print(f"[discover] cluster table -> {table_out}", file=sys.stderr)
print(f"[discover] every anchor candidate, its Q, distance and ranks -> {anchor_table_out}", file=sys.stderr)
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
N_OTH=0
if [[ "$WANT_OTH" == "1" ]]; then
    N_OTH=$(awk 'NF' "$OTH_ANCHORS" 2>/dev/null | wc -l)
    printf '[%s]   %-4s %4d  %s\n' "$(date +%T)" "OTH" "$N_OTH" "$(basename "$OTH_ANCHORS")"
    [[ "$N_OTH" -eq 0 ]] && EMPTY_LABELS="${EMPTY_LABELS} OTH"
fi
LABELS_SP="$SUPERPOPS"; [[ "$WANT_OTH" == "1" ]] && LABELS_SP="$SUPERPOPS OTH"
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
    for SP in $LABELS_SP; do
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
# the SUPERPOPS order, then OTH (only when K exceeds the named super-pops).
POP_SPEC="${ADM_DIR}/pop_spec.tsv"
: > "$POP_SPEC"
for SP in $SUPERPOPS; do printf '%s\t%s\n' "$SP" "$(anchors_for "$SP")" >> "$POP_SPEC"; done
[[ "$WANT_OTH" == "1" ]] && printf 'OTH\t%s\n' "$OTH_ANCHORS" >> "$POP_SPEC"
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
# 6. STAGE B -- supervised ADMIXTURE K, fitted on the unrelated set (its .pop
#    rows taken from the full .pop in fit order; every anchor is in the fit,
#    since discovery never picks a projected sample), relatives projected onto
#    the supervised allele frequencies, merged in .fam order.
# ---------------------------------------------------------------------------
SFIT="${ADM_DIR}/sup_fit"; mkdir -p "$SFIT"
python3 -I - "${PRUNED}.fam" "$POP_FILE" "${FIT_BED}.fam" "${SFIT}/lai_ref_fit.pop.new" <<'PYEOF'
import sys
fam, pop, ffam, out = sys.argv[1:5]
ids = [l.split()[1] for l in open(fam) if l.strip()]
lab = [l.strip() for l in open(pop)]
if len(ids) != len(lab):
    sys.exit(f"ERROR: {fam} ({len(ids)}) and {pop} ({len(lab)}) differ in length")
m = dict(zip(ids, lab))
fit = [l.split()[1] for l in open(ffam) if l.strip()]
with open(out, "w") as g:
    for s in fit:
        g.write(m[s] + "\n")
n_all = sum(1 for x in lab if x != "-"); n_fit = sum(1 for s in fit if m[s] != "-")
if n_fit != n_all:
    sys.exit(f"ERROR: {n_all - n_fit} anchor(s) are not in the ADMIXTURE fit set")
print(f"  supervised fit .pop: {len(fit)} rows, {n_fit} anchors", file=sys.stderr)
PYEOF
if [[ -s "${SFIT}/lai_ref_fit.pop" ]] && cmp -s "${SFIT}/lai_ref_fit.pop.new" "${SFIT}/lai_ref_fit.pop"; then
    rm -f "${SFIT}/lai_ref_fit.pop.new"
else
    mv "${SFIT}/lai_ref_fit.pop.new" "${SFIT}/lai_ref_fit.pop"
fi
for ext in bed bim fam; do ln -sf "${FIT_BED}.${ext}" "${SFIT}/lai_ref_fit.${ext}"; done
if newer_than "${SFIT}/lai_ref_fit.${K}.Q" "${SFIT}/lai_ref_fit.pop" "${FIT_BED}.bed" "${FIT_BED}.bim" "${FIT_BED}.fam"; then
    rm -f "${SFIT}/lai_ref_fit.${K}.Q" "${SFIT}/lai_ref_fit.${K}.P"
    echo "[$(date +%T)] admixture --supervised K=${K} on the ${N_FIT} fitted samples (slow)"
    ( cd "$SFIT" && admixture --supervised -j${THREADS} lai_ref_fit.bed "${K}" )
fi
[[ -s "${SFIT}/lai_ref_fit.${K}.Q" ]] || { echo "ERROR: supervised ADMIXTURE failed"; exit 1; }
project "${SFIT}/lai_ref_fit.${K}.P" "${ADM_DIR}/sup_proj"
merge_q "${SFIT}/lai_ref_fit.${K}.Q" "${ADM_DIR}/sup_proj/lai_ref_proj.${K}.Q" "${ADM_DIR}/lai_ref.${K}.Q"
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
#   scripts/build-almarri-manifest.sh <almarri_prepped.chrN.bcf> > reference_ids/almarri_manifest.tsv
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
# Donors' RELATIVES. The panel keeps relatives (only duplicates are removed,
# by 1e), but a relative pair is never split across donor and reference: a
# reference holding a donor's first-degree relative holds long stretches of
# that donor's haplotypes, which inflates accuracy exactly as the donor itself
# would. Same rule and threshold as split-homog-donors-reference.py
# (--role-kinship 0.177): a donor's relatives go where the donor is -- out of
# the homogeneous panel, which the HOMOG arms use as reference. ID checks
# (above, and 2b's) cannot see this.
ROLE_KINSHIP="${ROLE_KINSHIP:-0.177}"
KING_KIN0="${KING_KIN0:-${HOMOG_DIR}/king/king.kin0}"
DONOR_REL="${ADM_DIR}/donor_relatives.k${K}.tsv"
if [[ -s "$KING_KIN0" ]]; then
    python3 -I - "$KING_KIN0" "$ROLE_KINSHIP" "$DONOR_REL" "$REFS" <<'PYEOF'
import glob, os, sys
kin_f, thr, out_f, refs = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
donors = {}
for f in sorted(glob.glob(os.path.join(refs, "*_simulation.txt"))):
    for line in open(f):
        w = line.split()
        if w:
            donors.setdefault(w[0], []).append(os.path.basename(f)[:-len(".txt")])
with open(kin_f) as fh:
    head = fh.readline().lstrip("#").rstrip("\n").split("\t")
    idx = {n: i for i, n in enumerate(head)}
    miss = [c for c in ("IID1", "IID2", "KINSHIP") if c not in idx]
    if miss:
        sys.exit(f"ERROR: {kin_f} header {head} lacks {miss}")
    pairs = []
    for line in fh:
        f = line.rstrip("\n").split("\t")
        a, b, k = f[idx["IID1"]], f[idx["IID2"]], float(f[idx["KINSHIP"]])
        if k < thr or (a in donors) == (b in donors):
            continue
        d, r = (a, b) if a in donors else (b, a)
        pairs.append((r, d, k))
with open(out_f, "w") as g:
    g.write("relative\tdonor\tdonor_lists\tkinship\n")
    for r, d, k in sorted(pairs, key=lambda x: -x[2]):
        g.write(f"{r}\t{d}\t{','.join(donors[d])}\t{k:.4f}\n")
print(f"  donors' relatives at KING >= {thr}: {len({p[0] for p in pairs})} sample(s), "
      f"{len(pairs)} pair(s) -> {out_f}")
PYEOF
    awk -F'\t' 'NR > 1 {print $1}' "$DONOR_REL" >> "$SIM_DONORS"
    if [[ $(wc -l < "$DONOR_REL") -gt 1 ]]; then sed 's/^/    /' "$DONOR_REL"; fi
else
    echo "[$(date +%T)]   WARNING: no ${KING_KIN0} (1e not run): donors' relatives are NOT kept"
    echo "[$(date +%T)]   out of the homogeneous panel; 2b's kinship check cannot run either."
fi
sort -u -o "$SIM_DONORS" "$SIM_DONORS"
echo "[$(date +%T)] simulation donors (and their relatives) to exclude from the homog panel: $(wc -l < "$SIM_DONORS")"

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
          f"scripts/build-almarri-manifest.sh <almarri_prepped.chrN.bcf> > {almarri_meta_f}",
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

awk -F'\t' 'NR>1 && $NF==1 {print $1}' "$SAMPLES_TSV" > "${HOMOG_DIR}/homog_${K}pop_samples.txt.new"
# Replaced only on change: its timestamp decides whether the 22 homog BCFs are rebuilt.
if [[ -f "${HOMOG_DIR}/homog_${K}pop_samples.txt" ]] && cmp -s "${HOMOG_DIR}/homog_${K}pop_samples.txt.new" "${HOMOG_DIR}/homog_${K}pop_samples.txt"; then
    rm -f "${HOMOG_DIR}/homog_${K}pop_samples.txt.new"
else
    mv "${HOMOG_DIR}/homog_${K}pop_samples.txt.new" "${HOMOG_DIR}/homog_${K}pop_samples.txt"
fi

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
#
# AMR IS DIFFERENT (17 Nov): the 3b HOMOG arms are the 3-way NAT/EUR/AFR LAI
# pipeline, which keeps July's K=5 homogeneity calls and adds the AMR samples
# that are new since July and homogeneous now (scripts/build-3way-amr-homog.py,
# minus every donor, donor's relative, duplicate and exclusion). So
# amr_homog_rfmix.txt gets that list; amr_homog_rfmix.k${K}.txt keeps this
# run's own. HOMOG_3WAY_BASE= (empty) publishes this run's AMR list instead.
HOMOG_3WAY_BASE="${HOMOG_3WAY_BASE-${REFS}/lai_ref_panel_samples.tsv}"
for _SP in $SUPERPOPS; do
    _src="${HOMOG_DIR}/homog_${K}pop_by_superpop/${_SP}_homog_samples.txt"
    [[ -s "$_src" ]] || continue
    _lc=$(echo "$_SP" | tr 'A-Z' 'a-z')
    _ver="${REFS}/${_lc}_homog_rfmix.k${K}.txt"
    cp "$_src" "$_ver"
    _n=$(awk 'NF' "$_ver" | wc -l)
    echo "  ${_SP} homogeneous -> $(basename "$_ver") (${_n})"
    _pub="$_ver"; _what="the K=${K} set"
    if [[ "$_SP" == "AMR" && -n "$HOMOG_3WAY_BASE" ]]; then
        [[ -s "$HOMOG_3WAY_BASE" ]] || { echo "ERROR: HOMOG_3WAY_BASE ${HOMOG_3WAY_BASE} missing" >&2; exit 1; }
        python3 -I "${REPO_DIR}/scripts/build-3way-amr-homog.py" --base "$HOMOG_3WAY_BASE" \
            --table "$SAMPLES_TSV" --exclude "$SIM_DONORS" --out "${REFS}/amr_homog_rfmix.3way.txt" \
            | sed 's/^/    /'
        _pub="${REFS}/amr_homog_rfmix.3way.txt"; _what="the 3-way set (July K=5 + new since July)"
        _n=$(awk 'NF' "$_pub" | wc -l)
    fi

    _cur="${REFS}/${_lc}_homog_rfmix.txt"
    if [[ "${UPDATE_HOMOG_RFMIX:-1}" == "1" ]]; then
        if [[ -s "$_cur" ]] && ! cmp -s "$_cur" "$_pub"; then
            _bak="${_cur%.txt}.superseded.txt"
            cp "$_cur" "$_bak"
            echo "    previous $(basename "$_cur") ($(awk 'NF' "$_cur" | wc -l) ids)"
            echo "    preserved as $(basename "$_bak")"
        fi
        cp "$_pub" "$_cur"
        echo "    $(basename "$_cur") updated to ${_what} (${_n})"
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
    # Rebuilt when the sample list changed (it is only rewritten on change).
    if [[ -s "$DST" && -s "${DST}.csi" && "$DST" -nt "$HOMOG_KEEP" && "$DST" -nt "$SRC" ]]; then continue; fi
    bcftools view "$SRC" -S "$HOMOG_KEEP" --force-samples --threads "$THREADS" -Ob -o "${DST}.partial"
    bcftools index "${DST}.partial"
    mv "${DST}.partial.csi" "${DST}.csi"; mv "${DST}.partial" "$DST"
done

echo "[$(date +%T)] DONE."
echo "  samples TSV      : $SAMPLES_TSV"
echo "  discovered anchors (stage A, all ${K} labels):"
for SP in $SUPERPOPS; do
    F="$(anchors_for "$SP")"
    printf '    %-4s %4d  %s\n' "$SP" "$(awk 'NF' "$F" 2>/dev/null | wc -l)" "$F"
done
[[ "$WANT_OTH" == "1" ]] && printf '    %-4s %4d  %s\n' "OTH" "$N_OTH" "$OTH_ANCHORS"
echo "  curated *_rfmix.txt were NOT modified (they are the RFMix reference panels)"
echo "  homog BCFs       : ${HOMOG_SUB}/merged_chr*.homog_${K}pop.bcf"
