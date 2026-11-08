#!/usr/bin/env bash
# ============================================================================
# 6a_transfer_to_mab_clm.sh
#
# Rsync the four quadrants of the reference panel to mab.dldcc.bcm.edu:
#
#   /mount/genepi2/Joint_called_hgdp_1kg_plus_mxbiobank_wgs_ref_panel/
#     hg38/
#       all_samples/   # full QC'd phased panel (hg38)  + samples TSV
#       homog_5pop/    # 5-superpop homog subset (hg38) + samples TSV
#     hg19/
#       all_samples/   # full panel lifted to hg19      + samples TSV
#       homog_5pop/    # homog subset lifted to hg19    + samples TSV
#
# Builds the destination directory tree with `ssh ... mkdir -p` (rsync's
# --mkpath isn't available on all rsync versions; the explicit mkdir is
# more portable).
#
# Authentication: edit MAB_USER below to your username on mab. The script
# uses ssh/rsync's interactive password prompt -- no credentials are
# embedded in the script. Set up ssh-agent before running if you want to
# avoid re-typing the password for each of the four rsync calls:
#
#   eval "$(ssh-agent)"
#   ssh-add -t 4h   # cache cred for 4 hours
#
# Or use key-based auth (recommended):
#
#   ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519_mab
#   ssh-copy-id -i ~/.ssh/id_ed25519_mab.pub <MAB_USER>@mab.dldcc.bcm.edu
#
# Then set: export MAB_SSH_OPTS="-i ~/.ssh/id_ed25519_mab"
#
# This script does NOT run as a SLURM job; run it interactively from the
# head node (mhgcp-h00) where you have ssh outbound:
#
#   ./6a_transfer_to_mab_clm.sh           # transfer all four quadrants
#   ./6a_transfer_to_mab_clm.sh hg38      # only hg38 (both subdirs)
#   ./6a_transfer_to_mab_clm.sh hg19      # only hg19
#   ./6a_transfer_to_mab_clm.sh hg38/homog_5pop  # only one quadrant
# ============================================================================

# >>> EDIT THESE TWO LINES <<<
MAB_USER="${MAB_USER:-<USERNAME>}"
MAB_HOST="${MAB_HOST:-mab.dldcc.bcm.edu}"
# <<< ------------------------- >>>

MAB_BASE="${MAB_BASE:-/mount/genepi2/Joint_called_hgdp_1kg_plus_mxbiobank_wgs_ref_panel}"

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
HG19_DIR="${HG19_DIR:-${PROJECT_ROOT}/05_panel_hg19}"

# Repo reference_ids (anchor lists, SGDP candidate lists, manifests). 6a had
# no repo awareness at all -- it only ever shipped BCFs + two sample files.
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"

# K for the multi-K homogeneity panel (4a_build_homogeneity_panel_kN_clm.sh).
# The K=5 quadrants are kept for back-compat; the K-panel quadrants are new.
K="${K:-8}"

# The full panel is HGDP+1KG+MXB+SGDP once the publication rerun has run, so
# it ships under a name that says so rather than the bare "all_samples".
# Override to "all_samples" to keep the legacy remote path.
ALL_SAMPLES_SUBDIR="${ALL_SAMPLES_SUBDIR:-1KG-HGDP_MXB_SGDP_all_samples_combined}"
# The homogeneous panel carries exactly the 7 NAMED super-pops -- OTH (the
# leftover K=8 cluster) is deliberately excluded by 4a, since a mixture of
# leftovers is not a reference ancestry.
HOMOG_K_SUBDIR="${HOMOG_K_SUBDIR:-homog_${K}pop_7superpop}"

SAMPLES_TSV="${HOMOG_DIR}/lai_ref_panel_samples.tsv"
SAMPLES_TSV_K="${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv"
MAB_SSH_OPTS="${MAB_SSH_OPTS:-}"

set -euo pipefail

# The scopes that move nothing must not demand credentials: the whole point of
# "verify" is to check the split locally before deciding to ship, and asking
# for a username before printing the usage text is just unhelpful.
case "${1:-all}" in
    verify|help|--help|-h) NEEDS_CRED=0 ;;
    *)                     NEEDS_CRED=1 ;;
esac
if [[ "$MAB_USER" == "<USERNAME>" && "$NEEDS_CRED" == "1" ]]; then
    echo "ERROR: edit MAB_USER at the top of $0 (or export MAB_USER=...) before running."
    echo "       To check the two panels without transferring: $0 verify"
    echo "       For the list of scopes:                       $0 help"
    exit 1
fi

# Build local manifests (one per quadrant) with file sizes and md5s. The
# manifests get rsync'd alongside the data so you can verify on the mab side.
build_manifest () {
    local src="$1" manifest="$2"
    ( cd "$src" && \
      find . -maxdepth 2 -type f \( -name '*.bcf' -o -name '*.bcf.csi' \
                                    -o -name '*.tsv' -o -name '*.txt' \
                                    -o -name '*.vcf.gz' -o -name '*.vcf.gz.csi' \
                                    -o -name '*.vcf.gz.tbi' -o -name '*.md' \) \
        -exec stat -c '%s %n' {} + | sort
    ) > "$manifest"
}

# rsync one source dir -> one remote dir. Includes index files.
do_rsync () {
    local src="$1" remote_subdir="$2"
    local remote_full="${MAB_BASE}/${remote_subdir}"
    if [[ ! -d "$src" ]]; then
        echo "WARN: source $src missing, skipping $remote_subdir"
        return 0
    fi
    echo "==================================================================="
    echo "[$(date +%T)] $src  ->  ${MAB_HOST}:${remote_full}"
    echo "==================================================================="
    # Build local manifest first.
    local manifest="${src}/.manifest.txt"
    build_manifest "$src" "$manifest"

    # Ensure remote dir exists.
    ssh ${MAB_SSH_OPTS} "${MAB_USER}@${MAB_HOST}" "mkdir -p '${remote_full}'"

    # rsync. Include relevant data files only; exclude tmp / scratch / logs.
    # Excludes must precede the '*/' catch-all include -- rsync's first-match
    # wins, and '*/' would otherwise pull in tmp/ before the exclude rule fires.
    # -L / --copy-links: dereference symlinks and send the target file. The
    # hg38 handlers stage symlinks in .stage_*/ dirs, and without -L rsync
    # preserves those symlinks (pointing at source-cluster paths that don't
    # exist on mab), producing 188K of dead links instead of real BCFs.
    rsync -avzL --partial --progress \
        -e "ssh ${MAB_SSH_OPTS}" \
        --exclude='tmp/' --exclude='.stage_*' \
        --exclude='*.log' --exclude='*.tmp' \
        --exclude='*.rejected.vcf.gz' \
        --include='*/' \
        --include='*.bcf' --include='*.bcf.csi' \
        --include='*.vcf.gz' --include='*.vcf.gz.csi' --include='*.vcf.gz.tbi' \
        --include='*.tsv' --include='*.txt' \
        --include='.manifest.txt' \
        --include='PROVENANCE.txt' \
        `# the simulated cohort's own formats: .hanc is the truth tracts and` \
        `# .bp the recombination realisation that regenerates a unit, so both` \
        `# must travel -- without .bp the cohort cannot be reproduced.` \
        --include='*.hanc' --include='*.hanc2' --include='*.bp' \
        --include='*.phgeno' --include='*.haps' --include='*.sample' \
        --include='*.snp' --include='*.dat' \
        --exclude='*' \
        "$src/" \
        "${MAB_USER}@${MAB_HOST}:${remote_full}/"
}

# What to transfer. CLI arg filters; default = all four quadrants.
SCOPE="${1:-all}"

run_hg38_all () {
    # Stage a flat dir with the per-chr full-panel BCFs + samples TSV.
    local STAGE="${PROJECT_ROOT}/.stage_hg38_all_samples"
    mkdir -p "$STAGE"
    for f in "${PHASED_DIR}"/merged_chr*.shapeit5_phased.softunion_maf005.bcf{,.csi}; do
        ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    [[ -s "$SAMPLES_TSV" ]] && ln -sf "$SAMPLES_TSV" "$STAGE/lai_ref_panel_samples.tsv"
    do_rsync "$STAGE" "hg38/${ALL_SAMPLES_SUBDIR}"
}

run_hg38_homog () {
    local STAGE="${PROJECT_ROOT}/.stage_hg38_homog_5pop"
    mkdir -p "$STAGE"
    for f in "${HOMOG_DIR}/homog_5pop"/*.bcf{,.csi}; do
        [[ -e "$f" ]] || continue
        ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    [[ -s "$SAMPLES_TSV" ]] && ln -sf "$SAMPLES_TSV" "$STAGE/lai_ref_panel_samples.tsv"
    [[ -s "${HOMOG_DIR}/homog_5pop_samples.txt" ]] && \
        ln -sf "${HOMOG_DIR}/homog_5pop_samples.txt" "$STAGE/homog_5pop_samples.txt"
    do_rsync "$STAGE" "hg38/homog_5pop"
}

run_reference_lists () {
    # Everything the mab-side genotyping / ancestry-inference pipeline needs to
    # run SUPERVISED ADMIXTURE and FLARE against this panel, which previously
    # was never transferred at all:
    #   * <sp>_rfmix.txt      -- the supervised ADMIXTURE anchor sets (the 5
    #                            existing + OCE/MEN/OTH discovered by 4a)
    #   * homog_<K>pop_samples.txt and the per-super-pop ID lists
    #   * homog_composition.k<K>.tsv  -- per super_pop x population counts
    #   * population_cluster_meanQ.k<K>.tsv -- the refinement evidence table
    #   * sgdp_*_candidates.txt, sgdp_manifest.tsv -- SGDP provenance
    local STAGE="${PROJECT_ROOT}/.stage_reference_lists"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    for f in "${REFS}"/*_rfmix.txt "${REFS}"/sgdp_*_candidates.txt \
             "${REFS}/sgdp_manifest.tsv" "${REFS}/MXB50genomes_popinfo.tsv"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    for f in "$SAMPLES_TSV" "$SAMPLES_TSV_K" \
             "${HOMOG_DIR}/homog_${K}pop_samples.txt" \
             "${HOMOG_DIR}/homog_5pop_samples.txt" \
             "${HOMOG_DIR}/homog_composition.k${K}.tsv" \
             "${HOMOG_DIR}/population_cluster_meanQ.k${K}.tsv"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    if [[ -d "${HOMOG_DIR}/homog_${K}pop_by_superpop" ]]; then
        for f in "${HOMOG_DIR}/homog_${K}pop_by_superpop"/*.txt; do
            [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
        done
    fi
    do_rsync "$STAGE" "reference_lists"
}

run_hg38_homog_k () {
    local STAGE="${PROJECT_ROOT}/.stage_hg38_${HOMOG_K_SUBDIR}"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    local n=0
    for f in "${HOMOG_DIR}/homog_${K}pop"/*.bcf{,.csi}; do
        [[ -e "$f" ]] || continue
        ln -sf "$f" "$STAGE/$(basename "$f")"; n=$((n+1))
    done
    if [[ "$n" -eq 0 ]]; then
        echo "WARNING: no BCFs in ${HOMOG_DIR}/homog_${K}pop -- run 4a_build_homogeneity_panel_kN_clm.sh first; skipping"
        return 0
    fi
    [[ -s "$SAMPLES_TSV_K" ]] && ln -sf "$SAMPLES_TSV_K" "$STAGE/$(basename "$SAMPLES_TSV_K")"
    [[ -s "${HOMOG_DIR}/homog_${K}pop_samples.txt" ]] && \
        ln -sf "${HOMOG_DIR}/homog_${K}pop_samples.txt" "$STAGE/homog_${K}pop_samples.txt"
    [[ -s "${HOMOG_DIR}/homog_composition.k${K}.tsv" ]] && \
        ln -sf "${HOMOG_DIR}/homog_composition.k${K}.tsv" "$STAGE/homog_composition.k${K}.tsv"
    do_rsync "$STAGE" "hg38/${HOMOG_K_SUBDIR}"
}

run_hg19_homog_k () {
    if [[ ! -d "${HG19_DIR}/homog_${K}pop" ]]; then
        echo "WARNING: ${HG19_DIR}/homog_${K}pop missing -- run 5a to lift the K=${K} panel; skipping"
        return 0
    fi
    do_rsync "${HG19_DIR}/homog_${K}pop" "hg19/${HOMOG_K_SUBDIR}"
    for f in "$SAMPLES_TSV_K" "${HOMOG_DIR}/homog_${K}pop_samples.txt"; do
        [[ -s "$f" ]] && rsync -avz -e "ssh ${MAB_SSH_OPTS}" "$f" \
            "${MAB_USER}@${MAB_HOST}:${MAB_BASE}/hg19/${HOMOG_K_SUBDIR}/$(basename "$f")"
    done
}

run_hg19_all () {
    do_rsync "${HG19_DIR}/all_samples" "hg19/${ALL_SAMPLES_SUBDIR}"
    # Add the samples TSV alongside.
    if [[ -s "$SAMPLES_TSV" ]]; then
        rsync -avz -e "ssh ${MAB_SSH_OPTS}" "$SAMPLES_TSV" \
            "${MAB_USER}@${MAB_HOST}:${MAB_BASE}/hg19/${ALL_SAMPLES_SUBDIR}/lai_ref_panel_samples.tsv"
    fi
}

run_hg19_homog () {
    do_rsync "${HG19_DIR}/homog_5pop" "hg19/homog_5pop"
    if [[ -s "$SAMPLES_TSV" ]]; then
        rsync -avz -e "ssh ${MAB_SSH_OPTS}" "$SAMPLES_TSV" \
            "${MAB_USER}@${MAB_HOST}:${MAB_BASE}/hg19/homog_5pop/lai_ref_panel_samples.tsv"
    fi
    if [[ -s "${HOMOG_DIR}/homog_5pop_samples.txt" ]]; then
        rsync -avz -e "ssh ${MAB_SSH_OPTS}" "${HOMOG_DIR}/homog_5pop_samples.txt" \
            "${MAB_USER}@${MAB_HOST}:${MAB_BASE}/hg19/homog_5pop/homog_5pop_samples.txt"
    fi
}

# ===========================================================================
# THREE TREES: TWO REFERENCE PANELS AND THE COHORT ITSELF
#
# A reference panel and a cohort are different KINDS of thing -- one is what a
# method is given, the other is what it is tested on -- so they get separate
# trees and names that say which is which:
#
#   lai_reference_panel/        ADMIXTURE + LAI REFERENCE. Every homogeneous
#                               sample, donors INCLUDED. What the GWAS and
#                               FLARE run against. Nothing is simulated from
#                               these haplotypes, so there is no circularity
#                               to avoid and withholding donors would only
#                               shrink it.
#
# The two benchmark trees sit under one parent,
# ${HARMONIZATION_SUBDIR}/ (genotype_harmonization_pipeline_testing by
# default), so the harmonization project cannot be confused with the GWAS
# one next to it -- their donor policies are opposite:
#
#   <parent>/benchmark_reference_panel/
#                               ADMIXTURE + LAI REFERENCE. 1000 Genomes plus
#                               a population-stratified share of the non-1KG
#                               samples, donors EXCLUDED.
#
#   <parent>/benchmark_cohort/  THE COHORT. The test genomes themselves --
#                               real/ (real individuals drawn from HGDP, SGDP,
#                               MXB, Almarri) and simulated/ (the 2,160 with
#                               their truth tracts) -- plus the manifests that
#                               say who is in which chip x build cell.
#
# NEITHER PANEL IS AN IMPUTATION REFERENCE, and mistaking one for an
# imputation panel is the error this note exists to prevent. In the
# genotyping_pipeline flow, imputation runs against the EXTERNAL servers --
# TOPMed, All of Us and Michigan (1000 Genomes). The two trees here are used
# only for supervised ADMIXTURE (the global-ancestry Q) and for local-ancestry
# inference (FLARE / RFMix). They differ from each other only in donor policy.
#
# That split is what makes the 1KG rules look inconsistent until you see which
# role each is about:
#
#   * 1000 Genomes IS in benchmark_reference_panel/, and that is fine. This
#     panel is never used to impute, so 1KG here cannot leak a test genome's
#     haplotypes into an imputation answer; for ADMIXTURE and LAI it is simply
#     more reference data.
#
#   * 1000 Genomes is NEVER a donor, and never a real test sample. Those
#     genomes are what the external servers impute, and TOPMed, All of Us and
#     Michigan all contain 1KG -- so a 1KG test genome would be imputed
#     against its own haplotypes. Same reason sgdp_server_panel_ids.txt is
#     excluded from both roles: SGDP re-sequenced 1KG individuals under other
#     identifiers, and the servers hold those people whatever our cohort
#     column says.
#
#   * Donors are excluded from benchmark_reference_panel/ for a different
#     reason again: the simulated test genomes are built FROM donor
#     haplotypes, so a donor left in the LAI reference hands FLARE and RFMix
#     the exact haplotypes they are supposed to infer.
#
# Three exclusions, three distinct failures. The benchmark's imputation
# readings depend entirely on the test samples being absent from the external
# servers, which is what those exclusions guarantee -- nothing about these
# trees can establish it.
#
#   benchmark_cohort/           THE COHORT. The test genomes themselves --
#                               real/ (real individuals drawn from HGDP, SGDP,
#                               MXB, Almarri) and simulated/ (the 2,160 with
#                               their truth tracts) -- plus the manifests that
#                               say who is in which chip x build cell.
#
# An earlier revision shipped the benchmark REFERENCE lists under the name
# benchmark_cohort/ and shipped no cohort at all. That is the same failure as
# the suffix collision below: a name asserting one thing while holding
# another, with nothing inside to contradict it.
#
# Shipping the wrong reference panel is silent in both directions: a
# donor-inclusive reference inflates the benchmark, and a donor-excluded
# reference needlessly weakens the GWAS. Nothing in a BCF says which policy
# produced it.
#
# THE SUFFIX COLLISION. Both modes of split-homog-donors-reference.py write
# <sp>_reference_<suffix>.txt. Run them with the SAME suffix and the second
# overwrites the first -- one panel on disk, nothing in the filename saying
# which policy made it, and a 50/50 chance the benchmark ships a
# donor-inclusive reference. So the two modes get distinct suffixes, and
# verify_panel_distinction() refuses to ship if they end up identical anyway.
LAI_SUFFIX="${LAI_SUFFIX:-lai}"
BENCH_SUFFIX="${BENCH_SUFFIX:-bench}"
BENCH_OUT="${BENCH_OUT:-${REPO_DIR}/results}"
SUPERPOPS_T="${SUPERPOPS_T:-EUR EAS AFR AMR SAS MEN OCE}"
# The cohort's two halves on disk.
SIMSET_DIR="${SIMSET_DIR:-${PROJECT_ROOT}/06_benchmark_simset}"
REAL_DIR="${REAL_DIR:-${PROJECT_ROOT}/07_benchmark_real}"
# Everything for the genotype harmonization pipeline's testing lives under one
# parent on mab, kept apart from lai_reference_panel/ so the two projects
# cannot be picked up interchangeably. The donor policies are opposite and
# nothing inside a BCF says which one produced it.
HARMONIZATION_SUBDIR="${HARMONIZATION_SUBDIR:-genotype_harmonization_pipeline_testing}"

# split-homog-donors-reference.py LOWERCASES the super-pop in the filename
# (eur_reference_bench.txt). Looking only for the upper-case spelling finds
# nothing, and "nothing" here reads as "no lists were produced" -- which would
# make the distinction check pass vacuously on two empty sets. Resolve the
# spelling instead of assuming one.
role_file () {   # $1 = super-pop, $2 = reference|simulation, $3 = suffix
    local sp="$1" kind="$2" sfx="$3" f
    for f in "${REFS}/$(echo "$sp" | tr 'A-Z' 'a-z')_${kind}_${sfx}.txt" \
             "${REFS}/${sp}_${kind}_${sfx}.txt"; do
        if [[ -s "$f" ]]; then printf '%s\n' "$f"; return 0; fi
    done
    return 1
}

ids_union () {   # $1 = reference|simulation, $2 = suffix
    local kind="$1" sfx="$2" sp f
    for sp in $SUPERPOPS_T; do
        if f=$(role_file "$sp" "$kind" "$sfx"); then
            awk 'NF && $1 !~ /^#/ {print $1}' "$f"
        fi
    done | sort -u
}

verify_panel_distinction () {
    local t="${PROJECT_ROOT}/.verify_panels" rc=0 n ov
    rm -rf "$t"; mkdir -p "$t"
    ids_union reference  "$LAI_SUFFIX"   > "$t/lai_ref"
    ids_union reference  "$BENCH_SUFFIX" > "$t/bench_ref"
    ids_union simulation "$BENCH_SUFFIX" > "$t/bench_don"

    echo "-------------------------------------------------------------------"
    echo " panel distinction check"
    echo "-------------------------------------------------------------------"
    printf '  LAI reference   (suffix %-6s) %6d ids\n' "$LAI_SUFFIX"   "$(wc -l < "$t/lai_ref")"
    printf '  bench reference (suffix %-6s) %6d ids\n' "$BENCH_SUFFIX" "$(wc -l < "$t/bench_ref")"
    printf '  bench donors    (suffix %-6s) %6d ids\n' "$BENCH_SUFFIX" "$(wc -l < "$t/bench_don")"

    if [[ ! -s "$t/lai_ref" ]]; then
        echo "  ERROR: no LAI reference lists found. Run:"
        echo "    python3 -I scripts/split-homog-donors-reference.py \\"
        echo "        --panel ${SAMPLES_TSV_K} --mode lai --suffix ${LAI_SUFFIX}"
        rc=1
    fi
    if [[ ! -s "$t/bench_ref" ]]; then
        echo "  ERROR: no benchmark reference lists found. Run:"
        echo "    python3 -I scripts/split-homog-donors-reference.py \\"
        echo "        --panel ${SAMPLES_TSV_K} --mode benchmark --suffix ${BENCH_SUFFIX}"
        rc=1
    fi
    if [[ -s "$t/lai_ref" && -s "$t/bench_ref" ]] && cmp -s "$t/lai_ref" "$t/bench_ref"; then
        echo "  ERROR: the LAI and benchmark reference panels are IDENTICAL."
        echo "         Either both modes ran with the same suffix and one"
        echo "         overwrote the other, or --mode never took effect."
        echo "         Refusing to ship: one of the two would be wrong and"
        echo "         there is no way to tell which from the files alone."
        rc=1
    fi
    if [[ -s "$t/bench_ref" && -s "$t/bench_don" ]]; then
        ov=$(comm -12 "$t/bench_ref" "$t/bench_don" | wc -l)
        if [[ "$ov" -gt 0 ]]; then
            echo "  ERROR: ${ov} benchmark sample(s) are BOTH donor and reference."
            comm -12 "$t/bench_ref" "$t/bench_don" | head -5 | sed 's/^/           /'
            rc=1
        fi
    fi
    # No 1000 Genomes sample may be a benchmark donor: the imputation servers'
    # own panels contain 1KG, so a 1KG donor would be imputed against its own
    # haplotypes.
    if [[ -s "$SAMPLES_TSV_K" && -s "$t/bench_don" ]]; then
        n=$(awk -F'\t' 'NR==FNR{d[$1]=1; next}
              FNR==1{for(i=1;i<=NF;i++){if($i=="sample")si=i; if($i=="cohort")ci=i} next}
              (si && ci && ($si in d) && $ci=="1KG"){n++} END{print n+0}' \
              "$t/bench_don" "$SAMPLES_TSV_K")
        printf '  1KG samples among benchmark donors: %d (must be 0)\n' "$n"
        [[ "$n" -gt 0 ]] && rc=1
    fi
    rm -rf "$t"
    if [[ "$rc" -ne 0 ]]; then
        echo "-------------------------------------------------------------------"
        echo " REFUSING TO TRANSFER. Fix the above and re-run."
        echo "-------------------------------------------------------------------"
        return 1
    fi
    echo "  ok -- the two panels are distinct and the benchmark split is clean"
    echo "-------------------------------------------------------------------"
    return 0
}

write_provenance () {   # $1 = stage dir, $2 = title, $3 = donor policy, $4 = suffix
    local stage="$1"
    {
        echo "$2"
        printf '%s\n' "$(printf '=%.0s' $(seq 1 ${#2}))"
        echo
        echo "generated   $(date '+%Y-%m-%d %H:%M:%S %Z')"
        echo "host        $(hostname)"
        echo "source      ${PROJECT_ROOT}"
        echo "K           ${K}"
        echo "split mode  $4"
        echo "repo commit $(cd "$REPO_DIR" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
        echo "repo branch $(cd "$REPO_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
        echo
        echo "DONOR POLICY"
        echo "$3"
        echo
        echo "The HOMOG panel is HGDP + 1KG + MXB joint-called and filtered by"
        echo "supervised ADMIXTURE. It INCLUDES MXB samples -- it is not"
        echo "MXB-free."
    } > "${stage}/PROVENANCE.txt"
}

run_lai_refpanel () {
    local STAGE="${PROJECT_ROOT}/.stage_lai_reference_panel"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    local sp role f n=0
    for sp in $SUPERPOPS_T; do
        for role in reference simulation; do
            if f=$(role_file "$sp" "$role" "$LAI_SUFFIX"); then
                ln -sf "$f" "$STAGE/$(basename "$f")"; n=$((n+1))
            fi
        done
    done
    if [[ "$n" -eq 0 ]]; then
        echo "WARNING: no LAI split lists with suffix '${LAI_SUFFIX}'; skipping"
        return 0
    fi
    for f in "$SAMPLES_TSV_K" "${HOMOG_DIR}/homog_${K}pop_samples.txt" \
             "${HOMOG_DIR}/homog_composition.k${K}.tsv"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    write_provenance "$STAGE" \
        "LAI / GWAS reference panel (donors INCLUDED)" \
"  Donors are NOT withheld. Every homogeneous sample is in the reference.
  Use this panel for FLARE and for the GWAS, where no genome is simulated
  from these haplotypes, so there is no circularity to avoid and excluding
  donors would only shrink the panel.

  Do NOT use it to score the benchmarking cohort: those test genomes were
  simulated from donor haplotypes that are present here, which would make
  the method's own input its answer key." \
        "lai (suffix ${LAI_SUFFIX})"
    do_rsync "$STAGE" "lai_reference_panel"
}

run_harmonization_readme () {
    # A README at the parent, so the folder explains itself to whoever opens
    # it on mab rather than depending on this conversation being remembered.
    # The reference/imputation distinction is the thing most likely to be got
    # wrong, so it leads.
    local STAGE="${PROJECT_ROOT}/.stage_harmonization_readme"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    cat > "$STAGE/README.txt" <<EOF
GENOTYPE HARMONIZATION PIPELINE -- TESTING MATERIAL
===================================================

generated   $(date '+%Y-%m-%d %H:%M:%S %Z')
source      ${PROJECT_ROOT}
repo        $(cd "$REPO_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown) @ $(cd "$REPO_DIR" && git rev-parse --short HEAD 2>/dev/null || echo unknown)

WHAT IS IN HERE

  benchmark_reference_panel/   ADMIXTURE + LAI reference, donors EXCLUDED
  benchmark_cohort/            the test genomes
      real/                      real HGDP / SGDP / MXB / Almarri individuals
      simulated/                 the 2,160, per chromosome, with truth tracts

WHICH REFERENCE GOES WITH WHICH READING

  global ancestry Q (ADMIXTURE)        benchmark_reference_panel/
  local ancestry (FLARE / RFMix)       benchmark_reference_panel/
  imputation accuracy / calibration    TOPMed, All of Us, Michigan (1000
                                       Genomes) -- the EXTERNAL servers
  power and GWAS simulations           real genotypes + the ADMIXTURE Q as
                                       the ancestry covariate

NEITHER TREE IN HERE IS AN IMPUTATION PANEL. Imputation goes to the external
servers. The panel here is for ADMIXTURE and local-ancestry inference only.

WHAT NOT TO USE

  ../lai_reference_panel/ is the SAME KIND of panel with the OPPOSITE donor
  policy -- donors INCLUDED. It is for the GWAS and FLARE work, where nothing
  is simulated from those haplotypes. Scoring this cohort against it would
  hand FLARE and RFMix the exact haplotypes they are meant to infer, because
  these test genomes were simulated FROM those donors. Nothing inside a BCF
  says which policy produced it, which is why the two live in separate trees.

THE THREE EXCLUSIONS, AND WHY EACH IS DIFFERENT

  1KG is IN benchmark_reference_panel/ and that is correct: this panel never
  imputes, so 1KG here cannot leak a test genome's haplotypes into an
  imputation answer, and for ADMIXTURE and LAI it is simply more reference
  data.

  1KG is NEVER a donor and never a real test genome: the external servers all
  contain 1KG, so such a genome would be imputed against its own haplotypes.
  Same reason sgdp_server_panel_ids.txt is barred from both roles -- SGDP
  re-sequenced 1KG individuals under other identifiers, and the servers hold
  those people whatever the cohort column says.

  Donors are excluded from benchmark_reference_panel/ for a third reason: the
  simulated test genomes are built FROM donor haplotypes.

The imputation readings rest ENTIRELY on the test samples being absent from
TOPMed, All of Us and Michigan. Nothing about the panels in here can
establish that -- only those exclusions can.

SCORING RULES

  simulated individuals   local ancestry AND imputation. They have recorded
                          truth tracts (.hanc) and the .bp recombination
                          realisation that regenerates them.
  real individuals        imputation accuracy / calibration / retention, and
                          power / GWAS simulations. NEVER local ancestry:
                          nobody recorded which ancestry each segment came
                          from, and inferred tracts would be another method's
                          output, not truth.
  gen-0 controls          no admixture to get wrong, so their error is pure
                          method error, and their imputation accuracy is the
                          ceiling the admixed individuals of that component
                          are measured against.

Each subdirectory has its own PROVENANCE.txt with the donor policy, the split
mode and the repo commit. The HOMOG panel is HGDP + 1KG + MXB joint-called and
filtered by supervised ADMIXTURE; it INCLUDES MXB samples.
EOF
    do_rsync "$STAGE" "${HARMONIZATION_SUBDIR}"
}

run_benchmark_refpanel () {
    local STAGE="${PROJECT_ROOT}/.stage_benchmark_reference_panel"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    local sp role f n=0
    for sp in $SUPERPOPS_T; do
        for role in reference simulation; do
            if f=$(role_file "$sp" "$role" "$BENCH_SUFFIX"); then
                ln -sf "$f" "$STAGE/$(basename "$f")"; n=$((n+1))
            fi
        done
    done
    if [[ "$n" -eq 0 ]]; then
        echo "WARNING: no benchmark split lists with suffix '${BENCH_SUFFIX}'; skipping"
        return 0
    fi
    for f in "${REFS}/known_duplicate_pairs.tsv" \
             "${REFS}/sgdp_server_panel_ids.txt" \
             "${HOMOG_DIR}/homog_composition.k${K}.tsv" \
             "$SAMPLES_TSV_K"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    # Composition of the reference side, so "population balanced" is auditable
    # rather than asserted. Per super-pop: how many came from 1000 Genomes and
    # how many from everything else.
    if [[ -s "$SAMPLES_TSV_K" ]]; then
        {
            printf 'super_pop\tn_reference\tn_1KG\tn_other\n'
            for sp in $SUPERPOPS_T; do
                f=$(role_file "$sp" reference "$BENCH_SUFFIX") || continue
                awk -F'\t' -v sp="$sp" '
                  NR==FNR { if ($1 !~ /^#/ && NF) r[$1]=1; next }
                  FNR==1  { for(i=1;i<=NF;i++){if($i=="sample")si=i; if($i=="cohort")ci=i} next }
                  (si && ci && ($si in r)) { t++; if ($ci=="1KG") k++; else o++ }
                  END { printf "%s\t%d\t%d\t%d\n", sp, t+0, k+0, o+0 }' \
                  "$f" "$SAMPLES_TSV_K"
            done
        } > "$STAGE/reference_panel_composition_${BENCH_SUFFIX}.tsv"
    fi
    write_provenance "$STAGE" \
        "Benchmark REFERENCE PANEL (donors EXCLUDED, 1KG included)" \
"  NOT AN IMPUTATION PANEL. This is the reference for supervised ADMIXTURE
  (the global-ancestry Q) and for local-ancestry inference (FLARE / RFMix),
  and for nothing else. In the genotyping_pipeline flow imputation runs
  against the EXTERNAL servers -- TOPMed, All of Us and Michigan (1000
  Genomes) -- never against this.

  It is also not the cohort: the test genomes are in benchmark_cohort/.

  Composition: 1000 Genomes, plus a population-stratified share of the
  non-1KG samples so both sides carry a proportional slice of every
  population. See reference_panel_composition_${BENCH_SUFFIX}.tsv.

  WHY 1KG IS HERE BUT IS NEVER A DONOR OR A TEST GENOME. Because this panel
  never imputes, 1KG in it cannot leak a test genome's haplotypes into an
  imputation answer -- for ADMIXTURE and LAI it is simply more reference
  data. The external servers DO contain 1KG, so a 1KG donor or a 1KG test
  genome would be imputed against its own haplotypes. That is also why
  sgdp_server_panel_ids.txt is excluded from both roles: SGDP re-sequenced
  1KG individuals under other identifiers, and the servers hold those people
  whatever our cohort column says.

  Donors are excluded from this reference for a different reason again: the
  simulated test genomes are built FROM donor haplotypes, so a donor left
  here hands FLARE and RFMix the haplotypes they are meant to infer." \
        "benchmark (suffix ${BENCH_SUFFIX})"
    do_rsync "$STAGE" "${HARMONIZATION_SUBDIR}/benchmark_reference_panel"
}

run_benchmark_cohort () {
    # THE COHORT ITSELF: the test genomes, not the panel they are scored
    # against. real/ and simulated/ are separate subtrees because they are
    # scored differently -- see PROVENANCE.txt.
    local STAGE="${PROJECT_ROOT}/.stage_benchmark_cohort"
    rm -rf "$STAGE"; mkdir -p "$STAGE/real" "$STAGE/simulated"
    local f d n_real=0 n_sim=0

    for f in "${REAL_DIR}"/*.bcf "${REAL_DIR}"/*.bcf.csi; do
        [[ -e "$f" ]] || continue
        ln -sf "$f" "$STAGE/real/$(basename "$f")"; n_real=$((n_real+1))
    done
    for d in "${SIMSET_DIR}"/chr*; do
        [[ -d "$d" ]] || continue
        mkdir -p "$STAGE/simulated/$(basename "$d")"
        for f in "$d"/*.hanc "$d"/*.haps "$d"/*.sample "$d"/*.phgeno \
                 "$d"/*.snp "$d"/*.dat "$d"/*.bp; do
            [[ -e "$f" ]] || continue
            ln -sf "$f" "$STAGE/simulated/$(basename "$d")/$(basename "$f")"
            n_sim=$((n_sim+1))
        done
    done
    if [[ -d "${SIMSET_DIR}/provenance" ]]; then
        ln -sfn "${SIMSET_DIR}/provenance" "$STAGE/simulated/provenance"
    fi

    # The manifests belong with the cohort: they are what says who is in it.
    for f in "${BENCH_OUT}/benchmark_cohort_manifest.tsv" \
             "${BENCH_OUT}/benchmark_cohort_simset_individuals.tsv" \
             "${BENCH_OUT}/benchmark_cohort_simset_roster.tsv" \
             "${BENCH_OUT}/benchmark_cohort_counts.txt" \
             "${BENCH_OUT}/cohort_composition.txt" \
             "$SAMPLES_TSV_K"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done

    if [[ "$n_real" -eq 0 && "$n_sim" -eq 0 ]]; then
        echo "WARNING: no cohort data found."
        echo "         real genomes:  ${REAL_DIR}      (run 6c_extract_real_test_genomes_clm.sh)"
        echo "         simulated:     ${SIMSET_DIR}    (run 6b_simulate_benchmark_cohort_clm.sh)"
        echo "         Skipping benchmark_cohort -- shipping manifests with no"
        echo "         genomes would put an empty cohort on mab under a name"
        echo "         saying otherwise."
        return 0
    fi
    echo "[$(date +%T)] cohort: ${n_real} real files, ${n_sim} simulated files"

    write_provenance "$STAGE" \
        "Benchmark COHORT (the test genomes)" \
"  real/       real individuals from HGDP, SGDP, MXB and Almarri. Outside
              every imputation server panel, outside the donor set and
              outside both reference panels. Labelled real_homogeneous or
              real_admixed in benchmark_cohort_manifest.tsv.

              NEVER SCORED FOR LOCAL-ANCESTRY ACCURACY. A real genome has
              a WGS truth but nobody recorded which ancestry each segment
              came from, and inferred tracts would be another method's
              output, not truth.

              Used for (a) imputation accuracy, calibration and retention,
              and (b) POWER AND GWAS SIMULATIONS -- phenotypes simulated on
              real genotypes. (b) is why the real half cannot be replaced by
              more simulated individuals: power depends on the LD structure,
              allele-frequency spectrum and relatedness the genotypes
              actually carry, and simulated genomes drawn from a finite
              donor pool have compressed haplotype diversity plus
              donor-reuse cryptic relatedness.

              For (b) the ancestry covariate is the per-sample ADMIXTURE Q:
              max_Q_sup / argmax_pop_sup are in the manifest and the full Q
              table travels as lai_ref_panel_samples.k8.tsv.

  simulated/  the 2,160 (1,800 admixed + 360 single-ancestry controls),
              per chromosome, with .hanc truth tracts and the .bp
              recombination realisation that regenerates them. Scored for
              local ancestry AND imputation.

              The gen-0 controls have no admixture to get wrong, so their
              error is pure method error and their imputation accuracy is
              the ceiling the admixed individuals of that component are
              measured against.

  Every individual is on exactly one chip and one build; the cell is in
  benchmark_cohort_simset_individuals.tsv (simulated) and
  benchmark_cohort_manifest.tsv (real). Truth stays GRCh38 for everyone.

  WHICH REFERENCE GOES WITH WHICH READING:

    ADMIXTURE and LAI (FLARE / RFMix)
        -> benchmark_reference_panel/, NOT lai_reference_panel/. The latter
           contains the donors these genomes were simulated from, which
           would hand the method the haplotypes it is meant to infer.

    imputation
        -> the EXTERNAL servers: TOPMed, All of Us, Michigan (1000
           Genomes). Neither tree here is an imputation panel. The
           imputation readings are only valid because every genome in this
           cohort is absent from those servers -- no 1KG, and nothing in
           sgdp_server_panel_ids.txt, on either the donor or the test side." \
        "benchmark (suffix ${BENCH_SUFFIX})"
    do_rsync "$STAGE" "${HARMONIZATION_SUBDIR}/benchmark_cohort"
}

case "$SCOPE" in
    all)
        # The distinction check runs BEFORE any bytes move. Both panels ship
        # in this scope, so if they cannot be told apart the transfer must not
        # start -- a half-finished wrong transfer on mab is worse than none.
        verify_panel_distinction
        ssh ${MAB_SSH_OPTS} "${MAB_USER}@${MAB_HOST}" \
            "mkdir -p '${MAB_BASE}'/{hg38,hg19}/{'${ALL_SAMPLES_SUBDIR}',homog_5pop,'${HOMOG_K_SUBDIR}'} '${MAB_BASE}'/{reference_lists,lai_reference_panel} '${MAB_BASE}/${HARMONIZATION_SUBDIR}'/{benchmark_reference_panel,benchmark_cohort/{real,simulated}}"
        run_hg38_all
        run_hg38_homog
        run_hg38_homog_k
        run_hg19_all
        run_hg19_homog
        run_hg19_homog_k
        run_reference_lists
        run_lai_refpanel
        run_harmonization_readme
        run_benchmark_refpanel
        run_benchmark_cohort
        ;;
    hg38)               run_hg38_all; run_hg38_homog; run_hg38_homog_k ;;
    hg19)               run_hg19_all; run_hg19_homog; run_hg19_homog_k ;;
    hg38/all_samples)   run_hg38_all ;;
    hg38/homog_5pop)    run_hg38_homog ;;
    hg38/homog_k)       run_hg38_homog_k ;;
    hg19/all_samples)   run_hg19_all ;;
    hg19/homog_5pop)    run_hg19_homog ;;
    hg19/homog_k)       run_hg19_homog_k ;;
    reference_lists)    run_reference_lists ;;
    panels)             verify_panel_distinction; run_lai_refpanel; run_benchmark_refpanel ;;
    lai_reference_panel) verify_panel_distinction; run_lai_refpanel ;;
    benchmark_reference_panel) verify_panel_distinction; run_harmonization_readme; run_benchmark_refpanel ;;
    benchmark_cohort)   run_harmonization_readme; run_benchmark_cohort ;;
    benchmark)          verify_panel_distinction; run_harmonization_readme; run_benchmark_refpanel; run_benchmark_cohort ;;
    verify)             verify_panel_distinction; exit $? ;;
    help|--help|-h)     SCOPE_HELP=1 ;&
    *)
        [[ "${SCOPE_HELP:-0}" == "1" ]] || echo "Unknown scope: $SCOPE"
        echo "Usage: $0 <scope>"
        echo ""
        echo "  ADMIXTURE / LAI REFERENCE PANELS -- what a method is GIVEN."
        echo "  Neither is an imputation panel: imputation runs against the"
        echo "  external servers (TOPMed, All of Us, Michigan/1000 Genomes)."
        echo "    lai_reference_panel        donors INCLUDED; for FLARE and the GWAS"
        echo "    benchmark_reference_panel  1KG + population-stratified others,"
        echo "                               donors EXCLUDED; for scoring the"
        echo "                               benchmark cohort's ancestry calls"
        echo "    panels                     both of the above"
        echo ""
        echo "  Both benchmark trees ship under ${HARMONIZATION_SUBDIR}/"
        echo ""
        echo "  THE COHORT -- what a method is TESTED ON"
        echo "    benchmark_cohort           real/ + simulated/ test genomes, truth"
        echo "                               tracts and the cell manifests"
        echo "    benchmark                  benchmark_reference_panel + benchmark_cohort"
        echo ""
        echo "  OTHER"
        echo "    verify                     check the panels are distinct, move nothing"
        echo "    all | hg38 | hg19 | hg38/all_samples | hg38/homog_5pop | hg38/homog_k"
        echo "    hg19/all_samples | hg19/homog_5pop | hg19/homog_k | reference_lists"
        [[ "${SCOPE_HELP:-0}" == "1" ]] && exit 0
        exit 1
        ;;
esac

echo ""
echo "==================================================================="
echo "[$(date +%T)] transfer complete."
echo "On mab: ls -lh ${MAB_BASE}"
echo "==================================================================="
