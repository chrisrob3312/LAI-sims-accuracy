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

# The "verify" scope moves nothing, so it must not demand credentials: the
# whole point is to check the split locally before deciding to ship.
if [[ "$MAB_USER" == "<USERNAME>" && "${1:-all}" != "verify" ]]; then
    echo "ERROR: edit MAB_USER at the top of $0 (or export MAB_USER=...) before running."
    echo "       To check the two panels without transferring: $0 verify"
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
# THE TWO REFERENCE PANELS, AND WHY THEY MUST NOT BE CONFUSED
#
# Both go to mab, and they are different products of the same Q file:
#
#   lai_reference_panel/   every homogeneous sample, donors INCLUDED. This is
#                          what the GWAS and FLARE run against. Nothing is
#                          simulated there, so there are no donors to
#                          withhold, and withholding them would only shrink
#                          the panel for no benefit.
#
#   benchmark_cohort/      the donor-EXCLUDED reference, plus the cohort
#                          manifest. Here the same individuals' haplotypes
#                          were used to simulate the test genomes, so a donor
#                          left in the reference turns the method's own input
#                          into its answer key.
#
# Shipping the wrong one is silent in both directions: a donor-inclusive
# reference inflates the benchmark, and a donor-excluded reference needlessly
# weakens the GWAS. Nothing in a BCF says which policy produced it.
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

ids_union () {   # $1 = reference|simulation, $2 = suffix
    local kind="$1" sfx="$2" sp f
    for sp in $SUPERPOPS_T; do
        f="${REFS}/${sp}_${kind}_${sfx}.txt"
        if [[ -s "$f" ]]; then awk 'NF && $1 !~ /^#/ {print $1}' "$f"; fi
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
    local sp f n=0
    for sp in $SUPERPOPS_T; do
        for f in "${REFS}/${sp}_reference_${LAI_SUFFIX}.txt" \
                 "${REFS}/${sp}_simulation_${LAI_SUFFIX}.txt"; do
            [[ -s "$f" ]] && { ln -sf "$f" "$STAGE/$(basename "$f")"; n=$((n+1)); }
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

run_benchmark_cohort () {
    local STAGE="${PROJECT_ROOT}/.stage_benchmark_cohort"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    local sp f n=0
    for sp in $SUPERPOPS_T; do
        for f in "${REFS}/${sp}_reference_${BENCH_SUFFIX}.txt" \
                 "${REFS}/${sp}_simulation_${BENCH_SUFFIX}.txt"; do
            [[ -s "$f" ]] && { ln -sf "$f" "$STAGE/$(basename "$f")"; n=$((n+1)); }
        done
    done
    if [[ "$n" -eq 0 ]]; then
        echo "WARNING: no benchmark split lists with suffix '${BENCH_SUFFIX}'; skipping"
        return 0
    fi
    for f in "${BENCH_OUT}/benchmark_cohort_manifest.tsv" \
             "${BENCH_OUT}/benchmark_cohort_simset_roster.tsv" \
             "${BENCH_OUT}/benchmark_cohort_counts.txt" \
             "${BENCH_OUT}/cohort_composition.txt" \
             "${REFS}/known_duplicate_pairs.tsv" \
             "${REFS}/sgdp_server_panel_ids.txt" \
             "$SAMPLES_TSV_K"; do
        [[ -s "$f" ]] && ln -sf "$f" "$STAGE/$(basename "$f")"
    done
    write_provenance "$STAGE" \
        "Benchmarking cohort reference panel (donors EXCLUDED)" \
"  Donors are withheld from the reference, and no 1000 Genomes sample is a
  donor. Both exclusions exist to stop a method being scored against
  haplotypes it was given: the simulated test genomes are built from the
  donor half, and the imputation servers' own panels contain 1KG.

  benchmark_cohort_manifest.tsv carries the real test genomes with their
  sample_type (real_homogeneous / real_admixed) and chip x build cell.
  Real genomes have no recorded tracts, so they are scored for imputation
  accuracy, calibration and retention only -- never for local ancestry." \
        "benchmark (suffix ${BENCH_SUFFIX})"
    do_rsync "$STAGE" "benchmark_cohort"
}

case "$SCOPE" in
    all)
        # The distinction check runs BEFORE any bytes move. Both panels ship
        # in this scope, so if they cannot be told apart the transfer must not
        # start -- a half-finished wrong transfer on mab is worse than none.
        verify_panel_distinction
        ssh ${MAB_SSH_OPTS} "${MAB_USER}@${MAB_HOST}" \
            "mkdir -p '${MAB_BASE}'/{hg38,hg19}/{'${ALL_SAMPLES_SUBDIR}',homog_5pop,'${HOMOG_K_SUBDIR}'} '${MAB_BASE}'/{reference_lists,lai_reference_panel,benchmark_cohort}"
        run_hg38_all
        run_hg38_homog
        run_hg38_homog_k
        run_hg19_all
        run_hg19_homog
        run_hg19_homog_k
        run_reference_lists
        run_lai_refpanel
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
    panels)             verify_panel_distinction; run_lai_refpanel; run_benchmark_cohort ;;
    lai_reference_panel) verify_panel_distinction; run_lai_refpanel ;;
    benchmark_cohort)   verify_panel_distinction; run_benchmark_cohort ;;
    verify)             verify_panel_distinction; exit $? ;;
    *)
        echo "Unknown scope: $SCOPE"
        echo "Usage: $0 [all | hg38 | hg19 | hg38/all_samples | hg38/homog_5pop | hg38/homog_k | hg19/all_samples | hg19/homog_5pop | hg19/homog_k | reference_lists | panels | lai_reference_panel | benchmark_cohort | verify]"
        echo ""
        echo "  verify               check the two panels are distinct, move nothing"
        echo "  panels               both reference panels (LAI donor-inclusive + benchmark donor-excluded)"
        echo "  lai_reference_panel  donors INCLUDED -- for FLARE and the GWAS"
        echo "  benchmark_cohort     donors EXCLUDED -- for scoring the benchmarking cohort"
        exit 1
        ;;
esac

echo ""
echo "==================================================================="
echo "[$(date +%T)] transfer complete."
echo "On mab: ls -lh ${MAB_BASE}"
echo "==================================================================="
