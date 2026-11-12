#!/usr/bin/env bash
# =============================================================================
# verify-rebuild.sh -- where the rebuild is, and whether each finished stage is
# right. READ-ONLY. Run it any time: while the chain runs, or after.
#
# For each stage, in run order: its Slurm jobs (queued now, and finished since
# the rebuild epoch), then checks on what it wrote. Each check is PASS, WARN or
# FAIL with the numbers behind it; a stage with no output yet says "not yet".
#   0 SGDP download   1 SGDP prep   2 MXB lift   3 Almarri   4 1c merge
#   5 chunks   6 ligate   7 KING   8 4a K=8   9 simulation / RFMix / FLARE / scoring
#
# Run from the repo checkout inside the shapeit5 env:
#   bash scripts/verify-rebuild.sh
#   bash scripts/verify-rebuild.sh 2>&1 | tee verify_$(date +%Y%m%d_%H%M).txt
# When a stage FAILs, read that stage's own .out and .err in full before
# anything else (logs/ paths are printed per stage).
# =============================================================================
set -uo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
HOMOG_DIR="${HOMOG_DIR:-${PROJECT_ROOT}/04_homogeneity_panel}"
LOGDIR="${LOGDIR:-${PROJECT_ROOT}/logs}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"
FIXREF_MODE="${FIXREF_MODE:-flip-all}"
ADMIX_POP="${ADMIX_POP:-Brasa}"; GEN="${GEN:-12}"; K="${K:-8}"
EPOCH="${PROJECT_ROOT}/rebuild_epoch"
U="${USER:-$(id -un)}"
if [[ -f "$EPOCH" ]]; then SINCE=$(date -r "$EPOCH" +%Y-%m-%d); else SINCE=$(date -d '-7 days' +%Y-%m-%d); fi
cd "$REPO_DIR" || exit 99

NP=0; NW=0; NF=0
pass () { NP=$((NP + 1)); printf '  PASS  %s\n' "$*"; }
warn () { NW=$((NW + 1)); printf '  WARN  %s\n' "$*"; }
fail () { NF=$((NF + 1)); printf '  FAIL  %s\n' "$*"; }
info () { printf '  info  %s\n' "$*"; }
section () { printf '\n== %s ==\n' "$*"; }
# Slurm view of one job name: array elements queued now, then sacct states of
# the LATEST submission of that name since $SINCE. Earlier submissions (a
# cancelled duplicate, a run that was fixed and re-submitted) are listed as
# superseded and do not count as failures; a failure in the latest one does.
jobs_of () {
    local q rows latest cur old
    q=$(squeue -h -u "$U" -n "$1" -r -o '%T' | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')
    rows=$(sacct -n -X -u "$U" --name "$1" -S "$SINCE" -o JobID%30,State%20 \
           | awk 'NF {split($1, a, "_"); print a[1], $2}')
    latest=$(awk 'NF {print $1}' <<<"$rows" | sort -n | tail -1)
    cur=$(awk -v j="$latest" 'NF && $1 == j {print $2}' <<<"$rows" | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')
    old=$(awk -v j="$latest" 'NF && $1 != j {print $2}' <<<"$rows" | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')
    info "jobs '$1': queued [${q:-none}]  latest ${latest:-none} [${cur:-none}]${old:+  superseded [${old}]}  logs: ${LOGDIR}/${2}"
    if [[ "$cur" == *FAILED* || "$cur" == *TIMEOUT* || "$cur" == *OUT_OF_ME* || "$cur" == *NODE_FAIL* ]]; then
        fail "'$1' job ${latest} has failed tasks: read their .out/.err in full (sacct -X -j ${latest})"
    fi
}
nsamples () { bcftools query -l "$1" | wc -l; }
nrecs () { bcftools index -s "$1" | awk '{s += $3} END {print s + 0}'; }   # -n warns on .csi; -s does not
# Only what this rebuild wrote counts: older files are the November outputs
# that the running stages will replace.
is_new () { [[ -s "$1" && -f "$EPOCH" && "$1" -nt "$EPOCH" ]]; }
palindromic () {   # palindromic <bcf> [region] -> "n_snps pct"
    bcftools query ${2:+-r "$2"} -i 'TYPE="snp"' -f '%REF%ALT\n' "$1" \
      | awk '{n++} $1=="AT"||$1=="TA"||$1=="CG"||$1=="GC" {p++} END {printf "%d %.1f", n, (n ? 100 * p / n : 0)}'
}

echo "verify-rebuild  $(date '+%F %T')  commit $(git rev-parse --short HEAD)  epoch: $(cat "$EPOCH" 2>&1)"

section "0. SGDP download (scripts/fetch-sgdp-vcf-variants.sh)"
jobs_of sgdp_fetch 'sgdp_fetch_<job>.out'
if [[ -f "${SGDP_RAW_DIR}/vcfs/DONE" ]]; then
    pass "verified: $(tr '\n' ';' < "${SGDP_RAW_DIR}/vcfs/DONE")"
    grep -q 'mismatch 0' "${SGDP_RAW_DIR}/vcfs/DONE" && pass "every shipped md5 matched" || fail "md5 mismatches recorded in DONE"
else info "not yet (no ${SGDP_RAW_DIR}/vcfs/DONE)"; fi

section "1. SGDP prep (1b_sgdp_raw_prep_clm.sh)"
jobs_of 1b_sgdp_raw '1b_sgdp_raw_chr<N>_<job>.out/.err'
NB=0; NBAD=0; SN=""
for c in $(seq 1 22); do
    f="${PANEL_DIR}/sgdp_prepped.chr${c}.bcf"
    if [[ -s "$f" && -s "${f}.recipe" ]]; then
        NB=$((NB + 1)); r=$(cat "${f}.recipe")
        if [[ "$r" != sgdp_raw_v1* || "$r" != *"fixref=${FIXREF_MODE}"* ]]; then NBAD=$((NBAD + 1)); warn "chr${c}: recipe '${r}'"; fi
        SN+="$(nsamples "$f") "
    fi
done
if [[ "$NB" -eq 0 ]]; then info "not yet"
else
    [[ "$NB" -eq 22 ]] && pass "22 of 22 chromosomes built from the raw calls" || info "${NB} of 22 chromosomes built so far"
    [[ "$NBAD" -eq 0 ]] && pass "every stamp: sgdp_raw_v1 ... fixref=${FIXREF_MODE}" || fail "${NBAD} chromosome(s) with another recipe"
    US=$(tr ' ' '\n' <<<"$SN" | awk 'NF' | sort -u | tr '\n' ' ')
    [[ $(wc -w <<<"$US") -eq 1 ]] && pass "same sample count on every chromosome: ${US}" || fail "sample counts differ across chromosomes: ${US}"
    if [[ -s "${PANEL_DIR}/sgdp_prepped.chr21.bcf" ]]; then
        read -r n p <<<"$(palindromic "${PANEL_DIR}/sgdp_prepped.chr21.bcf")"
        awk -v p="$p" 'BEGIN {exit !(p > 5)}' && pass "chr21: ${n} SNPs, ${p}% A/T or C/G (kept)" || fail "chr21: ${n} SNPs, only ${p}% A/T or C/G -- palindromic SNPs dropped?"
    fi
fi

section "2. MXB lift (1a_prep-mxb-liftover_clm.sh)"
jobs_of 1a_prep_mxb_liftover '1a_prep_mxb_liftover_<job>.out/.err'
if [[ -s "${PANEL_DIR}/mxb_lifted.recipe" ]] && grep -q "fixref=${FIXREF_MODE}" "${PANEL_DIR}/mxb_lifted.recipe"; then
    pass "stamp: $(cat "${PANEL_DIR}/mxb_lifted.recipe")"
    f="${PANEL_DIR}/mxb_lifted.chr21.bcf"
    [[ "$(nsamples "$f")" -eq 50 ]] && pass "chr21: 50 samples" || fail "chr21: $(nsamples "$f") samples (50 expected)"
    read -r n p <<<"$(palindromic "$f")"
    awk -v p="$p" 'BEGIN {exit !(p > 5)}' && pass "chr21: ${n} SNPs, ${p}% A/T or C/G (kept)" || fail "chr21: ${n} SNPs, only ${p}% A/T or C/G"
elif [[ -s "${PANEL_DIR}/mxb_lifted.recipe" ]]; then info "not yet re-lifted with fixref=${FIXREF_MODE} (stamp: '$(cat "${PANEL_DIR}/mxb_lifted.recipe")')"
else info "not yet (no MXB lift stamp; the re-lift writes it when all 22 chromosomes are done)"; fi

section "3. Almarri (1d_almarri_prep_clm.sh)"
jobs_of 1d_almarri '1d_almarri_chr<N>_<job>.out/.err'
NV=0; NR=0
for c in $(seq 1 22); do
    f="${PANEL_DIR}/almarri_prepped.chr${c}.bcf"
    r=""; [[ -s "${f}.recipe" ]] && r=$(cat "${f}.recipe")
    if [[ "$r" == "v3 fixref=${FIXREF_MODE}" ]]; then
        NV=$((NV + 1)); bcftools view -h "$f" > /dev/null && NR=$((NR + 1))
    fi
done
if [[ "$NV" -eq 0 ]]; then info "not yet (no chromosome stamped v3 fixref=${FIXREF_MODE})"
else
    [[ "$NV" -eq 22 ]] && pass "22 of 22 chromosomes at v3 fixref=${FIXREF_MODE}" || info "${NV} of 22 chromosomes done so far"
    [[ "$NR" -eq "$NV" ]] && pass "all ${NR} read back" || fail "$((NV - NR)) stamped chromosome(s) do not read back"
    f="${PANEL_DIR}/almarri_prepped.chr21.bcf"
    [[ -s "$f" ]] && { [[ "$(nsamples "$f")" -eq 137 ]] && pass "chr21: 137 samples" || fail "chr21: $(nsamples "$f") samples (137 expected)"; }
fi

section "4. merge + site selection (1c_phasing-jointcall_clm.sh)"
jobs_of 1c_phase '1c_phase_chr<N>_<job>.out/.err'
NS=0; NSOLD=0; SC=""
for c in $(seq 1 22); do
    f="${PANEL_DIR}/merged_chr${c}.softunion.bcf"
    [[ -s "$f" ]] || continue
    is_new "$f" || { NSOLD=$((NSOLD + 1)); continue; }
    [[ "$NS" -eq 0 ]] && printf '  %-5s %-8s %-11s %-11s %-12s %s\n' chr samples softunion admix_sites fill_dropped stamp
    NS=$((NS + 1)); s=$(nsamples "$f"); SC+="${s} "
    fc="${PANEL_DIR}/provenance/chr${c}.fill_check.tsv"; fd="-"
    [[ -s "$fc" ]] && fd=$(awk 'NR == 1 {for (i = 1; i <= NF; i++) if ($i ~ /^positions_dropped=/) {sub(/.*=/, "", $i); print $i}}' "$fc")
    ab="${PANEL_DIR}/merged_chr${c}.admix_sites.bed"; na="-"; [[ -s "$ab" ]] && na=$(wc -l < "$ab")
    printf '  %-5s %-8s %-11s %-11s %-12s %s\n' "$c" "$s" "$(nrecs "$f")" "$na" "${fd:--}" \
        "$( [[ -s "${f}.recipe" ]] && echo yes || echo NO)"
done
if [[ "$NS" -eq 0 ]]; then
    if [[ "$NSOLD" -gt 0 ]]; then info "not yet (${NSOLD} soft-unions on disk are from before the rebuild; 1c replaces them)"; else info "not yet"; fi
else
    [[ "$NS" -eq 22 ]] && pass "22 of 22 soft-unions published" || info "${NS} of 22 published so far"
    US=$(tr ' ' '\n' <<<"$SC" | awk 'NF' | sort -u | tr '\n' ' ')
    [[ $(wc -w <<<"$US") -eq 1 ]] && pass "same sample count on every chromosome: ${US}" || fail "sample counts differ across chromosomes: ${US}"
    F21="${PANEL_DIR}/merged_chr21.softunion.bcf"
    if is_new "$F21"; then
        read -r n p <<<"$(palindromic "$F21")"
        awk -v p="$p" 'BEGIN {exit !(p > 5)}' && pass "chr21 soft-union: ${n} SNPs, ${p}% A/T or C/G" || fail "chr21 soft-union: only ${p}% A/T or C/G"
        bcftools query -l "$F21" | sort > "${TMPDIR:-/tmp}/verify_panel_ids.$$"
        for l in amr_rfmix afr_rfmix eur_rfmix pel_rfmix pel_eas_rfmix amr_hgdpmxb_rfmix mxb_rfmix mxb_simulation \
                 amr_simulation amr_hgdpmxb_simulation afr_simulation eur_simulation; do
            m=$(awk 'NF {print $1}' "reference_ids/${l}.txt" | sort -u | comm -23 - "${TMPDIR:-/tmp}/verify_panel_ids.$$" | wc -l)
            [[ "$m" -eq 0 ]] || warn "reference_ids/${l}.txt: ${m} IDs not in the panel (RFMix/simulation would lose them)"
        done
        rm -f "${TMPDIR:-/tmp}/verify_panel_ids.$$"
        if [[ -s "${PANEL_DIR}/provenance/chr21.qc.cohorts.tsv" ]]; then
            info "chr21 per-cohort provenance after QC (raw):"; sed 's/^/          /' "${PANEL_DIR}/provenance/chr21.qc.cohorts.tsv"
        fi
        [[ -s "${PANEL_DIR}/provenance/chr21.fill_check.tsv" ]] && { info "chr21 fill check (raw):"; sed 's/^/          /' "${PANEL_DIR}/provenance/chr21.fill_check.tsv"; }
    fi
fi

section "5. phasing chunks (1b_phase_chunks_clm.sh)"
jobs_of 1b_chunk '1b_chunk_<task>_<job>.out/.err'
if [[ -s "$CHUNKS_TSV" ]]; then
    NEXP=0; NFRESH=0; NOLD=0
    while IFS=$'\t' read -r cn cs ci rest; do
        [[ "$cn" == chr_num ]] && continue
        NEXP=$((NEXP + 1)); ch="${PANEL_DIR}/chunks/${cs}_chunk${ci}.phased.bcf"; su="${PANEL_DIR}/merged_${cs}.softunion.bcf"
        if is_new "$su" && [[ -s "$ch" && "$ch" -nt "$su" ]]; then NFRESH=$((NFRESH + 1)); elif [[ -s "$ch" ]]; then NOLD=$((NOLD + 1)); fi
    done < "$CHUNKS_TSV"
    if [[ "$NFRESH" -eq "$NEXP" ]]; then pass "${NFRESH} of ${NEXP} chunks phased from this rebuild's soft-union"
    else info "${NFRESH} of ${NEXP} chunks phased from this rebuild's soft-union (${NOLD} older chunks on disk will be re-phased)"; fi
else warn "no ${CHUNKS_TSV}"; fi

section "6. ligate (1b_phase_ligate_clm.sh)"
jobs_of 1b_ligate '1b_ligate_chr<N>_<job>.out/.err'
NL=0; NLBAD=0
for c in $(seq 1 22); do
    ph="${PANEL_DIR}/merged_chr${c}.shapeit5_phased.softunion_maf005.bcf"; rc="${ph%.bcf}.rechr.bcf"; su="${PANEL_DIR}/merged_chr${c}.softunion.bcf"
    if is_new "$su" && [[ -s "$ph" && -s "$rc" && "$ph" -nt "$su" ]]; then
        NL=$((NL + 1))
        [[ "$(nsamples "$ph")" -eq "$(nsamples "$su")" ]] || { NLBAD=$((NLBAD + 1)); warn "chr${c}: phased and soft-union sample counts differ"; }
    fi
done
if [[ "$NL" -eq 0 ]]; then info "not yet (no chromosome phased from this rebuild's soft-union)"
else
    [[ "$NL" -eq 22 ]] && pass "22 of 22 chromosomes ligated from the current soft-union (+ .rechr for RFMix)" || info "${NL} of 22 ligated so far"
    [[ "$NLBAD" -eq 0 ]] && pass "phased sample counts match the soft-union" || fail "${NLBAD} chromosome(s) with a different sample count"
fi

section "7. KING de-duplication (1e_king_dedup_clm.sh)"
jobs_of 1e_king_dedup '1e_king_dedup_<job>.out/.err'
DP="${HOMOG_DIR}/king/duplicate_pairs.tsv"; DROP="${REPO_DIR}/reference_ids/panel_duplicate_drop_ids.txt"
P22="${PANEL_DIR}/merged_chr22.shapeit5_phased.softunion_maf005.bcf"
if is_new "$P22" && [[ -s "$DP" && "$DP" -nt "$P22" ]]; then
    pass "duplicate pairs: $(( $(wc -l < "$DP") - 1 ))"
    if [[ -f "$DROP" && "$DROP" -nt "$P22" ]]; then pass "drop list for 4a: $(awk 'NF' "$DROP" | wc -l) IDs (${DROP})"
    else fail "1e wrote duplicate pairs but no current drop list at ${DROP} (4a reads it)"; fi
else info "not yet (no duplicate list from this rebuild's phased panel)"; fi

section "8. ADMIXTURE K=${K} (4a_build_homogeneity_panel_kN_clm.sh)"
jobs_of 4a_homog_kN '4a_homog_kN_<job>.out/.err'
ST="${HOMOG_DIR}/lai_ref_panel_samples.k${K}.tsv"; DT="${HOMOG_DIR}/discovery_clusters.k${K}.tsv"
if is_new "$P22" && [[ -s "$ST" && "$ST" -nt "$P22" ]]; then
    info "discovery clusters (raw):"; sed 's/^/          /' "$DT"
    shopt -s nullglob; ANCH=(reference_ids/*_anchors.k${K}.txt); shopt -u nullglob
    if [[ ${#ANCH[@]} -gt 0 ]]; then info "anchors per super-pop:"; wc -l "${ANCH[@]}" | sed 's/^/          /'
    else fail "4a wrote no reference_ids/*_anchors.k${K}.txt"; fi
    info "samples by super-pop: total / included (raw):"
    awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
                {t[$h["super_pop"]]++; if ($h["included"] == 1) n[$h["super_pop"]]++}
                END {for (s in t) printf "          %-6s %5d total %5d included\n", s, t[s], n[s] + 0}' "$ST" | sort
    NCH=$(grep -c 'CHMI' "$ST"); [[ "$NCH" -eq 0 ]] && pass "no CHMI sample in the table" || fail "CHMI present ${NCH} time(s)"
    NNA=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["super_pop"] == "NA"' "$ST" | wc -l)
    [[ "$NNA" -eq 0 ]] && pass "every sample has a super-pop label" || warn "${NNA} samples with super_pop NA (excluded by 4a)"
    for sp in afr eur eas sas amr men oce; do
        a="reference_ids/${sp}_anchors.k${K}.txt"
        [[ -s "$a" ]] || warn "no anchors for ${sp^^} (${a})"
    done
else info "not yet (no lai_ref_panel_samples.k${K}.tsv from this rebuild's phased panel)"; fi

section "9. simulation, RFMix, FLARE, scoring (publication stages 4-6)"
PROJECT_ROOT="$PROJECT_ROOT" bash scripts/reset-publication-outputs.sh --check | sed 's/^/          /'
[[ ${PIPESTATUS[0]} -eq 0 ]] && pass "publication outputs: nothing from before the rebuild epoch" || fail "publication outputs from before the rebuild are in place (above)"
for j in 2b_simu:'2b_simu_chr<N>_<job>' 3b_rfmix:'3b_rfmix_task<N>_<job>' 3b_homog:'3b_rfmix_task<N>_<job>' \
         3b_chip:'3b_rfmix_task<N>_<job>' 3b_chip_homog:'3b_rfmix_task<N>_<job>' 3b_flare:'3b_flare_task<N>_<job>' \
         finalize_pilot:'finalize_pilot_<job>'; do
    jobs_of "${j%%:*}" "${j#*:}.out/.err"
done
SIMD="${PROJECT_ROOT}/02_simulations/${ADMIX_POP}/gen${GEN}"
[[ -d "$SIMD" ]] && info "simulation outputs: $(find "$SIMD" -name '*.hanc' | wc -l) .hanc truth files in ${SIMD}" || info "simulation: not yet"
for d in "${PROJECT_ROOT}/03_rfmix/${ADMIX_POP}/gen${GEN}" "${PROJECT_ROOT}/03_rfmix/chip_gsa_jessica/${ADMIX_POP}/gen${GEN}"; do
    [[ -d "$d" ]] && info "RFMix: $(find "$d" -name '*.done' | wc -l) finished combos (.done) in ${d}" || info "RFMix: not yet (${d})"
done
FD="${PROJECT_ROOT}/03_flare/${ADMIX_POP}/gen${GEN}"
[[ -d "$FD" ]] && info "FLARE: $(find "$FD" -name '*.anc.vcf.gz' | wc -l) finished tasks in ${FD}" || info "FLARE: not yet"
for t in accuracy_pilot accuracy_pilot_chip_gsa; do
    f="results/${t}/accuracy_summary.tsv"
    if [[ -s "$f" && -f "$EPOCH" && "$f" -nt "$EPOCH" ]]; then pass "${f} written after the rebuild epoch ($(date -r "$f" '+%F %T'))"
    else info "${f}: not yet rewritten by this rebuild"; fi
done

printf '\n== summary: %d PASS, %d WARN, %d FAIL ==\n' "$NP" "$NW" "$NF"
exit $(( NF > 255 ? 255 : NF ))
