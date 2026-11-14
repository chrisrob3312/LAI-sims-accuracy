#!/usr/bin/env bash
# =============================================================================
# preflight-rebuild.sh -- does the cluster have everything the rebuild chain
# needs, before anything is queued? READ-ONLY.
#
# For every stage submit_panel_rebuild.sh queues -- SGDP prep, 1a, 1d, 1c,
# chunks, ligate, 1e, 4a, then 2b, the RFMix/FLARE arms and finalize -- this
# checks the code, the conda env and every tool, every input file, the SGDP
# download, the previous publication outputs, jobs already in the queue, and
# disk space. Each line is PASS / WARN / FAIL with its evidence; the exit
# status is the number of FAILs (0 = ready to submit).
#
# Run from the repo checkout, inside the shapeit5 env:
#   bash scripts/preflight-rebuild.sh
# Paths default to what the stage scripts use; the same env overrides apply.
# Nothing is written; no stderr is hidden.
# =============================================================================
set -uo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$(pwd)}"
PANEL_DIR="${PANEL_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
BRANCH="${BRANCH:-claude/setup-shapeit5-conda-xlcdl}"
SR=/storage/atkinson/shared_resources/reference
HGDP1KG_PATTERN="${HGDP1KG_PATTERN:-${SR}/ReferencePanels/TGP_HGDP_jointcall/original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chr__CHR__.filtered.SNV_INDEL.phased.shapeit5.bcf}"
# none (the rebuild default): 1c keeps relatives. A file: 1c drops its column 2.
OUTLIERS="${OUTLIERS:-none}"
GNOMAD_META="${GNOMAD_META:-${SR}/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
REF_FA="${REF_FA:-${SR}/reference_genomes/b38/Homo_sapiens_assembly38.fasta}"
CHAIN="${CHAIN:-${SR}/genetic_maps/liftover/hg19ToHg38.over.chain.gz}"
GMAP_DIR="${GMAP_DIR:-${SR}/genetic_maps/genetic_maps_shapeit4/genetic_maps_b38}"
MXB_HG19="${MXB_HG19:-${SR}/mexico_biobank/original_download/tosharemxb50wgs/mexican_50_autosomes.vcf.gz}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
CHUNKS_TSV="${CHUNKS_TSV:-${REPO_DIR}/resources/b38_chunks.tsv}"
ALMARRI_DIR="${ALMARRI_DIR:-${HOME}/almarri_vcfs}"
ALMARRI_VCF_DIR="${ALMARRI_VCF_DIR:-${ALMARRI_DIR}/vcfs}"
APPG_META="${ALMARRI_APPG_META:-${ALMARRI_DIR}/meta/appg.meta}"
SGDP_RAW_DIR="${SGDP_RAW_DIR:-${PROJECT_ROOT}/00_sgdp_raw/vcf_variants}"
SGDP_METADATA="${SGDP_METADATA:-${SGDP_RAW_DIR}/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt}"
SGDP_PROBE_META="${PROJECT_ROOT}/04_homogeneity_panel/diag_k8/sgdp_raw_probe/sgdp_meta.txt"
ADMIX_POP="${ADMIX_POP:-Brasa}"
GEN="${GEN:-12}"
DAT_FILE="${DAT_FILE:-${PROJECT_ROOT}/02_simulations/${ADMIX_POP}.dat}"
ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/simu-jointcall/admix-simu-master}"
RFMIX_DIR="${RFMIX_DIR:-/storage/atkinson/shared_resources/past_members/jessica_mauer/lai/RFMix_v1.5.4}"
FIXREF_MODE="${FIXREF_MODE:-flip-all}"
U="${USER:-$(id -un)}"

NP=0; NW=0; NF=0
pass () { NP=$((NP + 1)); printf '  PASS  %s\n' "$*"; }
warn () { NW=$((NW + 1)); printf '  WARN  %s\n' "$*"; }
fail () { NF=$((NF + 1)); printf '  FAIL  %s\n' "$*"; }
info () { printf '  info  %s\n' "$*"; }
section () { printf '\n== %s ==\n' "$*"; }
file_ok () {   # file_ok <label> <path>
    if [[ -s "$2" ]]; then pass "$1: $2"; else fail "$1 missing or empty: $2"; fi
}
data_rows () { awk 'NF && $1 !~ /^#/' "$1" | wc -l; }   # non-empty, non-comment lines

echo "preflight-rebuild  $(date '+%F %T')  host $(hostname)  user ${U}"

section "1. code"
cd "$REPO_DIR" || { fail "cannot cd to ${REPO_DIR}"; exit 99; }
HEAD_SHA=$(git rev-parse HEAD)
info "commit: $(git log -1 --format='%h %ci %s')"
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [[ "$CUR_BRANCH" == "$BRANCH" ]]; then pass "branch ${CUR_BRANCH}"; else fail "branch is ${CUR_BRANCH}, expected ${BRANCH}"; fi
MODS=$(git status --porcelain --untracked-files=no)
if [[ -z "$MODS" ]]; then pass "no local edits to tracked files (queued jobs run exactly the pushed code)"
else warn "local edits to tracked files -- queued jobs will run THESE:"; echo "$MODS" | sed 's/^/          /'; fi
REMOTE_SHA=$(timeout 60 git ls-remote origin "refs/heads/${BRANCH}" | awk 'NR == 1 {print $1}')
if [[ -z "$REMOTE_SHA" ]]; then warn "could not read origin/${BRANCH} (network?)"
elif [[ "$REMOTE_SHA" == "$HEAD_SHA" ]]; then pass "up to date with origin/${BRANCH}"
else fail "checkout ${HEAD_SHA:0:7} differs from origin/${BRANCH} ${REMOTE_SHA:0:7} -- git pull --ff-only first"; fi
MISSING=""
for s in submit_panel_rebuild.sh submit_publication_rerun.sh 1a_prep-mxb-liftover_clm.sh 1b_sgdp_raw_prep_clm.sh \
         1d_almarri_prep_clm.sh 1c_phasing-jointcall_clm.sh 1b_phase_chunks_clm.sh 1b_phase_ligate_clm.sh \
         1e_king_dedup_clm.sh 4a_build_homogeneity_panel_kN_clm.sh 2b_simulation_clm.sh \
         3b_wgs-rfmix-jointcall_clm.sh 3b_homog_wgs-rfmix-jointcall_clm.sh 3b_chip_wgs-rfmix-jointcall_clm.sh \
         3b_chip_homog_wgs-rfmix-jointcall_clm.sh 3b_flare_wgs-jointcall_clm.sh finalize_pilot.sbatch \
         run_accuracy_pilot.sh accuracy_v2.R scripts/fetch-sgdp-vcf-variants.sh scripts/sgdp-raw-sample-map.py \
         scripts/fill-contradictions.py scripts/softunion-sites.py scripts/build-missingness-groups.py \
         scripts/admix-site-filter.py scripts/site-provenance.py scripts/sbatch-when-ready.sh \
         scripts/reset-publication-outputs.sh; do
    [[ -f "$s" ]] || MISSING+=" $s"
done
if [[ -z "$MISSING" ]]; then pass "every stage script and helper is present"; else fail "missing:${MISSING}"; fi

section "2. environment and tools"
if [[ "${CONDA_DEFAULT_ENV:-}" == "shapeit5" ]]; then pass "conda env shapeit5 active"
else fail "conda env is '${CONDA_DEFAULT_ENV:-none}': run conda activate shapeit5 first"; fi
if [[ -n "${PYTHONHOME:-}${PYTHONPATH:-}" ]]; then warn "PYTHONHOME/PYTHONPATH set (unset PYTHONHOME PYTHONPATH)"; else pass "PYTHONHOME/PYTHONPATH unset"; fi
for t in bcftools bgzip tabix plink2 picard SHAPEIT5_phase_common SHAPEIT5_ligate admixture python3 Rscript java perl curl md5sum git; do
    p=$(command -v "$t")
    if [[ -n "$p" ]]; then pass "${t}: ${p}"; else fail "${t} not on PATH"; fi
done
info "$(bcftools --version | awk 'NR == 1')"
PLUGINS=$(bcftools plugin -l)
for pl in fill-tags fixref setGT; do
    if grep -qw -- "$pl" <<<"$PLUGINS"; then pass "bcftools plugin ${pl}"; else fail "bcftools plugin ${pl} missing"; fi
done
# R exactly as run_accuracy_pilot.sh runs it: with R_LIBS_USER on the shapeit5
# env's R library, because a system-wide /opt/biobuilds-2016.04 library
# otherwise shadows the env's tidyverse with a 2016 magrittr R 4 cannot load.
RLU="${R_LIBS_USER:-}"
if [[ -z "$RLU" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib/R/library" ]]; then RLU="${CONDA_PREFIX}/lib/R/library"
    elif [[ -d "$HOME/.conda/envs/shapeit5/lib/R/library" ]]; then RLU="$HOME/.conda/envs/shapeit5/lib/R/library"; fi
fi
RCHK=$(R_LIBS_USER="$RLU" Rscript -e 'for (p in c("dplyr","tidyr","readr","ggplot2","purrr")) suppressMessages(library(p, character.only = TRUE)); cat("OK\n")' 2>&1)
if [[ "$(tail -1 <<<"$RCHK")" == "OK" ]]; then pass "R packages for accuracy_v2.R load (dplyr tidyr readr ggplot2 purrr; R_LIBS_USER=${RLU:-unset})"
else fail "R packages for accuracy_v2.R do not load with R_LIBS_USER=${RLU:-unset}; R said:"; echo "$RCHK" | sed 's/^/          /'; fi
info "java: $(java -version 2>&1 | awk 'NR == 1')"

section "3. reference files"
file_ok "GRCh38 FASTA" "$REF_FA"
file_ok "FASTA index" "${REF_FA}.fai"
file_ok "Picard dictionary" "${REF_FA%.*}.dict"
file_ok "hg19->hg38 chain" "$CHAIN"
NG=0; for c in $(seq 1 22); do [[ -s "${GMAP_DIR}/chr${c}.b38.gmap.gz" ]] && NG=$((NG + 1)); done
if [[ "$NG" -eq 22 ]]; then pass "genetic maps: 22 of 22 in ${GMAP_DIR}"; else fail "genetic maps: ${NG} of 22 in ${GMAP_DIR}"; fi
if [[ -s "$CHUNKS_TSV" ]]; then
    NC=$(( $(wc -l < "$CHUNKS_TSV") - 1 ))
    if [[ "$NC" -eq 594 ]]; then pass "phasing chunks: ${NC} rows in ${CHUNKS_TSV}"; else warn "phasing chunks: ${NC} rows (594 expected) in ${CHUNKS_TSV}"; fi
else fail "missing ${CHUNKS_TSV} (python3 resources/build_b38_chunks.py)"; fi

section "4. HGDP+1KG (unfiltered Koenig v2)"
NH=0; NHI=0
for c in $(seq 1 22); do
    f="${HGDP1KG_PATTERN//__CHR__/$c}"
    [[ -s "$f" ]] && NH=$((NH + 1))
    [[ -s "${f}.csi" || -s "${f}.tbi" ]] && NHI=$((NHI + 1))
done
if [[ "$NH" -eq 22 && "$NHI" -eq 22 ]]; then pass "22 of 22 chromosome files, all indexed"; else fail "${NH} of 22 files, ${NHI} indexed (${HGDP1KG_PATTERN})"; fi
F22="${HGDP1KG_PATTERN//__CHR__/22}"
if [[ -s "$F22" ]]; then
    N22=$(bcftools query -l "$F22" | wc -l)
    if [[ "$N22" -eq 4091 ]]; then pass "chr22: ${N22} samples"; else warn "chr22: ${N22} samples (4,091 expected)"; fi
fi
if [[ "$OUTLIERS" == none ]]; then pass "relatedness drop: none (OUTLIERS=none -- relatives kept, 1e drops duplicates only)"
elif [[ -s "$OUTLIERS" ]]; then pass "relatedness drop list ${OUTLIERS}: $(wc -l < "$OUTLIERS") lines"
else fail "OUTLIERS=${OUTLIERS} is missing (none keeps every sample)"; fi
file_ok "gnomAD metadata (4a)" "$GNOMAD_META"
if [[ -s "$SAMPLE_GROUPS" ]]; then pass "sample_groups.tsv: $(wc -l < "$SAMPLE_GROUPS") samples"; else fail "missing ${SAMPLE_GROUPS} (make_sample_groups.sh)"; fi

section "5. MXB"
file_ok "MXB hg19 VCF" "$MXB_HG19"
file_ok "MXB index" "${MXB_HG19}.tbi"
STAMP_MXB="${PANEL_DIR}/mxb_lifted.recipe"
if [[ -s "$STAMP_MXB" ]] && grep -q "fixref=${FIXREF_MODE}" "$STAMP_MXB"; then info "MXB already lifted with fixref=${FIXREF_MODE}: 1a will be skipped"
elif [[ -s "$STAMP_MXB" ]]; then info "MXB lift stamp '$(cat "$STAMP_MXB")' -> the driver re-lifts MXB (1a)"
else info "no MXB lift stamp (lifted before stamps, with -m flip -d) -> the driver re-lifts MXB (1a)"; fi

section "6. Almarri"
if [[ -d "$ALMARRI_VCF_DIR" ]]; then
    NA=$(find "$ALMARRI_VCF_DIR" -maxdepth 1 -name '*.vcf.gz' | wc -l)
    NAI=$(find "$ALMARRI_VCF_DIR" -maxdepth 1 \( -name '*.vcf.gz.tbi' -o -name '*.vcf.gz.csi' \) | wc -l)
    if [[ "$NA" -eq 137 && "$NAI" -ge 137 ]]; then pass "${NA} per-sample VCFs, ${NAI} indexes in ${ALMARRI_VCF_DIR}"
    else fail "${NA} per-sample VCFs (137 expected), ${NAI} indexes in ${ALMARRI_VCF_DIR}"; fi
else fail "no Almarri VCF directory ${ALMARRI_VCF_DIR}"; fi
file_ok "appg.meta" "$APPG_META"
if [[ -s reference_ids/almarri_manifest.tsv ]]; then
    NM=$(( $(wc -l < reference_ids/almarri_manifest.tsv) - 1 ))
    if [[ "$NM" -eq 137 ]]; then pass "almarri_manifest.tsv: ${NM} samples"; else warn "almarri_manifest.tsv: ${NM} samples (137 expected)"; fi
else warn "no reference_ids/almarri_manifest.tsv: 4a labels Almarri without the country"; fi

section "7. SGDP (raw per-sample calls)"
if [[ -s reference_ids/sgdp_manifest.tsv ]]; then pass "sgdp_manifest.tsv: $(( $(wc -l < reference_ids/sgdp_manifest.tsv) - 1 )) samples"; else fail "missing reference_ids/sgdp_manifest.tsv"; fi
if [[ -f "${SGDP_RAW_DIR}/vcfs/DONE" ]]; then
    pass "download verified: $(tr '\n' ';' < "${SGDP_RAW_DIR}/vcfs/DONE")"
    META="$SGDP_METADATA"; [[ -s "$META" ]] || META="$SGDP_PROBE_META"
    if [[ -s "$META" ]]; then
        MAP=$(mktemp)
        MAPLOG=$(python3 -I scripts/sgdp-raw-sample-map.py --manifest reference_ids/sgdp_manifest.tsv \
                     --metadata "$META" --vcf-dir "${SGDP_RAW_DIR}/vcfs" --out "$MAP" 2>&1)
        NMAP=$(wc -l < "$MAP"); rm -f "$MAP"
        echo "$MAPLOG" | sed 's/^/          /'
        if [[ "$NMAP" -ge 277 ]]; then pass "manifest samples with a per-sample VCF: ${NMAP}"
        elif [[ "$NMAP" -gt 0 ]]; then warn "manifest samples with a per-sample VCF: ${NMAP} of 277 (listed above)"
        else fail "no manifest sample maps to a per-sample VCF"; fi
    else info "no SGDP metadata copy yet (the prep downloads it); sample mapping not checked here"; fi
else
    FJ=$(squeue -h -u "$U" -n sgdp_fetch -o '%i %T %M')
    if [[ -n "$FJ" ]]; then pass "download not finished; job ${FJ} -- the SGDP prep will wait on it"
    else
        LASTLOG=$(ls -t "${PROJECT_ROOT}"/logs/sgdp_fetch_*.out | awk 'NR == 1')
        fail "download not finished and no sgdp_fetch job queued. Newest log: ${LASTLOG:-none}"
        [[ -n "$LASTLOG" ]] && tail -5 "$LASTLOG" | sed 's/^/          /'
    fi
fi

section "8. reference lists in the repo"
if [[ -s reference_ids/panel_exclude_ids.txt ]]; then pass "panel_exclude_ids.txt: $(data_rows reference_ids/panel_exclude_ids.txt) IDs"; else fail "missing reference_ids/panel_exclude_ids.txt"; fi
for f in known_duplicate_pairs.tsv MXB50genomes_popinfo.tsv sgdp_server_panel_ids.txt; do file_ok "$f" "reference_ids/$f"; done
LMISS=""
for l in amr_rfmix afr_rfmix eur_rfmix pel_rfmix pel_eas_rfmix amr_hgdpmxb_rfmix mxb_rfmix mxb_simulation \
         amr_simulation amr_hgdpmxb_simulation afr_simulation eur_simulation; do
    [[ -s "reference_ids/${l}.txt" ]] || LMISS+=" ${l}.txt"
done
if [[ -z "$LMISS" ]]; then pass "12 RFMix panel / simulation donor lists present"; else fail "missing lists:${LMISS}"; fi
NCHIP=$(ls reference_ids/chip_gsa_jessica/ | grep -c '^chr[0-9]*\.snps$')
if [[ "$NCHIP" -eq 22 ]]; then pass "GSA chip site lists: 22 of 22"; else fail "GSA chip site lists: ${NCHIP} of 22"; fi

section "9. simulation, RFMix, FLARE, scoring"
file_ok "admix-simu model (.dat)" "$DAT_FILE"
for x in simu-mix.pl bp2anc.pl; do file_ok "admix-simu ${x}" "${ADMIXSIMU_DIR}/${x}"; done
file_ok "RFMix v1.5.4 RunRFMix.py" "${RFMIX_DIR}/RunRFMix.py"
FJAR=""
for p in "${FLARE_JAR:-}" "${PROJECT_ROOT}/software/flare.jar" "/storage/atkinson/shared_resources/software/flare/flare.jar" \
         "/storage/atkinson/shared_resources/software/flare.jar" "${HOME}/software/flare.jar" "${HOME}/flare.jar"; do
    if [[ -n "$p" && -s "$p" ]]; then FJAR="$p"; break; fi
done
if [[ -n "$FJAR" ]]; then pass "FLARE jar: ${FJAR}"; else fail "no flare.jar in the places 3b_flare looks (set FLARE_JAR)"; fi
if [[ -d results ]]; then pass "results/ present (finalize writes the accuracy TSVs there)"; else fail "results/ missing"; fi
RCHECK=$(PROJECT_ROOT="$PROJECT_ROOT" bash scripts/reset-publication-outputs.sh --check 2>&1)
if [[ $? -eq 0 ]]; then pass "previous simulation/RFMix/FLARE outputs are set aside (rebuild epoch present)"
else fail "previous simulation/RFMix/FLARE outputs still in place -- the arms would reuse them:"; echo "$RCHECK" | sed 's/^/          /'; fi
PUSH=$(timeout 60 git push --dry-run origin "HEAD:refs/heads/${BRANCH}" 2>&1)
if [[ $? -eq 0 ]]; then pass "git push access to origin/${BRANCH} (finalize pushes its TSVs)"
else warn "git push --dry-run failed -- finalize could not push its TSVs:"; echo "$PUSH" | sed 's/^/          /'; fi

section "10. queue and disk"
QJOBS=$(squeue -h -u "$U" -o '%i %j %T' \
    | awk '$2 ~ /^(1a_prep_mxb_liftover|1b_sgdp_raw|1d_almarri|1c_phase|1b_chunk|1b_ligate|1e_king_dedup|4a_homog_kN|2b_simu|3b_rfmix|3b_homog|3b_chip|3b_chip_homog|3b_flare|finalize_pilot)$/')
if [[ -z "$QJOBS" ]]; then pass "no rebuild-stage job already queued or running"
else fail "rebuild-stage jobs already queued/running -- a second chain would write the same files:"; echo "$QJOBS" | sed 's/^/          /'; fi
DF=$(df -Pk "$PROJECT_ROOT" | awk 'NR == 2 {print $4}')
info "$(df -h "$PROJECT_ROOT" | awk 'NR == 2 {print "disk: " $4 " free of " $2 " on " $6}')"
if [[ "$DF" -gt $((2 * 1024 * 1024 * 1024)) ]]; then pass "more than 2 TB free"; else warn "less than 2 TB free: 1c's per-chromosome intermediates are large"; fi

printf '\n== summary: %d PASS, %d WARN, %d FAIL ==\n' "$NP" "$NW" "$NF"
if [[ "$NF" -eq 0 ]]; then echo "ready: next, the dry run (DRY_RUN=1 bash submit_panel_rebuild.sh)"; else echo "NOT ready: fix every FAIL above first"; fi
exit $(( NF > 255 ? 255 : NF ))
