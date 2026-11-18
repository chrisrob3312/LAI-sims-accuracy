# Benchmark cohort & reference panel — build workflow

Reference for how the panel and the benchmarking cohort are built, what every
script does, and what each check catches. Companion to
`BENCHMARK_COHORT_PLAN.md` (who is in the cohort) and
`BENCHMARK_PANEL_DESIGN.md` (who may be a donor vs a reference).

Everything runs on **taco (mhgcp)**. Branch: `claude/setup-shapeit5-conda-xlcdl`.

---

## 0. Pipeline at a glance (copy-paste)

```text
LAI-sims-accuracy -- data path, script by script
One command submits A-F as one afterok chain:   bash submit_panel_rebuild.sh
Every stage skips work that is already current, so re-running it is safe.

A. SOURCES
  - HGDP+1KG  Koenig et al. 2024 phased_haplotypes_v2, UNFILTERED: 4,091 QC'd
              genomes, SNVs + indels, no MAF cut (GRCh38)
  - MXB       50 WGS, GATK HaplotypeCaller -> GenotypeGVCFs joint call, PASS (hg19)
  - SGDP      279 public per-sample WGS calls (vcf_variants): GATK
              UnifiedGenotyper, variant sites only, filter level FL 0-9 (hs37d5)
  - Almarri   137 per-sample Long Ranger WGS VCFs (GRCh38)

B. COHORT PREP (in parallel; FIXREF_MODE=flip-all in every prep)
  scripts/fetch-sgdp-vcf-variants.sh            SGDP download, once
    - resumable download of the 57 GB tar; exact byte count required
    - only bare-named regular files extracted, into a directory of their own
    - md5 of every file vs the .md5sum shipped beside it -> vcfs/DONE
  1b_sgdp_raw_prep_clm.sh                        SGDP, array 1-22
    - manifest sample -> SGDP_ID -> Illumina_ID -> per-sample VCF
    - bcftools merge -0 -m none -r N: no record in a sample = 0/0
    - variant call with FL < 1 -> ./. (SGDP's recommended level)
    - site where most SGDP carriers fail FL -> all of SGDP ./. there
    - GT only; rename to manifest IDs (S_/B_ ... .DG)
    - Picard LiftoverVcf hg19 -> hg38 (RECOVER_SWAPPED_REF_ALT)
    - norm -m- -f hg38 -> +fixref -m flip-all -d -> chrN -> sort
      -> sgdp_prepped.chrN.bcf + .recipe stamp
  1a_prep-mxb-liftover_clm.sh                    MXB, one job (only when its stamp is stale)
    - contigs 1..22 -> chr1..chr22; Picard LiftoverVcf hg19 -> hg38
    - norm -m- -f hg38 -> +fixref -m flip-all -d -> sort
      -> mxb_lifted.chrN.bcf + mxb_lifted.recipe
  1d_almarri_prep_clm.sh                         Almarri, array 1-22
    - per sample: low-quality call -> ./. (FILTER not PASS, GQ <= 20 or
      missing) -> norm -c w -m -snps -> biallelic SNPs -> +fixref -m flip-all -d
    - merge -m none -0 (absent = 0/0; a low-quality call stays ./.)
      -> site mask (more ./. than ALT: all of Almarri ./.) -> GT only
      -> almarri_prepped.chrN.bcf + .recipe (v4)

C. MERGE + SITE SELECTION   1c_phasing-jointcall_clm.sh (PHASE_MODE=chunked, array 1-22)
    - skip when the published soft-union is newer than every input and stamped
      with the same settings and code
    - HGDP+1KG: drop panel_exclude_ids.txt (Koenig's 60 QC removals, CHMI);
      relatives KEPT (OUTLIERS=none; July/November dropped related_outliers.txt);
      chr-prefix; norm -m -any -c w -> +fixref -m flip-all -d
    - add SGDP, Almarri (drop IDs already present or excluded)
    - refuse unless every cohort's prep stamp records the same fixref mode
    - split any multi-allelic record left in any input (norm -m -any)
    - site provenance per cohort (records, unique, absent / swap / flip / ALT)
    - 0/0 fill for HGDP1KG MXB SGDP ALMARRI: a record absent from a cohort's
      file is 0/0 for that cohort (union + placeholder sample + merge -0)
    - bcftools merge -m none (column-wise sample join)
    - plink2: biallelic ACGT SNPs, dedup
    - 3a fill check: drop a position where a filled cohort lacks an allele
      another cohort of the same super-pop carries (>= 5 copies,
      P(no copy) < 0.001, groups >= 20)
    - 3b: drop positions holding > 1 record (split multi-allelic sites)
    - soft-union: keep a site if MAF >= 0.005 (AN > 0) in ANY super-pop, every
      cohort's samples pooled by super-pop (groups >= 10)
    - ADMIXTURE site list: every cohort x super-pop group <= 15% missing and the
      panel <= 7% -> merged_chrN.admix_sites.bed
    - publish merged_chrN.softunion.bcf + .recipe stamp

D. JOINT PHASING
  1b_phase_chunks_clm.sh      594 chunks (8 cM + 2 cM overlap); 209, 302 in their own array
    - SHAPEIT5_phase_common on all cohorts together, --filter-maf 0.001;
      imputes any remaining ./.
    - a chunk is reused only if newer than the soft-union, same sample count
  1b_phase_ligate_clm.sh      array 1-22
    - SHAPEIT5_ligate -> merged_chrN.shapeit5_phased.softunion_maf005.bcf
      (+ .rechr.bcf with bare contigs for RFMix)

E. DE-DUPLICATION + ANCESTRY PANEL
  1e_king_dedup_clm.sh
    - concat -> 1c's called sites -> LD-prune -> plink2 KING-robust
    - duplicates (> 0.354) + known same-person pairs; keep MXB > ALMARRI >
      HGDP > 1KG > SGDP -> panel_duplicate_drop_ids.txt
  4a_build_homogeneity_panel_kN_clm.sh       K=7 (from 17 Nov)
    - concat -> 4a's OWN ADMIXTURE site list from 1c's pre-phasing soft-union
      (every super-pop <= 7% missing, every cohort overall <= 15%;
      scripts/admix-site-list.sh, 22 chromosomes in parallel, stamped)
    - fit set: unrelated (1e duplicate copies out, then one of each pair at
      KING >= 0.0884); relatives projected (admixture -P) in both runs below
    - LD prune 500 50 0.1 on the fit set only, after a global MAF >= 0.005
      floor on that set (18 Nov; sites were already MAF >= 0.005 in some
      super-pop); stops BEFORE the fits if the marker count implies more than
      ADMIX_MAX_HOURS (60) of ADMIXTURE
    - unsupervised ADMIXTURE K=7; PCA of the fit set (top 10 PCs)
    - data-first discovery: each cluster named from its top 50 at Q >= 0.94;
      anchors = 50 with the best rank of Q + distance to the group's centre
      in PC space, from candidates at Q >= 0.94 carrying the cluster's label
    - supervised K=7 (.pop) -> homogeneous at Q >= 0.94, every super-pop ->
      lai_ref_panel_samples.k7.tsv; amr_homog_rfmix.txt for the 3-way HOMOG
      arms = July's K=5 AMR panel + AMR samples new since July that pass now,
      minus donors and their relatives (scripts/build-3way-amr-homog.py)

F. LAI SIMULATION + ACCURACY   publication stages 4-6, chained to 4a
  2b_simulation_clm.sh        array 1-22
    - admix-simu donors from the phased panel (plink2 --export haps): NAT
      (HGDP) and NATMXB (HGDP + MXB) tracks; EUR = IBS, AFR = YRI
    - simulated admixed haplotypes + truth .bp / .hanc per track
  3b_{wgs,homog_wgs,chip_wgs,chip_homog_wgs}-rfmix + 3b_flare_wgs (x4)
    - RFMix v1.5.4 AND FLARE on the same grid: reference panels x 2 sim
      tracks (NAT = Brazilian-like, NATMXB = Mexican-like) x 22 chr; WGS and
      chip-density (GSA) arms, main and HOMOG panels -- 8 arrays (17 Nov:
      FLARE WGS HOMOG, chip main and chip HOMOG added; FLARE's HGDPMXB_FULL
      NAT list now built like 3b's, no longer from the panel5 keep file,
      which also held the IBS/YRI IDs and failed every such task)
    - references: 1KG IBS (EUR), 1KG YRI (AFR), NAT lists, amr_homog_rfmix
      (see reference_ids/rerun_jessicaH_3way_lai_sim/README.md)
  finalize_pilot.sbatch
    - checks no broken .classes, scores RFMix WGS + chip and FLARE WGS + chip
      with the same scorer (accuracy_v2.R; scripts/flare-to-lat3.py puts
      FLARE's .anc.vcf.gz in the .Lat3/.classes/.done layout, ancestry codes
      mapped by name from its ##ANCESTRY header), asserts 22/22 chromosomes in
      every table, commits + pushes results/accuracy_pilot{,_chip_gsa,_flare,
      _flare_chip_gsa}; slides/tables regenerate on the main-loop side

G. BENCHMARKING COHORT   after 4a, run by hand (section 13)
  scripts/split-homog-donors-reference.py   twice, distinct suffixes
    - --mode lai (suffix lai): any sample may be donor or reference
    - --mode benchmark (suffix bench): 1KG never a donor, OCE and SAS reference-only,
      pool = 4a's `homogeneous` (LAI donors included; `included` is the LAI panel),
      stratified by population
  scripts/assemble-benchmark-cohort.py
    - excludes donors, reference, 1KG, SGDP-resequenced 1KG, KING drops
    - sample_type from the ADMIXTURE argmax; design paired-builds
  6b_simulate_benchmark_cohort_clm.sh   13 models, 2,160 simulated genomes + truth
  6c_extract_real_test_genomes_clm.sh   real test genomes per chromosome
  6d_write_chip_cells_clm.sh            chip x build x chr cells (native hg19
                                        manifest; swap / flip handled)
  6e_cohort_kinship_clm.sh              genome-wide kinship -> N_eff
  6a_transfer_to_mab_clm.sh             verify, then copy panels + cohort to mab
```

**Running it** (all from the repo checkout, in the shapeit5 env):

| when | command | writes? |
|---|---|---|
| before submitting | `bash scripts/preflight-rebuild.sh` — every tool, input and setting the chain needs; exit = number of FAILs | no |
| once, before the first chained run | `bash scripts/reset-publication-outputs.sh` — moves the previous `02_simulations` / `03_rfmix` / `03_flare` trees of every cohort (`COHORTS`, default `Brasa:12:NAT MexLike:16:NATMXB`) aside and writes `rebuild_epoch`. The RFMix and FLARE arms reuse any output they find, whatever its age, so without this the old pilot would be skipped over as done and re-scored; the driver refuses to chain them until `--check` passes | moves (never deletes) |
| plan | `DRY_RUN=1 bash submit_panel_rebuild.sh` | no |
| submit / re-submit after any failure | `bash submit_panel_rebuild.sh` | queues jobs |
| after 1c | `sbatch --array=1-22 scripts/cohort-site-report.sh` (or `bash scripts/cohort-site-report.sh 21 22` in the env) — per soft-union record and cohort: where the cohort's own file stood (present / absent → filled 0/0 / swapped / flipped / other allele), its missing genotypes before phasing, its allele frequency; per-1-Mb bins; → `01_merged_phased_panel/provenance/chrN.site_report.{tsv.gz,bins.tsv}` | no |
| after 1c (again after the site report) | `bash scripts/merge-sites-collect.sh` — gathers 1c's per-chromosome site tables (per stage and cohort: records, unique, why a cohort lacks a panel record; which cohort combinations carry each record; the fill check by super-pop; soft-union and ADMIXTURE-site rules; per-chromosome totals) and, when the site report exists, per-cohort summaries of it (present / absent / allele mismatch, missingness, pairwise allele-frequency histograms, 1-Mb bins) into `results/merge_sites/*.tsv` + `MANIFEST.tsv`, to commit and plot off the cluster. It also writes `variant_accounting{,_genome,_loglines}.tsv` (scripts/variant-accounting.py): variants kept and lost at every step, per chromosome and genome-wide -- each cohort's release and prep (liftover rejects, masked sites), what each file brings to the merge, the union, plink2 site QC, the two position drops, the soft-union, the ADMIXTURE list, phased and RFMix input, 4a's candidates and LD-pruned SNPs -- every number with its source (a file's record count, a 1c table, or the stage's own log line). Seconds without the site report; `sbatch` it after the site report (~27M records) | `results/merge_sites/` only |
| after 4a | `python3 -I scripts/homog-by-population.py` — 4a's homogeneous panel per super-pop and population: n, in the panel, how many would pass at the cut given here (`--cut SAS=0.95` on a table built at 0.93), how many have no first-degree relative, why the rest are not in (donor or donor's relative / duplicate / strongest component another super-pop / below the cut), strongest-component counts, mean supervised Q per component, anchors per population; each cluster's anchors by cohort/population; for every population with simulation donors: relative-free count, donors with a relative, relatives they keep out | no |
| after 4a | `python3 -I scripts/homog-cut-preview.py [--k 8] [--cuts 0.90,0.92,0.94]` — per super-pop and cut: the benchmark pool split into donor-eligible / reference-only (1KG, server-panel ids), how many of them are LAI simulation donors or relatives (in the benchmark pool, out of the LAI panel), the mean `1 − share` of each donor pool, the samples each lower cut adds by population with their strongest other component, how many sit under the lowest cut, and how many MXB fall under each cut; the table's `homogeneous` count beside them as a check | no |
| after 4a | `REFS=$S RESULTS=$S ./run_benchmark_cohort.sh --to preflight` with `S` a new scratch folder — the benchmark's numbers from this table: the split (reference / donors per super-pop), real test genomes (homogeneous / admixed per component), the simulated roster, donor reuse, and the intake's rules; nothing submitted, nothing written in `reference_ids/` or `results/` | scratch only |
| after 4a | `python3 -I scripts/donor-switch-plan.py --table <lai_ref_panel_samples.k7.tsv, or .k8> --kin <king.kin0> --refs reference_ids --out-dir results/donor_plan` — per donor list and population, the same-size donor set from the same population that keeps the fewest relatives out of the reference (relative-free first, then both parents of a trio), never two related donors, never a curated-reference member or relative (a current donor who is one MUST be replaced: 2b refuses it), keeping as many current donors as possible; writes `<list>.proposed.txt`, `donor_plan_summary.tsv`, `donor_plan_changes.tsv`. Proposals only: copying one into `reference_ids/` is a separate step, and 4a must rerun before 2b | no |
| after 4a (4a runs it) | `python3 -I scripts/build-3way-amr-homog.py --base reference_ids/lai_ref_panel_samples.tsv --table <lai_ref_panel_samples.k7.tsv> --exclude <admixture_k7/simulation_donors.all.txt> --out <list>` -- the AMR reference of the 3-way HOMOG arms: July's K=5 AMR homogeneous samples still in the panel, plus AMR samples new since July (SGDP, restored HGDP, 1KG relatives) homogeneous in this run, minus donors, their relatives, duplicates and exclusions; reports what was kept, dropped and added | writes only the list |
| inside 4a | `bash scripts/admix-site-list.sh <chr> <out.bed>` -- one chromosome's ADMIXTURE site list under the 17 Nov rule (every super-pop <= 7%, every cohort <= 15%), from 1c's pre-phasing soft-union, groups rebuilt as 1c built them; `--groups` writes only the groups | writes only its outputs |
| finalize runs it | `python3 -I scripts/flare-to-lat3.py --flare-dir <03_flare[_gsa_jessica]/Brasa/gen12> --out-dir <same>/lat3 --gen 12 --jobs 4` -- FLARE's calls in the layout accuracy_v2.R reads (`.Lat3` per site and haplotype, `.classes` all 0, `.done`); ancestry codes mapped by name from the `##ANCESTRY` header; a missing header or an ancestry other than NAT/EUR/AFR stops it | writes only `lat3/` |
| any time after 1c | `bash scripts/admix-rule-compare.sh` — the ADMIXTURE site list under five versions of the missingness rule (`sp07_coh15`: every super-pop ≤ 7%, every cohort ≤ 15% -- the 17 Nov rule; `coded`: every cohort × super-pop ≤ 15%, panel ≤ 7%; `sp07_tot15`: every super-pop ≤ 7%, panel ≤ 15%; `grp07`: every cohort × super-pop ≤ 7%, panel ≤ 7%; `sp15_tot07`: every super-pop ≤ 15%, panel ≤ 7%), from the published pre-phasing soft-union with 1c's own groups; `coded` must come out IDENTICAL to 1c's list. `CHRS="$(seq 1 22)" RULES=<rule>` writes all 22 lists, which 4a reads with `ADMIX_SITES_TPL=<out>/<rule>/merged_chr%s.admix_sites.bed` — a rule change without re-running 1c | writes only its own folder |
| any time after a stage ends | `SWEEP=1 bash scripts/failure-report.sh` — the failure report plus section 6: every task's `.err` and `.out` of the current run, completed ones too, lines grouped across tasks with numbers made generic -- those with a warning/error word (htslib `[W::]`/`[E::]` included), every other `.err` line, and `.out` lines only a few tasks print; one raw example each. Screen shows `SWEEP_CAP` (15) groups per stage, the report file all | report file only |
| after any failure | `bash scripts/failure-report.sh` — per stage the current run's task states (a pending array range counted as its tasks; any job name no stage script sets, e.g. an interactive `bash`, marked "not a pipeline job"); every task that did not complete (state, exit code, run time vs limit, MaxRSS vs memory asked, node, log paths); failures grouped by stage/state/exit code with their nodes; the full `.err`/`.out` of the first failure of each kind; what is still queued and why. Whole report, logs uncut, in `logs/failure_report_<date>.txt` | report file only |
| any time | `bash scripts/verify-rebuild.sh` — per stage: jobs queued and finished, outputs, and the checks (stamps, sample counts, A/T + C/G kept, fill check, reference-list IDs in the panel and why any are missing, donor/reference relatives from KING, anchors, no CHMI, results rewritten) | no |

**Re-chaining after a stage fails.** Dependencies cannot be edited in place on
this cluster: its job-submit plugin refuses every `scontrol update` that does
not specify memory (`Memory not specified (--mem or --mem-per-cpu required)`),
including `Dependency=` (12 Nov 2026, job 3754278). So after a failed stage:
fix it, `scancel` the jobs still waiting behind it (they sit in
`DependencyNeverSatisfied`), and re-run `bash submit_panel_rebuild.sh`. Every
stage that already finished is reused from its stamp in seconds, so only the
failed stage and what follows it do real work.

---

## 1. Two products from one source

The same joint-called panel produces two deliverables with **opposite donor
policies**, and nothing inside a BCF says which policy produced it:

| product | donors | for |
|---|---|---|
| **LAI / GWAS reference panel** | **included** | FLARE, the GWAS, power simulations |
| **Benchmark reference panel** | **excluded** | scoring the benchmark cohort's ancestry calls |
| **Benchmark cohort** | — | the test genomes themselves (real + simulated) |

Neither panel is an **imputation** reference. Imputation runs against the
external servers — **TOPMed, All of Us, Michigan (1000 Genomes)**. The panels
here are for supervised ADMIXTURE and local-ancestry inference only.

---

## 2. Stage map, in run order

| stage | script | what it produces |
|---|---|---|
| 0 | `scripts/fetch-sgdp-vcf-variants.sh` | SGDP per-sample WGS calls, downloaded + md5-verified |
| 1a | `1a_prep-mxb-liftover_clm.sh` | MXB lifted to GRCh38 (+ `mxb_lifted.recipe`) |
| 1b | `1b_sgdp_raw_prep_clm.sh` | SGDP from the raw calls, lifted to GRCh38 (`1b_sgdp_liftover_clm.sh` = the imputed phased_data2021 release, superseded) |
| 1d | `1d_almarri_prep_clm.sh` | 137 Almarri per-chr BCFs |
| 1c | `PHASE_MODE=chunked` `1c_phasing-jointcall_clm.sh` | 4-cohort merge + soft-union sites (**3,921 samples** in November) |
| 1b | `1b_phase_chunks_clm.sh` → `1b_phase_ligate_clm.sh` | phased panel (594 chunks, 8 cM / 2 cM overlap) |

`submit_panel_rebuild.sh` submits the whole rebuild as one afterok chain:
[download job if still running] → SGDP prep, 1d, and 1a when MXB's stamp is
stale → 1c → chunks → ligate → 1e → 4a → publication stages 4–6
(`CHAIN_PUBLICATION=1`, default). Every stage skips work that is already
current (stamps and newer-than checks), so the same command can be re-submitted
after any failure. `scripts/check-sgdp-raw.sh` reports the SGDP download, the
sample mapping and the prepped chromosomes (read-only).

`1b_phase_prep_clm.sh` (and `submit_chunked_phasing.sh`, which submits it) is
the **2-cohort** July prep — HGDP+1KG + MXB only — and writes the same
`merged_chrN.softunion.bcf` path 1c publishes. It now refuses to run unless
`ALLOW_TWO_COHORT_PREP=1`.
| 1e | `1e_king_dedup_clm.sh` | duplicate drops (KING > 0.354 on called sites, + known same-person pairs) |
| 4a | `4a_build_homogeneity_panel_kN_clm.sh` | **ADMIXTURE K=7 Q table** (K=8 until 16 Nov) — everything downstream reads this |
| 5 | `scripts/split-homog-donors-reference.py` | donor / reference lists, per mode |
| 7+8 | `scripts/assemble-benchmark-cohort.py` | cohort manifests + chip/build assignment |
| 6 | `6b_simulate_benchmark_cohort_clm.sh` | the 2,160 simulated genomes + truth tracts |
| 7 | `6c_extract_real_test_genomes_clm.sh` | real test genomes, per-chr BCF |
| 9 | `scripts/build-chip-manifest-sites.py` | chip site lists, both builds |
| 10 | `6d_write_chip_cells_clm.sh` | chip × build × chr cell files |
| 10b | `6e_cohort_kinship_clm.sh` | genome-wide kinship → effective N |
| transfer | `6a_transfer_to_mab_clm.sh` | everything to mab |

Status helper: `bash scripts/pipeline-status.sh` — merges `sacct` + `squeue`
(deduplicated; a running array task is reported by both) and runs the
checkpoints.

---

## 3. Panel creation

### 3.0 Data path, script by script (what happens to the genotypes)

Top to bottom is the order the data moves. **NEW** marks each step the July
NAT accuracy run did not have. Numbers are chr22 from the November run
(1c job 3742123).

**Source files — one per cohort, each from a different calling pipeline**

| cohort | n | source | how it was called |
|---|---|---|---|
| HGDP+1KG | 3,398 (July) | **Default now:** `original_data_downloads/phased_haplotypes_v2/hgdp1kgp_chrN.filtered.SNV_INDEL.phased.shapeit5.bcf` — Koenig et al. 2024, **unfiltered** (4,091 QC'd samples; SNVs + indels; their site QC HWE / F_MISSING ≤ 0.1 / ExcHet; no MAF cut; chr21: 923,540 SNPs) | Joint call. July and November used `processed_data/phased_haplotypes_v2_filter1` = this release cut to biallelic SNPs at a **global** MAF ≥ 0.005 (chr21: 201,120 of 923,540 SNPs kept), which removed variants common in one super-pop but rare overall before the per-super-pop soft-union could see them. November's u235158 archive (4,151 samples) added Koenig's 60 QC removals; they are excluded by name |
| MXB | 50 | `mexican_50_autosomes.vcf.gz` | GATK HaplotypeCaller → GenotypeGVCFs **joint call**, PASS, singletons kept; hg19 → lifted by 1a. A site it lacks = no MXB carrier |
| SGDP | 276 | Reich lab `phased_data2021/chr.sgdp.pub.N.bcf` | **GLIMPSE-imputed against 1KG phase 3** from mpileup (README: "GT: Phased and imputed genotypes"), so only 1KG-p3 sites; hg19 → lifted by 1b_sgdp. **Being replaced** by the raw per-sample WGS calls, `vcf_variants/` (279 public samples, GATK UnifiedGenotyper, variant sites only, per-call filter level FL, hs37d5; one 57 GB tar): `scripts/fetch-sgdp-vcf-variants.sh` fetches and md5-checks it; the prep (FL filter → merge with absent = 0/0 like Almarri → liftover like 1a) follows. The rebuild driver refuses the imputed release unless `ALLOW_IMPUTED_SGDP=1` |
| Almarri | 137 | Long Ranger per-sample phased VCFs | 10x linked reads, per-sample GATK, GRCh38 — **NEW** |

**1a — `1a_prep-mxb-liftover_clm.sh`** (MXB; unchanged since July)
1. contigs `1..22` → `chr1..chr22` (the UCSC chain needs the prefix)
2. Picard `LiftoverVcf` hg19 → hg38, `RECOVER_SWAPPED_REF_ALT=true` (where
   hg38's base is the old ALT, swaps REF/ALT and the genotypes); unliftable →
   `mxb_lifted.rejected.vcf.gz`
3. `bcftools norm -f hg38 -m-` (split multi-allelics, left-align)
4. `bcftools +fixref -- -f hg38 -m ${FIXREF_MODE} -d` — **`flip-all` now;
   July/November used `-m flip -d`**, which discards every A/T and C/G SNP
   even when its REF matches GRCh38 (fixref.c: in `flip` mode "skip ambiguous
   pairs: A/T or C/G", and `-d` drops what it skipped). That is the
   array-data recipe; these are forward-strand WGS calls. `flip-all` keeps
   them: REF match → unchanged, ALT match → REF/ALT and GT swapped (tested).
5. sort, split → `01_merged_phased_panel/mxb_lifted.chrN.bcf`, then the stamp
   `mxb_lifted.recipe` (`fixref=<mode>`). The driver re-runs 1a only when the
   stamp does not match `FIXREF_MODE`.

This is the **liftover control**: the same recipe gave all 50 MXB samples
Q = 1.000 in July. **But in July and November MXB (and SGDP, same recipe) had
no palindromic SNP at all while HGDP+1KG and Almarri v1 did** — a whole-cohort
absence at every A/T and C/G site of the other cohorts' files, part of MXB's
61.9% and SGDP's 39.6% whole-cohort missingness. 1c now refuses to merge
cohorts whose stamps record different fixref modes.

**SGDP, raw — `scripts/fetch-sgdp-vcf-variants.sh` + `1b_sgdp_raw_prep_clm.sh`**
(**replaces** `1b_sgdp_liftover_clm.sh`, the GLIMPSE-imputed phased_data2021
release: only 1KG-phase-3 sites, so a site 1KG p3 lacked was `./.` for all of
SGDP — 39.6% of chr22's November soft-union.)
- *Source:* `vcf_variants/` — one GATK UnifiedGenotyper VCF per public sample
  (279), hs37d5, **variant sites only** (no `0/0` records: first sample
  SS6004478 = 2,280,169 `0/1`, 1,825,700 `1/1`, 3,860 `1/2`), SNPs only,
  FILTER always `.`, and a per-call filter level `FL` ("0-9 or no value
  (N,?) with zero being least reliable; to threshold at FL=n, use all levels
  n-9"; SS6004478: 0 = 9.9%, N = 23.9%, 1–9 = 66.1%). Sample = Illumina ID
  (SS6004478 = B_Australian-3 in the metadata).
- *Fetch* (one job): resumable download of the 57 GB tar (60,957,091,840
  bytes; 21 MB/s from mhgcp), exact byte count, members vetted (regular files,
  bare names) before extraction, md5 of every file vs the `.md5sum` shipped
  with it, `vcfs/DONE`.
- *Prep* (array 1–22):
  1. manifest sample → SGDP_ID → Illumina_ID (`sgdp-raw-sample-map.py`)
  2. `bcftools merge -0 -m none -r N` over the per-sample files: a sample with
     no record where another has one = 0/0 (as Almarri's merge -0)
  3. variant call with **FL < 1 → `./.`**, not 0/0. FL ≥ 1 is the SGDP
     authors' own: "At 'filter level 1' which we recommend for most
     analyses…", and their ADMIXTURE used "genotypes of at least filter
     level 1" (Mallick et al. 2016, doi:10.1038/nature18964). Per sample
     (`&` in the expression; `&&` matches across samples).
  4. **site mask:** where more SGDP carriers fail FL than pass, every SGDP
     sample is set `./.` there — a variant-only file has no FL for hom-ref
     calls, so leaving the others 0/0 would make SGDP look monomorphic
  5. GT only; rename to the manifest's IDs; Picard `LiftoverVcf` (1a recipe);
     `norm -m-` → `+fixref -m ${FIXREF_MODE} -d` → `chrN` → sort →
     `sgdp_prepped.chrN.bcf` + `.recipe`; the imputed output found there is
     moved to `sgdp_prepped.chrN.phased_data2021.bcf`, not deleted
- SGDP is then a complete callset and is filled like the others in 1c (1c
  checks the stamp: an imputed SGDP is never filled).

**1d — `1d_almarri_prep_clm.sh`** (Almarri — **NEW**; already GRCh38, no liftover)
1. per sample, per chr: low-quality calls → **`./.`** (v4; see below) →
   `norm -m -snps --fasta-ref` → biallelic SNPs → REF check
   - **November panel (recipe v1):** `norm -c ws` — a mismatched REF is
     *overwritten* from the FASTA while ALT and the genotypes are kept, and
     stderr was discarded, so how often that happened is not in the logs.
   - **Script now (v3), on any re-run:** `norm -c w` + `+fixref -- -f GRCh38
     -m ${FIXREF_MODE} -d` (flip-all) — the same check MXB and SGDP get;
     stderr kept; stamp `v3 fixref=<mode>`, per-sample cache tagged
     `v3-<mode>`, so v1/v2 files are not reused.
2. `bcftools merge -m none -0` (137 samples): a genotype ABSENT from a
   sample's file becomes 0/0; a call set missing in step 1 stays missing
3. site mask (`ALMARRI_SITE_MASK=1`): where more Almarri genotypes are missing
   than ALT, all of Almarri is set missing -- as for SGDP
4. strip to GT → `almarri_prepped.chrN.bcf`; refuses above `MAX_MISSING` 0.30
   (an unfilled merge of single-sample files is ~99% missing); logs mean hets
   and missing genotypes per sample (`bcftools stats -s`)

**v4 (12 Nov 2026).** v1–v3 followed the Almarri paper: keep PASS and
GQ > 20, then `missing2ref`. A dropped call then became a confident 0/0.
Measured on chr22 for 6 samples: **24–28% of each sample's ALT calls**,
~97% of them heterozygous, were non-PASS (21–24%) or GQ ≤ 20 (2–4%) -- so
Almarri's heterozygosity and ALT frequencies were cut by about a quarter, a
cohort signal (November's Almarri-only cluster) rather than ancestry. v4 makes
every such call **missing** (SHAPEIT5 fills it from LD, as for any gap):
FILTER not PASS, or GQ ≤ `MIN_GQ` (20) or missing. Long Ranger 2.2.2's four
filters on one sample's chr22 (16,629 non-PASS ALT calls, ~96% het): allele
fraction (AO < 2 or < 0.15), rescued molecule in high-divergence sequence,
phasing-inconsistent (the ALT does not segregate with the linked-read
haplotypes -- the false-het signature of paralogous sequence, high GQ because
GQ counts reads), and QUAL. All four are genotype filters, so
`ALMARRI_PASS_FILTERS` (for a filter about the phase only) stays empty. The
missingness guard runs on the merge before the site mask: masked, wholly
missing records (segmental duplications) would fail a correct merge, and after
the mask an unfilled one could pass. Two `setGT` passes, because a single expression mixing the site-level
FILTER with per-sample GQ applies only the GQ part (bcftools 1.19, tested).
Stamp `v4 fixref=<mode> lowq=missing minGQ=<n> pass=<filters> sitemask=<0|1>`;
per-sample cache tagged by the same settings.

**1c — `1c_phasing-jointcall_clm.sh` with `PHASE_MODE=chunked`** (merge + site
selection). In July this step was `1b_phase_prep_clm.sh`: the **same** merge →
QC → soft-union logic with **2 inputs** (HGDP+1KG + MXB). Now **4 inputs**.
- *Idempotent:* exits at once when `merged_chrN.softunion.bcf` is newer than
  every input and its `.recipe` stamp records the same settings and code (a
  hash of 1c and its helpers); otherwise rebuilds (`FORCE_1C=1` always does).
- *Fixref consistency:* reads each cohort's prep stamp and stops if any was
  made with a different `FIXREF_MODE` than 1c's own HGDP+1KG step.
1. HGDP+1KG: drop `panel_exclude_ids.txt` (`bcftools view -S ^drop`). Job
   3742123 resolved the exclusion list to Slurm's spool directory, found
   nothing, and kept `CHMI_CHMI3_WGS2`; 1c now stops instead, and 4a / the
   benchmark exclude CHMI by name. Then `norm -m -any -c w` + `+fixref -m
   ${FIXREF_MODE} -d` (flip-all): every multi-allelic record split (was
   `-m -snps`, which left SNP+indel records whole), REF checked against GRCh38.

   **Relatives are kept** (`OUTLIERS=none`, the default from 12 Nov 2026),
   the policy 1e and the benchmark already state: only duplicates / MZ twins
   are removed (1e, KING > 0.354); the GWAS models kinship with a GRM
   (GENESIS); and pruning relatives preferentially removes less-studied
   populations, where background relatedness is higher. The one rule: a pair
   at KING >= 0.177 is never split across donor / reference / test. It is
   applied by `split-homog-donors-reference.py` (benchmark), by 4a (a
   simulation donor's relatives are kept out of the homogeneous panel, listed
   in `admixture_k8/donor_relatives.k8.tsv`) and by 2b (refuses to simulate
   while a donor has such a relative in any reference list; its ID check
   cannot see relatives). `verify-rebuild.sh` section 7 shows the same pairs
   as soon as 1e has run, before 4a and 2b.

   *July and November dropped `related_outliers.txt` here* (693 IDs,
   relatives and outliers; 3,398 HGDP+1KG kept in July; November's u235158
   source still held Koenig's 60 QC removals). The reference and donor
   lists name 28 of them; no stage complained, because 2b and 3b select
   samples with `plink2 --keep` and FLARE with `bcftools view
   --force-samples`, which skip absent IDs. So those runs used fewer samples
   than the lists hold: `amr_rfmix` 20 of 31 (Maya 11, Pima 3, Karitiana 2,
   Surui 2, Colombian 2), `amr_hgdpmxb_rfmix` 45 of 56, `afr_rfmix` 48 of
   50, `eur_rfmix` 49 of 50, `pel_rfmix` 84 of 86, `pel_eas_rfmix` 576 of
   587, `amr_simulation` 28 of 31, `amr_hgdpmxb_simulation` 53 of 56. With
   relatives kept those 28 are back. `OUTLIERS=<2-col file>` restores
   a drop; 1c stamps which (`outliers=` in the soft-union's `.recipe`).

   *Five list IDs are Koenig QC removals*, absent from `phased_haplotypes_v2`
   and `_v2_filter1` and named in `panel_exclude_ids.txt`: HGDP00621 and
   HGDP01270 in `afr_rfmix`, HGDP00029 and HGDP00150 in `sas_rfmix` (the
   file's "no gnomAD hard filter" group: Koenig's ancestry outliers and
   duplicate drops), and NA18864 (YRI) in `afr_simulation`
   (`ambiguous_sex`). November's panel had all five. July's did not
   (`afr_rfmix` 46 of 50 with its two related outliers, `afr_simulation` 29
   of 30, `sas_rfmix` 48 of 50), and neither does this rebuild's (48 of 50,
   29 of 30, 48 of 50) until the lists are trimmed or the IDs replaced.
   `verify-rebuild.sh` section 4 names each missing ID's reason. 2b and 3b
   run without them; no arm reads `sas_rfmix` (4a loads it and no longer
   uses it). gnomAD's metadata labels HGDP00621 Bedouin and HGDP01270
   Mozabite (both MID), HGDP00029 Brahui and HGDP00150 Makrani (CSA).
   Resolved 14 Nov 2026: `afr_rfmix` takes HGDP00933 and HGDP00934 (Yoruba,
   in the panel, on no list, no recorded relationship), so 50 again;
   `sas_rfmix` drops its two (48). `afr_simulation` drops NA18864 (29, as
   July ran): the first replacement candidate, NA18487, has a parent-child
   relationship in gnomAD's inference, and 4a ran (job 3754843) before
   another was chosen. A donor added after 4a may already sit in its AFR
   homogeneous panel, which 2b refuses, so a 30th YRI donor needs 4a again.

   *17 Nov 2026: the 3-way EUR and AFR references are 1KG IBS and YRI again.*
   Commit 046baad (6 Jul) rebuilt `eur_rfmix`, `afr_rfmix`, `eas_rfmix`,
   `sas_rfmix` and `amr_rfmix` as Koenig-based ADMIXTURE anchor lists
   ("anchor rebalance"), and every 3b arm reads `eur_rfmix` / `afr_rfmix` as
   its EUR / AFR reference (`IBS_1KG` / `YRI_1KG` in the code). 2a1eb2c (1 Aug)
   restored `amr_rfmix` only, so from 6 Jul the EUR reference was HGDP French
   23, Orcadian 12, Sardinian 8, Russian 6 and the AFR reference HGDP Yoruba
   21, Mandenka 20, BantuKenya 6 -- not the design (Jessica's lists, described
   in 4a_build_homogeneity_panel_clm.sh: EUR = 77 1KG IBS, AFR = 77 1KG YRI,
   kept in `reference_ids/rerun_jessicaH_3way_lai_sim/`, named `backup_pre_koenig/` until 17 Nov). `eur_rfmix.txt` and
   `afr_rfmix.txt` now hold those lists exactly as July's panel used them
   (74 IBS, 72 YRI: the other 8 IDs were relatives July's panel did not
   contain); the 6 Jul lists (HGDP EUR/AFR/EAS/SAS, MXB AMR) are in `reference_ids/genotype_harmonization_benchmark_testing_anchor_list_20260706/` -- July's K=5 run used them as its supervised anchors, and the benchmark's pool is HGDP, Almarri, SGDP and MXB (never 1KG). The
   donors (IBS, YRI) stay disjoint from them and, with
   `scripts/donor-switch-plan.py`, unrelated to them. The four RFMix arms now
   key their cached reference haplotypes (`IBS_1KG_chrN.haps` ...) by the
   list's contents, so a changed list is re-extracted rather than reused.
   `eas_rfmix` / `sas_rfmix` stay the 6 Jul lists: no 3-way arm reads them.

   *2b's donor/reference check* reads every `reference_ids/*_rfmix.txt`
   (with 4a's `amr_homog_rfmix.txt`) and the 3-way simulation's own
   homogeneous panel, `04_homogeneity_panel/lai3way_k5/homog_5pop_samples.txt`
   (4a `HOMOG_SET=lai3way`, from 17 Nov; `HOMOG_K` / `HOMOG_LIST_DIR`). Job 3754844
   (16 Nov 2026) failed 22/22 because it also read July's K=5
   `homog_5pop_samples.txt`, built before 4a excluded donors: 102 of the 116
   donors are in it, although no arm reads that list. `verify-rebuild.sh`
   section 7 runs the same check on the same sets before 2b starts.
   `python3 -I scripts/related-outliers-breakdown.py` (read-only) shows what
   that list holds by the gnomAD metadata flags (relatives vs PCA outliers,
   per population) and, where the metadata records relationships, every pair
   among the curated lists with the side each member is on.
2. add SGDP, then Almarri: chr-prefix contigs, drop IDs already present or
   excluded → **3,921** in November (relatives dropped; more with them kept). Then **every input's remaining
   multi-allelic records are split** (`norm -m -any`): `merge -m none` joins
   records only on identical REF+ALT, so a SNP allele inside a cohort's
   multi-allelic or SNP+indel record never met another cohort's biallelic
   record, and with the fill below that cohort would read 0/0 on it.
3. Site provenance logged (from these files, before any fill) -- COUNTS per
   cohort and class (absent / swapped / flipped / other ALT / other REF), at
   three stages: `provenance/chrN.{merged,qc,softunion}.{cohorts,patterns}.tsv`
   and `[provenance ...]` lines in the 1c log. The POSITIONS, with each
   cohort's missingness and AF, come from `scripts/cohort-site-report.sh`
   after 1c (`provenance/chrN.site_report.tsv.gz`). Then the
   **0/0 fill** (`FILL_REF_COHORTS`, rebuild default **"HGDP1KG MXB
   ALMARRI"**; 1c on its own: off): for a cohort whose file is a complete
   variant callset, a record absent from it is 0/0 for that cohort instead of
   `./.`. Done by merging each such file with `-0` against a union of every
   input's records that carries one placeholder sample (GT `./.`), then
   removing the placeholder. The earlier sites-only union could not work:
   `-0` only fills records at which some input has FORMAT/GT, so the records
   a cohort lacked came out with no genotype at all (tested on bcftools 1.19;
   the fill had never run on the cluster). `./.` calls inside a file stay
   missing. Never for the imputed SGDP release (the driver refuses).
   Then `bcftools merge -m none` of the four files (1,153,425 records in
   November).
4. plink2 site QC: `--max-alleles 2 --snps-only just-acgt --rm-dup exclude-all
   --set-missing-var-ids @:#[b38]`. Then:
   - **3a — groups and fill check.** Cohort × super-pop groups
     (`build-missingness-groups.py`; built once, used by 3a, 5 and 6). For
     every record a filled cohort's file lacked, `fill-contradictions.py`
     asks whether another cohort of the same super-pop makes "nobody here
     carries it" implausible: ≥ 5 copies there, and P(no copy in this
     cohort's haplotypes | that frequency) = (1 − AF)^AN < 1e-3, both groups
     ≥ 20 samples. If so the source's own QC dropped a real site, and the
     position is dropped rather than entering as a cohort-wide 0/0.
     `provenance/chrN.fill_check.tsv`.
   - **3b** — every position still holding more than one record (split
     multi-allelic sites: cohorts saw different ALTs) and every 3a position
     is dropped. Logged with counts; provenance logged again.
5. **Soft-union site gate**, one `+fill-tags -S` pass
   (`softunion-sites.py`; same rule, read once instead of once per group).
   Per group with ≥ 10 samples keep `AN>0 && MAF ≥ 0.005`, union of
   positions. `SOFTUNION_GROUPS` decides whose samples form the groups:
   - **`superpop` (rebuild default):** every cohort's samples pooled by
     super-pop (HGDP1KG_AMR + MXB_AMR + SGDP_AMR → AMR; Almarri → MEN). With
     the fill, AN includes the filled 0/0s, so a cohort-unique variant's MAF
     is its frequency in the whole super-pop.
   - **`sample_groups`** = July/November: `sample_groups.tsv` (gnomAD
     HGDP+1KG metadata + MXB) only; **SGDP and Almarri are in no group**
     (groups summed to 3,507 of 3,921), so their private variants could not
     enter the panel. November: AFR 768 → 221,495 · AMR 448 → 231,038 · CSA
     691 → 164,673 · EAS 736 → 140,477 · EUR 682 → 159,766 · MEN 155 (HGDP
     Middle East, not Almarri) → 192,260 · OCE 27 → 109,295 ⇒ **301,309
     positions, 301,913 records** (604 positions with two records — now
     removed in 3b).
   Group sizes and passing counts: `provenance/chrN.softunion_groups.tsv`.
6. **ADMIXTURE site list** (`merged_chrN.admix_sites.bed` + `.summary.tsv`):
   the positions where every cohort × super-pop group with ≥ 5 samples is
   ≤ 15% missing and the whole panel ≤ 7%, on these PRE-phasing genotypes.
   With the fill, filled cohorts have no source-file gaps left, so this rule
   acts on unfilled cohorts and on genuinely uncalled genotypes. Only 4a and
   1e use it; the soft-union, and so the LAI panel, keeps every site.
7. publish `01_merged_phased_panel/merged_chrN.softunion.bcf` and stop

> `REPO_DIR` must resolve correctly or `EXCLUDE_IDS` points at Slurm's spool
> directory. That is how CHMI got into the panel: 3,921 *looked* like
> confirmation but was the arithmetic of the exclusion not happening.

**1b_phase_chunks — `1b_phase_chunks_clm.sh`** (joint phasing; **same method
as July** — the chunked path dates from June 15 — only the input changed: 3,921
samples from 4 cohorts instead of ~3,449 from 2)
- `SHAPEIT5_phase_common --input merged_chrN.softunion.bcf --map
  chrN.b38.gmap.gz --region <8 cM chunk + 2 cM overlap> --filter-maf 0.001`,
  594 chunks
- **All 3,921 samples from all four cohorts are phased together** in every
  chunk.
- `phase_common` **imputes every `./.`** from the haplotypes of samples called
  at that site. This is the only way phasing changes a genotype. Phase itself
  — switch errors, chunk boundaries — cannot affect ADMIXTURE, which reads
  unphased genotypes.
- `--filter-maf 0.001` drops sites under 0.1% across all samples.
- Idempotent with a staleness check: an existing chunk is reused only if its
  sample count equals the soft-union's.

**1b_phase_ligate — `1b_phase_ligate_clm.sh`**
- `SHAPEIT5_ligate` per chromosome →
  `merged_chrN.shapeit5_phased.softunion_maf005.bcf` (+ `.rechr.bcf` with bare
  contigs for RFMix) — same command and filename as July.

**1e — `1e_king_dedup_clm.sh`** (**NEW**) — concat → restrict to 1c's called
sites → LD-prune → `plink2 --make-king-table` → `king/king.kin0`,
`king/duplicate_pairs.tsv`, `reference_ids/panel_duplicate_drop_ids.txt`.
Pairs in `known_duplicate_pairs.tsv` (an SGDP sample whose alias is in the
panel) are duplicates whatever KING says. November, on all sites: 48 KING
duplicates, but 26 known same-person pairs scored only 0.30–0.35 because
SGDP's filled-in genotypes blurred identity — those people were in the
panel twice. **The HGDP/1KG copy is kept** (keep-priority MXB > ALMARRI >
HGDP > 1KG > SGDP): it is joint-called and complete, the SGDP copy lacks ~40%
of sites. `KING_INPUT=softunion` reads the pre-phasing genotypes instead.

**4a — `4a_build_homogeneity_panel_kN_clm.sh`** (**K=7 from 17 Nov; K=8 16 Nov; K=5 in July**)
1. concat the 22 phased chromosomes → `admixture_k8/all_chr.phased.bcf`
2. restrict to 1c's ADMIXTURE site list (`--extract bed1`), then
   `plink2 --indep-pairwise 500 50 0.1` over all samples at once →
   `--max-alleles 2 --snps-only just-acgt --make-bed` → `admixture_k8/lai_ref`.
   **Global MAF floor 0.005 from 18 Nov** (`ADMIX_PRUNE_MAF`; `none` turns it
   off), measured on the unrelated fit set in the same plink2 call as the
   prune. Every site already passed the soft-union rule (MAF ≥ 0.005 in at
   least one super-pop), so a global floor removes the variants common in one
   super-pop and rare overall. Where it has been:

   | date | global floor | pruned markers | per fit |
   |---|---|---|---|
   | July, until 16 Nov | `--maf 0.01` | 528,876 (Nov, no site list); 352,550 on 16 Nov (K=8), when it removed 11,285,750 of 18,880,026 sites | both fits 8 h 23 min |
   | 17 Nov | none | 3,196,141 (K=7), 70% of them MAF < 0.005 overall | ~33 h (K=7), ~22 h (K=5) |
   | **18 Nov** | **`--maf 0.005`** | fewer than the 958,622 of the 17 Nov set at MAF ≥ 0.005: the 500-*variant* window spans more genome once rare sites are gone (0.01 kept 352,550 where the unfloored set held 506,129 at ≥ 0.01) | roughly 10 h or less; `PRUNE_ONLY=1` gives the number |

   0.005 is between the two: a variant needs about 37 copies among ~3,700
   fitted samples, so a variant private to OCE (~55 samples) needs roughly a
   third of OCE's alleles. Rare-overall variants are also where call-set
   differences between the merged sources show most. A count floor is not the
   lever: no marker of the 17 Nov K=7 set had a minor allele count under 6, so
   `--mac 3` would remove nothing. The rule is stamped in `lai_ref.recipe`; a
   pruned set made under another rule is rebuilt. `lai_ref.mac_bins.tsv`
   counts the markers by minor allele count over every sample of the run --
   with the floor, bins 1-2 must be empty and bin 3 nearly so (the floor is
   measured on the unrelated fit set). `PRUNE_ONLY=1` stops after this step
   and prints the count and the expected ADMIXTURE hours; `ADMIX_MAX_SNPS=<n>`
   thins at random (seed 20261116) if the count is too slow to fit.
3. **fit set** (from 12 Nov 2026, with relatives kept in the panel):
   ADMIXTURE assumes unrelated individuals, so both runs are fitted on an
   unrelated subset -- 1e's duplicate copies out, then greedily the sample
   with the most relatives at KING ≥ `FIT_KINSHIP` (0.0884, second degree)
   until no such pair is left -- and every relative is **projected** onto the
   fitted allele frequencies (`admixture -P`). Relatives still get a Q, can
   still be homogeneous and enter the panel, but do not shape the clusters,
   and cannot be anchors. The merged Q keeps every sample in `.fam` order, so
   discovery and stage C read it as before. Record:
   `admixture_k8/admixture_fit_set.k8.tsv` (sample, fit / projected, why,
   related_to, kinship). `FIT_UNRELATED=0` fits on everyone (July/November).
   Then unsupervised ADMIXTURE K=7 on the fit set
   (`admixture_k8/unsup/fit/`), relatives projected (`unsup/proj/`), merged
   → `admixture_k8/unsup/lai_ref.8.Q` (~10 h; reused while its inputs are
   unchanged)
4. **data-first discovery:** per cluster, rank every eligible sample, name
   the cluster from its top 50 at Q ≥ 0.93, then **choose the anchors —
   highest Q and central** (from 16 Nov; see §3, 4a): PCA of the fit set
   (`admixture_k8/fit_pca.*`), and of the samples at Q ≥ 0.93 carrying the
   cluster's label, the 50 with the smallest sum of their Q rank and their
   rank of distance to the group's median → `anchor_selection.k8.tsv`
5. supervised ADMIXTURE K=7 on those anchors, fitted on the same set
   (`sup_fit/`), relatives projected onto its allele frequencies
   (`sup_proj/`), merged → `admixture_k8/lai_ref.8.Q` → labels → tiered
   homogeneity → `lai_ref_panel_samples.k8.tsv`

July had no step 4: the curated `*_rfmix.txt` lists were the supervised labels
directly, and the unsupervised run only reported proportions.

#### Where a cohort difference can enter the genotypes

| step | cohorts | mechanism | what tests it |
|---|---|---|---|
| source calling | all four differ | caller, filters, coverage | same person twice (duplicate pairs): genotype discordance, Q difference |
| liftover (1a, 1b_sgdp) | MXB, SGDP | strand / allele errors | duplicate pairs: 0/0 ↔ 1/1 swaps; MXB is the working control |
| `+fixref -m flip -d` | July/Nov: MXB, SGDP (not HGDP+1KG, Almarri v1) | every A/T and C/G SNP discarded → whole-cohort absence at palindromic sites | stamps checked by 1c; count A/T + C/G records per cohort file |
| `GQ>20` + `missing2ref` (1d v1–v3) | Almarri | uncertain ALT calls → 0/0: 24–28% of ALT calls per sample (chr22) | fixed in v4 (→ missing); 1d logs mean hets per sample |
| merge without `-0` (1c) | whoever lacks a site | whole-cohort `./.` | pre-phasing missingness; sites with AN = 0 per cohort |
| SHAPEIT5 imputation | the same | `./.` filled from other cohorts' haplotypes | het / F on the ADMIXTURE input vs pre-phasing |
| split records (1c) | 604 positions on chr22 | each cohort called on a different record | listing the duplicated positions |

`scripts/diagnose-panel-batch.py` runs each test (modes `q`, `psc`, `het`,
`pair-ids`, `pairs-gt`).

### 1e — KING de-duplication
Genome-wide KING-robust between the merge and ADMIXTURE. Drops **only**
pairs > 0.354 (duplicate/MZ). Keep-priority: non-1KG → MXB/Almarri →
lexicographic. Catches the same person under two IDs across HGDP / 1KG / MXB /
Almarri, which `1c`'s ID-level de-duplication cannot.

### 4a — ADMIXTURE K=7 (K=8 until 16 Nov)
Three stages: anchor discovery across all 7 pools → supervised ADMIXTURE at
K=7 → homogeneity filter.

**Data-first anchor discovery** (stage A). For each of the K unsupervised
clusters: rank **every** eligible panel sample by its loading on it, take the
top 50 at Q ≥ 0.93, and **name the cluster from who they are** — the
plurality metadata super-pop among them, holding ≥ 50% of the set. If two
clusters claim one super-pop, the one with more supporting samples keeps it; a
single unnamed cluster becomes OTH. Then the anchors are chosen from the
cluster's samples at Q ≥ 0.93 by Q and centrality (see "Which anchors" below).
Eligible = everyone except named exclusions, 1e's duplicate drops and the
projected relatives. No candidate pools and no a-priori naming: the curated
`*_rfmix.txt` lists are not read (their overlap with the discovered anchors is
reported).

Every cluster's composition goes to `discovery_clusters.k8.tsv` with the
reason it was or was not named. Two checks print in the log:
- **COHORT-DOMINATED** (warning): one minority cohort supplies ≥ 80% of a
  cluster's anchors *and* they span several super-pops — the batch signature
  (November's c7: 187 SGDP across EAS/EUR/SAS/OCE).
- **single-cohort cluster** (note): one cohort, one super-pop — reported with
  how strongly that super-pop's samples from *other* cohorts load on it
  (November's Almarri cluster: HGDP Middle Easterners 0.09).

> **History.** The pool design named clusters a priori by curated label lists.
> In job 3747721 `seed_for()` read `$sp` in the same `local` statement that
> assigned it, the lists came back empty, and the pools fell back to SGDP
> alone — 276 of 3,921 samples were the only candidates for every super-pop.
> Ranking everyone within the cluster is the method July settled on; in July's
> K=5 Q it gives SAS 279 samples at ≥ 0.93 (top 50 = 1KG ITU/STU/BEB).

**Markers.** Two switches decide what ADMIXTURE reads:

| switch | values | effect |
|---|---|---|
| `ADMIX_SITE_FILTER` | **1** / 0 | restrict to 1c's per-cohort × super-pop missingness site list |
| `ADMIX_INPUT` | **phased** / softunion | phased panel (gaps filled by SHAPEIT5) or the pre-phasing soft-union (gaps stay missing; ADMIXTURE skips missing genotypes) |

Within a cohort, missingness is all-or-nothing per site (it comes from the
cohort's source file lacking the site: every SGDP sample misses exactly 39.55%
of chr22, every MXB sample 61.9%), so **any per-group threshold selects the same
sites — those every cohort has, i.e. the intersection**, which drops the
MXB-specific variants. `ADMIX_INPUT=softunion ADMIX_SITE_FILTER=0` keeps every
site without using a single imputed genotype. `scripts/site-presence-check.py`
measures how many "absent" sites are really allele mismatches (REF/ALT swapped,
strand flipped, different ALT) that harmonisation would recover.

**The rule (17 Nov):** every **super-pop** (its cohorts pooled) ≤ 7% missing,
and every **cohort** (HGDP1KG, MXB, SGDP, ALMARRI; its super-pops pooled) ≤
15%; a margin under 5 samples is not judged alone. 4a measures it itself
(`ADMIX_SITE_SOURCE=4a`, the default) on 1c's published pre-phasing
soft-union, per chromosome with `scripts/admix-site-list.sh` (groups rebuilt
exactly as 1c built them) and `scripts/admix-site-margins.py`, into
`04_homogeneity_panel/admix_sites_4a/`, each list stamped with the rule.
1c is not re-run: the list only decides which markers ADMIXTURE reads.

1c still writes its own list under the 12 Nov reading -- each cohort ×
super-pop group ≤ 15%, panel ≤ 7%, which put the two numbers the wrong way
round; `ADMIX_SITE_SOURCE=1c` uses it. On chr22 1c's rule kept 279,155 of
367,930 soft-union sites. `scripts/admix-rule-compare.sh` counts the sites each
version keeps (its `coded` row must reproduce 1c's list byte for byte, which
checks the group reconstruction).

The site list itself: In November, without it, every SGDP sample had no call
at 38% of chr22's ADMIXTURE markers and every MXB sample at 62%; SHAPEIT5 had
filled those in, and three of the eight clusters tracked cohorts (c0 ≈ 0.4 in
every 1KG/HGDP sample, c4 Almarri-only, c7 SGDP-only).

**Rebuilds.** Concat, LD-pruned bed, unsupervised Q and supervised Q are each
reused only if newer than all of their inputs; the `.pop` is replaced only if
the anchors changed, so a discovery-only rerun repeats the supervised fit only
when it has to.

Writes `lai_ref_panel_samples.k8.tsv`, the Q table everything downstream
reads: `sample, cohort, pop, super_pop, Q_sup_*, max_Q_sup, argmax_pop_sup,
Q_unsup_*, sup_unsup_agree, included`.

**Homogeneity threshold: 0.94 for every super-pop (17 Nov),** and the same
0.94 for anchor candidates. Before: 0.95 for EUR/EAS/AFR (SAS too from 16
Nov) and 0.93 for AMR/OCE/MEN. One number is the compromise: the
under-represented groups lose the few samples between 0.93 and 0.94, and every
super-pop meets the same standard. The cut acts in stage C only, so changing
it later re-runs minutes of work, not the fits (`HOMOG_STRICT_SUPERPOPS` +
`HOMOG_THRESHOLD_STRICT` restore tiers).

**K=7 (17 Nov), not 8.** At K=8 the eighth cluster (OTH) had 2 anchors
(SGDP Ulchi, HGDP Oroqen) and the supervised fit turned it into a second
African component: San, Ju'hoan and Khomani 1.00, Mbuti 0.39-0.77, Masai 0.32,
Dinka 0.31, LWK 0.29 (0 of 97 LWK homogeneous). At K=7 there is one component
per named super-pop and no OTH; 4a brings OTH back only when K exceeds the
named super-pops.

**Two homogeneity panels, one script (17 Nov): `HOMOG_SET`.** The 3-way
NAT/EUR/AFR LAI simulation and the genotype-harmonization benchmark need
different panels. The 3-way one is judged at K=5 on July's five super-pops
(never by the K=7 run), with its anchors found by its own unsupervised fit on
the expanded panel rather than taken by name as July's were:

| | `HOMOG_SET=global` (default) | `HOMOG_SET=lai3way` |
|---|---|---|
| for | the benchmark, the global panel | the 3-way simulation's HOMOG arms |
| K, super-pops | 7: EUR AFR AMR EAS SAS OCE MEN | 5: EUR AFR AMR EAS SAS |
| samples | every panel sample | MEN, OCE and no-super-pop samples left out of the missingness rule, the LD prune and both fits |
| anchors | discovered: central and high-Q (≥ 0.94) | discovered the same way at K=5, among its own samples (`reference_ids/<sp>_anchors.k5.txt`); the log counts how many are July's by-name anchors again (the 6 Jul lists: HGDP EUR 49, AFR 46, EAS 50, SAS 48; MXB 50) |
| homogeneous | Q ≥ 0.94 and argmax = label | the same |
| writes | `04_homogeneity_panel/*.k7.*`, `reference_ids/<sp>_homog_rfmix.k7.txt`, `homog_7pop/` BCFs | `04_homogeneity_panel/lai3way_k5/`, and `reference_ids/amr_homog_rfmix.txt` = July's AMR panel + every AMR sample homogeneous here (`HOMOG_3WAY_MODE=union`; `keep-july` adds only samples July never saw, `current` uses this run alone), minus donors and their relatives |

**One ADMIXTURE fit per job (17 Nov).** With no global MAF floor the LD-pruned
set was ~3.2 million markers (K=7: 3,196,141 of 19,164,602 sites after the
missingness rule; 70% of them rare overall but ≥ 0.005 in some super-pop), so
the 17 Nov PRUNE_ONLY runs estimated ~33 h per fit at K=7 and ~22 h at K=5 —
two fits in one job (67 h at K=7) would not fit the 72 h limit. The drivers
therefore run each set as two jobs: `ADMIX_STOP_AFTER=unsup` (the unsupervised
fit and its projection), then the full run, which reuses it. The time guard
counts only the fits the job itself will run. With the 0.005 floor (18 Nov)
both fits would probably fit in one job; the split is kept because it costs
nothing and keeps each job short.

Both use the same rule (each super-pop ≤ 7% missing, each cohort ≤ 15%,
measured on the set's own samples; global MAF ≥ 0.005 on the set's unrelated
fit set, 18 Nov; LD on the unrelated fit set;
relatives projected). lai3way also drops sites with no minor allele among its
fitted samples (`--mac 1`, not a frequency floor): 1c's soft-union keeps OCE
and MEN singletons, which are invariant without those samples. The two jobs
share 1e's KING table and one concatenated input
(`04_homogeneity_panel/admixture_input/`, built under a lock; a current concat
left in `admixture_k*/` is adopted by hard link). The driver submits both;
stages 4-6 wait on lai3way, the benchmark on global.

**Pruning on the unrelated set (17 Nov).** The LD prune used every sample
until 17 Nov, so 638 relatives and 153 duplicate copies entered every r2; it
now measures LD on the ADMIXTURE fit set and makes the bed for every
sample.

**Which anchors — highest Q and central (from 16 Nov).** The top 50 by
unsupervised Q are each ancestry's most *drifted* members, because a sample
loads 0.999 on a cluster exactly when it is furthest from every other cluster.
On 16 Nov that made the anchors: EUR all northern (GBR 20, Orcadian 11, CEU 10,
FIN 4, Estonian 2, Icelandic 2, Norwegian 1); AFR 24 of 50 rainforest
hunter-gatherers (Biaka 16, Mbuti 8); EAS 45 CDX; SAS 49 ITU/STU. The
supervised components, which the anchors define, centred on those ends.

`ANCHOR_SELECT=central` (default): candidates are the samples at
Q ≥ `ANCHOR_MIN_Q` on the cluster that carry its metadata super-pop. Each is
ranked by Q (highest first) and by Euclidean distance to the candidates'
median in the top `ANCHOR_PCS` (10) principal components of the unrelated fit
set (each PC scaled by √eigenvalue; nearest first); the 50 with the smallest
rank sum are the anchors. The centre is the median of the *candidates*, not of
every sample with the label: for AMR that median sits among admixed Latin
Americans, and the nearest unadmixed samples would be the ones with the most
European ancestry. `ANCHOR_CENTER_MIN_Q` lowers the Q a sample needs to count
toward the centre. PCs use the markers' common subset (`ANCHOR_PCA_MAF` 0.05)
because PCA up-weights rare variants by 1/√(p(1−p)); that floor shapes the PCA
only. `anchor_selection.k8.tsv` lists every candidate with its Q, distance,
both ranks and whether it was chosen. `ANCHOR_SELECT=top_q` restores the old
rule.

What this is not expected to fix: at K=8 the MEN component takes 23–38% of
every southern European population's supervised Q (TSI 0.35, Sardinian 0.38,
IBS 0.23), so none of them can be an EUR anchor at Q ≥ 0.93, and more central
anchors move their EUR share up only as far as the model's EUR/MEN split
allows. July's K=5 had no MEN component and kept every European population.

**Anchor floors — two, not one.** Stage A takes the candidates with cluster Q
**≥ `ANCHOR_MIN_Q` (0.93)**, **capped at `ANCHOR_CAP` (50)** — the "top 25–50
≥ 0.93" rule from the NAT/Mexican accuracy pipeline. Two separate floors
matter, and they fail differently:

| floor | what it stops | failure if absent |
|---|---|---|
| `ANCHOR_MIN_Q` = 0.93 | an **admixed** sample becoming an anchor | that sample *defines* the ancestry in the prior, and every Q estimate for the component inherits it |
| `ANCHOR_FLOOR` = 10 | a super-pop's prior resting on **2–3 samples** | the label is populated so ADMIXTURE runs happily, but the component is estimated from a handful |

The second is the quiet one. 4a previously failed only on a label with **zero**
anchors — because that is what makes `ADMIXTURE --supervised` refuse to run —
so a super-pop with two passed silently and then fed the homogeneity call, the
donor/reference split, the test genomes' `sample_type` and the GWAS ancestry
covariate. **OCE, MEN and SAS are where this bites.**

A thin set is a loud **warning**, not an error — a thin super-pop may be a fact
for the methods rather than a reason to stop. `ANCHOR_FLOOR_FATAL=1` makes it
stop. Counts, floors and the cap are written to
`04_homogeneity_panel/anchor_summary.k8.tsv`, so "how many exemplars defined
this ancestry" is answerable later without re-reading logs.

> **Do not lower `ANCHOR_MIN_Q` to fill a thin set.** That admits admixed
> samples as the definition of a pure ancestry, which is worse than a small set
> of clean ones.

`ANCHOR_MIN_Q_STRICT` optionally tiers the Q floor (0.95 for EUR/EAS/AFR, as
the homogeneity threshold does). **Off by default** — 0.93 uniform is what the
previous pipeline used, and changing it silently would move every anchor set.

**`included == 1` requires `max_Q ≥ threshold` AND `argmax == assigned
super_pop`**, and is 0 for donors, `super_pop = NA`, OTH, named exclusions and
KING duplicate drops. It is *not* a general-purpose eligibility flag — see §5.

`super_pop_of()` harmonises the vocabulary: **MXB → AMR**, SGDP via its own
table, **CSA → SAS**, **MID → MEN**. Anything else → `NA`. It aborts if its two
vocabularies don't overlap.

> The HOMOG panel is HGDP + 1KG + MXB joint-called and Q-filtered. **It
> includes MXB samples** — it is not MXB-free.

---

## 4. Step 5 — the donor / reference split

`scripts/split-homog-donors-reference.py`

- `reference = clamp(round(0.5 × N), 15, 65)`, donors take the rest.
- Stratified by population (largest-remainder) so **both sides carry a
  proportional share of every sub-population** — otherwise a difference
  between sides is a difference between datasets read as an ancestry effect.
- `--mode benchmark`: **1000 Genomes is never a donor**, because the external
  imputation servers contain 1KG. Reference = 1KG + a stratified share of the
  non-1KG samples.
- `--mode lai`: no source rule; every sample eligible for either role.
- `--no-donor-superpops OCE,SAS` (SAS from 18 Nov): reference only, **no
  donors**, remainder left *unallocated* so the assembler takes it as test
  genomes. The reference is filled from 1000 Genomes first (spread over its
  populations, skipping a 1KG sample related to a test-bound genome) and from
  the other sources only for what 1KG cannot cover. OCE has no 1KG population
  and no model: **~15 reference + ~5 real test, 0 donors**. SAS has ~10
  homogeneous genomes outside 1KG, too few to simulate a truth set from, so
  `reference_ids/models.tsv` has no SAS model and SAS is **65 1KG reference,
  0 donors, and its ~10 homogeneous genomes become real_homogeneous test
  genomes** beside the many real_admixed SAS ones. No SAS local-ancestry
  truth set until more homogeneous SAS genomes arrive (DUAs).
- **The pool (18 Nov): `--pool-column homogeneous` for the benchmark.** 4a
  writes two flags: `homogeneous`, the pool both projects choose from (the
  homogeneity rule, 1e duplicates and named exclusions out, no project's
  roles), and `included`, the LAI project's panel — `homogeneous` minus ITS
  simulation donors and their relatives. The benchmark split reads
  `homogeneous`: until 18 Nov it read `included`, so 56 Native American
  genomes (31 HGDP/SGDP + 25 MXB) and their relatives were missing from the
  benchmark for a rule of the other project. `--mode lai` keeps `included`.
  A table from before the column stops the benchmark split with an error
  rather than silently splitting the old pool.

### Fill order: 1000 Genomes first

**1KG can only ever be reference** — barred from the donor side because the
imputation servers contain it. A non-1KG sample can be either. So the
single-use resource is spent first and the dual-use one only on what 1KG
cannot cover. That is what `--ref-share default=0` expresses, and it explains
the contract's share values rather than leaving them arbitrary. The 1KG part of
the reference is spread over the 1KG populations by the same largest-remainder
rule as the rest (until 18 Nov it was the first IDs alphabetically: EUR got
GBR and FIN only, AFR ACB and GWD, EAS CHS -- and with no YRI in the reference
the coverage rule then pulled HGDP Yoruba donors to cover them). A 1KG sample
related to a donor is left out of it.

| super-pop | 1KG coverage | share |
|---|---|---|
| AFR, EUR, EAS, SAS | good | `0` — 1KG fills the reference |
| **MEN** | none — 1KG has no Middle Eastern population | `0.5` (0.33 until 18 Nov: 18 reference genomes of 59), and `--ref-min 25` |
| **NAT** | only *admixed* American (PEL/MXL/CLM/PUR), mostly below 0.93 | `NAT@MXB=0.5` |
| **OCE** | none | all reference is non-1KG |

### Sub-population coverage, on ancestry GROUPS not labels

`--ref-min-per-pop default=3,EUR=0` guarantees three reference samples per
uncovered population, because 1KG carries none of the deeply divergent
lineages: no Biaka, Mbuti, San or Mandenka in AFR, none of the HGDP minority
populations in EAS. It never takes a super-pop's donors below the donor floor
(15), and never moves a sample whose relative is a donor.

**EUR is exempt (18 Nov).** With 3, the K=8 preview moved 3 of EUR's 18
non-1KG genomes into the reference for Basque/Sardinian/Orcadian-type
coverage and left EUR exactly at the floor of 15 donors; EUR donors are scarce
enough that they stay donors. The split still names the EUR lineages the
reference then lacks, under "decided already", so the methods can state them.
A super-pop named in the terms is a decision and its gaps do not fail the
step; super-pops on `default` still go under "NEEDS A DECISION".

`--ref-min 25` (18 Nov) is the smallest reference for a group whose reference
comes from its own pool -- MEN above all -- when the pool allows it with the
donor floor kept; never past half a pool under 50 (`plan_simset.R --ref-min`'s
rule). OCE, reference-only, keeps its own sizing so it still has test genomes.

> **Coverage is counted on ancestry groups, not raw labels.** The sources name
> the same people differently — 1000 Genomes `YRI` and HGDP `Yoruba` are one
> population sampled twice. Comparing raw strings reports Yoruba as uncovered
> while YRI sits in the reference, and the floor then pulls donor-eligible
> samples to fill a gap that does not exist. On a test panel that wasted 3 of
> 9 pulls. `reference_ids/population_aliases.tsv` holds the groupings with a
> `confidence` column: only `same` applies by default, `close` is opt-in via
> `--alias-confidence close`, because "near enough for a reference panel" is a
> judgement about the analysis rather than a fact about the data.
>
> Basque, Sardinian, Kalash, San, Mbuti and Biaka are **deliberately not
> aliased** — treating an isolate as covered by a neighbour defeats the check.
> The 1000 Genomes American populations are listed as **not** equivalent to
> Indigenous American ones.

### How the cut sizes the donor pool (`scripts/homog-cut-preview.py`, 18 Nov)

Read-only. For each super-pop and each cut (default 0.90, 0.92, 0.94) it
counts who 4a would call homogeneous and splits them into donor-eligible and
reference-only (1KG or a server-panel id), and says how many of them are LAI
simulation donors or their relatives — out of the LAI panel (`included=0`),
in the benchmark pool (`homogeneous=1`, 18 Nov). It lists the samples each
lower cut adds, by population, with their strongest other component, and the
mean `1 − share` of each pool: the part of a donor's genome the simulation
labels with the wrong ancestry. At the table's own cut, donor + reference must
equal the table's `homogeneous` count (the check on the script itself; on a
table from before 18 Nov, `included` = donor + reference − LAI donors).

**What the K=8 table showed (18 Nov).** SAS: 10 donor-eligible at 0.94 and only
11 at 0.90 — most HGDP Pakistani genomes sit under 0.70 SAS, so no cut helps;
hence SAS reference-only (§4). EUR: 19 at 0.94, 23 at 0.92, 44 at 0.90 (the
0.90 step is 16 HGDP Russians at ~0.906 EUR with ~0.06 EAS). The cut stays 0.94.
MXB: 4a counts every MXB sample homogeneous whatever its Q (July's curated
set); the preview prints how many fall under each cut by their own AMR share.

0.94 is one number on both sides — anchors, the homogeneous panel, the split,
`plan_simset.R --threshold` and the contract's "donors are homogeneous: ≥ 0.94"
(`donors.tsv` has no `top_fraction`, so intake cannot check it). A looser
donor cut for a thin group would be a documented exception made on both
sides, and the samples it adds are today `real_admixed` test genomes, so
they would leave the test pool.

### The suffix trap
Both modes write `<sp>_reference_<suffix>.txt`. **Run them with the same
suffix and the second silently overwrites the first** — one panel on disk,
nothing in the filename saying which policy made it, and an even chance the
benchmark is scored against a donor-inclusive reference. Use distinct
suffixes: `lai` and `bench`.

Filenames are **lowercase** (`eur_reference_bench.txt`). Consumers resolve
either case; looking only for uppercase finds nothing, and "nothing" reads as
"no lists were produced".

---

## 5. Steps 7 + 8 — cohort assembly

`scripts/assemble-benchmark-cohort.py` — needs only the Q file and the split,
so it runs within minutes of 4a landing.

### The four exclusions (test eligibility)
Each closes a different way a method could already have seen the haplotypes:

1. **donor** — its haplotypes built the simulated test genomes
2. **reference panel** — it's in the ADMIXTURE/LAI reference
3. **1000 Genomes** — the external servers contain 1KG
4. **`sgdp_server_panel_ids.txt`** — SGDP re-sequenced 41 1KG individuals under
   other IDs; a cohort column reading "SGDP" does *not* establish absence from
   a server panel
Plus named exclusions and KING duplicate drops. (**Supervised anchors are
NOT excluded** — see below.)

### Supervised anchors and role assignment

The `assigned` column comes from ADMIXTURE's **`.pop` prior**, so it records
which samples were supervised anchors and for which super-pop.

**Anchors are eligible for every role.** The supervised ADMIXTURE run is not a
method under evaluation — it is the upstream classifier that identifies the
pools (who clears the homogeneity bar, who is admixed), and its Q is never an
*outcome* of the benchmark. What the benchmark scores is imputation against
WGS truth, and local ancestry on the **simulated** individuals against their
`.hanc` tracts; real genomes are never scored for local ancestry at all. So an
anchor is not a sample the benchmark was trained on, and it may be a donor, a
reference member or a real test genome under the same rules as everyone else —
including **no 1KG donor and no 1KG real test genome**.

**But an anchor's homogeneity is read from the `max_Q_unsup` column.** Its
supervised Q is ~1.0 by construction, because that sample *is* the prior, so
reading homogeneity off it would report the assumption back. The unsupervised
Q is both the honest reading and the one that selected it in the first place:
4a's stage A takes the top-N by *unsupervised* Q per super-pop, so anchor
homogeneity is an empirical finding — just recorded in a different column. The
manifest's `q_source` column says which Q was used per sample, and `is_anchor`
flags them.

A worked consequence: an anchor with supervised Q 1.000 but unsupervised Q 0.80
is classified **`real_admixed`**, not `real_homogeneous`.

`--anchors-to-reference` on the splitter remains available as a *preference* —
it concentrates the purest samples in the reference panel at the cost of a
smaller donor pool — not as a correctness fix. It is a **move, not an
addition**, and each anchor moved out of the donor half is swapped for a
donor-side sample from the same population so the proportional sub-population
matching survives.

### Sample types

| `sample_type` | rule |
|---|---|
| `real_homogeneous` | `max_Q ≥ threshold(argmax)` |
| `real_admixed` | below that threshold |
| `real_unassigned` | argmax is **OTH** — no reference ancestry fits |

**The threshold is keyed on the ADMIXTURE argmax, not the metadata label.**
"Is this genome homogeneous?" is about the ancestry the Q actually found, so a
sample at Q=0.97 on EUR is measured against EUR's 0.95 whatever its label
says. Label agreement is recorded separately in `label_concordant`; folding it
into the type would turn a metadata problem into a false biological claim.

`reference_support` = `full` / `thin` (<25 reference samples) / `none`. Kept
separate from `sample_type` because whether a matching reference panel exists
is a property of the **panel**, not the person — that separation is what lets
"imputation with a thin reference" be asked apart from "imputation on someone
no reference fits". **OCE and OTH are in the cohort precisely because they
can't be pooled into a stratified analysis.**

### The three designs

| `--design` | imputation runs | build effect | chip effect | sim per cell |
|---|---|---|---|---|
| `one-cell` | 1× ≈ 2,660 | between-individual | between-individual | 360 |
| **`paired-builds`** (default) | 2× ≈ 5,320 | **within-individual** | between-individual | **720** |
| `full` | 6× ≈ 15,960 | within-individual | **within-individual** | 2,160 |

The two readings are not symmetric. The **build effect** is "did the
liftover change *this person's* answer?" — a within-individual question, so
every individual goes on both builds. The **chip effect** is marker density, a
property of the array, so an ancestry-balanced between-individual comparison
is adequate — and splitting chips is what keeps the cost at 2× not 6×.

`chips` and `builds` are **comma-separated lists**; 6d tests membership, not
equality (an equality test against `"hg19,hg38"` matches nothing and would
yield empty cells rather than an error).

The deal **rotates each group's starting unit**. A fixed start sends every
group's remainder to unit 0, and a cell-size imbalance is indistinguishable
from a chip effect.

### Outputs
`benchmark_cohort_manifest.tsv` (real), `..._simset_individuals.tsv` (2,160
rows, one per simulated individual), `..._simset_roster.tsv` (counts),
`..._counts.txt` (the human-readable report).

> Per-individual rows matter: step 6 simulates and step 10 masks, so both must
> agree which individual is in which cell. Counts alone can't express that,
> and independent deals would disagree.

---

## 6. Step 6 — simulation

`6b_simulate_benchmark_cohort_clm.sh` — one task per chromosome, all 13 models
inside. Components are extracted **once** per chromosome and reused, since
`plink2 --export haps` over 3,921 samples dominates the runtime (22 tasks, not
286).

Models, generations and proportions come from the roster, not from the script,
so Table 1 lives in exactly one place.

**Gen-0 controls bypass admix-simu.** A single-ancestry individual has no
admixture to recombine, so donor haplotype *pairs* are sampled directly and a
constant-ancestry `.hanc` written. That is what makes a control's error pure
method error.

### Reproducibility is by preservation, not seeding
`mixer` draws the recombination breakpoints, takes no seed and ships without
source. So each unit archives its `.bp` — which **is** the recombination
realisation — plus `SIMU_SEED`. `FROM_BP` replays exactly; `FROM_BP` set with
a missing `.bp` is a hard error, not a silent fresh draw.

### The one assumption that can't fail loudly
The `.dat` **generation lines**. `simu-mix.pl` validates the header's
population labels, but the generation lines go to the `mixer` binary, which
accepts a differently-shaped history **without complaining** and simulates the
wrong demography. `--dat-only` prints all thirteen `.dat` files and exits,
needing no panel, conda env or Slurm context.

The format itself is settled (admix-simu README, "Dat file", read 17 Nov):
the first generation line is a pulse at generation 1 and must have
p_admixed = 0; the population then mixes with itself up to the generation on
the LAST line. So Jessica's one-line `Brasa.dat` (`12 0 0.15 0.60 0.25`) and
6b's two-line form (`1 0 p...` then `G 1 0 0 0`) describe the same
G-generation pulse model, and the README's own example says so. What stays
worth checking is the output: `validate-simset-tracts.py --hanc2` reads one
2b truth file and reports the generations its tract lengths imply
(`verify-rebuild.sh` runs it on chr22 of each cohort).

### Did it simulate what we planned?

Three things must agree, and `cohort-ancestry-summary.py` names the failing
link rather than just reporting a mismatch:

```
roster (the plan)  →  .dat (what write_dat wrote)  →  .hanc (what mixer did)
```

| verdict | meaning |
|---|---|
| **LINK 1** | the archived `.dat` ≠ the roster. `write_dat` is wrong; `mixer` simulated its input faithfully. Re-simulating without fixing `write_dat` reproduces the error |
| **PERMUTED** | the `.dat` matched the plan but the realized ancestry is a *reordering* of it — the `.dat` header order and the `-POP` argument order disagree, so `.hanc`'s ancestry indices point at the wrong components. Not a demography error |
| **LINK 2** | `.dat` matched the plan, realized means don't match the `.dat` — `mixer` didn't produce the proportions in its own input |
| **spread only** | proportions right, SD wrong = the demography is wrong while the pulse survived (it comes off the `.dat`'s *first* line) |

Two bounds worth knowing:

- **The mean's bound comes from the OBSERVED spread (`SD/√n`), not from the
  intended generations.** Deriving it from `gen` is circular — whether the sim
  ran at `gen` is the question — and data generated at the wrong generation
  count is *more* variable than `gen` predicts, so a gen-derived bound is too
  tight and reports a mean failure for a demography failure.
- **A permutation is only claimed when the identity assignment actually fails
  and some reordering fits clearly better** (identity misfit > bound, best
  reordering < bound/2). LAT3's 0.50/0.45/0.05 sit closer together than the
  noise at modest n, so sorted-value agreement alone proves nothing.

---

## 7. Step 7 (data) — real test genomes

`6c_extract_real_test_genomes_clm.sh` — subsets the phased panel to the
manifest's test-eligible individuals, one BCF per chromosome.

**It applies no policy of its own.** The four exclusions live in the
assembler; a second implementation here would drift, and the drifting copy
would put a donor in the test set. It reads column 1 and hard-fails if any
listed sample is absent from the panel rather than letting `--force-samples`
ship a smaller cohort than the manifest claims.

---

## 8. Step 10 — chip × build × chr cells

`6d_write_chip_cells_clm.sh` + `scripts/lift-cell-to-hg19.py`

Real individuals (from 6c) and simulated ones (from 6b's haplotypes) are
subset to the chip's sites, merged, and written per chip × build × chr.

**hg19 cells use the chip's NATIVE hg19 manifest**, never a down-lift of the
hg38 sites. A down-lift contains only probes that map back cleanly, by
construction, so it would hide every probe the pipeline's own CrossMap drops,
moves or strand-flips — which is part of what the benchmark measures.

### The allele-orientation problem
**GT is an index into REF/ALT.** Between builds a site can:

| class | | GT |
|---|---|---|
| `same` | alleles identical | unchanged |
| `flip` | alleles complemented (strand flip) | **unchanged** — complementing both leaves indices alone |
| `swap` | REF and ALT exchanged | **FLIPPED 0↔1** |
| `flip_swap` | complemented *and* exchanged | **FLIPPED** |
| `unresolvable` | no correspondence | dropped |

Missing `swap` writes a perfectly valid VCF with **inverted alt-allele dosage**
at exactly the sites where the builds disagree about the reference allele. An
imputation benchmark on it reports inflated error that reads as a *pipeline*
problem, not a coordinate bug.

**Palindromic sites are irreducibly ambiguous** — at A/T or C/G, complement
*equals* swap. Those are resolved from the vendor strand column alone and
**counted separately** in the report rather than absorbed. GSA: 6,422 of
596,591 (1.08%).

### Guards
- `bcftools reheader -s` renames **positionally**, so each model's `.haps`
  individual count is checked against the roster — a mismatch would shift
  every ID and file genotypes under the wrong person.
- Each finished cell's sample count is checked against the number assigned.
  A cell that silently loses individuals makes a chip effect and an ancestry
  effect inseparable.

---

## 9. Step 10b — kinship and effective N

`6e_cohort_kinship_clm.sh` + `scripts/cohort-power-summary.py`

**Both halves carry power and GWAS simulations.** ~500 real genomes alone is
too few; ~2,660 together is workable. But simulated individuals from a finite
donor pool are cryptic relatives **by construction** — the thinnest donor
pool is the binding case (SAS was once projected at ~24; that projection is
stale, and the split prints the real counts).

That need not bias anything: the GWAS uses **GENESIS with a GRM**, which models
relatedness rather than pruning it. What relatedness costs is **effective**
sample size:

```
N_eff = N² / Σᵢⱼ Rᵢⱼ        R = the GRM (2 × kinship off-diagonal)
```

**Power calculations must use N_eff, not the roster count.** On a synthetic
thin-donor fixture, 60 SASEUR individuals came out behaving like 5.

Not an array — kinship is genome-wide and per-chromosome estimates can't be
averaged. LD-pruned first; kinship on correlated markers is inflated. Drops
nobody: this is measurement, not de-duplication (that's 1e).

**It also catches a leak nothing else can:** a real test genome related to a
simulated one at second degree or closer means a donor, or a donor's close
relative, is in the test set.

---

## 10. Step 9 — chip manifests

`scripts/build-chip-manifest-sites.py` joins vendor **GRCh37 + GRCh38**
manifests on probe ID and resolves REF/ALT against each FASTA (half the array
is minus-strand, so the vendor allele column can't be trusted).

GSA v1.0: 618,540 probes → **596,591 sites**, 6,422 palindromic.

`scripts/inventory-chip-manifests.py` reports what's on disk. Two
**independent** blockers per chip:
- **both native builds** — one alone isn't buildable, and a liftover is not a
  substitute
- **a readable format** — `read_annotated()` requires tab-separated
  `Name/Chr/MapInfo/Alleles`, so an Affymetrix `annot.csv` can't be read at all

---

## 11. Transfer — `6a_transfer_to_mab_clm.sh`

```
/mount/genepi2/Joint_called_hgdp_1kg_plus_mxbiobank_wgs_ref_panel/
├── lai_reference_panel/                        donors INCLUDED
│   ├── lai_full_panel_samples.txt              ← FOR FLARE + THE GWAS
│   ├── <sp>_reference_lai.txt / _simulation_lai.txt   ← LAI accuracy paper
│   └── lai_ref_panel_samples.k8.tsv            full Q table
└── genotype_harmonization_pipeline_testing/
    ├── README.txt
    ├── benchmark_reference_panel/               donors EXCLUDED
    └── benchmark_cohort/{real,simulated}/       the test genomes
```

**Three sample lists in the LAI tree are not interchangeable.**
`<sp>_reference_lai.txt` is about **half** its super-pop (the split's reference
half). `lai_full_panel_samples.txt` is the union — every homogeneous sample,
donors and 1KG included. Reading the wrong one costs half the panel with no
error.

`homog_8pop_samples.txt` is 4a's `included == 1`, and 4a drops anything in an
**unversioned** `reference_ids/*_simulation.txt` — a glob the suffixed lists
don't match. So whether donors are in that file depends on leftover files; the
transfer log reports when that happened.

Scopes: `help`, `verify` (no credentials needed), `panels`, `benchmark`, `all`.

---

## 12. Every check — what it catches, and what it can't

| check | where | catches | blind to |
|---|---|---|---|
| **panel distinction** | `6a verify` | the suffix collision; donor-in-reference; 1KG donors. Runs **before any bytes move** | — |
| **label vocabulary** | assembler | a cohort's super_pop outside the canonical set — which would silently take the relaxed threshold, read reference size 0, and never match the argmax | a *canonical* but wrong label (reported as discordance instead) |
| **label concordance** | assembler | a super-pop discordant for most of its samples = unharmonised cohort | — |
| **soft-union sample count** | `pipeline-status.sh` | chunks phasing the old 3,448 panel — available ~10 h before ligate | — |
| **chunk staleness** | `1b` | a chunk left from a previous panel | — |
| **tract lengths** | `validate-simset-tracts.py` | the `.dat` terminal line ignored → every model runs as one pulse | **gen-2 models** (gen 2 *is* the failure mode) |
| **ancestry mean + SD** | `cohort-ancestry-summary.py` | same failure via a **different statistic**. The **means stay correct** under it — the between-individual SD is what carries the signal | gen-2 models |
| **roster → .dat → .hanc chain** | `cohort-ancestry-summary.py` | *which* link broke: `.dat` ≠ roster blames `write_dat`; realized ≠ `.dat` blames `mixer`; a reordering blames the component order | nothing, when the `.dat` was never archived — it says so |
| **allele orientation** | `lift-cell-to-hg19.py` | `swap` → inverted dosage on hg19 cells | palindromic sites (strand-resolved, counted separately) |
| **positional reheader** | 6d, 6e | a model's `.haps` count ≠ roster → every ID shifted | — |
| **cell sample count** | 6d | a cell that lost or gained individuals | — |
| **real↔simulated kinship** | `6e` | a donor or close relative in the test set | — |
| **N_eff** | `cohort-power-summary.py` | donor reuse inflating nominal N | — |
| **manifest ⊆ panel** | 6c | manifest built against a different panel | — |
| **anchor sufficiency** | 4a | a super-pop whose supervised prior rests on fewer than `ANCHOR_FLOOR` exemplars — ADMIXTURE runs regardless, so nothing else would say so | nothing; it writes `anchor_summary.k8.tsv` either way |
| **supervised anchors** | assembler | a 1KG anchor used as a donor (anchor status is no exemption from the no-1KG rule); reports anchor counts and roles | an empty `assigned` column, which it reports rather than assumes |

Two checks deliberately use **different statistics over the same truth**
(switches per Morgan vs fractions per individual) so a bug in one doesn't
reproduce in the other. Both are blind to gen-2 models, which is the argument
for keeping gen 7/16/30 in Table 1.

---

## 13. Run order after 4a lands

```bash
# verify 4a
grep -E 'named exclusions|PRESENT in the panel' $P/logs/4a_homog*.out
grep -c CHMI $P/04_homogeneity_panel/lai_ref_panel_samples.k8.tsv   # must be 0

# split (BOTH modes, DISTINCT suffixes)
python3 -I scripts/split-homog-donors-reference.py --panel $T --mode lai       --suffix lai
python3 -I scripts/split-homog-donors-reference.py --panel $T --mode benchmark --suffix bench \
        --no-donor-superpops OCE --anchors-to-reference

# assemble + report
python3 -I scripts/assemble-benchmark-cohort.py --panel $T --refs reference_ids \
        --suffix bench --design paired-builds --out results/benchmark_cohort
python3 -I scripts/cohort-composition-report.py --panel $T --refs reference_ids --suffix bench

# build
sbatch 6b_simulate_benchmark_cohort_clm.sh     # 22 tasks, simulate
sbatch 6c_extract_real_test_genomes_clm.sh     # 22 tasks, real
sbatch 6d_write_chip_cells_clm.sh              # 22 tasks, cells
sbatch 6e_cohort_kinship_clm.sh                # 1 job, genome-wide

# transfer
./6a_transfer_to_mab_clm.sh verify
./6a_transfer_to_mab_clm.sh all
```

---

## 14. Known gaps

| gap | state |
|---|---|
| **Step 11** — tables, checksums, provenance | not built |
| **Affy 6.0 / 500K site lists** | need a GRCh38 vendor manifest **and** an Affymetrix reader |
| **GSA v2.0** | needs its vendor GRCh38 A1 manifest |
| **`.dat` generation-line format** | settled from the admix-simu README (17 Nov); still checked from output |
| **6d / 6e bcftools pipelines** | untested against real data (orientation logic *is* unit-tested) |
| **FLARE invocation** | untested against the real jar |
| **SAS decision** | deferred to real counts: reduce n, accept reuse and report N_eff, or control-only |
| **`SIMSET_CONTRACT.md`** | the copy in hand is **v1** (4 models, one build, no real genomes); Doc 20's appendix is v2 |
