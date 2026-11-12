# Benchmark cohort & reference panel — build workflow

Reference for how the panel and the benchmarking cohort are built, what every
script does, and what each check catches. Companion to
`BENCHMARK_COHORT_PLAN.md` (who is in the cohort) and
`BENCHMARK_PANEL_DESIGN.md` (who may be a donor vs a reference).

Everything runs on **taco (mhgcp)**. Branch: `claude/setup-shapeit5-conda-xlcdl`.

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
| 1a | `1a_prep-mxb-liftover_clm.sh` | MXB lifted to GRCh38 |
| 1b | `1b_sgdp_liftover_clm.sh` | SGDP lifted to GRCh38 |
| 1d | `1d_almarri_prep_clm.sh` | 137 Almarri per-chr BCFs |
| 1c | `PHASE_MODE=chunked` `1c_phasing-jointcall_clm.sh` | 4-cohort merge + soft-union sites (**3,921 samples**) |
| 1b | `1b_phase_chunks_clm.sh` → `1b_phase_ligate_clm.sh` | phased panel (594 chunks, 8 cM / 2 cM overlap) |

`submit_panel_rebuild.sh` submits 1d → 1c → chunks → ligate → 1e → 4a as one
afterok chain (SGDP prep and, optionally, the MXB liftover ahead of 1c), and
prints the command that chains the publication stages 4–6 to its 4a job.

`1b_phase_prep_clm.sh` (and `submit_chunked_phasing.sh`, which submits it) is
the **2-cohort** July prep — HGDP+1KG + MXB only — and writes the same
`merged_chrN.softunion.bcf` path 1c publishes. It now refuses to run unless
`ALLOW_TWO_COHORT_PREP=1`.
| 1e | `1e_king_dedup_clm.sh` | duplicate drops (KING > 0.354 on called sites, + known same-person pairs) |
| 4a | `4a_build_homogeneity_panel_kN_clm.sh` | **ADMIXTURE K=8 Q table** — everything downstream reads this |
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
4. `bcftools +fixref -- -f hg38 -m flip -d` (flip strand where needed, drop
   what cannot be fixed)
5. sort, split → `01_merged_phased_panel/mxb_lifted.chrN.bcf`

This is the **liftover control**: the same recipe gave all 50 MXB samples
Q = 1.000 in July.

**1b_sgdp — `1b_sgdp_liftover_clm.sh`** (SGDP — **NEW**)
1. download `chr.sgdp.pub.N.bcf` (hg19, phased) to `00_sgdp_raw/` (cached;
   size and header checked)
2. `bcftools view -S sgdp_combined_candidates.txt --force-samples` (277
   requested, 276 present) → `bcftools norm -m -any`
3. contigs → chr-prefixed
4. Picard `LiftoverVcf` — **identical recipe to 1a**; rejects →
   `sgdp_prepped.chrN.rejected.vcf.gz`, reject % in the log
5. `bcftools norm -f hg38 -m-` → `+fixref -m flip -d` → keep `chrN` → sort →
   `sgdp_prepped.chrN.bcf`

Different from MXB in the INPUT, not the recipe: a released, already-phased
dataset from another caller. A site this release lacks is `./.` for all 276
after the merge.

**1d — `1d_almarri_prep_clm.sh`** (Almarri — **NEW**; already GRCh38, no liftover)
1. per sample, per chr: `bcftools view -f PASS -i 'FMT/GQ>20'` → `norm -m -snps
   --fasta-ref` → biallelic SNPs → REF check
   - **November panel (recipe v1):** `norm -c ws` — a mismatched REF is
     *overwritten* from the FASTA while ALT and the genotypes are kept, and
     stderr was discarded, so how often that happened is not in the logs.
   - **Script now (v2), on any re-run:** `norm -c w` + `+fixref -- -f GRCh38
     -m flip -d` — the same check MXB and SGDP get; stderr kept; per-sample
     cache tagged `v2` so v1 files are not reused.
2. `bcftools merge -m none` (137 samples) → **`+missing2ref`**: every absent
   genotype becomes 0/0
3. strip to GT → `almarri_prepped.chrN.bcf` (393,138 records); asserts ~0
   missingness

A call removed by `GQ>20`, or a site with no coverage, becomes **hom-ref, not
missing**: an uncertain het is recorded as 0/0. That can lower Almarri
heterozygosity relative to cohorts whose uncertain calls stay missing.

**1c — `1c_phasing-jointcall_clm.sh` with `PHASE_MODE=chunked`** (merge + site
selection). In July this step was `1b_phase_prep_clm.sh`: the **same** merge →
QC → soft-union logic with **2 inputs** (HGDP+1KG + MXB). Now **4 inputs**.
1. HGDP+1KG: drop `related_outliers.txt` + `panel_exclude_ids.txt`
   (`bcftools view -S ^drop`) → 3,398 + MXB. Job 3742123 resolved the
   exclusion list to Slurm's spool directory, found nothing, and kept
   `CHMI_CHMI3_WGS2`; 1c now stops instead, and 4a / the benchmark exclude
   CHMI by name. Then `norm -m -any -c w` + `+fixref -m flip -d`: every
   multi-allelic record split (was `-m -snps`, which left SNP+indel records
   whole), REF checked against GRCh38.
2. add SGDP, then Almarri: chr-prefix contigs, drop IDs already present or
   excluded → **3,921** in November. Then **every input's remaining
   multi-allelic records are split** (`norm -m -any`): `merge -m none` joins
   records only on identical REF+ALT, so a SNP allele inside a cohort's
   multi-allelic or SNP+indel record never met another cohort's biallelic
   record, and with the fill below that cohort would read 0/0 on it.
3. Site provenance logged (from these files, before any fill). Then the
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

**4a — `4a_build_homogeneity_panel_kN_clm.sh`** (**K=8 now, K=5 in July**)
1. concat the 22 phased chromosomes → `admixture_k8/all_chr.phased.bcf`
2. restrict to 1c's ADMIXTURE site list (`--extract bed1`), then
   `plink2 --maf 0.01 --indep-pairwise 500 50 0.1` over all samples at once
   (same as July) → `--max-alleles 2 --snps-only just-acgt --make-bed` →
   `admixture_k8/lai_ref` (528,876 SNPs in November, without the site list)
3. unsupervised ADMIXTURE K=8 → `admixture_k8/unsup/lai_ref.8.Q` (~10 h;
   reused on a rerun unless deleted)
4. **NEW — data-first discovery:** per cluster, rank every eligible sample,
   keep the top 50 at Q ≥ 0.93, and name the cluster from those anchors'
   metadata (see §3, 4a)
5. supervised ADMIXTURE K=8 on those anchors → labels → tiered homogeneity →
   `lai_ref_panel_samples.k8.tsv`

July had no step 4: the curated `*_rfmix.txt` lists were the supervised labels
directly, and the unsupervised run only reported proportions.

#### Where a cohort difference can enter the genotypes

| step | cohorts | mechanism | what tests it |
|---|---|---|---|
| source calling | all four differ | caller, filters, coverage | same person twice (duplicate pairs): genotype discordance, Q difference |
| liftover (1a, 1b_sgdp) | MXB, SGDP | strand / allele errors | duplicate pairs: 0/0 ↔ 1/1 swaps; MXB is the working control |
| `GQ>20` + `missing2ref` (1d) | Almarri | uncertain ALT calls → 0/0 | pre-phasing het rate by cohort |
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

### 4a — ADMIXTURE K=8
Three stages: anchor discovery across all 7 pools → supervised ADMIXTURE at
K=8 → homogeneity filter.

**Data-first anchor discovery** (stage A). For each of the K unsupervised
clusters: rank **every** eligible panel sample by its loading on it, keep the
top 50 at Q ≥ 0.93, and **name the cluster from who those anchors are** — the
plurality metadata super-pop among them, holding ≥ 50% of the set. If two
clusters claim one super-pop, the one with more supporting samples keeps it; a
single unnamed cluster becomes OTH. Eligible = everyone except named
exclusions and 1e's duplicate drops. No candidate pools and no a-priori
naming: the curated `*_rfmix.txt` lists are not read (their overlap with the
discovered anchors is reported).

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

**Tiered homogeneity thresholds:** 0.95 for EUR/EAS/AFR, 0.93 for
AMR/OCE/MEN/SAS.

**Anchor selection — two floors, not one.** Stage A takes each super-pop's
draw-pool members with cluster Q **≥ `ANCHOR_MIN_Q` (0.93)**, sorted
descending, **capped at `ANCHOR_CAP` (50)** — the "top 25–50 ≥ 0.93" rule from
the NAT/Mexican accuracy pipeline. Two separate floors matter, and they fail
differently:

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
- `--no-donor-superpops OCE`: reference only, **no donors**, remainder left
  *unallocated* so the assembler takes it as test genomes. Table 1 has no OCE
  model, so an OCE donor is withheld for nothing. OCE ends up **~15 reference
  + ~5 real test, 0 donors** — never simulated.

### Fill order: 1000 Genomes first

**1KG can only ever be reference** — barred from the donor side because the
imputation servers contain it. A non-1KG sample can be either. So the
single-use resource is spent first and the dual-use one only on what 1KG
cannot cover. That is what `--ref-share default=0` expresses, and it explains
the contract's share values rather than leaving them arbitrary:

| super-pop | 1KG coverage | share |
|---|---|---|
| AFR, EUR, EAS, SAS | good | `0` — 1KG fills the reference |
| **MEN** | none — 1KG has no Middle Eastern population | `0.33` |
| **NAT** | only *admixed* American (PEL/MXL/CLM/PUR), mostly below 0.93 | `NAT@MXB=0.5` |
| **OCE** | none | all reference is non-1KG |

### Sub-population coverage, on ancestry GROUPS not labels

`--ref-min-per-pop 3` guarantees three reference samples per uncovered
population, because 1KG carries none of the deeply divergent lineages: no
Biaka, Mbuti, San or Mandenka in AFR, no Basque, Sardinian or Orcadian in EUR,
none of the HGDP minority populations in EAS.

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
donor pool are cryptic relatives **by construction** — SAS at ~24 donors for
200 SASEUR individuals is the binding case.

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
| **`.dat` generation-line format** | assumption; validated from output, not from the format spec |
| **6d / 6e bcftools pipelines** | untested against real data (orientation logic *is* unit-tested) |
| **FLARE invocation** | untested against the real jar |
| **SAS decision** | deferred to real counts: reduce n, accept reuse and report N_eff, or control-only |
| **`SIMSET_CONTRACT.md`** | the copy in hand is **v1** (4 models, one build, no real genomes); Doc 20's appendix is v2 |
