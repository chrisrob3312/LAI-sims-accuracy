# Building the benchmark cohort: who is in it, and how they are made

Companion to `BENCHMARK_PANEL_DESIGN.md`. That file says who may be a donor
and who may be a reference. **This file says what is actually built from
them** — the admixture models, the real genomes that join them, and the chip ×
build cells.

Source of truth: Document 19 (Table 1, the models) and Document 20 (the
thirteen build steps), contract v2.

> **Version warning.** The `SIMSET_CONTRACT.md` currently in hand is **v1**:
> four models, one build, no real test genomes. Document 20's appendix
> describes **v2**. Building against v1 silently omits CAR3, AFAM2, MENEUR,
> all six controls, every hg19 cell and every real genome. Confirm which
> version `ingest_simset.sh` enforces before building.

## 1. The simulated individuals — 2,160 of them

| model | components | proportions | gen | n | per chip × build cell | focal | priority |
|---|---|---|---|---|---|---|---|
| **LAT3** | NAT, EUR, AFR | 0.50 / 0.45 / 0.05 | 16 | **600** | 100 | NAT | primary |
| **CAR3** | NAT, EUR, AFR | 0.10 / 0.50 / 0.40 | 16 | 200 | 33–34 | AFR | secondary |
| **AFAM2** | AFR, EUR | 0.80 / 0.20 | 7 | 200 | 33–34 | AFR | secondary |
| **MENA3** | EUR, MEN, AFR | 0.30 / 0.50 / 0.20 | 30 | 200 | 33–34 | MEN | secondary |
| **MENEUR** | MEN, EUR | 0.50 / 0.50 | 2 | 200 | 33–34 | MEN | secondary |
| **EASEUR** | EAS, EUR | 0.50 / 0.50 | 2 | 200 | 33–34 | EAS | secondary |
| **SASEUR** | SAS, EUR | 0.50 / 0.50 | 2 | 200 | 33–34 | SAS | secondary |
| NAT1, EUR1, AFR1, EAS1, SAS1, MEN1 | one each | 1.0 | 0 | 60 each | 10 | that one | control |

1,800 admixed + 360 single-ancestry controls = **2,160**, exactly 360 per
chip × build cell.

**LAT3 is the Mexican-like cohort** and the primary model: half Indigenous
American, with the small African fraction a Mexican cohort really carries.
**CAR3 is the Caribbean / Dominican-like one**: little NAT, heavy AFR and EUR.
The two share components and generations and differ only in proportions,
which is deliberate — it isolates *proportion* as the variable, so a
difference between them is not confounded with a difference in admixture
timing.

**The controls are not filler.** A single-ancestry individual with zero
generations has no admixture to get wrong, so any ancestry error on it is pure
method error, and its imputation accuracy is the ceiling the admixed
individuals of that component are measured against.

**The generations differ on purpose**: 2 for the recent two-way models, 7 for
AFAM2, 16 for the Latino models, 30 for MENA3. Tract length falls as
generations rise, so this spans the range over which local-ancestry callers
degrade — a method that looks fine at gen 2 can fail at gen 30.

## 1b. Which reference goes with which reading

The two reference panels this repository ships are for **supervised ADMIXTURE
and local-ancestry inference only.** Imputation in the `genotyping_pipeline`
flow runs against the **external servers — TOPMed, All of Us and Michigan
(1000 Genomes)** — never against a panel built here.

| reading | reference used |
|---|---|
| global ancestry Q (ADMIXTURE) | `benchmark_reference_panel/` |
| local ancestry (FLARE / RFMix) | `benchmark_reference_panel/` |
| imputation accuracy, calibration, retention | TOPMed / All of Us / Michigan |
| power and GWAS simulations | real genotypes + the ADMIXTURE Q as covariate |

This is what makes the 1000 Genomes rules look inconsistent until you see
which role each is about:

- **1KG is in `benchmark_reference_panel/`, and that is correct.** The panel
  never imputes, so 1KG in it cannot leak a test genome's haplotypes into an
  imputation answer. For ADMIXTURE and LAI it is simply more reference data.
- **1KG is never a donor and never a real test genome.** The external servers
  all contain 1KG, so such a genome would be imputed against its own
  haplotypes. Same reason `sgdp_server_panel_ids.txt` is barred from both
  roles: SGDP re-sequenced 1KG individuals under other identifiers, and the
  servers hold those people whatever our cohort column says.
- **Donors are excluded from `benchmark_reference_panel/`** for a third,
  separate reason: the simulated test genomes are built *from* donor
  haplotypes, so a donor left in the LAI reference hands FLARE and RFMix the
  haplotypes they are supposed to infer.

Three exclusions, three distinct failures. And note that the imputation
readings rest **entirely** on the test samples being absent from the external
servers — nothing about the panels in this repository can establish that, only
the exclusions can.

`lai_reference_panel/` is the same kind of panel with the opposite donor
policy: donors included, for the GWAS and FLARE work where nothing is
simulated from these haplotypes.

## 2. Where each component's haplotypes come from

Donors are the `--mode benchmark` donor half from
`split-homog-donors-reference.py`, so **no 1000 Genomes sample is ever a
donor** (the imputation servers' panels contain 1000 Genomes; a 1KG donor
would be imputed against its own haplotypes).

| component | donor populations | projected donors |
|---|---|---|
| **NAT** | MXB (Nahua, Zapotec/Mazatec, Maya, Totonac, Tzotzil, Tarahumara, Huichol) + HGDP (Karitiana, Surui, Pima, Colombian) + SGDP | ~46 |
| **EUR** | HGDP (French, Sardinian, Basque, Tuscan, Orcadian, Russian, Adygei) | ~59 |
| **AFR** | HGDP (Yoruba, Mandenka, Biaka, Mbuti, Bantu, San) | ~47 |
| **MEN** | Almarri (Emirati, Omani, Iraqi, Saudi, Jordanian…) + HGDP (Bedouin, Druze, Palestinian, Mozabite) + SGDP | ~251 |
| **EAS** | HGDP (Han, Japanese, Dai, Tujia, Yi, Miao, She, Lahu, Naxi) | ~106 |
| **SAS** | HGDP (Balochi, Brahui, Sindhi, Pathan, Burusho, Hazara, Makrani, Kalash) | **~24** |

MXB sits on **both** sides for NAT — half donors, half reference — so the
Mexican-like cohort is simulated from Mexican haplotypes *and* tested against
a reference panel that represents them. The 25/25 split is already built and
stratified by genetic cluster (`split_mxb_samples.R`).

**SAS at ~24 donors is the binding constraint.** 48 haplotypes for 200 SASEUR
individuals is heavy reuse, which makes simulated SAS individuals cryptic
relatives of each other — the intake's `7.sim_kinship` check is what detects
it. See §5.

## 3. The real genomes — test samples, not donors

Real individuals from HGDP, SGDP, MXB and Almarri join the cohort as test
samples. They are **outside every imputation server panel**, so no server has
seen them.

| sample_type | who | what it adds |
|---|---|---|
| `real_homogeneous` | top Q ≥ threshold **and** argmax == assigned super-pop | real rare variation and real haplotype structure that simulation cannot produce |
| `real_admixed` | below the threshold, or argmax disagrees with the label | real admixture, with the correlation structure between ancestry and genotype that a simulation idealises away |
| `real_unassigned` | argmax is **OTH**, the empirical leftover cluster at K=8 | the multiracial-like individuals a real diverse cohort carries. No reference ancestry fits them, so they are neither homogeneous nor admixed between named components |

Thresholds are **tiered**, matching what 4a applies to the panel: 0.95 for
EUR/EAS/AFR, 0.93 for AMR/OCE/MEN/SAS. Classifying test genomes on a flat
0.93 while the panel uses a tiered bar would admit samples to
`real_homogeneous` that the panel itself rejects.

**Eligibility is NOT `included == 1`.** 4a sets that only when
`max_Q >= threshold` AND the argmax matches the assigned super-pop, so gating
on it admits only homogeneous samples — which makes `real_admixed` empty by
construction and drops every OTH sample, since 4a excludes the leftover
cluster from the homogeneous panel on purpose. Eligibility is decided from the
role lists instead (not a donor, not in either reference panel, not 1KG, not
in `sgdp_server_panel_ids.txt`, not a named exclusion, not a KING duplicate
drop), and `included` is carried as a reported column.

**Include every OCE and OTH sample not in a reference panel.** These are the
groups that cannot be pooled into a stratified analysis, so they have to be
tested directly rather than inferred from a neighbouring super-pop. A separate
`reference_support` column records `full` / `thin` (< 25 reference samples) /
`none`, because whether a matching reference panel exists is a property of the
PANEL, not of the person — keeping them apart is what lets "imputation with a
thin reference" be asked separately from "imputation on someone no reference
fits".

**They have no tracts.** A real genome has a WGS truth, but nobody recorded
which ancestry each of its segments came from. Inferred tracts would be
another method's output, not truth. So real genomes are **never scored for
local-ancestry accuracy.**

What they ARE used for:

1. **Imputation accuracy, calibration and retention.**
2. **Power and GWAS simulations** — phenotypes simulated on real genotypes.

The second is not a side use, and it is a reason the real half cannot be
replaced by more simulated individuals. Power depends on the LD structure, the
allele-frequency spectrum and the relatedness structure that the genotypes
actually carry. Simulated genomes drawn from a finite donor pool have
compressed haplotype diversity and carry donor-reuse-induced cryptic
relatedness — see §5, where SAS reuse is the binding constraint — so power
estimated on them would be optimistic in a way that is hard to quantify and
impossible to correct after the fact. Real genomes are where a power estimate
can be believed.

That has a practical consequence for the transfer: the real half has to ship
as **genotypes plus the per-sample ADMIXTURE Q**, since Q is the ancestry
covariate a GWAS simulation conditions on. `benchmark_cohort_manifest.tsv`
carries `max_Q_sup`, `argmax_pop_sup` and `sample_type`, and the full Q table
travels as `lai_ref_panel_samples.k8.tsv`.

This is what absorbs the thin-donor problem: SAS and OCE can lean on real
genomes for the imputation readings, and only their *local-ancestry* precision
is limited by donor count.

## 4. Chip × build assignment — six cells, balanced

Three chips (`snp6`, `affy500k`, `gsa`) × two builds (hg19, hg38) = six cells.

- Every individual is on **exactly one chip and one build**. No one is on two.
- Every model, and every real `sample_type` × top-component group, is split
  **equally across all six cells**, within 1.
- hg19 cells use the chip's **native hg19 manifest**, never a down-lift of the
  hg38 sites — that is what exercises the pipeline's own liftover. A down-lift
  would contain only the sites that map back cleanly, by construction, and
  would hide every probe the pipeline's CrossMap drops, moves or strand-flips.
- The truth stays GRCh38 for everyone; `grch38_site` in the hg19 manifest maps
  each hg19 probe back to the truth site it stands for.

Balance matters because the arms harmonize separate platforms. If a group
concentrated on one chip or one build, a chip effect would read as an ancestry
effect, and the reverse.

## 5. The two decisions the real counts will force

**SAS.** ~24 donors for 200 SASEUR individuals. Options, in preference order:
reduce SASEUR's n; accept the reuse and report the kinship distribution; or
drop SASEUR to control-only and rely on real admixed SAS for the imputation
readings. The LAI reading is the one that suffers either way.

**OCE.** **Decided: no OCE-containing admixture model.** Doc 19 already assumes
it, and Table 1 has no OCE model, so every OCE sample withheld as a donor was
withheld for nothing. Run the split with `--no-donor-superpops OCE`: the
reference takes its planned share and the remainder is left *unallocated*, so
the assembler picks it up as test genomes.

Note the sizing trap this exposed — `min(n, ref_target)` would take all ~20
OCE as reference and leave **zero** OCE test samples, which is the gap putting
OCE in the cohort was meant to close. The reference-only path uses the same
`plan_sizes` policy as a split group with no donor floor, so ~20 OCE gives
15 reference / 5 test. That is thin either way and should be reported as a
limitation rather than smoothed over: the OCE imputation reading rests on a
handful of genomes.

**SAS is not defaulted.** It is decided when stage 3 prints the counts.

## 6. Build order

Steps 1–4 are what this repository produces. Steps 5–13 run in
`genotyping_pipeline/benchmarking/simulate/`.

| step | what | where | state |
|---|---|---|---|
| 1 | assemble the real-genome pool on GRCh38 | 1a / 1b / 1d | SGDP + Almarri done, MXB done |
| 2 | de-duplicate people across sources (KING ≥ 0.354) | `1e_king_dedup_clm.sh` | **built**, chained before 4a |
| 3 | ADMIXTURE K = 8 → the Q file everything reads | 4a | pending stage 3 |
| 4 | phase the donor panel jointly | 1c | pending stage 2 |
| 5 | plan donor roles and model sizes | `split-homog-donors-reference.py --mode benchmark` | **built** |
| 6 | simulate admixed + control individuals | admix-simu | **not built** |
| 7 | add real test genomes to the truth | `assemble-benchmark-cohort.py` | **built** |
| 8 | deal everyone into the six chip × build cells | `assemble-benchmark-cohort.py` | **built** |
| 9 | chip site lists in both builds | `build-chip-manifest-sites.py` | **GSA done** (596,591 sites); Affy 6.0 and 500K pending |
| 10 | write the chip files per chip × build × chr | — | **not built** |
| 11 | tables, checksums, provenance | — | **not built** |
| 12 | `ingest_simset.sh` intake check | — | — |
| 13 | hand over | — | — |

Steps 2, 5, 7 and 8 need only the Q file, so all four run within minutes of
stage 3 landing. Steps 6, 10 and 11 are the remaining gap, and 6 and 10 belong
in `genotyping_pipeline/benchmarking/simulate/` rather than here.

**Step 2 is built.** SGDP re-sequenced individuals who are also in HGDP, and
the same person can appear under two IDs across HGDP, 1KG, MXB and Almarri. If
one copy is a donor and the other a test sample, the method has seen the test
sample's haplotypes. `1c` deduplicates by *sample ID* at the merge, which
catches identical identifiers but not the same person under two different ones.
`1e_king_dedup_clm.sh` runs KING-robust genome-wide over the merged panel
between the merge and ADMIXTURE, drops only pairs above 0.354, and keeps the
non-1KG copy so the surviving record stays donor-eligible.

**The same four exclusions apply to test samples as to donors.** A real genome
is test-eligible only if it is not a donor, not in the reference panel, not
1000 Genomes, and not in `sgdp_server_panel_ids.txt` — an SGDP cohort label
does not establish absence from an imputation server's panel, which is the
donor-leak class that file exists to close. `assemble-benchmark-cohort.py`
enforces all four and refuses to run if the split is missing, because without
it every panel sample looks test-eligible.
