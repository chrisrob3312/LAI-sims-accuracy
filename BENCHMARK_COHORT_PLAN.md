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
| `real_homogeneous` | top ADMIXTURE component ≥ 0.93 | real rare variation and real haplotype structure that simulation cannot produce; includes OCE and OTH, which have no simulated model at all |
| `real_admixed` | top component < 0.93 — naturally admixed | real admixture, with the correlation structure between ancestry and genotype that a simulation idealises away |

**They have no tracts.** A real genome has a WGS truth, but nobody recorded
which ancestry each of its segments came from. Inferred tracts would be
another method's output, not truth. So real genomes are scored for
**imputation accuracy, calibration and retention only — never for
local-ancestry accuracy.**

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

**OCE.** ~20 donors and ~20 reference. Likely answer: **no OCE-containing
admixture model at all** — OCE enters as `real_homogeneous` test samples and as
an ADMIXTURE/LAI reference group, and is never simulated. Doc 19 already
assumes this: OCE and OTHER have no reference panel and no model.

Neither is defaulted. Both are decided when stage 3 prints the counts.

## 6. Build order

Steps 1–4 are what this repository produces. Steps 5–13 run in
`genotyping_pipeline/benchmarking/simulate/`.

| step | what | where | state |
|---|---|---|---|
| 1 | assemble the real-genome pool on GRCh38 | 1a / 1b / 1d | SGDP + Almarri done, MXB done |
| 2 | de-duplicate people across sources (KING ≥ 0.354) | — | **not built** |
| 3 | ADMIXTURE K = 8 → the Q file everything reads | 4a | pending stage 3 |
| 4 | phase the donor panel jointly | 1c | pending stage 2 |
| 5 | plan donor roles and model sizes | `plan_simset.R` + our splitter | splitter ready |
| 6 | simulate admixed + control individuals | admix-simu | — |
| 7 | add real test genomes to the truth | — | — |
| 8 | deal everyone into the six chip × build cells | — | — |
| 9 | chip site lists in both builds | `build-chip-manifest-sites.py` | **GSA done** (596,591 sites); Affy 6.0 and 500K pending |
| 10 | write the chip files per chip × build × chr | — | — |
| 11 | tables, checksums, provenance | — | — |
| 12 | `ingest_simset.sh` intake check | — | — |
| 13 | hand over | — | — |

**Step 2 has no owner yet.** SGDP re-sequenced individuals who are also in
HGDP, and the same person can appear under two IDs across HGDP, 1KG, MXB and
Almarri. If one copy is a donor and the other a test sample, the method has
seen the test sample's haplotypes. `1c` deduplicates by *sample ID* at the
merge, which catches identical identifiers but not the same person under two
different ones. A KING run over the merged panel is what catches that, and it
is not written.
