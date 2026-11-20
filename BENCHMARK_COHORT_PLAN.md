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
| ~~SASEUR~~ | SAS, EUR | 0.50 / 0.50 | 2 | 200 | — | SAS | dropped 18 Nov (§5) |
| NAT1, EUR1, AFR1, EAS1, MEN1 (~~SAS1~~ dropped 18 Nov) | one each | 1.0 | 0 | 60 each | 10 | that one | control |

**Agreed 18 Nov: about 2,350 genomes in all** — 1,600 simulated
(`reference_ids/models.tsv`: LAT3 300, no SAS models) plus roughly 750 real
test genomes (every eligible genome outside 1000 Genomes: 1,391 of the 4,557 in
the K=8 table are HGDP/SGDP/Almarri/MXB, less duplicates, server-panel ids,
donors and their reference). Both halves carry every reading — imputation
accuracy, retention, power — and the power analysis uses the donor-based
N_eff. The original plan, in the table above: 1,800 admixed + 360
single-ancestry controls = **2,160** simulated, exactly 360 per
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

## 1a. Where it lands on mab

Both benchmark trees ship under one parent, kept apart from the GWAS panel
beside it, because their donor policies are opposite and nothing inside a BCF
says which one produced it:

```
/mount/genepi2/Joint_called_hgdp_1kg_plus_mxbiobank_wgs_ref_panel/
├── lai_reference_panel/                        donors INCLUDED — GWAS + FLARE
└── genotype_harmonization_pipeline_testing/
    ├── README.txt                              scoring rules, the 3 exclusions
    ├── benchmark_reference_panel/               donors EXCLUDED — ADMIXTURE/LAI
    └── benchmark_cohort/
        ├── real/                                real test genomes
        └── simulated/                           the 2,160 + truth tracts
```

The parent carries a `README.txt` so the folder explains itself to whoever
opens it, rather than depending on this being remembered. Override the parent
name with `HARMONIZATION_SUBDIR`.

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

## 1c. Two independent checks on the simulated demography

The `.dat` generation lines are the one assumption in step 6 that cannot fail
loudly: `mixer` simulates whatever history it parses and says nothing if that
is not the history intended. Two checks run at the end of every chromosome in
6b, using **different statistics over the same truth** so a bug in one does
not reproduce in the other.

| script | statistic | catches |
|---|---|---|
| `validate-simset-tracts.py` | ancestry **switches** per Morgan | tract lengths inconsistent with intended generations |
| `cohort-ancestry-summary.py` | ancestry **fractions** per individual | mean vs target, and between-individual spread |

**The ancestry means are NOT a sufficient check, and this is the trap.** The
pulse proportions come off the `.dat`'s first line and survive even when the
terminal generation line is ignored. So a cohort where every model silently ran
a single generation shows means matching target to within sampling noise. The
signal is in the **between-individual SD**:

$$\mathrm{SD} \approx \sqrt{\frac{p(1-p)}{2(g-1)L}}$$

Per chromosome (L ≈ 1.4 M) that is ≈ 0.30 at g = 2 against ≈ 0.077 at g = 16 —
a factor of four. Verified on synthetic fixtures: a cohort that ran one pulse
instead of its intended generations shows **correct means with SDs 4–5× too
wide**, flagged uniformly across every high-generation model.

Both checks are blind to the gen-2 models (MENEUR, EASEUR, SASEUR), since gen
2 *is* the failure mode. The gen 7/16/30 models are what make it detectable.

The real half has no truth tracts, so its ancestry comes from the supervised
ADMIXTURE Q and is reported for **inspection, not pass/fail** — whether each
group looks like the ancestry it claims, and whether the cohort spans the
space the models were chosen to cover.

## 2. Where each component's haplotypes come from

Donors are the `--mode benchmark` donor half from
`split-homog-donors-reference.py`, so **no 1000 Genomes sample is ever a
donor** (the imputation servers' panels contain 1000 Genomes; a 1KG donor
would be imputed against its own haplotypes).

| component | donor populations | projected donors (STALE -- see below) |
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

**These projected counts are stale (17 Nov 2026) and should not be quoted.**
They were projected before the current panel: before SAS was judged beyond
July's 48 by-name anchors (July's table labelled every other Central/South
Asian sample NA -- 504 1KG SAS and 131 HGDP CSA never judged), before SGDP and
Almarri were in the fit, and before the 0.94 cut and K=7. The 16 Nov K=8 run
already had far more homogeneous SAS and fewer homogeneous MEN than this
table. The real counts come from the split on the K=7 table:
`run_benchmark_cohort.sh` prints them, and `split-homog-donors-reference.py
--mode benchmark ... --dry-run` previews them on any table without writing.
Benchmark donors are the NON-1KG homogeneous samples (AFR/EUR/EAS/SAS take
their reference from 1KG), so the SAS donor pool is the homogeneous HGDP CSA
and SGDP SAS outside the imputation panels.

If SAS turns out thin after all, the reuse is what the intake's
`7.sim_kinship` check detects. See §5.

## 3. The real genomes — test samples, not donors

Real individuals from HGDP, SGDP, MXB and Almarri join the cohort as test
samples. They are **outside every imputation server panel**, so no server has
seen them.

| sample_type | who | what it adds |
|---|---|---|
| `real_homogeneous` | top Q ≥ threshold for the **argmax** super-pop | real rare variation and real haplotype structure that simulation cannot produce |
| `real_admixed` | top Q below that threshold | real admixture, with the correlation structure between ancestry and genotype that a simulation idealises away |
| `real_unassigned` | argmax is **OTH**, the leftover cluster (only when 4a runs above K=7; none at K=7) | the multiracial-like individuals a real diverse cohort carries. No reference ancestry fits them, so they are neither homogeneous nor admixed between named components |

The threshold is **0.94 for every super-pop** (17 Nov 2026), the same number 4a
applies to anchors and the homogeneous panel, and the contract's flat threshold
(`run_benchmark_cohort.sh FLAT_THRESHOLD`, `check-simset-contract.py`, and on
the Nextflow side `plan_simset.R --threshold` and intake 2.real_samples). Until
then the panel was tiered (0.95 EUR/EAS/AFR, 0.93 AMR/OCE/MEN/SAS) while the
contract was a flat 0.93, so the two could disagree for samples between 0.93
and 0.95; one number for both removes that.

**The threshold is keyed on the ADMIXTURE argmax, not on the metadata label,
and homogeneity does not depend on the two agreeing.** "Is this genome
homogeneous?" is a question about the ancestry the supervised Q actually
found, so a sample at Q = 0.97 on EUR is measured against EUR's 0.95 whatever
its label says. Agreement is recorded separately in `label_concordant`.

This matters because the argmax and the metadata super-pop are two
independently built vocabularies, and three different lookups are keyed on a
super-pop name — the tiered threshold, the reference-panel size, and the
comparison itself. Every one of them **fails silently** on a label this
vocabulary does not know: the threshold falls back to the relaxed 0.93, the
reference size reads 0 so support becomes `none`, and the comparison is never
equal so the sample is called admixed. A cohort whose metadata gives a
geographic region rather than a super-pop would be misclassified wholesale
with nothing in the output to show it.

4a does harmonise on the way in — `super_pop_of()` maps MXB → AMR, SGDP via
its own table, CSA → SAS and MID → MEN, and aborts if its two vocabularies do
not overlap — so MXB samples carry `AMR`, not a Mexican geographic label. The
assembler **re-checks rather than trusting it**, prints both censuses, and
refuses to run on any non-canonical label, because a newly added cohort is
exactly what would introduce one. It also reports per-label concordance and
flags any label that disagrees with the argmax for most of its samples, since
that is a harmonisation signal rather than a biological finding.

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
2. **Power and GWAS simulations** — phenotypes simulated on real genotypes,
   **alongside the simulated individuals.**

**Both halves are used for power, because the real half alone is too small.**
A few hundred genomes cannot carry a power calculation; ~500 real plus 2,160
simulated is about 2,660, which can.

The cost is that simulated individuals drawn from a finite donor pool are
cryptic relatives *by construction* — two individuals built from the same
donor haplotypes share long stretches, and the thinner the pool the worse it
is (SAS was projected as the binding case; see §2 for why that projection is stale). That
does **not** have to bias anything: the GWAS uses GENESIS with a GRM, which
models relatedness rather than pruning it, and power simulations run through
the same machinery inherit that. What relatedness costs is **effective**
sample size, and that is a number rather than a caveat:

$$N_{\text{eff}} = \frac{N^2}{\sum_{ij} R_{ij}}$$

with $R$ the GRM. `6e_cohort_kinship_clm.sh` computes genome-wide kinship over
the combined cohort and `cohort-power-summary.py` turns it into $N_{\text{eff}}$
per group. **Power calculations must use $N_{\text{eff}}$, not the roster
count** — on a synthetic thin-donor fixture, 60 SASEUR individuals came out
behaving like 5.

That also makes §5's SAS options decidable: accepting the reuse is defensible
*provided* the reuse is reported and the power calculation uses
$N_{\text{eff}}$. It is now reportable.

**How related, at the 18 Nov donor counts (20 Nov).** Reuse is not kinship:
simu-mix gives every breakpoint a fresh donor haplotype and never lets two
haplotypes of one run copy the same donor at a marker, so two simulated
individuals from different runs have expected kinship
$\sum_P p_P^2 / (2 \times \text{donors}_P)$ and two from one run share
nothing. A Monte Carlo of the two programs' code found **no MENEUR or EASEUR
pair at 3rd degree** (closest 0.042 and 0.029; KING's 3rd-degree bin starts at
0.0442), so the two-generation models stay. The gen-0 controls are the
related ones, because each chromosome comes whole from a small pool: EUR1
(32 donors) put 59 of 1,770 pairs at 3rd degree, MEN1 (39) 33, none at 2nd.
The GRM models that; it is reported, not pruned.

**Batch size (20 Nov).** simu-mix dies when every donor haplotype of a
component is in use at a marker. 6b used to size each run so the AVERAGE use
was 80% of the pool; the use at a marker is binomial and varies along the
chromosome, so it crossed 100% on most chromosomes (the same Monte Carlo: a
chromosome task died ~90% of the time for LAT3, ~78% CAR3, ~67% MENA3, ~15%
MENEUR/EASEUR). Each run is now sized so the mean sits 5 SD below the pool
(`FEAS_Z=5`): LAT3 44 per run, CAR3/MENEUR/EASEUR 40, MENA3 51, AFAM2 85 —
0 deaths in 4,400 chromosome runs per model.

The real half still contributes what simulation cannot — real LD structure,
real allele-frequency spectrum, real rare variation — so it is not
interchangeable with more simulated individuals even though both are used.

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

**Which design, and why builds are paired by default.** The two readings are
not symmetric:

- The **build effect** is "did the pipeline's liftover change *this person's*
  answer?" That is a within-individual question. Comparing hg19 in one group
  of people against hg38 in a different group makes the liftover effect
  compete with between-person variance, which is much larger. So every
  individual goes on **both** builds.
- The **chip effect** is marker density — a property of the array, not an
  interaction with the person — so an ancestry-balanced between-individual
  comparison is adequate, and splitting individuals across chips is what keeps
  the cost at 2× rather than 6×.

Imputation runs on TOPMed / All of Us / Michigan are the binding constraint,
not disk, which is what makes that asymmetry worth exploiting.

| `--design` | imputation runs | build effect | chip effect | sim per cell |
|---|---|---|---|---|
| `one-cell` | 1× ≈ 2,660 | between-individual | between-individual | 360 |
| **`paired-builds`** (default) | 2× ≈ 5,320 | **within-individual** | between-individual | **720** |
| `full` | 6× ≈ 15,960 | within-individual | **within-individual** | 2,160 |

- Under the default every individual is on **one chip and both builds**. The
  manifests carry `chips` and `builds` as comma-separated lists, and 6d tests
  membership rather than equality — an equality test against `"hg19,hg38"`
  matches nothing and would yield empty cells rather than an error.
- Every model, and every real `sample_type` × top-component group, is split
  **equally across the chips** (or cells, under `one-cell`), within 1.
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

**SAS.** **Decided 18 Nov: no SAS simulation; SAS is read from real
genomes.** The K=8 counts settled it: about 10 homogeneous SAS genomes outside
1000 Genomes (all SGDP South Indian), and loosening the cut does not help --
11 at 0.90, because the HGDP Pakistani populations sit far below it (most
under 0.70 SAS), not just under the line. 10 donors are 20 haplotypes for
260 SASEUR + SAS1 individuals. So SASEUR and SAS1 are out of
`reference_ids/models.tsv` (kept there as comments, with how to restore them),
and SAS is reference-only like OCE (`NO_DONOR_SUPERPOPS=OCE,SAS`): its
reference is 1000 Genomes only, so those ~10 genomes are real_homogeneous test
genomes beside the many real_admixed SAS ones. Real genomes are never scored
for local ancestry, so **SAS has no LAI truth set** until more homogeneous SAS
genomes arrive (DUA-gated data are being pursued); its imputation reading is
unaffected.

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
handful of genomes. (K=7, 18 Nov: 30 homogeneous OCE, about 15 / 15. The
reference takes whole relative clusters since 20 Nov, so it can land a
genome or two under 15.)

**Real homogeneous test genomes (18 Nov).** Donor-first left the real half
~96% admixed (~27 homogeneous: OCE 15, SAS 10, EAS 1, EUR 1). The split now
takes a few homogeneous genomes back from the donors as `real_homogeneous`
test genomes — AFR 15, EAS 15, EUR 10 (`--test-holdout`; EUR at the low end
of 10–15 because its 42 donors are the thinnest pool most models use). The
13 populations with 3 or fewer donors stay donors: a reference that does
not cover every lineage is what real data look like, and their gaps are
reported for the methods.

**One pool, two projects (18 Nov).** The benchmark splits 4a's `homogeneous`
flag, not `included`: `included` is the LAI project's panel and leaves out
THAT project's simulation donors and their relatives (56 Native American
genomes plus relatives at K=8), which is no reason to withhold them here.

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
| 10 | write the chip files per chip × build × chr | `6d_write_chip_cells_clm.sh` + `lift-cell-to-hg19.py` | **built**; orientation logic unit-tested, bcftools pipeline untried |
| 10b | genome-wide kinship over the combined cohort | `6e_cohort_kinship_clm.sh` + `cohort-power-summary.py` | **built** |
| 11 | tables, checksums, provenance | — | **not built** |
| 12 | `ingest_simset.sh` intake check | — | — |
| 13 | hand over | — | — |

Steps 2, 5, 7 and 8 need only the Q file, so all four run within minutes of
stage 3 landing. Step 11 is the remaining gap. Steps 6 and 10 run on taco
(mhgcp) alongside the rest, not in `genotyping_pipeline`.

**Step 10's hg19 path is where a silent wrong answer would come from.** GT is
an index into REF/ALT, so a genotype means nothing apart from the allele pair
written beside it. Between GRCh37 and GRCh38 a site can keep its alleles, have
them complemented (a strand flip — indices unchanged), or have REF and ALT
*exchanged* (the reference allele changed between builds — indices must be
flipped). Missing that last case writes a perfectly valid VCF with inverted
alternate-allele dosage at exactly the sites where the builds disagree, and an
imputation benchmark run on it would report inflated error that reads as a
pipeline problem. `lift-cell-to-hg19.py` classifies all four cases plus the
unresolvable one and is unit-tested on each.

Palindromic sites (A/T, C/G) are irreducibly ambiguous: complement equals swap
there, so orientation comes from the vendor strand column alone and is counted
separately in the per-chromosome report. GSA has 6,422 of them, 1.08%.

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
