# benchmark_truth — simulated truth set for the genotyping-pipeline benchmark

Purpose: a **non-circular** truth set for benchmarking harmonization + imputation
(TOPMed, All of Us AnVIL, Michigan 1000G, Michigan HRC) and downstream simulated-GWAS
power by ancestry. Every genome here is spliced from **HGDP + MX Biobank** haplotypes,
which sit in none of those reference panels, so no arm can copy its own truth.
1KG samples are never donors (they are inside the 1000G panel and inside HRC).

Design summary (edit `cohorts.tsv`, everything else reads it):

| cohort | type | n | sources (props) | gen | build | platforms |
|---|---|---|---|---|---|---|
| HOMOG_{AFR,EUR,EAS,CSA,MEN} | homogeneous | 100 | one superpop, spliced *within* population, pooled | 10 | 38 / 37 mixed | SNP6 + GSA |
| HOMOG_{AMR,OCE} | homogeneous | 50 | capped by donor count (64 / 28 donors) | 10 | 38 | GSA (+SNP6 for AMR) |
| ADM_AFR_EUR | admixed | 100 | AFR 0.80 / EUR 0.20 | 8 | 38 | SNP6 + GSA |
| ADM_MEX, ADM_CAR | admixed | 100 | EUR/AMR/AFR two profiles | 12 | 38 / 37 | SNP6 + GSA |
| ADM_EUR_CSA, ADM_EAS_EUR | admixed | 100 | 0.5 / 0.5 | 10 | 38 | as listed |
| ADM_4WAY, ADM_5WAY | admixed | 100 | 4- and 5-way | 12 / 15 | 38 / 37 | SNP6 + GSA |

`build` decides which cohorts are *delivered* on GRCh37 (arrays lifted hg38→hg19 in 5e);
the truth always stays on GRCh38, so the pipeline's harmonization/liftover is exercised
and scored against one reference.

## Stages (all under `$PROJECT_ROOT/06_benchmark_truth/`)

| stage | script | in → out | run as |
|---|---|---|---|
| 5a | `5a_build_donor_lists.sh` | `sample_groups.tsv`, gnomAD meta, MXB popinfo, (optional) `lai_ref_panel_samples.tsv` → `donors/donors.tsv`, per-population keep-lists, `donor_counts.tsv` | one-shot, seconds |
| 5b | `5b_simulate_cohorts_clm.sh` | phased panel + donors + `cohorts.tsv` → `sim/<COHORT>/*.haps/.sample` + local-ancestry truth `.bp/.hanc2` | sbatch array 1-22 |
| 5c | `5c_truth_vcfs_clm.sh` | haps → `truth_hg38/<COHORT>.chr*.truth.vcf.gz`, `ALL.chr*.truth.vcf.gz` | sbatch array 1-22 |
| 5d | `5d_insilico_arrays_clm.sh` (+ `insilico_array.py`) | truth → `arrays_hg38/{snp6,gsa}/<COHORT>.*.vcf.gz/.bed` + `.changes.tsv` | sbatch array 1-22 |
| 5e | `5e_liftover_arrays_hg19_clm.sh` | build-37 cohorts' arrays → `arrays_hg19/…` | sbatch array 1-22 |
| 5f | `5f_qc_truth_set_clm.sh` (+ `qc_truth_set.R`) | counts, PCA by cohort, `qc/MANIFEST.tsv` | sbatch one-shot |

Prerequisites already produced by the LAI pipeline: `1b` merged phased panel
(`01_merged_phased_panel/merged_chr*.…rechr.bcf`), `sample_groups.tsv`, `resources/gmap_hg38/`,
and (recommended) `4a` → `lai_ref_panel_samples.tsv` so homogeneous donors are Q-filtered.
SNP6 site lists: build once with `scripts/prepare-chip-snps.sh --name snp6` from the
Affymetrix GenomeWideSNP_6 hg38 annotation (see the comment block in 5d).

## What the in-silico arrays add (and log)

`insilico_array.py` restricts to the platform's sites and then injects, per platform
(rates set at the top of 5d, defaults below):

| platform | minus-strand reporting (non-palindromic) | palindromic A/T, C/G | genotype error | missingness |
|---|---|---|---|---|
| GSA | 50 % of sites | random strand | 0.2 % | 1 % |
| SNP6 | 35 % of sites | random strand | 0.4 % | 2 % |

Every modification is written to `<file>.changes.tsv` and tagged in `INFO/INSILICO`, so the
harmonization step can be **scored** (flips recovered, palindromes resolved, errors passed
through) rather than eyeballed. Genotypes are unphased in the array files.

## Run order (read Block A before each Block B)

Block A — verify inputs exist:
```bash
ls $PROJECT_ROOT/01_merged_phased_panel/merged_chr22.shapeit5_phased.softunion_maf005.rechr.bcf
ls $PROJECT_ROOT/sample_groups.tsv $PROJECT_ROOT/04_homogeneity_panel/lai_ref_panel_samples.tsv
ls reference_ids/chip_gsa_jessica/chr22.snps reference_ids/chip_snp6/chr22.snps
```
Expect: all five paths listed. If `lai_ref_panel_samples.tsv` is missing, 5a still runs
(no Q filter) and says so.

Block B — only if Block A listed every file:
```bash
bash benchmark_truth/5a_build_donor_lists.sh          # then READ donors/donor_counts.tsv
```
Check `donor_counts.tsv`: superpop totals must be ≥ the `n` in `cohorts.tsv` for every
homogeneous cohort (OCE/AMR are capped at 50 for this reason). Adjust `cohorts.tsv` if not.

Block C — only after the counts check:
```bash
sbatch benchmark_truth/5b_simulate_cohorts_clm.sh
# when all 22 finish:
sbatch benchmark_truth/5c_truth_vcfs_clm.sh
sbatch benchmark_truth/5d_insilico_arrays_clm.sh
sbatch benchmark_truth/5e_liftover_arrays_hg19_clm.sh
sbatch benchmark_truth/5f_qc_truth_set_clm.sh
```

## Hand-off to the genotyping pipeline

Give the benchmark the **array** files (`arrays_hg38/`, `arrays_hg19/`) as inputs and the
**truth** files (`truth_hg38/ALL.chr*.truth.vcf.gz`) as the comparison set; `qc/MANIFEST.tsv`
lists everything with md5s; `truth_hg38/ALL.chr*.samples.tsv` maps sample → cohort for
per-ancestry metrics; `sim/<COHORT>/*.hanc2` is the local-ancestry truth for
power-by-haplotype-ancestry. The design table (`cohorts.tsv`) is the methods table.

## Notes / caveats to state in the methods

- Spliced individuals from a small donor pool (OCE 28, AMR 64 + 50 MXB) share long IBD
  tracts; fine for imputation accuracy, but do not treat them as independent for anything
  that assumes unrelated samples without LD-aware relatedness handling.
- Within-population splicing avoids manufacturing within-superpop admixture; superpop
  cohorts are pools of population-level splices (sizes ∝ donor counts).
- `admix-simu` has no seed; `SEED` in 5b is recorded, not applied. Re-runs produce
  different genomes — keep the delivered set under version control via `qc/MANIFEST.tsv`.
- Donors that are also RFMix references (`reference_ids/*_rfmix.txt`) are excluded by
  default (`EXCLUDE_RFMIX_REFS=1`) so the pipeline's LAI module cannot copy its own donors;
  set to 0 to widen the pool and say so.
- The `.dat` header written by `make_dat.py` (`<n> ADMIX <POP…>` / `<gen> <props…>`) follows
  admix-simu's layout; compare with an existing `02_simulations/*.dat` on first use.
