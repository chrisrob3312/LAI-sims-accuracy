# LAI-sims-accuracy
*Code for determining local ancestry inference accuracy*

simulation.sh - general code for simulating admixed individuals with admix-simu https://github.com/williamslab/admix-simu

phasing-jointcall.sh - Filtering sites and phasing of the 1kg-hgdp joint call dataset (legacy SHAPEIT4 pipeline).

envs/shapeit5.yml - conda environment for the merge pipeline: SHAPEIT5, bcftools/htslib, Picard, plink2.
Create with `conda env create -f envs/shapeit5.yml`, activate with `conda activate shapeit5`.

1a_prep-mxb-liftover_clm.sh - one-shot SLURM job: liftover MXB hg19 -> hg38 (Picard with
`RECOVER_SWAPPED_REF_ALT=true`), normalize, fix REF/ALT, split per chromosome.

1c_phasing-jointcall_clm.sh - SLURM array (1-22) merging HGDP+1KG (relatives kept, `OUTLIERS=none`;
earlier panels used the post-outlier set) with the lifted MXB, SGDP and Almarri,
applying a soft-union per-superpop MAF filter (>=0.005 in any of AFR/AMR/EUR/EAS/SAS/CSA/OCE/MEN),
plink2 site QC, and SHAPEIT5_phase_common joint phasing. Imputation will be done with
TOPMed/All-of-Us so rare-variant phasing is OFF by default; set `RUN_PHASE_RARE=1` at submit
time to also run `SHAPEIT5_phase_rare` and emit a full common+rare phased panel.

2b_simulation_clm.sh - modernized replacement for simulation.sh. Builds admix-simu donor
.phgeno files for both NAT donor configurations (HGDP-NAT-only and HGDP-NAT+MXB) directly
from the merged phased panel via plink2 (no SHAPEIT2 dependency), runs admix-simu, and
emits simulated admixed haps + truth files. SLURM array 1-22.

3b_wgs-rfmix-jointcall_clm.sh - modernized replacement for wgs-simulation-rfmix-jointcall.sh.
Builds RFMix v1 inputs and runs RFMix across the 4 reference panels x 2 sim-tracks grid:
panels {NAT_HGDP, NAT_PEL, NAT_PEL_EAS, NAT_HGDPMXB} x tracks {NAT, NATMXB}. SLURM array 1-22.

make_sample_groups.sh - emits the 2-column sample->superpop TSV that drives the per-superpop MAF
filter (consumed by `bcftools +fill-tags -S`).

build-panel-keep-files.sh - composes RFMix-reference and SIMU-donor keep-files for the 4-panel
comparison (Panel 1 = HGDP-NAT + IBS + YRI; Panel 4 = HGDP-NAT + MXB + IBS + YRI; SIMU donor
tracks for HGDP-NAT-only and HGDP-NAT+MXB) from the lists in `reference_ids/`.

build-pel-panels.sh - derives `pel_rfmix.txt` and `pel_eas_rfmix.txt` from
`gnomad_meta_updated.tsv` so panels 2 (PEL) and 3 (PEL+EAS) run alongside 1 and 4.
Pass the metadata path: `./build-pel-panels.sh /path/to/gnomad_meta_updated.tsv`.

All outputs are anchored under `$PROJECT_ROOT`, default
`/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/`.
Override at submit time if running elsewhere: `PROJECT_ROOT=/path/to/proj sbatch ...`.

Project layout:
```
$PROJECT_ROOT/
  01_merged_phased_panel/   1a + 1b outputs (lifted MXB BCFs, merged + phased panel)
  02_simulations/           2b outputs (donor .phgeno, admix-simu, sim haps); place
                              ${ADMIX_POP}.dat and ${ADMIX_POP}.sample.txt here
  03_rfmix/                 3b outputs (RFMix runs, Lat3 accuracy-prep)
  04_accuracy/              accuracy.R outputs
  panel_keep_files/         build-panel-keep-files.sh output
  logs/                     SLURM stdout/stderr
  sample_groups.tsv         make_sample_groups.sh output
```

Run order:
1. `conda env create -f envs/shapeit5.yml && conda activate shapeit5`
2. `./init-project-tree.sh` (creates the dirs above + builds sample_groups.tsv,
   pel_rfmix.txt, pel_eas_rfmix.txt, and panel keep-files in one shot)
3. `sbatch 1a_prep-mxb-liftover_clm.sh` (one-shot)
4. `sbatch 1c_phasing-jointcall_clm.sh` (array 1-22)
5. Place `${ADMIX_POP}.dat` and `${ADMIX_POP}.sample.txt` in `02_simulations/`
6. `ADMIX_POP=Brasa GEN=12 sbatch 2b_simulation_clm.sh` (array 1-22; per (admix-pop, gen))
7. `ADMIX_POP=Brasa GEN=12 sbatch 3b_wgs-rfmix-jointcall_clm.sh` (array 1-22)
8. accuracy.R / accuracy_error_plots.R unchanged (consume RFMix Lat3 + truth files)

wgs-simulation-rfmix-jointcall.sh - code used for generating all simulated models and RFMix v1 runs, and preparing files for the accuracy calculation.

accuracy.R - code for calculating true positive rates of RFMix calls of simulations, getting counts of miscalls per error mode between ancestry groups, and getting the positions with highest number of miscalls.

accuracy_error_plots.R - code used to generate the accuracy and error modes manuscript figures.

viterbi2msp.R - code to convert RFMix v1 .viterbi output to RFMix v2 .msp output file formats for use with Tractor.

*Which sample sets serve which project (17 Nov 2026).* The projects need
different sets, so the lists are kept apart:

| project | references | donors / test genomes | rule set |
|---|---|---|---|
| 3-way NAT/EUR/AFR LAI accuracy simulation (rerun of Jessica Honorato-Mauer's design; RFMix + FLARE, WGS + GSA chip, 22 chromosomes; two cohorts, below) | top-level `*_rfmix.txt`: 1KG IBS (EUR), 1KG YRI (AFR), the NAT lists, and `amr_homog_rfmix.txt` (July's K=5 AMR panel + AMR samples new since July that pass a K=5 run at Q >= 0.94: 4a `HOMOG_SET=lai3way`, July's anchors, MEN/OCE left out) | `*_simulation.txt` (1KG IBS, 1KG YRI, HGDP NAT, MXB) | references and donors disjoint, never first-degree relatives; 1KG allowed. Originals in `rerun_jessicaH_3way_lai_sim/` |
| genotyping-harmonization / imputation benchmark (Nextflow) and the global panel | the K=7 panel (4a `HOMOG_SET=global`: `lai_ref_panel_samples.k7.tsv`, every sample, Q >= 0.94, central high-Q anchors), split by `run_benchmark_cohort.sh` | donors and real test genomes from HGDP, Almarri, SGDP and MXB; never 1KG (1KG is in the imputation panels) | contract threshold 0.94 (`plan_simset.R` must match). The 6 Jul anchor lists (July's K=5 anchors) are in `genotype_harmonization_benchmark_testing_anchor_list_20260706/` |

*The 3-way simulation's two cohorts (17 Nov).* As the August pilot paired
them -- but August ran both off `Brasa.dat`, so its "Mexican-like" cohort had
Brazilian proportions and 12 generations:

| cohort | `ADMIX_POP` | model (`.dat`) | gen | NAT donors (track) | results |
|---|---|---|---|---|---|
| Brazilian-like | `Brasa` | NAT .15 / EUR .60 / AFR .25 (Jessica's, `02_simulations/Brasa.dat` on the cluster) | 12 | HGDP NAT (`NAT`) | `results/accuracy_pilot*/` |
| Mexican-like | `MexLike` | NAT .50 / EUR .45 / AFR .05 (`simulation_models/MexLike.dat`) | 16 | HGDP NAT + 25 MXB (`NATMXB`) | `results/accuracy_pilot*_MexLike/` |

Both: 60 haplotypes (30 individuals), panels NAT_HGDP, +PEL, +PEL+EAS,
+25 MXB, +50 MXB (Brazilian-like only: the Mexican-like cohort's MXB donors
would be in the reference) and the AMR homogeneous panel. The driver's
`COHORTS` (default `Brasa:12:NAT MexLike:16:NATMXB`) sets them; 2b checks each
`.dat` against its `GEN` and writes the `.notref` 3b reads.

*IDs of reference samples used to generate admixed haplotypes and RFMix Reference:*

amr_simulation.txt - HGDP AMR population IDs used for simulation 

eur_simulation.txt - 1KG IBS population IDs used for simulation

afr_simulation.txt - 1KG YRI population IDs used for simulation


amr_rfmix.txt - HGDP AMR population IDs used as RFMix reference

eur_rfmix.txt - 1KG IBS IDs (74) used as the EUR RFMix / FLARE reference in every
3-way arm. Restored on 17 Nov 2026 to the design (Jessica's 77-ID IBS list in
rerun_jessicaH_3way_lai_sim/, as July's panel used it: the other 3 were relatives July's
panel did not contain). From 6 Jul to 17 Nov it held HGDP French, Orcadian,
Sardinian and Russian anchors (genotype_harmonization_benchmark_testing_anchor_list_20260706/).

afr_rfmix.txt - 1KG YRI IDs (72) used as the AFR RFMix / FLARE reference in every
3-way arm. Restored on 17 Nov 2026 the same way (77-ID YRI list; 5 relatives not
in July's panel). From 6 Jul to 17 Nov it held HGDP Yoruba, Mandenka and Bantu
anchors (genotype_harmonization_benchmark_testing_anchor_list_20260706/).

sas_rfmix.txt - HGDP Central/South Asian IDs (Balochi, Brahui, Makrani, Sindhi).
HGDP00029 (Brahui) and HGDP00150 (Makrani), Koenig removals, were dropped on
14 Nov 2026 (48 IDs). No simulation or RFMix arm reads this list.

