# LAI-sims-accuracy
*Code for determining local ancestry inference accuracy*

simulation.sh - general code for simulating admixed individuals with admix-simu https://github.com/williamslab/admix-simu

phasing-jointcall.sh - Filtering sites and phasing of the 1kg-hgdp joint call dataset (legacy SHAPEIT4 pipeline).

envs/shapeit5.yml - conda environment for the merge pipeline: SHAPEIT5, bcftools/htslib, Picard, plink2.
Create with `conda env create -f envs/shapeit5.yml`, activate with `conda activate shapeit5`.

1a_prep-mxb-liftover_clm.sh - one-shot SLURM job: liftover MXB hg19 -> hg38 (Picard with
`RECOVER_SWAPPED_REF_ALT=true`), normalize, fix REF/ALT, split per chromosome.

1b_phasing-jointcall_clm.sh - SLURM array (1-22) merging HGDP+1KG postoutlier with the lifted MXB,
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

Run order:
1. `conda env create -f envs/shapeit5.yml && conda activate shapeit5`
2. `./make_sample_groups.sh <gnomad_meta_updated.tsv> reference_ids/MXB50genomes_popinfo.tsv > sample_groups.tsv`
3. `sbatch 1a_prep-mxb-liftover_clm.sh` (one-shot)
4. `sbatch 1b_phasing-jointcall_clm.sh` (array 1-22)
5. `./build-panel-keep-files.sh` (one-shot, optional)
6. `ADMIX_POP=Brasa GEN=12 sbatch 2b_simulation_clm.sh` (array 1-22; per (admix-pop, gen))
7. `ADMIX_POP=Brasa GEN=12 sbatch 3b_wgs-rfmix-jointcall_clm.sh` (array 1-22)
8. accuracy.R / accuracy_error_plots.R unchanged (consume RFMix Lat3 + truth files)

wgs-simulation-rfmix-jointcall.sh - code used for generating all simulated models and RFMix v1 runs, and preparing files for the accuracy calculation.

accuracy.R - code for calculating true positive rates of RFMix calls of simulations, getting counts of miscalls per error mode between ancestry groups, and getting the positions with highest number of miscalls.

accuracy_error_plots.R - code used to generate the accuracy and error modes manuscript figures.

viterbi2msp.R - code to convert RFMix v1 .viterbi output to RFMix v2 .msp output file formats for use with Tractor.

*IDs of reference samples used to generate admixed haplotypes and RFMix Reference:*

amr_simulation.txt - HGDP AMR population IDs used for simulation 

eur_simulation.txt - 1KG IBS population IDs used for simulation

afr_simulation.txt - 1KG YRI population IDs used for simulation


amr_rfmix.txt - HGDP AMR population IDs used as RFMix reference

eur_rfmix.txt - 1KG IBS population IDs used as RFMix reference

afr_rfmix.txt - 1KG YRI population IDs used as RFMix reference

