# LAI-sims-accuracy
*Code for determining local ancestry inference accuracy*

simulation.sh - general code for simulating admixed individuals with admix-simu https://github.com/williamslab/admix-simu

phasing-jointcall.sh - Filtering sites and phasing of the 1kg-hgdp joint call dataset (legacy SHAPEIT4 pipeline).

envs/shapeit5.yml - conda environment for the merge pipeline: SHAPEIT5, bcftools/htslib, Picard, plink2.
Create with `conda env create -f envs/shapeit5.yml`, activate with `conda activate shapeit5`.

prep-mxb-liftover.sh - one-shot SLURM job: liftover MXB hg19 -> hg38 (Picard with
`RECOVER_SWAPPED_REF_ALT=true`), normalize, fix REF/ALT, split per chromosome.

merge-mxb-hgdp1kg.sh - SLURM array (1-22) merging HGDP+1KG postoutlier with the lifted MXB,
applying a soft-union per-superpop MAF filter (>=0.005 in any of AFR/AMR/EUR/EAS/SAS/CSA/OCE/MEN),
plink2 site QC, and SHAPEIT5 joint phasing. Run after `prep-mxb-liftover.sh` and after building
`sample_groups.tsv` via `make_sample_groups.sh`.

make_sample_groups.sh - emits the 2-column sample->superpop TSV that drives the per-superpop MAF
filter (consumed by `bcftools +fill-tags -S`).

build-panel-keep-files.sh - composes RFMix-reference and SIMU-donor keep-files for the 4-panel
comparison (Panel 1 = HGDP-NAT + IBS + YRI; Panel 4 = HGDP-NAT + MXB + IBS + YRI; SIMU donor
tracks for HGDP-NAT-only and HGDP-NAT+MXB) from the lists in `reference_ids/`.

Run order:
1. `conda env create -f envs/shapeit5.yml && conda activate shapeit5`
2. `./make_sample_groups.sh <gnomad_meta_updated.tsv> reference_ids/MXB50genomes_popinfo.tsv > sample_groups.tsv`
3. `sbatch prep-mxb-liftover.sh` (one-shot)
4. `sbatch merge-mxb-hgdp1kg.sh` (array 1-22)
5. `./build-panel-keep-files.sh` (one-shot)
6. Existing `wgs-simulation-rfmix-jointcall.sh` updated to use the new merged phased
   panel as VCF input and the panel/simu keep-files from step 5.

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

