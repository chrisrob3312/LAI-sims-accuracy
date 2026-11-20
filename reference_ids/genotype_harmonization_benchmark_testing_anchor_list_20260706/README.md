# genotype_harmonization_benchmark_testing_anchor_list_20260706

The 6 Jul 2026 anchor lists (commit 046baad, "anchor rebalance": Koenig et al.
metadata, cap 50): one list per super-pop, none of them 1000 Genomes.

| file | IDs | who |
|---|---|---|
| `eur_rfmix.txt` | 50 | HGDP French 23, Orcadian 12, Sardinian 8, Russian 6, +1 |
| `afr_rfmix.txt` | 50 | HGDP Yoruba 21, Mandenka 20, BantuKenya 6, BantuSouthAfrica 1, +2 (as of 14 Nov: HGDP00933/HGDP00934 replaced two Koenig QC removals, HGDP00621 and HGDP01270) |
| `eas_rfmix.txt` | 50 | HGDP Japanese, Han, Yakut |
| `sas_rfmix.txt` | 48 | HGDP Sindhi, Makrani, Balochi, Brahui (as of 14 Nov: Koenig removals HGDP00029 and HGDP00150 dropped) |
| `amr_rfmix.txt` | 50 | MXB (the 6 Jul MXB-only AMR list; the top-level `amr_rfmix.txt` went back to 31 HGDP NAT on 1 Aug) |

Until 17 Nov 2026 the EUR and AFR lists sat at the top level of
`reference_ids/`, where every 3-way RFMix / FLARE arm read them; they were
then restored to 1KG IBS and YRI (see `../rerun_jessicaH_3way_lai_sim/`).
This folder was `backup_koenig_jul06/`, then `benchmark_hgdp_anchors_jul06/`,
on 17 Nov; it is named for the benchmark because its pool is not HGDP alone:
the genotype-harmonization benchmark draws donors and real test genomes from
HGDP, Almarri, SGDP and MXB, never from 1000 Genomes (1KG is inside the
imputation servers' panels).

## What reads them

- **July's K=5 run used them as its supervised anchors.** The `assigned`
  column of `../lai_ref_panel_samples.tsv` is these lists less the samples
  July's panel did not contain: EUR 49, AFR 46, EAS 50, SAS 48, AMR 50.
- **4a `HOMOG_SET=lai3way`** (K=5, the 3-way LAI simulation's AMR panel) does
  NOT reuse them: it finds its anchors from its own unsupervised fit on the
  expanded panel (central, Q >= 0.94), as the K=7 run does, and reports how
  many of these by-name anchors it picks again. `ANCHOR_SOURCE=table`
  reproduces July's by-name anchors for a comparison run.
- **The benchmark's own anchors are not these.** 4a `HOMOG_SET=global` (K=7,
  every sample) discovers central, high-Q anchors from the data and writes
  them to `../<sp>_anchors.k7.txt`; the benchmark's donor / reference split
  (`run_benchmark_cohort.sh`) reads the K=7 table.

No pipeline stage reads the files in this folder directly. `2b` checks
donors against the top-level `*_rfmix.txt` only.
