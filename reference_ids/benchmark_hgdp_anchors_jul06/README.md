# benchmark_hgdp_anchors_jul06 — HGDP (non-1KG) EUR and AFR anchor lists

The `eur_rfmix.txt` and `afr_rfmix.txt` that sat at the top level of
`reference_ids/` from commit 046baad (6 Jul 2026, "anchor rebalance": Koenig
et al. metadata, cap 50) until 17 Nov 2026. This folder was
`backup_koenig_jul06/` for a few hours on 17 Nov.

| file | IDs |
|---|---|
| `eur_rfmix.txt` | 50 HGDP: French 23, Orcadian 12, Sardinian 8, Russian 6, +1 |
| `afr_rfmix.txt` | 50 HGDP: Yoruba 21 (incl. HGDP00933/HGDP00934, swapped in on 14 Nov for two Koenig QC removals), Mandenka 20, BantuKenya 6, BantuSouthAfrica 1, +2 |

Why they exist: the genotyping-harmonization / imputation benchmark cannot use
1000 Genomes samples as simulation donors or as real test genomes, because
1KG is inside the imputation reference panels — so it needs non-1KG
(HGDP, SGDP, MXB, Almarri) material for those roles. These are HGDP lists of
that kind.

Not read by any pipeline stage now. The 3-way LAI simulation's EUR and AFR
references are the 1KG IBS and YRI lists again (see
`../rerun_jessicaH_3way_lai_sim/README.md`); because every 3b arm reads the
top-level `eur_rfmix.txt` / `afr_rfmix.txt`, the August 2026 accuracy pilot
most likely ran with these HGDP lists as its EUR and AFR references. The
benchmark's own donor / reference split comes from
`scripts/split-homog-donors-reference.py` on the K=7 panel table (via
`run_benchmark_cohort.sh`), not from these files.

`eas_rfmix.txt` and `sas_rfmix.txt` at the top level are the same 6 Jul
Koenig lists for EAS and SAS; no 3-way arm reads them.
