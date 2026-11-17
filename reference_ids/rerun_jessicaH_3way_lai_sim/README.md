# rerun_jessicaH_3way_lai_sim — the 3-way (NAT / EUR / AFR) LAI simulation lists

Jessica Honorato-Mauer's original lists for the 3-way LAI accuracy simulation,
as they stood before commit 046baad (6 Jul 2026) replaced several of the
top-level copies with Koenig-based ADMIXTURE anchor lists. This folder was
`backup_pre_koenig/` until 17 Nov 2026.

The rerun of the 3-way simulation (2b → 3b RFMix and FLARE, WGS and GSA chip)
reads the **top-level** `reference_ids/*_rfmix.txt` and `*_simulation.txt`.
Those relate to the files here as follows:

| file here | IDs | top-level list the rerun reads |
|---|---|---|
| `eur_rfmix.txt` | 77 1KG IBS | `eur_rfmix.txt`: the 74 of them July's panel contained (the other 3 were relatives July's panel dropped) |
| `afr_rfmix.txt` | 77 1KG YRI | `afr_rfmix.txt`: the 72 July's panel contained (5 relatives dropped) |
| `amr_rfmix.txt` | 31 HGDP NAT + legacy MXB | identical (restored on 1 Aug, commit 2a1eb2c) |
| `amr_hgdpmxb_rfmix.txt` | 56 | identical |
| `mxb_rfmix.txt` | 25 MXB | identical |
| `pel_rfmix.txt` | 86 (84 PEL) | identical |
| `pel_eas_rfmix.txt` | 587 PEL + 1KG EAS | identical |
| `eas_rfmix.txt` | 501 1KG EAS | not used by the 3-way arms (top level holds the 6 Jul Koenig list) |
| `sas_rfmix.txt` | 507 1KG SAS | not used by the 3-way arms (top level holds the 6 Jul Koenig list) |

The NAT_HOMOG arms add `amr_homog_rfmix.txt`, which 4a writes: July's K=5
AMR homogeneous samples plus the AMR samples new since July that pass in the
current run, minus simulation donors and their relatives
(`scripts/build-3way-amr-homog.py`).

Rules for this project: reference and donor samples are disjoint and never
first-degree relatives (KING ≥ 0.177; 2b refuses otherwise); 1KG samples may
be references and donors here, unlike in the genotyping-harmonization
benchmark.
