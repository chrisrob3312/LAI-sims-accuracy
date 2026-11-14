# Reference and donor panels: who may be a reference, who may be a donor

Decided 2026-11-04/05. There are **three** panel sets built from the same
homogeneous pool, under different rules, for different purposes. Confusing
them is the main way this goes wrong, so they are named here and kept apart
everywhere else.

| # | set | built by | rule |
|---|---|---|---|
| 1 | LAI / ADMIXTURE reference panels (project A) | `split-homog-donors-reference.py --mode lai` | even stratified split, every super-pop, no source rule |
| 2 | Benchmarking panels (project B) | `... --mode benchmark` | same split, **plus** 1000 Genomes may never be a donor |
| 3 | The full homogeneous panel | `4a` stage C, unsplit | everything that clears its tier. Never used as a benchmark reference |

## The one rule all three serve

**A haplotype used as a simulation donor must never appear in the reference
panel the method is given.** Otherwise the method is scored on haplotypes it
was handed, and the result measures memorisation.

`4a` enforces this for the homogeneous panel by excluding the union of every
`reference_ids/*_simulation.txt`, for every ancestry. `2b` has the matching
gate on the donor side. The splitter adds two more checks, both fatal: a
sample on both sides within a super-pop, and a sample that is a reference for
one super-pop and a donor for another — same person, same haplotypes.

## The split, for sets 1 and 2

One formula, applied per super-pop to the number of homogeneous samples `N`:

```
reference = clamp(round(0.5 × N), 15, 65)
donors    = N − reference
```

| N | ref | donors | % donors |
|---|---|---|---|
| 400 | 65 | 335 | 84% |
| 316 | 65 | 251 | 79% |
| 200 | 65 | 135 | 68% |
| 130 | 65 | 65 | 50% |
| 100 | 50 | 50 | 50% |
| 40 | 20 | 20 | 50% |
| 30 | 15 | 15 | 50% |
| <30 | — | — | refuses to split |

Even until the cap bites at N = 130, then donors take the remainder — which
clears 2/3 donors above N ≈ 195 without that being a separate rule. The cap of
65 gives 130 haplotypes per ancestry, comfortably above the curated 50-sample
panels. Below 30 the script writes nothing and exits non-zero: 10 reference
individuals is not a panel, and the real options there are decisions, not
defaults.

### Both sides are stratified by population

Neither side is a biased sample of its super-pop. An unstratified split would
let one population land mostly on the donor side and another mostly on the
reference side, so the benchmark would measure cross-population transfer while
reporting it as within-super-pop accuracy.

Largest-remainder allocation, deterministic throughout (populations by name,
samples by ID), so it reproduces with no seeded RNG. Flooring cannot erase a
population from the reference side — verified on a 40-sample OCE whose
smallest strata are 2 samples.

Measured on the projected panel: MEN matches within **1.0 percentage point**
across 8 populations, OCE within **5.0** across 6, AMR within **10.0**.

`split_mxb_samples.R` is the precedent: 50 MXB genomes, 25/25, stratified by
inferred cluster (Nahua 9/9, Zapotec/Mazatec 5/5, Maya 3/4 …), deterministic
within each stratum. Verified 25 and 25, zero overlap.

## What set 2 adds: 1000 Genomes may never be a donor

Not a general tidiness rule — a specific one. **The imputation servers'
reference panels (TOPMed, Michigan, All of Us) all contain 1000 Genomes.** A
1000 Genomes simulation donor therefore has its own haplotypes sitting inside
the server's panel, and the imputation is scored on data the server has
memorised. No donor/reference split repairs that, because the leak is into a
panel we do not control. So the exclusion is **by source**, and it applies to
every super-pop, not only the ones with enough 1000 Genomes to build a
reference from.

Set 1 has no such rule. The LAI and ADMIXTURE panels are not imputed against
those servers, so the circularity does not exist there and withholding 1000
Genomes from donors would only shrink the pool for nothing.

### How set 2 keeps the compositions matched anyway

The non-1000G samples are split evenly and stratified, so both sides carry the
same mix; the reference is then topped up from 1000 Genomes to its target.
Taking the reference whole from one source and donors from another would make
any difference between them a difference between *datasets*, which the
benchmark would read as an ancestry effect.

**Two consequences, both reported by the script rather than discovered later.**

A population that exists only in 1000 Genomes can only ever be a reference.
The script names every population that lands on one side only. PEL is the AMR
case and accounts for that group's entire 10-point gap; GIH is the SAS case.
These belong in the methods; they are not defects to tune away.

Donor pools for the 1000G-rich super-pops are bounded by their non-1000G
supply. On the projection: **SAS 24 donors, AFR 47, EUR 59** — against 335,
569 and 555 in set 1. SAS at 24 is thin.

## The SAS consequence, and what absorbs it

A thin SAS donor pool limits how many simulated SAS-containing individuals can
be made before they become cryptic relatives of each other. That cost does not
fall evenly across the benchmark's three readings:

- **Imputation accuracy, calibration and retention** are scored on real test
  genomes as well as simulated ones, so SAS can lean on the **real admixed SAS
  genomes** here.
- **Local-ancestry accuracy** cannot. Real genomes have a WGS truth but nobody
  recorded which ancestry each of their segments came from, so they carry no
  tracts and are never scored for LAI. Only simulated individuals have
  recorded truth tracts.

So a thin SAS donor pool costs **SAS local-ancestry precision specifically**,
while the imputation metrics stay well supplied. That is the right trade to
make if one has to be made, and it is worth saying in the methods rather than
reporting a wide SAS LAI interval without explaining it.

The same reasoning applies to OCE, harder: with ~40 homogeneous samples it
splits 20/20, and 20 donors is few. The decision to make once the real count
lands is whether to simulate an OCE-containing admixture model at all, or to
use OCE as a reference/control group only. Surfaced rather than buried in a
ratio.

## Set 3: the full homogeneous panel

Everything clearing its tier — 0.95 for EUR/EAS/AFR, 0.93 for AMR/OCE/MEN/SAS
(`HOMOG_STRICT_SUPERPOPS`) — kept whole and unsplit, as a resource for later
work.

**It must never be used as a benchmark reference panel.** It contains the
donors. Keeping it intact and keeping it out of project B are the same
decision, and the reason it is written down here is that it is otherwise the
obvious thing for a later reader to reach for.

## Build order

Everything waits on the K = 8 supervised Q file from stage 3. None of it can
be built from the candidate lists, because membership of the homogeneous pool
is exactly what stage 3 decides.

1. Stage 3 writes `lai_ref_panel_samples.k8.tsv` and `homog_8pop_samples.txt`.
2. `split-homog-donors-reference.py --mode lai --suffix a` → set 1.
3. `split-homog-donors-reference.py --mode benchmark --suffix b2` → set 2.
4. Read the printed composition table and the per-super-pop counts before
   simulating anything. The donor counts are what bound the cohort size.
5. `scripts/check-donor-panel-disjointness.sh` must report zero overlap.

## Project A's three comparison arms are unchanged

| arm | reference panel |
|---|---|
| hgdp+pel | the curated `reference_ids/*_rfmix.txt` lists |
| homog | the homogeneous panel, tiered Q |
| mxb | HGDP-NAT + MXB (`amr_hgdpmxb_rfmix.txt`) |

This is why `4a` writes discovered anchors to `<sp>_anchors.k8.txt` and never
to `<sp>_rfmix.txt`: overwriting the curated lists would destroy the hgdp+pel
arm, which is one of the three things the paper compares.
