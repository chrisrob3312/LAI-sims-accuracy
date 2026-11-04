# Reference and donor panels: who may be a reference, who may be a donor

Decided 2026-11-04. This governs the **genotyping-pipeline benchmark set**
(project B). The LAI simulation and accuracy paper (project A) keeps its own
three comparison arms and is described at the bottom.

## The one rule everything else serves

**A haplotype used as a simulation donor must never appear in the reference
panel the ancestry method is given.** If it does, the method is scored on
haplotypes it was handed as reference, and the benchmark is circular — it
measures memorisation, not inference.

`4a_build_homogeneity_panel_kN_clm.sh` enforces this for the homogeneous panel
by unioning every `reference_ids/*_simulation.txt` and excluding that union,
for every ancestry. `2b_simulation_clm.sh` has the matching gate on the donor
side. Neither is advisory.

## The policy, per super-pop

Two tiers, because the two situations are genuinely different.

### EUR, EAS, AFR — reference is 1000 Genomes only

These have enough 1000 Genomes samples to build a reference panel without
touching anything else. So:

| role | drawn from |
|---|---|
| LAI reference | homogeneous **1000 Genomes** samples only, Q ≥ 0.95 |
| simulation donors | HGDP, SGDP, Almarri, MXB — never 1000 Genomes |

Donor/reference disjointness here is guaranteed **by source**, not by a split:
no 1000 Genomes sample is ever a donor, no non-1000G sample is ever in these
reference panels. That is the strongest form of the guarantee available and it
costs nothing, because the samples are plentiful.

Q ≥ 0.95 is the strict tier (`HOMOG_STRICT_SUPERPOPS`): these groups are
richly sampled, so the stricter cut costs them little and buys cleaner
reference haplotypes.

### SAS, AMR, MEN, OCE — reference and donors both come from the homogeneous pool

There are not enough 1000 Genomes samples to build a reference panel for these
from 1000G alone. So both roles draw from the same homogeneous pool
(HGDP + 1KG + MXB + SGDP + Almarri), and disjointness has to be **constructed**
rather than inherited:

| role | drawn from |
|---|---|
| LAI reference | homogeneous pool, Q ≥ 0.93, **disjoint half** |
| simulation donors | homogeneous pool, Q ≥ 0.93, **the other half** |

Q ≥ 0.93 is the floor tier. At 0.95 these groups lose members the panel cannot
spare, and a thin reference panel for exactly the ancestries whose local
ancestry is hardest to call is the thing this project exists to measure.

**AMR specifically includes MXB in the reference panel**, so the Mexican-like
cohort is tested against a reference that actually represents it rather than
against HGDP Indigenous American alone.

#### The cost of this tier, stated plainly

Splitting one limited pool into two halves makes **both** halves smaller. The
reference panel for SAS/AMR/MEN/OCE is therefore thinner than it would be if
donors came from elsewhere, and the simulated cohorts draw from fewer distinct
donor haplotypes — which raises the chance that simulated individuals become
cryptic relatives of each other. That is a real cost, accepted deliberately
because the alternative (reusing samples in both roles) is not a benchmark at
all. The per-super-pop counts on each side must be reported, not assumed
adequate, and the kinship among simulated individuals must be checked.

## The precedent: MXB is already split this way

`split_mxb_samples.R` splits the 50 MXB WGS genomes 25 donors / 25 reference,
stratified by inferred genetic cluster so neither half is a biased sample of
the whole:

| cluster | n | donors | reference |
|---|---|---|---|
| Nahua | 18 | 9 | 9 |
| Zapotec/Mazatec | 10 | 5 | 5 |
| Maya | 7 | 3 | 4 |
| Totonac | 6 | 3 | 3 |
| Tzotzil | 4 | 2 | 2 |
| Tarahumara | 3 | 2 | 1 |
| Huichol | 2 | 1 | 1 |

The split within a cluster is deterministic (lexicographic sample ID, first
*n* to donors) so it reproduces without seeded RNG. Verified: 25 and 25, zero
overlap.

**The SAS / MEN / OCE splits should follow this shape**: stratify by
population within the super-pop, split deterministically inside each stratum.
Splitting a super-pop at random instead would let one population land mostly
on the donor side and another mostly on the reference side, so the benchmark
would be measuring cross-population transfer while reporting it as
within-super-pop accuracy.

## Build order

Everything below waits on the K = 8 supervised Q file from stage 3. Nothing
here can be built from the candidate lists, because membership of the
homogeneous pool is exactly what stage 3 decides.

1. Stage 3 writes `lai_ref_panel_samples.k8.tsv` and `homog_8pop_samples.txt`.
2. Partition the homogeneous pool by super-pop and by source (1000 Genomes vs
   everything else).
3. EUR/EAS/AFR: reference = the 1000 Genomes slice. No split needed.
4. SAS/AMR/MEN/OCE: stratified deterministic split into donor and reference
   halves, following the MXB recipe.
5. Write `*_simulation_b2.txt` and `*_reference_b2.txt` per super-pop.
6. Run `scripts/check-donor-panel-disjointness.sh` over the result. It must
   report zero overlap before anything is simulated.

## Project A (the LAI accuracy paper) is unchanged

Project A compares three reference-panel arms and needs all three to coexist:

| arm | reference panel |
|---|---|
| hgdp+pel | the curated `reference_ids/*_rfmix.txt` lists |
| homog | the homogeneous panel, tiered Q (0.95 EUR/EAS/AFR, 0.93 others) |
| mxb | HGDP-NAT + MXB (`amr_hgdpmxb_rfmix.txt`) |

This is why `4a` writes discovered anchors to `<sp>_anchors.k8.txt` and never
to `<sp>_rfmix.txt`: overwriting the curated lists would destroy the hgdp+pel
arm, which is one of the three things the paper compares.
