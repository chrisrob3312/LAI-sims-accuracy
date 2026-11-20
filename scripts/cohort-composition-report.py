#!/usr/bin/env python3
"""
cohort-composition-report.py

Write ONE text file describing the finished cohort: who is in the panel, who
became a donor, who became a reference, how hard each donor pool is reused by
each model, and what was dropped and why.

It is a reporting step, not a build step: it reads what the pipeline produced
and changes nothing. Run it at the end and keep the output -- it is the
record the methods section is written from, and most of what it holds cannot
be reconstructed later (the drop lists get overwritten, the counts move every
time a threshold changes).

WHY RELATEDNESS IS COMPUTED HERE AND NOT GUESSED. A finite donor pool makes
simulated individuals cryptic relatives of each other, and reuse (draws per
donor haplotype) does not say how much. What does is how 6b draws: simu-mix
gives every breakpoint a fresh donor haplotype and never lets two haplotypes
of one run copy the same donor at a marker, and the gen-0 controls take each
chromosome's two haplotypes independently from the pool. So two individuals
in DIFFERENT runs copy the same donor haplotype at a marker with probability
sum_P p_P^2 / (2 x donors_P), which IS their expected kinship, and two in the
same run share nothing. Until 20 Nov this report tagged gen-2 models HIGH
from reuse x generations; a Monte Carlo of the two programs' code (mixer's
pulse and random mating, simu-mix's draw rule) at the 18 Nov donor counts
found no MENEUR or EASEUR pair at 3rd degree (closest 0.04), while EUR1 and
MEN1 -- the gen-0 controls on 32 and 39 donors -- put 2-3% of pairs there,
and none anywhere reached 2nd degree. The tag is calibrated on that run.

Usage:
  python3 -I scripts/cohort-composition-report.py \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k7.tsv \\
      --refs  reference_ids --suffix b2 \\
      --out   results/cohort_composition.txt
"""
import argparse
import os
import sys
from collections import defaultdict, OrderedDict

# Document 19, Table 1. Overridable with --models (a TSV with the same
# columns) so the report follows the set actually built.
DEFAULT_MODELS = [
    # name,    components,            proportions,        gen, n,   priority
    ("LAT3",   ("NAT", "EUR", "AFR"), (0.50, 0.45, 0.05),  16, 600, "primary"),
    ("CAR3",   ("NAT", "EUR", "AFR"), (0.10, 0.50, 0.40),  16, 200, "secondary"),
    ("AFAM2",  ("AFR", "EUR"),        (0.80, 0.20),         7, 200, "secondary"),
    ("MENA3",  ("EUR", "MEN", "AFR"), (0.30, 0.50, 0.20),  30, 200, "secondary"),
    ("MENEUR", ("MEN", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    ("EASEUR", ("EAS", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    # SASEUR and SAS1 dropped 18 Nov: too few homogeneous SAS donors; SAS is
    # read from real genomes (reference_ids/models.tsv says how to restore).
    ("NAT1",   ("NAT",),              (1.0,),               0,  60, "control"),
    ("EUR1",   ("EUR",),              (1.0,),               0,  60, "control"),
    ("AFR1",   ("AFR",),              (1.0,),               0,  60, "control"),
    ("EAS1",   ("EAS",),              (1.0,),               0,  60, "control"),
    ("MEN1",   ("MEN",),              (1.0,),               0,  60, "control"),
]

# The simulated components map onto super-pop donor pools. NAT is the AMR pool.
COMPONENT_POOL = {"NAT": "AMR", "EUR": "EUR", "AFR": "AFR",
                  "MEN": "MEN", "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}

def expected_kinship(comps, props, n_donors):
    """Mean kinship of two simulated individuals of one model that are not
    in the same simu-mix run (any two gen-0 controls): at a marker, two
    haplotypes copy the same donor haplotype with probability
    sum_P p_P^2 / (2 x donors_P). Pairs inside one run share nothing."""
    return sum(p * p / (2.0 * n_donors[c])
               for c, p in zip(comps, props) if p > 0 and n_donors.get(c))


def kinship_tag(phi, gen):
    """Calibrated on the 20 Nov Monte Carlo (one full run of each model at
    the 18 Nov donor counts; pairs at >= 0.0442, KING's 3rd degree):
      MENEUR gen 2  expected 0.0071   0 of 19,900   (closest 0.042)
      EASEUR gen 2  expected 0.0050   0 of 19,900   (closest 0.029)
      EUR1   gen 0  expected 0.0156  59 of 1,770    (closest 0.072)
      MEN1   gen 0  expected 0.0125  33 of 1,770
      NAT1 / AFR1 / EAS1  0.0066 / 0.0061 / 0.0044   2 / 4 / 3 of 1,770
    No pair at 2nd degree (0.0884) anywhere. Gen-0 controls take whole
    chromosomes, so their tail is wider at the same mean."""
    if phi >= 0.0221:
        return "HIGH", "the average pair is 4th-degree; expect 2nd-degree pairs"
    if gen == 0 and phi >= 0.004:
        return "watch", ("whole chromosomes are shared, so some pairs reach "
                         "3rd degree; none closer")
    if phi >= 0.008:
        return "watch", "some pairs may reach 3rd degree; none closer"
    return "ok", "no pair expected at 3rd degree"


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            ff = line.rstrip("\n").split("\t")
            if len(ff) < len(header):
                ff += [""] * (len(header) - len(ff))
            rows.append(dict(zip(header, ff)))
    return rows


def read_list(path):
    try:
        with open(path) as fh:
            return [l.strip() for l in fh
                    if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--suffix", default="b2")
    ap.add_argument("--models", default=None,
                    help="TSV: model components proportions generations "
                         "n_individuals priority (comma lists). Defaults to "
                         "Document 19 Table 1.")
    ap.add_argument("--out", default="results/cohort_composition.txt")
    ap.add_argument("--kg-cohort", default="1KG")
    args = ap.parse_args()

    models = DEFAULT_MODELS
    if args.models:
        # reference_ids/models.tsv opens with a '#' comment block, so the
        # header is the first line that is not a comment (read_tsv would take
        # the first comment line as the header).
        with open(args.models) as fh:
            body = [l.rstrip("\n") for l in fh
                    if l.strip() and not l.lstrip().startswith("#")]
        if not body or not body[0].startswith("model\t"):
            sys.exit(f"ERROR: {args.models}: no 'model' header row after the comments")
        mhead = body[0].split("\t")
        models = []
        for r in (dict(zip(mhead, l.split("\t"))) for l in body[1:]):
            models.append((
                r["model"],
                tuple(x.strip() for x in r["components"].split(",")),
                tuple(float(x) for x in r["proportions"].split(",")),
                int(float(r["generations"])), int(float(r["n_individuals"])),
                r.get("priority", "")))

    panel = read_tsv(args.panel)
    # The benchmark's pool is 4a's `homogeneous` flag (18 Nov); `included` is
    # the LAI project's panel, which also leaves out that project's donors.
    inc_col = next((c for c in ("homogeneous", "included") if c in panel[0]), None)
    homog = [r for r in panel if not inc_col or r.get(inc_col) == "1"]

    L = []
    def w(s=""):
        L.append(s)

    w("=" * 78)
    w("COHORT COMPOSITION")
    w("=" * 78)
    w(f"panel table : {args.panel}")
    w(f"lists       : {args.refs}/*_{{reference,simulation}}_{args.suffix}.txt")
    w("")

    # ---- 1. the panel ----------------------------------------------------
    w("-" * 78)
    w("1. THE MERGED PANEL")
    w("-" * 78)
    w(f"  samples in panel      : {len(panel)}")
    w(f"  homogeneous ({inc_col + '=1' if inc_col else 'no flag column'}): {len(homog)}")
    w(f"  NOT homogeneous       : {len(panel) - len(homog)}")
    w("")
    w("  A sample that is not homogeneous is not necessarily admixed: at a K")
    w("  with no cluster for its ancestry it simply has nowhere to be")
    w("  homogeneous. Read the by-super-pop table before calling these admixed.")
    w("")
    by_cohort = defaultdict(lambda: [0, 0])
    for r in panel:
        is_h = (not inc_col) or r.get(inc_col) == "1"
        by_cohort[r.get("cohort", "NA")][0 if is_h else 1] += 1
    w(f"  {'cohort':<10}{'homogeneous':>13}{'other':>9}{'total':>9}")
    for c in sorted(by_cohort):
        h, o = by_cohort[c]
        w(f"  {c:<10}{h:>13}{o:>9}{h+o:>9}")
    w("")

    by_sp = defaultdict(lambda: [0, 0])
    for r in panel:
        is_h = (not inc_col) or r.get(inc_col) == "1"
        by_sp[r.get("super_pop", "NA")][0 if is_h else 1] += 1
    w(f"  {'super_pop':<10}{'homogeneous':>13}{'other':>9}{'total':>9}")
    for s in sorted(by_sp):
        h, o = by_sp[s]
        w(f"  {s:<10}{h:>13}{o:>9}{h+o:>9}")
    w("")

    # ---- 2. donors and reference ----------------------------------------
    w("-" * 78)
    w("2. DONOR AND REFERENCE SPLIT")
    w("-" * 78)
    donors, refs = {}, {}
    for sp in sorted({r.get("super_pop", "NA") for r in panel}):
        d = read_list(os.path.join(args.refs, f"{sp.lower()}_simulation_{args.suffix}.txt"))
        f = read_list(os.path.join(args.refs, f"{sp.lower()}_reference_{args.suffix}.txt"))
        if d is not None:
            donors[sp] = d
        if f is not None:
            refs[sp] = f
    if not donors:
        w("  NO SPLIT LISTS FOUND. Run split-homog-donors-reference.py first;")
        w("  everything below that depends on them is omitted.")
        w("")
    else:
        w(f"  {'super_pop':<10}{'homog':>8}{'reference':>11}{'donors':>8}"
          f"{'%don':>7}{'unused':>8}")
        for sp in sorted(set(donors) | set(refs)):
            nh = by_sp[sp][0]
            nr, nd = len(refs.get(sp, [])), len(donors.get(sp, []))
            pct = 100.0 * nd / (nr + nd) if (nr + nd) else 0
            w(f"  {sp:<10}{nh:>8}{nr:>11}{nd:>8}{pct:>6.0f}%{nh-nr-nd:>8}")
        w("")
        w("  'unused' are homogeneous samples in neither role. For the")
        w("  benchmark set they are the real_homogeneous TEST pool, minus any")
        w("  that are 1000 Genomes (a 1000 Genomes test sample would be")
        w("  imputed against its own haplotypes in the server's panel).")
        w("")
        kg_ids = {r["sample"] for r in panel
                  if r.get("cohort", "").upper() == args.kg_cohort.upper()}
        used = {s for v in donors.values() for s in v} | \
               {s for v in refs.values() for s in v}
        w(f"  {'super_pop':<10}{'unused':>8}{'of which 1KG':>14}{'TEST pool':>11}")
        for sp in sorted(by_sp):
            unused = [r["sample"] for r in homog
                      if r.get("super_pop") == sp and r["sample"] not in used]
            n_kg = sum(1 for s in unused if s in kg_ids)
            w(f"  {sp:<10}{len(unused):>8}{n_kg:>14}{len(unused)-n_kg:>11}")
        w("")

    # ---- 3. models and reuse --------------------------------------------
    w("-" * 78)
    w("3. MODELS, GENERATIONS AND DONOR REUSE")
    w("-" * 78)
    w("  reuse = haplotype-draws needed / donor haplotypes available")
    w("        = (n_individuals x 2 x proportion) / (2 x n_donors)")
    w("")
    w(f"  {'model':<8}{'gen':>5}{'n':>6}  {'component':<10}{'prop':>6}"
      f"{'donors':>8}{'reuse':>8}")
    tot = 0
    kin_rows = []
    for name, comps, props, gen, n, _prio in models:
        tot += n
        nd_of = {}
        for comp, p in zip(comps, props):
            pool = COMPONENT_POOL.get(comp, comp)
            nd = len(donors.get(pool, []))
            if nd == 0:
                w(f"  {name:<8}{gen:>5}{n:>6}  {comp:<10}{p:>6.2f}"
                  f"{'--':>8}{'--':>8}  no donor list for {pool}")
                continue
            nd_of[comp] = nd
            reuse = (n * 2 * p) / (2 * nd)
            w(f"  {name:<8}{gen:>5}{n:>6}  {comp:<10}{p:>6.2f}"
              f"{nd:>8}{reuse:>7.1f}x")
        if len(nd_of) == len(comps):
            kin_rows.append((name, gen, n, expected_kinship(comps, props, nd_of)))
    w("")
    w(f"  total simulated individuals: {tot}")
    w("")
    w("  Reuse alone does not say how related the simulated individuals")
    w("  are. Expected kinship does, given how 6b draws (see the header):")
    w("  two individuals from different simu-mix runs (any two gen-0")
    w("  controls) copy the same donor haplotype at a marker with")
    w("  probability sum_P p_P^2 / (2 x donors_P); inside one run, never.")
    w("")
    w(f"  {'model':<8}{'gen':>5}{'n':>6}{'kinship':>10}  tag")
    for name, gen, n, phi in kin_rows:
        tag, why = kinship_tag(phi, gen)
        w(f"  {name:<8}{gen:>5}{n:>6}{phi:>10.4f}  {tag} ({why})")
    w("")
    w("  KING bins: 0.0442 3rd degree, 0.0884 2nd, 0.177 1st. Tags are")
    w("  calibrated on a Monte Carlo of mixer + simu-mix (20 Nov); 6e")
    w("  measures the realised distribution once the cohort is built.")
    w("")

    # ---- 4. what was removed --------------------------------------------
    w("-" * 78)
    w("4. REMOVED, AND WHY")
    w("-" * 78)
    # These three live in the repo's reference_ids/ whatever --refs says: a
    # scratch run (REFS=$S) put the split lists elsewhere, and the report then
    # printed "step not run" for steps that had run.
    repo_refs = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "reference_ids")
    def fixed_list(name):
        for d in (args.refs, repo_refs):
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
        return os.path.join(args.refs, name)
    for label, path, note in [
        ("duplicates dropped (KING)",
         fixed_list("panel_duplicate_drop_ids.txt"),
         "same person twice under two ids; one copy kept"),
        ("excluded by name",
         fixed_list("panel_exclude_ids.txt"),
         "e.g. CHMI_CHMI3_WGS2, a hydatidiform-mole cell line"),
        ("1000G-via-SGDP (donor-ineligible)",
         fixed_list("sgdp_server_panel_ids.txt"),
         "inside TOPMed/Michigan panels under an NA/HG id"),
    ]:
        ids = read_list(path)
        if ids is None:
            w(f"  {label:<36} (no file -- step not run)")
        else:
            w(f"  {label:<36} {len(ids):>5}   {note}")
    dups = None
    dp = os.path.join(args.refs, "known_duplicate_pairs.tsv")
    if os.path.exists(dp):
        rows = read_tsv(dp)
        dups = sum(1 for r in rows if r.get("alias_in_panel") == "1")
        w(f"  {'duplicate pairs named by metadata':<36} {dups:>5}   "
          f"cross-check: KING should find at least these")
    w("")

    # ---- 5. relatedness policy ------------------------------------------
    w("-" * 78)
    w("5. RELATEDNESS POLICY")
    w("-" * 78)
    w("  Only duplicates and MZ twins (KING kinship > 0.354) are removed.")
    w("  First, second and third degree relatives are KEPT and reported.")
    w("")
    w("  This is deliberate and matches the downstream analysis: the GWAS uses")
    w("  GENESIS with a GRM, which models relatedness rather than requiring")
    w("  unrelated samples. Excluding everything that looks related would")
    w("  preferentially remove individuals from less-studied populations,")
    w("  where background relatedness is higher and reference panels are")
    w("  thinner -- shrinking exactly the groups this project exists to")
    w("  measure. The relatedness distribution belongs in the methods.")
    w("")
    w("=" * 78)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as g:
        g.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\n[report] written to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
