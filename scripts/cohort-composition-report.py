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

WHY DONOR REUSE IS COMPUTED HERE AND NOT GUESSED. The number of distinct donor
haplotypes bounds how many independent individuals can be simulated before
they become cryptic relatives of each other. Reuse alone does not say how bad
that is: it interacts with GENERATIONS. At gen 2 a haplotype carries ~70 long
segments, so two individuals drawing on the same donor share long tracts; at
gen 16 it carries ~560 short ones and the same reuse fragments into near
independence. So the report prints reuse AND generations AND a severity that
accounts for both, rather than a ratio that looks alarming for LAT3 and
harmless for SASEUR when the truth is the reverse.

Usage:
  python3 -I scripts/cohort-composition-report.py \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k8.tsv \\
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
    ("SASEUR", ("SAS", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    ("NAT1",   ("NAT",),              (1.0,),               0,  60, "control"),
    ("EUR1",   ("EUR",),              (1.0,),               0,  60, "control"),
    ("AFR1",   ("AFR",),              (1.0,),               0,  60, "control"),
    ("EAS1",   ("EAS",),              (1.0,),               0,  60, "control"),
    ("SAS1",   ("SAS",),              (1.0,),               0,  60, "control"),
    ("MEN1",   ("MEN",),              (1.0,),               0,  60, "control"),
]

# The simulated components map onto super-pop donor pools. NAT is the AMR pool.
COMPONENT_POOL = {"NAT": "AMR", "EUR": "EUR", "AFR": "AFR",
                  "MEN": "MEN", "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}

# Rough genetic length of the autosomes, Morgans. Used only to turn
# generations into an expected segment count, which is what decides whether a
# given reuse level matters.
GENOME_MORGANS = 35.0


def severity(reuse, gen):
    """How much a reuse level matters AT THIS NUMBER OF GENERATIONS."""
    if reuse <= 1.0:
        return "ok", "every individual can have its own donor haplotype"
    if gen == 0:
        # A zero-generation control is one donor's haplotypes with no
        # recombination at all, so there is no mosaic to fragment and the
        # segment argument below does not apply. Reuse above 1 means two
        # controls share a WHOLE haplotype -- kinship ~0.25, the
        # half-sib/avuncular range, not subtle at all.
        #
        # The fix is not fewer donors' worth of controls: each control is
        # DIPLOID, so n donors give C(n,2) distinct haplotype PAIRS (46 donors
        # -> 1,035). Drawing each control's two haplotypes from two different
        # donors, and no pair twice, removes the problem without reducing n.
        return "HIGH", (f"gen 0: controls would share whole haplotypes "
                        f"(kinship ~0.25). Draw each control's 2 haplotypes "
                        f"from 2 DIFFERENT donors, no pair reused")
    # Expected independent segments per simulated haplotype.
    seg = max(1.0, gen * GENOME_MORGANS)
    # Reuse is tolerable when the mosaic is fine enough that two individuals
    # sharing a donor share only short, scattered pieces of it.
    ratio = seg / reuse
    if ratio >= 100:
        return "ok", f"~{seg:.0f} segments/hap fragments the reuse"
    if ratio >= 30:
        return "watch", f"~{seg:.0f} segments/hap; check sim-kinship"
    return "HIGH", (f"~{seg:.0f} segments/hap against {reuse:.1f}x reuse -- "
                    f"expect cryptic relatedness")


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
        models = []
        for r in read_tsv(args.models):
            models.append((
                r["model"],
                tuple(x.strip() for x in r["components"].split(",")),
                tuple(float(x) for x in r["proportions"].split(",")),
                int(float(r["generations"])), int(float(r["n_individuals"])),
                r.get("priority", "")))

    panel = read_tsv(args.panel)
    inc_col = "included" if "included" in panel[0] else None
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
    w(f"  homogeneous (included): {len(homog)}")
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
      f"{'donors':>8}{'reuse':>8}  severity")
    tot = 0
    for name, comps, props, gen, n, _prio in models:
        tot += n
        for comp, p in zip(comps, props):
            pool = COMPONENT_POOL.get(comp, comp)
            nd = len(donors.get(pool, []))
            if nd == 0:
                w(f"  {name:<8}{gen:>5}{n:>6}  {comp:<10}{p:>6.2f}"
                  f"{'--':>8}{'--':>8}  no donor list for {pool}")
                continue
            reuse = (n * 2 * p) / (2 * nd)
            tag, why = severity(reuse, gen)
            w(f"  {name:<8}{gen:>5}{n:>6}  {comp:<10}{p:>6.2f}"
              f"{nd:>8}{reuse:>7.1f}x  {tag} ({why})")
    w("")
    w(f"  total simulated individuals: {tot}")
    w("")
    w("  Severity accounts for GENERATIONS, not reuse alone. A haplotype after")
    w("  g generations carries roughly g x 35 independent segments, so the")
    w("  same reuse that produces cryptic relatives at gen 2 fragments into")
    w("  near-independence at gen 16. A high reuse on a low-generation model")
    w("  is the one to act on; the simulated-individual kinship check is what")
    w("  confirms it.")
    w("")

    # ---- 4. what was removed --------------------------------------------
    w("-" * 78)
    w("4. REMOVED, AND WHY")
    w("-" * 78)
    for label, path, note in [
        ("duplicates dropped (KING)",
         os.path.join(args.refs, "panel_duplicate_drop_ids.txt"),
         "same person twice under two ids; one copy kept"),
        ("excluded by name",
         os.path.join(args.refs, "panel_exclude_ids.txt"),
         "e.g. CHMI_CHMI3_WGS2, a hydatidiform-mole cell line"),
        ("1000G-via-SGDP (donor-ineligible)",
         os.path.join(args.refs, "sgdp_server_panel_ids.txt"),
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
