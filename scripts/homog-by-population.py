#!/usr/bin/env python3
"""
homog-by-population.py -- 4a's homogeneous panel, population by population.
READ-ONLY. For each super-pop and each population in lai_ref_panel_samples.kK.tsv:

  n | in the homogeneous panel | would pass at the cut used HERE (so a table
      built at 0.93 shows what 0.95 would keep) | with no first-degree
      relative in the panel
  | why the others are not, counted in this order:
      donor or a donor's KING >= 0.177 relative / 1e duplicate drop /
      strongest supervised component is ANOTHER super-pop / below the cut
  | strongest-component counts | mean supervised Q of every component >= 0.02
  | how many of the population are anchors, and for which cluster

then each cluster's anchors by cohort/population, and for EVERY population
that has simulation donors: its size, how many have no first-degree relative,
the donors, how many of them have one, and how many relatives they keep out of
the reference. scripts/donor-switch-plan.py finds the donor sets that keep out
the fewest.

Reads what 4a and 1e wrote: the Q table, reference_ids/*_anchors.kK.txt,
reference_ids/*_simulation.txt, reference_ids/panel_duplicate_drop_ids.txt,
admixture_kK/donor_relatives.kK.tsv and king/king.kin0. The cut is 4a's
default, 0.94 for every super-pop from 17 Nov (before: 0.95 for EUR/EAS/AFR,
and SAS from 16 Nov; 0.93 for the rest); --cut SP=Q matches a run made with
others, e.g. --k 8 --cut EUR=0.95 --cut AMR=0.93 for the 16 Nov K=8 table.

  python3 -I scripts/homog-by-population.py
  python3 -I scripts/homog-by-population.py --table <tsv> --refs reference_ids --homog <04_homogeneity_panel>
"""
import argparse
import csv
import glob
import os
from collections import Counter, defaultdict

ROOT = os.environ.get("PROJECT_ROOT",
                      "/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument("--homog", default=os.path.join(ROOT, "04_homogeneity_panel"))
    ap.add_argument("--table", default="", help="default: <homog>/lai_ref_panel_samples.k<K>.tsv")
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--cut", action="append", default=[], help="SP=Q, overriding 4a's default cut for SP")
    a = ap.parse_args()
    tsv = a.table or os.path.join(a.homog, f"lai_ref_panel_samples.k{a.k}.tsv")

    cut = {}                                                      # 4a's cut: 0.94 for every super-pop
    for c in a.cut:
        sp, q = c.split("=", 1)
        cut[sp.upper()] = float(q)
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    comps = [c[6:] for c in rows[0] if c.startswith("Q_sup_")]
    print(f"{tsv}: {len(rows)} samples; supervised components: {' '.join(comps)}")

    donors = set()
    for f in glob.glob(os.path.join(a.refs, "*_simulation.txt")):
        donors |= set(open(f).read().split())
    rel_f = os.path.join(a.homog, f"admixture_k{a.k}", f"donor_relatives.k{a.k}.tsv")
    rel = ({l.split("\t")[0] for l in open(rel_f) if l.strip() and not l.startswith("relative")}
           if os.path.exists(rel_f) else set())
    dup_f = os.path.join(a.refs, "panel_duplicate_drop_ids.txt")
    dup = ({l.split()[0] for l in open(dup_f) if l.strip() and not l.startswith("#")}
           if os.path.exists(dup_f) else set())
    anchors = {}
    for f in sorted(glob.glob(os.path.join(a.refs, f"*_anchors.k{a.k}.txt"))):
        for s in open(f).read().split():
            anchors[s] = os.path.basename(f).split("_")[0].upper()
    kin_f = os.path.join(a.homog, "king", "king.kin0")
    has_rel = set()
    kin_of = defaultdict(set)
    if os.path.exists(kin_f):
        with open(kin_f) as fh:
            h = fh.readline().lstrip("#").split()
            ix = {c: i for i, c in enumerate(h)}
            for line in fh:
                f = line.split()
                if float(f[ix["KINSHIP"]]) >= 0.177:
                    x, y = f[ix["IID1"]], f[ix["IID2"]]
                    has_rel.update((x, y))
                    kin_of[x].add(y); kin_of[y].add(x)
    else:
        print(f"WARNING: no {kin_f}: relative counts below are all zero")
    print(f"donors {len(donors)}, donor relatives kept out {len(rel)} ({rel_f}), duplicate drops {len(dup)},"
          f" anchors {len(anchors)}, samples with a KING >= 0.177 relative {len(has_rel)}")

    in_panel = {r["sample"] for r in rows} - dup
    relfree = {s for s in in_panel if not (kin_of[s] & in_panel)}

    def passes(r, t):
        s = r["sample"]
        if s in donors or s in rel or s in dup:
            return False
        if r["cohort"] == "MXB":                    # 4a includes every non-donor MXB sample
            return True
        return r["argmax_pop_sup"] == r["super_pop"] and float(r["max_Q_sup"]) >= t

    by = defaultdict(list)
    for r in rows:
        by[(r["super_pop"], r["pop"])].append(r)
    print("\nper population: n | in the homogeneous panel | would pass at the cut used here | with no"
          " first-degree relative | not in the panel, because: donor or donor's relative / duplicate /"
          " strongest component is another super-pop / below the cut | strongest-component counts |"
          " mean supervised Q (components >= 0.02) | anchors")
    order = ["EUR", "AFR", "EAS", "SAS", "MEN", "AMR", "OCE"]
    for sp in order + sorted({s for s, _ in by} - set(order)):
        pops = sorted(p for (s, p) in by if s == sp)
        if not pops:
            continue
        t = cut.get(sp, 0.94)
        allr = [r for p in pops for r in by[(sp, p)]]
        print(f"\n== {sp} (cut used here {t}): {len(allr)} samples,"
              f" {sum(r['included'] == '1' for r in allr)} in the panel,"
              f" {sum(passes(r, t) for r in allr)} would pass at {t}")
        for p in pops:
            rs = by[(sp, p)]
            n, inc = len(rs), sum(r["included"] == "1" for r in rs)
            npass = sum(passes(r, t) for r in rs)
            nfree = sum(r["sample"] in relfree for r in rs)
            dr = du = wrong = below = 0
            for r in rs:
                if r["included"] == "1":
                    continue
                s = r["sample"]
                if s in donors or s in rel:
                    dr += 1
                elif s in dup:
                    du += 1
                elif r["argmax_pop_sup"] != sp:
                    wrong += 1
                elif float(r["max_Q_sup"]) < t:
                    below += 1
            am = " ".join(f"{k}:{v}" for k, v in Counter(r["argmax_pop_sup"] for r in rs).most_common(3))
            mq = {c: sum(float(r[f"Q_sup_{c}"]) for r in rs) / n for c in comps}
            mqs = " ".join(f"{c}:{v:.2f}" for c, v in sorted(mq.items(), key=lambda kv: -kv[1]) if v >= 0.02)
            na = Counter(anchors[r["sample"]] for r in rs if r["sample"] in anchors)
            print(f"  {p:18s} n={n:<4d} in={inc:<4d} pass@{t}={npass:<4d} relfree={nfree:<4d}"
                  f" not: {dr}/{du}/{wrong}/{below}  strongest {am}  meanQ {mqs}"
                  f"  anchors {' '.join(f'{k}:{v}' for k, v in na.items()) or '-'}")

    print("\n== each cluster's anchors, by cohort/population")
    pop_of = {r["sample"]: f"{r['cohort']}/{r['pop']}" for r in rows}
    for sp in sorted(set(anchors.values())):
        c = Counter(pop_of.get(s, "NOT-IN-TABLE") for s, x in anchors.items() if x == sp)
        print(f"  {sp} ({sum(c.values())}): " + ", ".join(f"{k} {v}" for k, v in c.most_common()))

    print("\n== every population with simulation donors: relatives at KING >= 0.177 (duplicate drops not counted)")
    print(f"  {'population':18s} {'n':>4s} {'relfree':>7s} {'donors':>6s} {'w/ relative':>11s} {'relatives kept out':>18s}")
    pop_of_s = {r["sample"]: r["pop"] for r in rows}
    dpops = Counter(pop_of_s.get(s, "NOT-IN-TABLE") for s in donors)
    for p, nd in sorted(dpops.items()):
        ids_p = [s for s, q in pop_of_s.items() if q == p]
        cur = [s for s in donors if pop_of_s.get(s) == p]
        out = set().union(*(kin_of[s] & in_panel for s in cur)) - donors if cur else set()
        print(f"  {p:18s} {len(ids_p):4d} {sum(s in relfree for s in ids_p):7d} {nd:6d} "
              f"{sum(bool(kin_of[s] & in_panel) for s in cur):11d} {len(out):18d}")
    print("  -> scripts/donor-switch-plan.py proposes, for every list, the donors that keep out the fewest")


if __name__ == "__main__":
    main()
