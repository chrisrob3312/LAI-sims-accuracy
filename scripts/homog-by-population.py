#!/usr/bin/env python3
"""
homog-by-population.py -- 4a's homogeneous panel, population by population.
READ-ONLY. For each super-pop and each population in lai_ref_panel_samples.kK.tsv:

  n | in the homogeneous panel | why the others are not, counted in this order:
      donor or a donor's KING >= 0.177 relative / 1e duplicate drop /
      strongest supervised component is ANOTHER super-pop / below the cut
  | strongest-component counts | mean supervised Q of every component >= 0.02
  | how many of the population are anchors, and for which cluster

then each cluster's anchors by cohort/population, and how many IBS and YRI
in the panel have no first-degree relative (candidates for relative-free
simulation donors).

Reads what 4a and 1e wrote: the Q table, reference_ids/*_anchors.kK.txt,
reference_ids/*_simulation.txt, reference_ids/panel_duplicate_drop_ids.txt,
admixture_kK/donor_relatives.kK.tsv and king/king.kin0. The cuts are 4a's
defaults (HOMOG_THRESHOLD 0.93, HOMOG_THRESHOLD_STRICT 0.95 for EUR, EAS, AFR);
pass --cut SP=Q to match a run made with others.

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
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--homog", default=os.path.join(ROOT, "04_homogeneity_panel"))
    ap.add_argument("--table", default="", help="default: <homog>/lai_ref_panel_samples.k<K>.tsv")
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--cut", action="append", default=[], help="SP=Q, overriding 4a's default cut for SP")
    a = ap.parse_args()
    tsv = a.table or os.path.join(a.homog, f"lai_ref_panel_samples.k{a.k}.tsv")

    cut = {"EUR": 0.95, "EAS": 0.95, "AFR": 0.95}          # 4a's tiers; every other super-pop 0.93
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
    if os.path.exists(kin_f):
        with open(kin_f) as fh:
            h = fh.readline().lstrip("#").split()
            ix = {c: i for i, c in enumerate(h)}
            for line in fh:
                f = line.split()
                if float(f[ix["KINSHIP"]]) >= 0.177:
                    has_rel.update((f[ix["IID1"]], f[ix["IID2"]]))
    print(f"donors {len(donors)}, donor relatives kept out {len(rel)} ({rel_f}), duplicate drops {len(dup)},"
          f" anchors {len(anchors)}, samples with a KING >= 0.177 relative {len(has_rel)}")

    by = defaultdict(list)
    for r in rows:
        by[(r["super_pop"], r["pop"])].append(r)
    print("\nper population: n | in the homogeneous panel | not, because: donor or donor's relative / duplicate /"
          " strongest component is another super-pop / below the cut | strongest-component counts |"
          " mean supervised Q (components >= 0.02) | anchors")
    order = ["EUR", "AFR", "EAS", "SAS", "MEN", "AMR", "OCE"]
    for sp in order + sorted({s for s, _ in by} - set(order)):
        pops = sorted(p for (s, p) in by if s == sp)
        if not pops:
            continue
        t = cut.get(sp, 0.93)
        print(f"\n== {sp} (cut {t}): {sum(len(by[(sp, p)]) for p in pops)} samples,"
              f" {sum(r['included'] == '1' for p in pops for r in by[(sp, p)])} in the panel")
        for p in pops:
            rs = by[(sp, p)]
            n, inc = len(rs), sum(r["included"] == "1" for r in rs)
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
            print(f"  {p:18s} n={n:<4d} in={inc:<4d} not: {dr}/{du}/{wrong}/{below}  strongest {am}  meanQ {mqs}"
                  f"  anchors {' '.join(f'{k}:{v}' for k, v in na.items()) or '-'}")

    print("\n== each cluster's anchors, by cohort/population")
    pop_of = {r["sample"]: f"{r['cohort']}/{r['pop']}" for r in rows}
    for sp in sorted(set(anchors.values())):
        c = Counter(pop_of.get(s, "NOT-IN-TABLE") for s, x in anchors.items() if x == sp)
        print(f"  {sp} ({sum(c.values())}): " + ", ".join(f"{k} {v}" for k, v in c.most_common()))

    print("\n== IBS and YRI in the panel with no KING >= 0.177 relative (candidates for relative-free donors)")
    for p in ("IBS", "YRI"):
        ids = [r["sample"] for r in rows if r["pop"] == p]
        cur = [s for s in ids if s in donors]
        print(f"  {p}: {len(ids)} in the panel, {sum(s not in has_rel for s in ids)} with no first-degree relative;"
              f" current donors {len(cur)}, of which {sum(s in has_rel for s in cur)} have one")


if __name__ == "__main__":
    main()
