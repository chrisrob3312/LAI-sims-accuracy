#!/usr/bin/env python3
"""
homog-cut-preview.py -- READ-ONLY. How the homogeneous pool, and the part of it
the benchmark can use as DONORS, changes with the homogeneity cut.

For each super-pop and each cut (default 0.90, 0.92, 0.94) it counts the
samples 4a would call homogeneous at that cut -- strongest supervised
component is the sample's own super-pop, at a share >= the cut; MXB always,
as 4a does -- and splits them into

  donor-eligible  not 1000 Genomes and not on sgdp_server_panel_ids.txt
  reference-only  1000 Genomes or a server-panel id: the benchmark never makes
                  them donors (their haplotypes are in the imputation panels)

and counts, among them, the LAI simulation donors (reference_ids/
*_simulation.txt) and their relatives. Those are out of the LAI project's panel
(4a's `included`) but IN the benchmark's pool (4a's `homogeneous`, 18 Nov): the
two projects choose from one pool, each under its own rules. 1e's duplicate
drops and the named exclusions (panel_exclude_ids.txt) are out of both.

Then, for every cut below the highest, the donor-eligible samples it ADDS, by
population: how many, their mean share of their own super-pop's component,
and their strongest OTHER component. A donor at share q hands the simulation
1 - q of its genome labelled with the wrong ancestry, so the mean of 1 - q
over each cut's donor-eligible pool is printed beside its size.

The table's own pool count is printed beside the recomputed numbers as a check
on this script: at the table's own cut (0.94 for every super-pop from 17 Nov)
"donor-eligible + reference-only" must equal the table's `homogeneous` count
(tables from before 18 Nov have only `included`, which also leaves out the LAI
donors, so there it equals donor + reference - lai). Writes nothing.

  python3 -I scripts/homog-cut-preview.py                       # K=7 table
  python3 -I scripts/homog-cut-preview.py --k 8 --cuts 0.88,0.90,0.92,0.94,0.95
  python3 -I scripts/homog-cut-preview.py --table <tsv> --refs reference_ids --homog <04_homogeneity_panel>
"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter, defaultdict

ROOT = os.environ.get("PROJECT_ROOT",
                      "/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank")
ORDER = ["EUR", "SAS", "AFR", "EAS", "MEN", "AMR", "OCE"]


def fc(t):
    """0.9 -> '0.90', 0.945 -> '0.945'."""
    return f"{t:.2f}" if round(t, 2) == t else f"{t:g}"


def ids_from(path, first_col_header=None):
    """First whitespace field of every non-empty, non-# line. Missing file =
    empty set, said so by the caller."""
    out = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()[0]
            if first_col_header and f == first_col_header:
                continue
            out.add(f)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k", type=int, default=7)
    ap.add_argument("--homog", default=os.path.join(ROOT, "04_homogeneity_panel"))
    ap.add_argument("--table", default="",
                    help="default: <homog>/lai_ref_panel_samples.k<K>.tsv (the global set's)")
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--cuts", default="0.90,0.92,0.94")
    ap.add_argument("--kg-cohort", default="1KG")
    ap.add_argument("--server-panel-ids", default="",
                    help="default: <refs>/sgdp_server_panel_ids.txt")
    ap.add_argument("--show-pops", type=int, default=12,
                    help="populations listed per super-pop for the samples a lower cut adds")
    a = ap.parse_args()

    tsv = a.table or os.path.join(a.homog, f"lai_ref_panel_samples.k{a.k}.tsv")
    try:
        cuts = sorted({float(x) for x in a.cuts.split(",") if x.strip()})
    except ValueError:
        sys.exit(f"ERROR: --cuts '{a.cuts}' is not a comma-separated list of numbers")
    if not cuts or not all(0 < c <= 1 for c in cuts):
        sys.exit(f"ERROR: --cuts must be numbers in (0, 1], got '{a.cuts}'")

    with open(tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        rows = list(rd)
    need = ["sample", "cohort", "pop", "super_pop", "max_Q_sup", "argmax_pop_sup", "included"]
    miss = [c for c in need if c not in cols]
    comps = [c[len("Q_sup_"):] for c in cols if c.startswith("Q_sup_")]
    if miss or not comps:
        sys.exit(f"ERROR: {tsv} lacks columns {miss or 'Q_sup_*'}; header: {cols}")
    print(f"[table] {tsv}: {len(rows)} samples")
    print(f"[table] header (raw): {' '.join(cols)}")
    for r in rows[:2]:
        print("[table] row (raw, key fields): "
              + " ".join(f"{c}={r[c]}" for c in need)
              + " " + " ".join(f"Q_sup_{c}={r['Q_sup_' + c]}" for c in comps))

    # Who 4a leaves out, read from the same files 4a reads.
    def opt_ids(path, label, header=None):
        if os.path.exists(path):
            s = ids_from(path, header)
            print(f"[input] {label}: {len(s)} ids ({path})")
            return s
        print(f"[input] {label}: NONE -- {path} not found; counted as empty")
        return set()

    dup = opt_ids(os.path.join(a.refs, "panel_duplicate_drop_ids.txt"), "1e duplicate drops")
    named = opt_ids(os.path.join(a.refs, "panel_exclude_ids.txt"), "named exclusions")
    lai = set()
    for f in sorted(glob.glob(os.path.join(a.refs, "*_simulation.txt"))):
        s = ids_from(f)
        print(f"[input] LAI simulation donors: {len(s)} ids ({f})")
        lai |= s
    lai_rel = opt_ids(os.path.join(a.homog, f"admixture_k{a.k}", f"donor_relatives.k{a.k}.tsv"),
                      "relatives of LAI donors (4a)", header="relative")
    held = lai | lai_rel
    spi = opt_ids(a.server_panel_ids or os.path.join(a.refs, "sgdp_server_panel_ids.txt"),
                  "server-panel ids (donor-ineligible)")

    def own_q(r):
        try:
            return float(r.get(f"Q_sup_{r['super_pop']}", "nan"))
        except ValueError:
            return float("nan")

    def other(r):
        sp = r["super_pop"]
        best = max(((float(r[f"Q_sup_{c}"]), c) for c in comps if c != sp), default=(0.0, "-"))
        return best[1], best[0]

    def homog_at(r, t):
        if r["cohort"] == "MXB":                      # 4a includes every non-donor MXB sample
            return True
        return r["argmax_pop_sup"] == r["super_pop"] and float(r["max_Q_sup"]) >= t

    def role(r):
        s = r["sample"]
        if s in dup or s in named:
            return "dropped"
        if r["cohort"].strip().upper() == a.kg_cohort.upper() or s in spi:
            return "refonly"
        return "donor"

    by_sp = defaultdict(list)
    for r in rows:
        by_sp[r["super_pop"]].append(r)
    sps = [s for s in ORDER if s in by_sp] + sorted(set(by_sp) - set(ORDER))

    pool_col = "homogeneous" if "homogeneous" in cols else "included"
    print("\n== the benchmark pool at each cut (strongest component = own super-pop, share >= cut; MXB always)")
    print("   donor = donor-eligible (not 1KG, not a server-panel id); ref = reference-only (1KG or "
          "server-panel id);\n   lai = of donor + ref, LAI simulation donors and their relatives (in "
          "the benchmark pool, out of the LAI panel);\n   wrong = mean (1 - own-component share) over "
          f"the donor-eligible pool; table = the table's {pool_col}=1 count")
    W = 22                                   # one cut's block: donor ref lai wrong
    print(f"  {'':<4} {'':>5} {'':>5} | " + " | ".join(f"{'at ' + fc(t):^{W}}" for t in cuts))
    print(f"  {'sp':<4} {'n':>5} {'table':>5} | " + " | ".join(
        f"{'donor':>6} {'ref':>4} {'lai':>4} {'wrong':>5}" for _ in cuts))
    added = {}
    for sp in sps:
        rs = by_sp[sp]
        n_tab = sum(r[pool_col].strip() in ("1", "true", "True") for r in rs)
        cells = []
        for t in cuts:
            hom = [r for r in rs if homog_at(r, t) and role(r) != "dropped"]
            c = Counter(role(r) for r in hom)
            n_lai = sum(r["sample"] in held for r in hom)
            dq = [own_q(r) for r in hom if role(r) == "donor" and r["cohort"] != "MXB"]
            wrong = f"{sum(1 - q for q in dq) / len(dq):.3f}" if dq else "-"
            cells.append(f"{c['donor']:>6} {c['refonly']:>4} {n_lai:>4} {wrong:>5}")
        print(f"  {sp:<4} {len(rs):>5} {n_tab:>5} | " + " | ".join(cells))
        # samples each lower cut adds to the donor-eligible pool
        for lo, hi in zip(cuts, cuts[1:]):
            add = [r for r in rs if r["cohort"] != "MXB" and role(r) == "donor"
                   and r["argmax_pop_sup"] == sp and lo <= float(r["max_Q_sup"]) < hi]
            if add:
                added.setdefault(sp, []).append((lo, hi, add))
    mxb = [r for r in rows if r["cohort"] == "MXB" and role(r) != "dropped"]
    if mxb:
        # 4a counts every MXB sample homogeneous whatever its Q (July's curated
        # set); the contract's rule is >= the cut of one component. How many
        # that rule would leave out, by their own AMR share:
        print(f"  (MXB: {len(mxb)} samples, all counted homogeneous as 4a does, "
              f"{sum(r['sample'] in held for r in mxb)} of them LAI donors or relatives; under "
              "the cut by their own AMR share: "
              + ", ".join(f"{fc(t)}: {sum(own_q(r) < t for r in mxb)}" for t in cuts) + ")")

    print("\n== what each lower cut ADDS to the donor-eligible pool, by population")
    print("   n | mean own-component share | strongest other component (mean share)")
    if not added:
        print("  nothing: no donor-eligible sample sits between the cuts")
    for sp in sps:
        for lo, hi, add in added.get(sp, []):
            print(f"  {sp} [{fc(lo)}, {fc(hi)}): +{len(add)}")
            byp = defaultdict(list)
            for r in add:
                byp[f"{r['cohort']}/{r['pop']}"].append(r)
            for p, prs in sorted(byp.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:a.show_pops]:
                oc = Counter(other(r)[0] for r in prs).most_common(1)[0][0]
                om = sum(other(r)[1] for r in prs) / len(prs)
                mq = sum(own_q(r) for r in prs) / len(prs)
                print(f"      {p:<28} {len(prs):>3}   {sp} {mq:.3f}   {oc} {om:.3f}")
            if len(byp) > a.show_pops:
                print(f"      ... {len(byp) - a.show_pops} more population(s)")

    # Below the lowest cut but still strongest in their own component: how
    # much further the pool could grow, and with what.
    lo = cuts[0]
    print(f"\n== donor-eligible, strongest component = own super-pop, but under {fc(lo)} (not admitted at any cut above)")
    for sp in sps:
        rest = [r for r in by_sp[sp] if r["cohort"] != "MXB" and role(r) == "donor"
                and r["argmax_pop_sup"] == sp and float(r["max_Q_sup"]) < lo]
        if not rest:
            continue
        bins = Counter("0.85+" if own_q(r) >= 0.85 else "0.70-0.85" if own_q(r) >= 0.70 else "<0.70"
                       for r in rest)
        print(f"  {sp}: {len(rest)}  ({', '.join(f'{k} {v}' for k, v in sorted(bins.items(), reverse=True))})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
