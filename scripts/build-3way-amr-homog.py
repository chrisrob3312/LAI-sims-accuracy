#!/usr/bin/env python3
"""
build-3way-amr-homog.py -- the AMR reference for the 3-way (NAT/EUR/AFR) LAI
HOMOG arms (3b_homog, 3b_chip_homog, FLARE homog), which read
reference_ids/amr_homog_rfmix.txt.

Decided 17 Nov 2026: the 3-way pipeline keeps July's K=5 homogeneity calls
and adds the AMR samples that are NEW since July -- SGDP, the restored HGDP
samples, the 1KG relatives now kept in the panel -- where the current run
calls them homogeneous. So:

  base  = July's K=5 table (reference_ids/lai_ref_panel_samples.tsv),
          super_pop AMR and included == 1 (88 samples: MXB force-included,
          HGDP NAT, 9 PEL), still in the panel
  new   = the current table's AMR samples with included == 1 that are not in
          July's table at all
  out   = base + new, minus every simulation donor, donor's relative
          (KING >= 0.177), 1e duplicate drop and named exclusion -- the set 4a
          writes to admixture_k<K>/simulation_donors.all.txt

A sample that WAS in July's table keeps July's call; the current run is used
only for samples July never saw.

  python3 -I scripts/build-3way-amr-homog.py --base reference_ids/lai_ref_panel_samples.tsv \\
      --table <04_homogeneity_panel>/lai_ref_panel_samples.k7.tsv \\
      --exclude <04_homogeneity_panel>/admixture_k7/simulation_donors.all.txt \\
      --out reference_ids/amr_homog_rfmix.3way.txt
"""
import argparse
import csv
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="July's K=5 lai_ref_panel_samples.tsv")
    ap.add_argument("--table", required=True, help="the current lai_ref_panel_samples.k<K>.tsv")
    ap.add_argument("--exclude", required=True, help="donors + relatives + duplicates + exclusions")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-q", type=float, default=None,
                    help="re-cut the NEW samples at this supervised Q instead of the table's "
                         "included flag (argmax AMR and max_Q_sup >= Q; MXB kept as 4a keeps it) "
                         "-- e.g. 0.94 on the 16 Nov K=8 table, cut at 0.93, as a preview")
    a = ap.parse_args()

    base = list(csv.DictReader(open(a.base), delimiter="\t"))
    cur = list(csv.DictReader(open(a.table), delimiter="\t"))
    for name, rows in (("base", base), ("table", cur)):
        miss = {"sample", "super_pop", "included"} - set(rows[0]) if rows else {"rows"}
        if miss:
            sys.exit(f"ERROR: --{name} lacks columns {sorted(miss)}")
    excl = {l.split()[0] for l in open(a.exclude) if l.strip() and not l.startswith("#")}
    panel = {r["sample"]: r for r in cur}
    base_ids = {r["sample"] for r in base}

    base_amr = [r["sample"] for r in base if r["super_pop"] == "AMR" and r["included"] == "1"]
    absent = [s for s in base_amr if s not in panel]
    dropped = [s for s in base_amr if s in panel and s in excl]
    keep_base = [s for s in base_amr if s in panel and s not in excl]
    def passes(r):
        if a.min_q is None:
            return r["included"] == "1"
        if r.get("cohort") == "MXB":                 # 4a includes every non-donor MXB sample
            return True
        return r.get("argmax_pop_sup") == "AMR" and float(r.get("max_Q_sup", 0)) >= a.min_q
    new = [r["sample"] for r in cur if r["super_pop"] == "AMR" and passes(r)
           and r["sample"] not in base_ids and r["sample"] not in excl]
    out = keep_base + new
    if len(set(out)) != len(out):
        sys.exit("ERROR: a sample is listed twice")
    with open(a.out, "w") as g:
        g.write("\n".join(out) + "\n")

    def comp(ids):
        return ", ".join(f"{k} {v}" for k, v in Counter(
            f"{panel[s].get('cohort', '?')}/{panel[s].get('pop', '?')}" for s in ids).most_common()) or "-"
    print(f"3-way AMR reference -> {a.out}: {len(out)} samples")
    print(f"  July K=5 AMR homogeneous: {len(base_amr)}; kept {len(keep_base)}; "
          f"{len(dropped)} now donors, donors' relatives, duplicate drops or excluded; "
          f"{len(absent)} not in the current panel{(' (' + ', '.join(absent[:5]) + ')') if absent else ''}")
    print(f"    kept:    {comp(keep_base)}")
    if dropped:
        print(f"    dropped: {comp(dropped)}")
    print(f"  new since July and homogeneous now"
          f"{f' (supervised AMR Q >= {a.min_q})' if a.min_q is not None else ''}: {len(new)}")
    print(f"    added:   {comp(new)}")


if __name__ == "__main__":
    main()
