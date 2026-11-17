#!/usr/bin/env python3
"""
build-3way-amr-homog.py -- the AMR reference for the 3-way (NAT/EUR/AFR) LAI
HOMOG arms (3b_homog, 3b_chip_homog, FLARE homog), which read
reference_ids/amr_homog_rfmix.txt.

Two sources of homogeneity calls meet here:

  July   reference_ids/lai_ref_panel_samples.tsv: the K=5 run the August pilot
         used (anchors taken BY NAME from the 6 Jul lists; 88 AMR included: the
         50 MXB force-included, HGDP NAT, 9 PEL)
  now    4a HOMOG_SET=lai3way: K=5 again, on the expanded panel (HGDP, 1KG, MXB,
         SGDP), MEN / OCE / no-super-pop samples left out, anchors FOUND by the
         unsupervised fit (central, Q >= 0.94), homogeneous at Q >= 0.94. Never
         the K=7 global run: "should not judge new amr homog by k=7, that's not
         consistent with earlier work" (17 Nov 2026).

--mode decides how they combine (4a passes HOMOG_3WAY_MODE):
  union      (default) July's AMR panel kept, plus EVERY AMR sample homogeneous
             now -- new since July, or judged otherwise in July. "Keep that run
             plus the additional samples that broaden the homog panel", with the
             rerun able to make the panel fuller.
  keep-july  July's AMR panel kept, plus only samples July never saw (SGDP, the
             restored HGDP, 1KG relatives). July's verdict stands for everyone
             July judged.
  current    this run's AMR homogeneous samples only.
Every mode then removes each simulation donor, donor's relative (KING >=
0.177), 1e duplicate drop and named exclusion: --exclude, the set 4a writes to
admixture_k<K>/simulation_donors.all.txt. A sample July included that fails
now is kept under union and keep-july and reported.

  python3 -I scripts/build-3way-amr-homog.py --base reference_ids/lai_ref_panel_samples.tsv \\
      --table <04_homogeneity_panel>/lai3way_k5/lai_ref_panel_samples.k5.tsv \\
      --exclude <04_homogeneity_panel>/lai3way_k5/admixture_k5/simulation_donors.all.txt \\
      --out reference_ids/amr_homog_rfmix.3way.txt [--mode union]
"""
import argparse
import csv
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="July's K=5 lai_ref_panel_samples.tsv")
    ap.add_argument("--table", required=True, help="this run's lai_ref_panel_samples.k<K>.tsv")
    ap.add_argument("--exclude", required=True, help="donors + relatives + duplicates + exclusions")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=("union", "keep-july", "current"), default="union")
    ap.add_argument("--min-q", type=float, default=None,
                    help="judge this run's samples at this supervised Q instead of the table's "
                         "included flag (argmax AMR and max_Q_sup >= Q; MXB kept as 4a keeps it) "
                         "-- for a preview on another run's table")
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

    def passes(r):
        if a.min_q is None:
            return r["included"] == "1"
        if r.get("cohort") == "MXB":                 # 4a includes every non-donor MXB sample
            return True
        return r.get("argmax_pop_sup") == "AMR" and float(r.get("max_Q_sup", 0)) >= a.min_q

    july = [r["sample"] for r in base if r["super_pop"] == "AMR" and r["included"] == "1"]
    absent = [s for s in july if s not in panel]
    july_out = [s for s in july if s in panel and s in excl]
    july_in = [s for s in july if s in panel and s not in excl]
    now = [r["sample"] for r in cur if r["super_pop"] == "AMR" and passes(r) and r["sample"] not in excl]
    now_set = set(now)
    new_unseen = [s for s in now if s not in base_ids]                      # July never judged them
    new_rejudged = [s for s in now if s in base_ids and s not in set(july)]  # July: not homogeneous
    july_fail_now = [s for s in july_in if s not in now_set]                 # July: yes; now: no

    if a.mode == "union":
        out = july_in + [s for s in now if s not in set(july_in)]
    elif a.mode == "keep-july":
        out = july_in + new_unseen
    else:
        out = list(now)
    if len(set(out)) != len(out):
        sys.exit("ERROR: a sample is listed twice")
    with open(a.out, "w") as g:
        g.write("\n".join(out) + "\n")

    def comp(ids):
        return ", ".join(f"{k} {v}" for k, v in Counter(
            f"{panel[s].get('cohort', '?')}/{panel[s].get('pop', '?')}" for s in ids).most_common()) or "-"
    q = f" (supervised AMR Q >= {a.min_q})" if a.min_q is not None else ""
    print(f"3-way AMR reference -> {a.out}: {len(out)} samples (mode {a.mode})")
    print(f"  July K=5 AMR homogeneous: {len(july)}; {len(july_in)} still eligible, "
          f"{len(july_out)} now donors, donors' relatives, duplicate drops or excluded, "
          f"{len(absent)} not in this run{(' (' + ', '.join(absent[:5]) + ')') if absent else ''}")
    print(f"    eligible: {comp(july_in)}")
    if july_out:
        print(f"    out:      {comp(july_out)}")
    print(f"  homogeneous in this run{q}, eligible: {len(now)}")
    print(f"    new since July ({len(new_unseen)}): {comp(new_unseen)}")
    print(f"    July judged NOT homogeneous ({len(new_rejudged)}): {comp(new_rejudged)}"
          f"{'' if a.mode == 'union' else '  [not added: --mode union adds them]'}")
    print(f"  July homogeneous but not in this run ({len(july_fail_now)}): {comp(july_fail_now)}"
          f"{'  [kept: July' + chr(39) + 's call]' if a.mode != 'current' else '  [dropped: --mode current]'}")


if __name__ == "__main__":
    main()
