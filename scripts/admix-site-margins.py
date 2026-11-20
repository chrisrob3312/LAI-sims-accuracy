#!/usr/bin/env python3
"""admix-site-margins.py -- which soft-union sites ADMIXTURE may use, judged on
each SUPER-POP and each COHORT as a whole (the rule set on 17 Nov 2026):

  * every super-pop, all of its cohorts pooled, <= --max-superpop-miss (0.07)
  * every cohort, all of its super-pops pooled, <= --max-cohort-miss  (0.15)

Optional extra caps, off by default (1.0): --max-group-miss on each cohort x
super-pop group, --max-total-miss on the whole panel. A margin (super-pop or
cohort) smaller than --min-n samples is not judged on its own; its samples
still count toward the other margin. Missingness = 1 - AN / (2n), on the
PRE-phasing genotypes. A position is kept only if every record at it is kept.

Groups are COHORT_SUPERPOP (scripts/build-missingness-groups.py): the cohort is
the text before the first "_", the super-pop the text after it.

stdin: `bcftools query -f '%CHROM\\t%POS\\t%INFO/AN_<G1>...\\n'` in the order of
--groups' distinct values (--print-format prints the format string).
Writes --out (CHROM<TAB>POS<TAB>POS, plink2 --extract bed1) and --summary.
"""
import argparse
import sys
from collections import Counter, OrderedDict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", required=True, help="sample<TAB>COHORT_SUPERPOP")
    ap.add_argument("--print-format", action="store_true")
    ap.add_argument("--max-superpop-miss", type=float, default=0.07)
    ap.add_argument("--max-cohort-miss", type=float, default=0.15)
    ap.add_argument("--max-group-miss", type=float, default=1.0)
    ap.add_argument("--max-total-miss", type=float, default=1.0)
    ap.add_argument("--min-n", type=int, default=5)
    ap.add_argument("--out")
    ap.add_argument("--summary")
    a = ap.parse_args()

    sizes = OrderedDict()
    with open(a.groups) as fh:
        for l in fh:
            f = l.split()
            if len(f) >= 2:
                sizes[f[1]] = sizes.get(f[1], 0) + 1
    G = list(sizes)
    if a.print_format:
        print("%CHROM\\t%POS" + "".join(f"\\t%INFO/AN_{g}" for g in G) + "\\n")
        return
    if not (a.out and a.summary):
        sys.exit("ERROR: --out and --summary are required unless --print-format")
    bad_names = [g for g in G if "_" not in g]
    if bad_names:
        sys.exit(f"ERROR: group(s) {bad_names} are not COHORT_SUPERPOP")

    co_of = {g: g.split("_", 1)[0] for g in G}
    sp_of = {g: g.split("_", 1)[1] for g in G}
    margins = {}                                  # name -> (kind, member group indices, n)
    for kind, key in (("superpop", sp_of), ("cohort", co_of)):
        for m in sorted(set(key.values())):
            idx = [i for i, g in enumerate(G) if key[g] == m]
            margins[f"{kind}:{m}"] = (kind, idx, sum(sizes[G[i]] for i in idx))
    cap = {"superpop": a.max_superpop_miss, "cohort": a.max_cohort_miss}
    judged = {k: v for k, v in margins.items() if v[2] >= a.min_n}
    n_total = sum(sizes.values())

    fail = Counter(); only = Counter()
    n_rec = n_keep = n_tot_fail = 0
    pos_ok = OrderedDict()
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if len(f) != 2 + len(G):
            sys.exit(f"ERROR: expected {2 + len(G)} columns, got {len(f)}: {line[:200]!r}")
        an = [0 if x == "." else int(x) for x in f[2:]]
        n_rec += 1
        bad = []
        for name, (kind, idx, n) in judged.items():
            if 1 - sum(an[i] for i in idx) / (2 * n) > cap[kind] + 1e-12:
                bad.append(name)
        if a.max_group_miss < 1.0:
            for i, g in enumerate(G):
                if sizes[g] >= a.min_n and 1 - an[i] / (2 * sizes[g]) > a.max_group_miss + 1e-12:
                    bad.append(f"group:{g}")
        tot_bad = 1 - sum(an) / (2 * n_total) > a.max_total_miss + 1e-12
        for b in bad:
            fail[b] += 1
        if len(bad) == 1 and not tot_bad:
            only[bad[0]] += 1
        n_tot_fail += tot_bad
        ok = not bad and not tot_bad
        n_keep += ok
        key = (f[0], f[1])
        pos_ok[key] = pos_ok.get(key, True) and ok

    kept = [k for k, v in pos_ok.items() if v]
    with open(a.out, "w") as g:
        for chrom, pos in kept:
            g.write(f"{chrom}\t{pos}\t{pos}\n")
    with open(a.summary, "w") as g:
        g.write(f"# max_superpop_miss={a.max_superpop_miss} max_cohort_miss={a.max_cohort_miss} "
                f"max_group_miss={a.max_group_miss} max_total_miss={a.max_total_miss} min_n={a.min_n}\n")
        g.write(f"# records={n_rec} records_kept={n_keep} positions={len(pos_ok)} "
                f"positions_kept={len(kept)} records_failing_total={n_tot_fail}\n")
        g.write("margin\tn_samples\tjudged\trecords_failing\trecords_failing_only_this\n")
        for name, (kind, idx, n) in margins.items():
            g.write(f"{name}\t{n}\t{int(name in judged)}\t{fail[name]}\t{only[name]}\n")
        for gname in G:
            if fail[f"group:{gname}"]:
                g.write(f"group:{gname}\t{sizes[gname]}\t1\t{fail['group:' + gname]}\t{only['group:' + gname]}\n")
    pct = 100 * len(kept) / max(len(pos_ok), 1)
    print(f"[admix-margins] {n_rec:,} records at {len(pos_ok):,} positions; kept {len(kept):,} ({pct:.1f}%) "
          f"with every super-pop <= {a.max_superpop_miss:.0%} and every cohort <= {a.max_cohort_miss:.0%} missing",
          file=sys.stderr)
    for name, n in fail.most_common(8):
        print(f"[admix-margins]   {name:20s} (n={margins[name][2] if name in margins else sizes.get(name[6:], '?')}) "
              f"over its limit at {n:,} records ({only[name]:,} of them for no other reason)", file=sys.stderr)
    small = [f"{k}({v[2]})" for k, v in margins.items() if k not in judged]
    if small:
        print(f"[admix-margins]   not judged alone (n < {a.min_n}): {', '.join(small)}", file=sys.stderr)


if __name__ == "__main__":
    main()
