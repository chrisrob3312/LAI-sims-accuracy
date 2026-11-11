#!/usr/bin/env python3
"""admix-site-filter.py -- which soft-union sites ADMIXTURE may use.

A site is kept only if, on the PRE-phasing genotypes:
  * every cohort x super-pop group with >= --min-group-n samples has
    missingness <= --max-group-miss, and
  * the whole panel has missingness <= --max-total-miss.
A position is kept only if every record at it is kept (a position that holds
two records -- the sources disagreed on the alleles -- is dropped if either is).

Why both rules: missingness here is mostly a whole cohort absent from a site
(its source file had no record), and a small cohort can be entirely absent
while the panel-wide fraction stays tiny (Almarri is 3.5% of samples). SHAPEIT5
then fills every such genotype, so the phased panel has no missingness left to
filter and ADMIXTURE treats the filled-in values as observed.

stdin: `bcftools query -f '%CHROM\\t%POS[...]\\t%INFO/AN_<G1>\\t%INFO/AN_<G2>...\\n'`
       with the groups in the order of --groups' distinct values (see --print-format).
Writes --out (CHROM<TAB>POS<TAB>POS, plink2 --extract bed1) and --summary (TSV).
"""
import argparse
import sys
from collections import Counter, OrderedDict


def load_groups(path):
    sizes = OrderedDict()
    with open(path) as fh:
        for l in fh:
            f = l.split()
            if len(f) >= 2:
                sizes[f[1]] = sizes.get(f[1], 0) + 1
    return sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", required=True, help="sample<TAB>GROUP (build-missingness-groups.py)")
    ap.add_argument("--print-format", action="store_true",
                    help="print the bcftools query format string for --groups and exit")
    ap.add_argument("--max-group-miss", type=float, default=0.10)
    ap.add_argument("--max-total-miss", type=float, default=0.07)
    ap.add_argument("--min-group-n", type=int, default=5,
                    help="groups smaller than this are not judged on their own (they still count overall)")
    ap.add_argument("--out")
    ap.add_argument("--summary")
    a = ap.parse_args()

    sizes = load_groups(a.groups)
    G = list(sizes)
    if a.print_format:
        print("%CHROM\\t%POS" + "".join(f"\\t%INFO/AN_{g}" for g in G) + "\\n")
        return
    if not (a.out and a.summary):
        sys.exit("ERROR: --out and --summary are required unless --print-format")
    judged = [g for g in G if sizes[g] >= a.min_group_n]
    n_total = sum(sizes.values())

    fail_group = Counter()      # records where this group broke its rule
    only_group = Counter()      # records where it was the ONLY reason
    n_rec = n_rec_keep = n_total_fail = 0
    pos_ok = OrderedDict()      # (chrom, pos) -> all records ok so far
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if len(f) != 2 + len(G):
            sys.exit(f"ERROR: expected {2 + len(G)} columns, got {len(f)}: {line[:200]!r}")
        an = [0 if x == "." else int(x) for x in f[2:]]
        n_rec += 1
        bad = []
        for g, v in zip(G, an):
            if sizes[g] >= a.min_group_n and 1 - v / (2 * sizes[g]) > a.max_group_miss + 1e-12:
                bad.append(g)
        tot_miss = 1 - sum(an) / (2 * n_total)
        tot_bad = tot_miss > a.max_total_miss + 1e-12
        for g in bad:
            fail_group[g] += 1
        if len(bad) == 1 and not tot_bad:
            only_group[bad[0]] += 1
        if tot_bad:
            n_total_fail += 1
        ok = not bad and not tot_bad
        n_rec_keep += ok
        key = (f[0], f[1])
        pos_ok[key] = pos_ok.get(key, True) and ok

    kept = [k for k, v in pos_ok.items() if v]
    with open(a.out, "w") as g:
        for chrom, pos in kept:
            g.write(f"{chrom}\t{pos}\t{pos}\n")
    with open(a.summary, "w") as g:
        g.write(f"# max_group_miss={a.max_group_miss} max_total_miss={a.max_total_miss} "
                f"min_group_n={a.min_group_n}\n")
        g.write(f"# records={n_rec} records_kept={n_rec_keep} positions={len(pos_ok)} "
                f"positions_kept={len(kept)} records_failing_overall={n_total_fail}\n")
        g.write("group\tn_samples\tjudged\trecords_failing\trecords_failing_only_this_group\n")
        for grp in G:
            g.write(f"{grp}\t{sizes[grp]}\t{int(sizes[grp] >= a.min_group_n)}\t"
                    f"{fail_group[grp]}\t{only_group[grp]}\n")
    pct = 100 * len(kept) / max(len(pos_ok), 1)
    print(f"[admix-sites] {n_rec:,} records at {len(pos_ok):,} positions; kept {len(kept):,} "
          f"positions ({pct:.1f}%) with every judged group <= {a.max_group_miss:.0%} missing and "
          f"overall <= {a.max_total_miss:.0%}", file=sys.stderr)
    for grp, n in fail_group.most_common(6):
        print(f"[admix-sites]   {grp:16s} (n={sizes[grp]}) over the limit at {n:,} records "
              f"({only_group[grp]:,} of them for no other reason)", file=sys.stderr)
    small = [x for x in G if x not in judged]
    if small:
        print(f"[admix-sites]   not judged alone (n < {a.min_group_n}): "
              f"{', '.join(f'{x}({sizes[x]})' for x in small)}", file=sys.stderr)


if __name__ == "__main__":
    main()
