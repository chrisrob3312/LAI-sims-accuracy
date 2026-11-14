#!/usr/bin/env python3
"""softunion-sites.py -- the per-super-pop MAF soft-union, in one pass.

A site is kept if, in ANY group with >= --min-n samples, AN > 0 and
MAF >= --maf, written as AC >= maf*AN and (AN-AC) >= maf*AN, exactly the rule
1c applied before with one `bcftools view -S` subset per group. Here every
group's AC/AN comes from one `bcftools +fill-tags -- -S <groups> -t AC,AN`,
so the panel is read once instead of once per group.

stdin   bcftools query: %CHROM %POS then AC_<g> AN_<g> per judged group, in
        the order --print-format writes them
--groups   sample<TAB>GROUP
--out      passing positions, CHROM<TAB>POS, in input order
--counts   TSV: group, samples, judged (yes/no), passing records
"""
import argparse
import sys
from collections import OrderedDict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--maf", type=float, default=0.005)
    ap.add_argument("--print-format", action="store_true")
    ap.add_argument("--chrom", default="")
    ap.add_argument("--out")
    ap.add_argument("--counts")
    a = ap.parse_args()

    sizes = OrderedDict()
    with open(a.groups) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[1]] = sizes.get(f[1], 0) + 1
    judged = [g for g in sizes if sizes[g] >= a.min_n]
    if a.print_format:
        print("%CHROM\\t%POS" + "".join(f"\\t%INFO/AC_{g}\\t%INFO/AN_{g}" for g in judged) + "\\n")
        return
    if not (a.out and a.counts):
        sys.exit("ERROR: --out and --counts are required unless --print-format")
    if not judged:
        sys.exit(f"ERROR: no group has >= {a.min_n} samples: {dict(sizes)}")

    passing = [0] * len(judged)
    n = kept = 0
    width = 2 + 2 * len(judged)
    m = a.maf
    with open(a.out, "w") as out:
        for line in sys.stdin:
            f = line.rstrip("\n").split("\t")
            if len(f) != width:
                sys.exit(f"ERROR: expected {width} columns, got {len(f)}: {line[:160]!r}")
            n += 1
            ok = False
            for k in range(len(judged)):
                acs, ans = f[2 + 2 * k], f[3 + 2 * k]
                if ans in (".", "0", "") or acs in (".", "") or "," in acs:
                    continue
                ac, an = int(acs), int(ans)
                if ac >= m * an and (an - ac) >= m * an:
                    passing[k] += 1
                    ok = True
            if ok:
                kept += 1
                out.write(f"{f[0]}\t{f[1]}\n")
    tag = f"[softunion {a.chrom}]" if a.chrom else "[softunion]"
    with open(a.counts, "w") as g:
        g.write("group\tsamples\tjudged\tpassing_records\n")
        for grp, ns in sizes.items():
            k = judged.index(grp) if grp in judged else None
            g.write(f"{grp}\t{ns}\t{'yes' if k is not None else 'no'}\t{passing[k] if k is not None else ''}\n")
            if k is None:
                print(f"{tag}   skip {grp}: {ns} samples < {a.min_n}")
            else:
                print(f"{tag}   {grp}: {ns} samples, {passing[k]} passing sites")
    print(f"{tag} {n:,} records in; {kept:,} kept (MAF >= {m:g} in at least one group)")


if __name__ == "__main__":
    main()
