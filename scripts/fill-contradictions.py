#!/usr/bin/env python3
"""fill-contradictions.py -- sites where filling "absent" as 0/0 contradicts the data.

1c fills a record that a cohort's file lacks as 0/0 for that cohort
(FILL_REF_COHORTS), on the grounds that the file is a complete variant
callset: no record means nobody in the cohort carries the allele. That is
right where the allele is truly absent from the cohort, and wrong where the
file dropped a real site for another reason -- the source's own site QC
(HWE, missingness, VQSR), a mask, a liftover reject. Filled there, the whole
cohort reads 0/0 at a site its own ancestry carries: a cohort-batch signal,
the thing the fill exists to remove.

The other cohorts tell the two cases apart. For a filled cohort C lacking a
record, and another cohort D with a group in the SAME super-pop:

    p = (1 - AF_D) ** AN_C

is the probability that C's AN_C haplotypes show no copy at all if the allele
were as common in C as it is in D. When p < --max-p (and D carries at least
--min-ac copies, both groups having >= --min-n samples), the absence is not
believable as "nobody carries it", so the POSITION is dropped from the panel
instead of filled. Frequency differences between populations inside a
super-pop make the rule conservative (it drops a little more than it must),
never permissive. Groups whose super-pop is NA are never compared.

stdin   bcftools query on the merged, filled panel after
        `bcftools +fill-tags -- -S <groups> -t AC,AN`: %CHROM %POS %REF %ALT,
        then AC_<g> AN_<g> for every group (--print-format writes the string)
--groups        sample<TAB>COHORT_SP (build-missingness-groups.py)
--filled        the filled cohorts, e.g. "HGDP1KG MXB ALMARRI"
--cohort-sites  LABEL=file, one per filled cohort: the records its ORIGINAL
                (pre-fill, split) file had, as POS<TAB>REF<TAB>ALT
--out           positions to drop (CHROM<TAB>POS); --summary: TSV
"""
import argparse
import sys
from collections import Counter, OrderedDict, defaultdict


def load_groups(path):
    sizes = OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[1]] = sizes.get(f[1], 0) + 1
    return sizes


def split_group(g):
    co, _, sp = g.rpartition("_")
    return co, sp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--print-format", action="store_true")
    ap.add_argument("--filled", default="")
    ap.add_argument("--cohort-sites", action="append", default=[])
    ap.add_argument("--min-n", type=int, default=20, help="samples in BOTH groups compared")
    ap.add_argument("--min-ac", type=int, default=5, help="copies the comparator group must carry")
    ap.add_argument("--max-p", type=float, default=1e-3)
    ap.add_argument("--chrom", default="")
    ap.add_argument("--out")
    ap.add_argument("--summary")
    a = ap.parse_args()

    sizes = load_groups(a.groups)
    G = list(sizes)
    if a.print_format:
        print("%CHROM\\t%POS\\t%REF\\t%ALT" + "".join(f"\\t%INFO/AC_{g}\\t%INFO/AN_{g}" for g in G) + "\\n")
        return
    if not (a.out and a.summary):
        sys.exit("ERROR: --out and --summary are required unless --print-format")
    filled = a.filled.split()
    had = {}
    for spec in a.cohort_sites:
        lab, path = spec.split("=", 1)
        s = set()
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    for alt in f[2].split(","):
                        s.add(f"{f[0]}\t{f[1]}\t{alt}")
        had[lab] = s
    missing_lists = [c for c in filled if c not in had]
    if missing_lists:
        sys.exit(f"ERROR: no --cohort-sites for filled cohort(s) {missing_lists}")

    # super-pop -> [(cohort, column index)] for groups big enough to judge
    by_sp = defaultdict(list)
    for i, g in enumerate(G):
        co, sp = split_group(g)
        if sp != "NA" and sizes[g] >= a.min_n:
            by_sp[sp].append((co, i))
    # what each filled cohort is compared on: (its column, [comparator columns])
    plan = {}
    for c in filled:
        plan[c] = [(sp, i, [(d, j) for d, j in members if d != c])
                   for sp, members in by_sp.items() for co, i in members if co == c]
        plan[c] = [p for p in plan[c] if p[2]]

    lacked, contra = Counter(), Counter()
    by_sp_count = defaultdict(Counter)
    examples = defaultdict(list)
    drop = OrderedDict()
    n = 0
    width = 4 + 2 * len(G)
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if len(f) != width:
            sys.exit(f"ERROR: expected {width} columns, got {len(f)}: {line[:160]!r}")
        n += 1
        rec = f"{f[1]}\t{f[2]}\t{f[3]}"
        ac = f[4::2]
        an = f[5::2]
        for c in filled:
            if rec in had[c]:
                continue
            lacked[c] += 1
            worst = None
            for sp, i, comps in plan[c]:
                an_c = 0 if an[i] in (".", "") else int(an[i])
                if an_c == 0:
                    continue
                for d, j in comps:
                    if ac[j] in (".", "") or "," in ac[j] or an[j] in (".", "0", ""):
                        continue
                    ac_d, an_d = int(ac[j]), int(an[j])
                    if ac_d < a.min_ac:
                        continue
                    p = (1.0 - ac_d / an_d) ** an_c
                    if p < a.max_p and (worst is None or p < worst[0]):
                        worst = (p, sp, d, ac_d / an_d, ac_d, an_c)
            if worst:
                contra[c] += 1
                by_sp_count[c][worst[1]] += 1
                drop[(f[0], f[1])] = True
                if len(examples[c]) < 3:
                    p, sp, d, af, ac_d, an_c = worst
                    examples[c].append(f"{f[0]}:{f[1]} {f[2]}>{f[3]}: not in {c}'s file; {d} {sp} carries "
                                       f"{ac_d} copies (AF {af:.3f}); P(0 of {an_c} {c} {sp} haplotypes) = {p:.1e}")

    with open(a.out, "w") as g:
        for chrom, pos in drop:
            g.write(f"{chrom}\t{pos}\n")
    with open(a.summary, "w") as g:
        g.write(f"# min_n={a.min_n} min_ac={a.min_ac} max_p={a.max_p} records={n} positions_dropped={len(drop)}\n")
        g.write("cohort\tcompared_on\trecords_absent_from_file\tcontradicted\tcontradicted_by_superpop\n")
        for c in filled:
            on = ",".join(sp for sp, _, _ in plan[c]) or "-"
            bysp = ",".join(f"{k}:{v}" for k, v in sorted(by_sp_count[c].items())) or "-"
            g.write(f"{c}\t{on}\t{lacked[c]}\t{contra[c]}\t{bysp}\n")
    tag = f"[fill-check {a.chrom}]" if a.chrom else "[fill-check]"
    print(f"{tag} {n:,} records; {len(drop):,} positions dropped: absent from a filled cohort's file "
          f"although another cohort of the same super-pop makes that implausible "
          f"(>= {a.min_ac} copies, P(no copy) < {a.max_p:g}, groups of >= {a.min_n})")
    for c in filled:
        on = ",".join(sp for sp, _, _ in plan[c]) or "nothing (no comparable group)"
        print(f"{tag}   {c:8s} absent at {lacked[c]:>9,} records; contradicted at {contra[c]:>7,}; compared on {on}")
        for e in examples[c]:
            print(f"{tag}       e.g. {e}")


if __name__ == "__main__":
    main()
