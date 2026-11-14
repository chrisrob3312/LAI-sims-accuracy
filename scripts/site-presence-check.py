#!/usr/bin/env python3
"""site-presence-check.py -- when a cohort has NO call at a panel site, is the
site really missing from that cohort's file, or is it there with alleles the
merge could not join (REF/ALT swapped, strand flipped, a different ALT)?

The merge in 1c joins records on position AND alleles. If a cohort's prepped
file carries a site with its alleles swapped or on the other strand, or saw a
different ALT, the merge keeps two records, and on the panel's record that
cohort is entirely missing -- which looks exactly like the cohort never having
the site. This tells the two apart, by POSITION, without dropping anything.

--absent        POS<TAB>REF<TAB>ALT of panel records where the cohort has AN==0
                (bcftools view -S <cohort> | +fill-tags -t AN | query -i 'AN==0')
--cohort-sites  POS<TAB>REF<TAB>ALT of EVERY record in the cohort's prepped file
                (bcftools query -f '%POS\\t%REF\\t%ALT\\n' <cohort.bcf>)

Categories, first match wins:
  no record at this position         the cohort's file lacks the site entirely
  same alleles                       present and identical: the merge should have
                                     joined it (contig/normalisation problem)
  REF/ALT swapped                    fixable: swap alleles and genotypes
  strand-flipped                     fixable: complement both alleles
  flipped + swapped                  fixable: both
  same REF, different ALT            a different variant at the same position
  different REF                      the files disagree on the reference base
A/T and C/G sites are reported separately: for them a swap and a flip cannot be
told apart from the alleles alone.
"""
import argparse
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
ORDER = ["no record at this position", "same alleles", "REF/ALT swapped", "strand-flipped",
         "flipped + swapped", "same REF, different ALT", "different REF"]


def classify(ref, alt, recs):
    if not recs:
        return ORDER[0]
    rank = {c: i for i, c in enumerate(ORDER)}
    best = None
    for r, a in recs:
        if (r, a) == (ref, alt):
            c = "same alleles"
        elif (r, a) == (alt, ref):
            c = "REF/ALT swapped"
        elif (r.translate(COMP), a.translate(COMP)) == (ref, alt):
            c = "strand-flipped"
        elif (r.translate(COMP), a.translate(COMP)) == (alt, ref):
            c = "flipped + swapped"
        elif r == ref:
            c = "same REF, different ALT"
        else:
            c = "different REF"
        if best is None or rank[c] < rank[best]:
            best = c
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--absent", required=True)
    ap.add_argument("--cohort-sites", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--examples", type=int, default=3)
    a = ap.parse_args()

    cohort = defaultdict(list)
    n_coh = 0
    with open(a.cohort_sites) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3: continue
            n_coh += 1
            for alt in f[2].split(","):
                cohort[f[0]].append((f[1], alt))
    cats, pal, ex = Counter(), Counter(), defaultdict(list)
    n = 0
    with open(a.absent) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3: continue
            n += 1
            c = classify(f[1], f[2], cohort.get(f[0], []))
            cats[c] += 1
            if (f[1], f[2]) in PAL: pal[c] += 1
            if len(ex[c]) < a.examples:
                there = ",".join(f"{r}>{x}" for r, x in cohort.get(f[0], [])) or "-"
                ex[c].append(f"{f[0]} panel {f[1]}>{f[2]} | {a.label} {there}")
    print(f"[{a.label}] cohort file: {n_coh:,} records at {len(cohort):,} positions")
    print(f"[{a.label}] panel records with NO {a.label} call: {n:,}")
    for c in ORDER:
        if cats[c]:
            print(f"  {c:28s} {cats[c]:9,d}  {100 * cats[c] / n:5.1f}%"
                  f"{f'   (A/T or C/G: {pal[c]:,})' if pal[c] else ''}")
            for e in ex[c]:
                print(f"      e.g. {e}")
    fixable = cats["REF/ALT swapped"] + cats["strand-flipped"] + cats["flipped + swapped"]
    print(f"[{a.label}] recoverable by allele harmonisation: {fixable:,} ({100 * fixable / max(n, 1):.1f}%);"
          f" truly absent: {cats[ORDER[0]]:,} ({100 * cats[ORDER[0]] / max(n, 1):.1f}%)")


if __name__ == "__main__":
    main()
