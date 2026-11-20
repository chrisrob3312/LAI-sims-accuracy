#!/usr/bin/env python3
"""site-provenance.py -- where every panel record came from, per chromosome.

Run by 1c on the merged panel (before the soft-union filter) and again on the
soft-union, with the site list of each cohort file that went into the merge.
It answers, per cohort:
  * how many records its file had, and how many reached the panel;
  * how many are UNIQUE to that study (no other cohort file has the record);
  * for every panel record the cohort has NO call at, whether its file lacks
    the position entirely, or has it with alleles the merge could not join
    (REF/ALT swapped, strand flipped, both, a different ALT, a different REF,
    or identical alleles that still did not join);
and, over the whole panel, which COMBINATIONS of cohorts carry each record.

Inputs are `bcftools query -f '%POS\\t%REF\\t%ALT\\n'` lists (multi-ALT records
are split on commas). Nothing here reads genotypes.

  --cohort NAME=sites.tsv   one per merged file, in merge order
  --panel  sites.tsv        the panel records at this stage
  --stage  label            e.g. merged | softunion
  --out-cohorts / --out-patterns   TSVs; a summary goes to stdout
"""
import argparse
import sys
from collections import Counter, OrderedDict, defaultdict

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
CLASSES = ["absent (no record at this position)", "same alleles, not joined", "REF/ALT swapped",
           "strand-flipped", "flipped + swapped", "same REF, different ALT", "different REF"]


def read_sites(path):
    recs, by_pos = set(), defaultdict(list)
    n = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3: continue
            n += 1
            for alt in f[2].split(","):
                recs.add((f[0], f[1], alt))
                by_pos[f[0]].append((f[1], alt))
    return n, recs, by_pos


def classify(ref, alt, there):
    if not there:
        return CLASSES[0]
    rank = {c: i for i, c in enumerate(CLASSES)}
    best = None
    for r, a in there:
        if (r, a) == (ref, alt): c = CLASSES[1]
        elif (r, a) == (alt, ref): c = CLASSES[2]
        elif (r.translate(COMP), a.translate(COMP)) == (ref, alt): c = CLASSES[3]
        elif (r.translate(COMP), a.translate(COMP)) == (alt, ref): c = CLASSES[4]
        elif r == ref: c = CLASSES[5]
        else: c = CLASSES[6]
        if best is None or rank[c] < rank[best]: best = c
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", action="append", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--chrom", default="")
    ap.add_argument("--out-cohorts", required=True)
    ap.add_argument("--out-patterns", required=True)
    a = ap.parse_args()

    co = OrderedDict()
    for spec in a.cohort:
        name, path = spec.split("=", 1)
        co[name] = read_sites(path)
    n_panel, panel, _ = read_sites(a.panel)
    names = list(co)

    owners = Counter()                 # record -> how many cohort files carry it
    for nm in names:
        for r in co[nm][1]: owners[r] += 1

    patterns = Counter()
    lacks = {nm: Counter() for nm in names}
    for r in panel:
        have = [nm for nm in names if r in co[nm][1]]
        patterns["+".join(have) if have else "(none)"] += 1
        for nm in names:
            if r not in co[nm][1]:
                lacks[nm][classify(r[1], r[2], co[nm][2].get(r[0], []))] += 1

    tag = f"[provenance {a.chrom} {a.stage}]".replace("  ", " ")
    with open(a.out_cohorts, "w") as g:
        g.write("stage\tcohort\tfile_records\tfile_alleles\tin_panel\tunique_to_cohort\tunique_in_panel\t"
                "panel_records_lacking\t" + "\t".join(c.split(" (")[0].replace(" ", "_").replace(",", "")
                                                    for c in CLASSES) + "\n")
        print(f"{tag} panel: {n_panel:,} records ({len(panel):,} alleles)")
        for nm in names:
            n_file, recs, _ = co[nm]
            uniq = {r for r in recs if owners[r] == 1}
            row = [a.stage, nm, n_file, len(recs), len(recs & panel), len(uniq), len(uniq & panel),
                   sum(lacks[nm].values())] + [lacks[nm][c] for c in CLASSES]
            g.write("\t".join(str(x) for x in row) + "\n")
            lk = sum(lacks[nm].values())
            fix = lacks[nm][CLASSES[2]] + lacks[nm][CLASSES[3]] + lacks[nm][CLASSES[4]]
            print(f"{tag} {nm:8s} file {n_file:>9,} | in panel {len(recs & panel):>9,} | unique to {nm} "
                  f"{len(uniq):>9,} ({len(uniq & panel):,} in panel) | panel records it lacks {lk:>9,}: "
                  f"absent {lacks[nm][CLASSES[0]]:,}, swap/flip {fix:,}, other ALT {lacks[nm][CLASSES[5]]:,}, "
                  f"other REF {lacks[nm][CLASSES[6]]:,}, same-but-unjoined {lacks[nm][CLASSES[1]]:,}")
    with open(a.out_patterns, "w") as g:
        g.write("stage\tcohorts_carrying_record\tn_records\tpct\n")
        for p, n in patterns.most_common():
            g.write(f"{a.stage}\t{p}\t{n}\t{100 * n / max(len(panel), 1):.2f}\n")
    print(f"{tag} which cohorts carry each panel record (top 8):")
    for p, n in patterns.most_common(8):
        print(f"{tag}   {p:32s} {n:>9,}  {100 * n / max(len(panel), 1):5.1f}%")


if __name__ == "__main__":
    main()
