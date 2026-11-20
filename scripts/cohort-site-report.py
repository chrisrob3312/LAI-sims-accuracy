#!/usr/bin/env python3
"""
cohort-site-report.py -- per site, per cohort, over a published soft-union.

For every record of 1c's PRE-phasing soft-union (merged_chrN.softunion.bcf):
  src_<C>   where that record stood in cohort C's OWN file (the file 1c merged):
              P   present
              A   absent: no record at the position -- 1c filled it 0/0 for a
                  complete callset (HGDP1KG, MXB, SGDP, ALMARRI)
              S   REF/ALT swapped        F   strand-flipped    FS  both
              X   same REF, other ALT    R   other REF
              U   same alleles, yet not joined (should not happen)
  miss_<C>  fraction of C's genotypes missing (./.) at the record in the
            soft-union -- what SHAPEIT5 will impute for C there
  af_<C>    C's ALT allele frequency there (AC/AN), NA when AN = 0
  admix     1 if the record's position is on 1c's ADMIXTURE site list
plus a per-bin table (default 1 Mb): per cohort, records it lacked in its own
file (A), records it had with other alleles (S/F/FS/X/R/U), and its mean
missingness -- so where a cohort is absent or missing shows by region.

Reads `bcftools query` text only (no genotypes): --su is
  POS REF ALT [AN_<C> AC_<C>]... for C in --labels order
(fill-tags -S <sample TAB cohort> -t AN,AC on the soft-union), --cohort
C=<POS REF ALT list of C's file>. Run by scripts/cohort-site-report.sh.
"""
import argparse
import gzip
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
CODES = ["P", "A", "U", "S", "F", "FS", "X", "R"]
NAMES = {"P": "present", "A": "absent (filled 0/0)", "U": "same alleles, not joined",
         "S": "REF/ALT swapped", "F": "strand-flipped", "FS": "flipped + swapped",
         "X": "same REF, other ALT", "R": "other REF"}


def classify(ref, alt, there):
    """Best (closest) relation between the panel record and the cohort's records at its position."""
    if not there:
        return "A"
    order = {c: i for i, c in enumerate(CODES)}
    best = None
    for r, a in there:
        if (r, a) == (ref, alt): c = "U"
        elif (r, a) == (alt, ref): c = "S"
        elif (r.translate(COMP), a.translate(COMP)) == (ref, alt): c = "F"
        elif (r.translate(COMP), a.translate(COMP)) == (alt, ref): c = "FS"
        elif r == ref: c = "X"
        else: c = "R"
        if best is None or order[c] < order[best]:
            best = c
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--su", required=True)
    ap.add_argument("--labels", required=True, help="cohorts in the --su column order, space-separated")
    ap.add_argument("--n", action="append", default=[], help="C=number of C's samples in the soft-union")
    ap.add_argument("--cohort", action="append", default=[], help="C=POS REF ALT list of C's own file")
    ap.add_argument("--admix-bed", default="")
    ap.add_argument("--bin-mb", type=float, default=1.0)
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--out-bins", required=True)
    a = ap.parse_args()

    labels = a.labels.split()
    nsamp = {k: int(v) for k, v in (x.split("=", 1) for x in a.n)}
    files = dict(x.split("=", 1) for x in a.cohort)
    for c in labels:
        if c not in nsamp or nsamp[c] <= 0:
            sys.exit(f"ERROR: no --n for cohort {c}")

    # 1. the soft-union records (positions kept to filter the cohort lists)
    su = []
    with open(a.su) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 + 2 * len(labels):
                sys.exit(f"ERROR: {a.su}: {len(f)} columns, expected {3 + 2 * len(labels)} -- "
                         f"labels and query format disagree")
            su.append(f)
    pos_needed = {f[0] for f in su}

    # 2. each cohort's own records, only at soft-union positions
    have, by_pos, n_file = {}, {}, {}
    for c in labels:
        recs, bp, n = set(), defaultdict(list), 0
        if c in files:
            with open(files[c]) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 3:
                        continue
                    n += 1
                    if f[0] not in pos_needed:
                        continue
                    for alt in f[2].split(","):
                        recs.add((f[0], f[1], alt))
                        bp[f[0]].append((f[1], alt))
        have[c], by_pos[c], n_file[c] = recs, bp, (n if c in files else None)

    admix = set()
    if a.admix_bed:
        with open(a.admix_bed) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    admix.add(f[1])

    # 3. per record
    binw = int(a.bin_mb * 1_000_000)
    tally = {c: defaultdict(int) for c in labels}
    miss_sum = defaultdict(float); miss_gt0 = defaultdict(int); miss_gt15 = defaultdict(int)
    admix_miss = defaultdict(float); n_admix = 0
    bins = defaultdict(lambda: {"n": 0, "admix": 0, **{f"A_{c}": 0 for c in labels},
                                **{f"O_{c}": 0 for c in labels}, **{f"m_{c}": 0.0 for c in labels}})
    hdr = ["pos", "ref", "alt", "admix"] + [x for c in labels for x in (f"src_{c}", f"miss_{c}", f"af_{c}")]
    with gzip.open(a.out_sites, "wt") as g:
        g.write("\t".join(hdr) + "\n")
        for f in su:
            pos, ref, alt = f[0], f[1], f[2]
            ad = pos in admix
            n_admix += ad
            b = bins[int(pos) // binw]
            b["n"] += 1; b["admix"] += ad
            row = [pos, ref, alt, "1" if ad else "0"]
            for i, c in enumerate(labels):
                an_s, ac_s = f[3 + 2 * i], f[4 + 2 * i]
                an = int(an_s) if an_s not in (".", "") else 0
                ac = int(ac_s) if ac_s not in (".", "") else 0
                src = "P" if (pos, ref, alt) in have[c] else classify(ref, alt, by_pos[c].get(pos, []))
                miss = 1.0 - an / (2.0 * nsamp[c])
                af = f"{ac / an:.4f}" if an > 0 else "NA"
                tally[c][src] += 1
                miss_sum[c] += miss; miss_gt0[c] += miss > 0; miss_gt15[c] += miss > 0.15
                if ad: admix_miss[c] += miss
                if src == "A": b[f"A_{c}"] += 1
                elif src != "P": b[f"O_{c}"] += 1
                b[f"m_{c}"] += miss
                row += [src, f"{miss:.4f}", af]
            g.write("\t".join(row) + "\n")

    with open(a.out_bins, "w") as g:
        g.write("\t".join(["bin_start", "bin_end", "records", "admix_records"]
                          + [x for c in labels for x in (f"absent_{c}", f"other_alleles_{c}", f"mean_miss_{c}")]) + "\n")
        for k in sorted(bins):
            b = bins[k]
            g.write("\t".join([str(k * binw), str((k + 1) * binw), str(b["n"]), str(b["admix"])]
                              + [x for c in labels for x in (str(b[f"A_{c}"]), str(b[f"O_{c}"]),
                                                              f"{b[f'm_{c}'] / b['n']:.4f}")]) + "\n")

    # 4. summary
    n = len(su)
    tag = f"[site-report {a.chrom}]"
    print(f"{tag} soft-union records: {n:,}; on the ADMIXTURE site list: {n_admix:,}")
    for c in labels:
        t = tally[c]
        parts = ", ".join(f"{NAMES[k]} {t[k]:,} ({100 * t[k] / max(n, 1):.1f}%)" for k in CODES if t[k])
        nf = f"{n_file[c]:,} records in its file" if n_file[c] is not None else "no file given"
        print(f"{tag} {c:8s} n={nsamp[c]:<5d} {nf}; at soft-union records: {parts}")
        print(f"{tag} {c:8s} missing genotypes: mean {miss_sum[c] / max(n, 1):.4f}; records with any "
              f"{miss_gt0[c]:,}, with > 15% {miss_gt15[c]:,}; mean on ADMIXTURE records "
              f"{(admix_miss[c] / n_admix) if n_admix else 0:.4f}")
    print(f"{tag} where each cohort lacked the most records in its own file (top 5 {a.bin_mb:g}-Mb bins, "
          f"absent / records):")
    for c in labels:
        top = sorted(bins.items(), key=lambda kv: -kv[1][f"A_{c}"] / max(kv[1]["n"], 1))[:5]
        cells = "; ".join(f"{k * binw / 1e6:.0f}-{(k + 1) * binw / 1e6:.0f} Mb {v[f'A_{c}']}/{v['n']}" for k, v in top)
        print(f"{tag}   {c:8s} {cells}")
    print(f"{tag} per record: {a.out_sites}")
    print(f"{tag} per bin:    {a.out_bins}")


if __name__ == "__main__":
    main()
