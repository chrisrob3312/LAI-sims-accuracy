#!/usr/bin/env python3
"""Turn a truth VCF into an in-silico genotyping-array VCF for one platform.

What a real array file has that a truth VCF does not, and what we inject so the
harmonization step has real work to do (every change is logged so it can be scored):

  1. site restriction   : only the platform's sites (chr:pos list, hg38)
  2. strand reporting   : a fraction of non-palindromic SNPs reported on the
                          minus strand (REF/ALT complemented) -- "--flip-frac"
  3. palindromic SNPs   : A/T and C/G sites are reported on a RANDOM strand,
                          which is undetectable without frequency information
  4. genotype error     : each called genotype replaced by another with prob e
  5. missingness        : each genotype set to ./. with prob m
  6. unphased           : arrays are unphased; '|' becomes '/'

Outputs: <out>.vcf.gz (bgzipped, unindexed) and <out>.changes.tsv with columns
  site  change  detail   (flip / palindrome_flip / palindrome_kept / error / missing counts per site)

Usage:
  insilico_array.py --vcf ALL.chr22.truth.vcf.gz --sites chr22.snps --out gsa.chr22 \
     --flip-frac 0.5 --error 0.002 --missing 0.01 --seed 22
"""
import argparse, gzip, random, sys

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALIN = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

p = argparse.ArgumentParser()
p.add_argument("--vcf", required=True)
p.add_argument("--sites", required=True, help="chr:pos per line (hg38); 'chr' prefix optional")
p.add_argument("--out", required=True)
p.add_argument("--flip-frac", type=float, default=0.5)
p.add_argument("--error", type=float, default=0.002)
p.add_argument("--missing", type=float, default=0.01)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
rng = random.Random(a.seed)

sites = set()
for line in open(a.sites):
    s = line.strip()
    if not s: continue
    c, pos = s.split(":")
    sites.add((c.replace("chr", ""), pos))

GTS = ["0/0", "0/1", "1/1"]
n_in = n_out = 0
with gzip.open(a.vcf, "rt") as f, gzip.open(a.out + ".vcf.gz", "wt") as g, open(a.out + ".changes.tsv", "w") as log:
    log.write("site\tchange\tdetail\n")
    for line in f:
        if line.startswith("##"):
            g.write(line); continue
        if line.startswith("#CHROM"):
            g.write('##INFO=<ID=INSILICO,Number=1,Type=String,Description="in-silico array modification: flip|palflip|palkept|none">\n')
            g.write(line); continue
        n_in += 1
        cols = line.rstrip("\n").split("\t")
        chrom, pos, vid, ref, alt = cols[0].replace("chr", ""), cols[1], cols[2], cols[3], cols[4]
        if (chrom, pos) not in sites: continue
        if len(ref) != 1 or len(alt) != 1 or ref not in COMP or alt not in COMP: continue  # arrays: biallelic SNPs only
        # strand reporting
        tag = "none"
        if (ref, alt) in PALIN:
            if rng.random() < 0.5:
                ref, alt, tag = COMP[ref], COMP[alt], "palflip"
            else:
                tag = "palkept"
        elif rng.random() < a.flip_frac:
            ref, alt, tag = COMP[ref], COMP[alt], "flip"
        # genotypes: unphase, error, missingness
        fmt = cols[8].split(":"); gti = fmt.index("GT")
        n_err = n_mis = 0
        newg = []
        for s in cols[9:]:
            parts = s.split(":"); gt = parts[gti].replace("|", "/")
            if gt in ("./.", "."):
                newg.append("./."); continue
            al = sorted(gt.split("/")); gt = "/".join(al)
            r = rng.random()
            if r < a.missing:
                gt = "./."; n_mis += 1
            elif r < a.missing + a.error:
                gt = rng.choice([x for x in GTS if x != gt]); n_err += 1
            newg.append(gt)
        info = ("INSILICO=" + tag) if cols[7] in (".", "") else (cols[7] + ";INSILICO=" + tag)
        g.write("\t".join([cols[0], pos, vid, ref, alt, cols[5], cols[6], info, "GT"] + newg) + "\n")
        n_out += 1
        if tag != "none" or n_err or n_mis:
            log.write(f"{chrom}:{pos}\t{tag}\terr={n_err};mis={n_mis}\n")
print(f"[insilico_array] {a.out}: {n_out} of {n_in} truth sites kept; flip_frac={a.flip_frac} error={a.error} missing={a.missing} seed={a.seed}", file=sys.stderr)
