#!/usr/bin/env python3
"""
flare-to-lat3.py -- FLARE's local-ancestry calls in the layout the accuracy
scorer reads, so RFMix and FLARE are scored by the same code.

accuracy_v2.R reads, per (track, panel, chromosome), from its WORKDIR:
  <track>.<panel>.gen<G>_chr<C>.Lat3     one row per site: position, snp id,
                                         then one ancestry code per haplotype
  <track>.<panel>.gen<G>_chr<C>.classes  one integer per haplotype; 0 = the
                                         admixed haplotypes it scores
  <track>.<panel>.gen<G>_chr<C>.done     the marker that the combo finished
with ancestry codes 0 = NAT, 1 = EUR, 2 = AFR. 3b_flare_wgs-jointcall_clm.sh
writes <track>.<panel>.gen<G>_chr<C>.anc.vcf.gz: FORMAT/AN1 and AN2 per
simulated individual per site, with the codes named in an ##ANCESTRY header.

The codes are mapped BY NAME from that header (never assumed), so a FLARE run
that numbered the ancestries differently is still scored correctly; a name
other than NAT, EUR or AFR, or a missing header, stops the conversion.
FLARE's target VCF holds only the simulated individuals, in 2b's .sample
order, so every haplotype is admixed (.classes all 0) and the haplotype order
is the truth (.hanc) order: sample 1 hap 1, sample 1 hap 2, ...

A combo is converted again only when its .anc.vcf.gz is newer than its .done.

  python3 -I scripts/flare-to-lat3.py --flare-dir <03_flare/Brasa/gen12> \\
      --out-dir <03_flare/Brasa/gen12/lat3> --gen 12 --jobs 4
"""
import argparse
import glob
import os
import re
import subprocess
import sys
from multiprocessing import Pool

SCORER = {"NAT": "0", "EUR": "1", "AFR": "2"}
NAME_RE = re.compile(r"^(?P<track>[^.]+)\.(?P<panel>[^.]+)\.gen(?P<gen>\d+)_chr(?P<chr>\d+)\.anc\.vcf\.gz$")


def ancestry_map(vcf):
    hdr = subprocess.run(["bcftools", "view", "-h", vcf], check=True, capture_output=True, text=True).stdout
    line = next((l for l in hdr.splitlines() if l.startswith("##ANCESTRY=")), None)
    if line is None:
        raise ValueError(f"{vcf}: no ##ANCESTRY header")
    body = line[len("##ANCESTRY="):].strip().strip("<>")
    pairs = dict(p.split("=", 1) for p in body.split(",") if "=" in p)
    unknown = sorted(set(pairs) - set(SCORER))
    if unknown or set(pairs) != set(SCORER):
        raise ValueError(f"{vcf}: ancestries {sorted(pairs)} are not exactly NAT, EUR, AFR ({line})")
    codes = {v: SCORER[k] for k, v in pairs.items()}
    if any(len(c) != 1 for c in codes):
        raise ValueError(f"{vcf}: multi-digit ancestry codes in {line}")
    return str.maketrans({c: s for c, s in codes.items()}), line


def convert(job):
    vcf, prefix = job
    done = prefix + ".done"
    if os.path.exists(done) and os.path.getmtime(done) >= os.path.getmtime(vcf):
        return (vcf, "current", 0, 0)
    for ext in (".done", ".Lat3", ".classes"):
        if os.path.exists(prefix + ext):
            os.remove(prefix + ext)
    table, hdr = ancestry_map(vcf)
    samples = subprocess.run(["bcftools", "query", "-l", vcf], check=True, capture_output=True,
                             text=True).stdout.split()
    nhap = 2 * len(samples)
    if nhap == 0:
        raise ValueError(f"{vcf}: no samples")
    n = 0
    q = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS[\t%AN1\t%AN2]\n", vcf],
                         stdout=subprocess.PIPE, text=True)
    with open(prefix + ".Lat3.partial", "w") as g:
        for line in q.stdout:
            chrom, pos, rest = line.rstrip("\n").split("\t", 2)
            calls = rest.translate(table).split("\t")
            if len(calls) != nhap or any(c not in ("0", "1", "2") for c in calls):
                q.kill()
                raise ValueError(f"{vcf}: site {chrom}:{pos} has {len(calls)} calls or a code outside "
                                 f"NAT/EUR/AFR (expected {nhap})")
            g.write(f"{pos} {chrom}:{pos} " + " ".join(calls) + "\n")
            n += 1
    if q.wait() != 0:
        raise ValueError(f"{vcf}: bcftools query failed")
    if n == 0:
        raise ValueError(f"{vcf}: no sites")
    with open(prefix + ".classes.partial", "w") as g:
        g.write(" ".join(["0"] * nhap) + "\n")
    os.replace(prefix + ".Lat3.partial", prefix + ".Lat3")
    os.replace(prefix + ".classes.partial", prefix + ".classes")
    with open(prefix + ".done", "w") as g:
        g.write(f"from {vcf}\n{hdr}\n{len(samples)} individuals, {n} sites\n")
    return (vcf, "converted", len(samples), n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flare-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gen", default="12")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = []
    for vcf in sorted(glob.glob(os.path.join(a.flare_dir, "*.anc.vcf.gz"))):
        m = NAME_RE.match(os.path.basename(vcf))
        if not m or m["gen"] != str(a.gen):
            continue
        jobs.append((vcf, os.path.join(a.out_dir, os.path.basename(vcf)[:-len(".anc.vcf.gz")])))
    if not jobs:
        sys.exit(f"ERROR: no <track>.<panel>.gen{a.gen}_chr<N>.anc.vcf.gz in {a.flare_dir}")
    failed = 0
    with Pool(max(1, a.jobs)) as pool:
        for res in pool.imap_unordered(_safe, jobs):
            if isinstance(res, str):
                print(f"ERROR: {res}", file=sys.stderr); failed += 1
            elif res[1] == "converted":
                print(f"  {os.path.basename(res[0])}: {res[2]} individuals, {res[3]} sites")
    n_cur = sum(1 for _, p in jobs if os.path.exists(p + ".done"))
    print(f"FLARE -> scorer layout: {n_cur} of {len(jobs)} combos ready in {a.out_dir}; {failed} failed")
    sys.exit(1 if failed else 0)


def _safe(job):
    try:
        return convert(job)
    except Exception as e:          # report every failure, keep converting the rest
        msg = str(e)
        return msg if job[0] in msg else f"{job[0]}: {msg}"


if __name__ == "__main__":
    main()
