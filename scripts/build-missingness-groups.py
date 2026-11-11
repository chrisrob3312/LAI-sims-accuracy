#!/usr/bin/env python3
"""build-missingness-groups.py -- sample -> COHORT_SUPERPOP group, for the
per-group missingness rule that decides which sites ADMIXTURE may use.

Cohort comes from MEMBERSHIP, never from the shape of an ID: 1c passes the
sample list of each file it merged (--cohort NAME=FILE). Samples in none of
them get --default-cohort (the HGDP+1KG source). Super-pop: SGDP manifest for
SGDP, MEN for ALMARRI, AMR for MXB, otherwise sample_groups.tsv (CSA->SAS,
MID->MEN); a sample with none of those is NA and still gets a group
(e.g. HGDP1KG_NA) so it is counted, never silently dropped.

Output (stdout): sample <TAB> GROUP, in the order of --samples. GROUP is a
valid INFO-tag suffix, so `bcftools +fill-tags -- -S <this> -t AN` writes one
AN_<GROUP> per group. Group sizes go to stderr.
"""
import argparse
import csv
import re
import sys
from collections import Counter


def read_ids(path):
    with open(path) as fh:
        return [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="every sample, in BCF order (bcftools query -l)")
    ap.add_argument("--cohort", action="append", default=[], help="NAME=FILE (repeatable)")
    ap.add_argument("--default-cohort", default="HGDP1KG")
    ap.add_argument("--groups", required=True, help="sample_groups.tsv (sample<TAB>superpop)")
    ap.add_argument("--sgdp-manifest", help="reference_ids/sgdp_manifest.tsv")
    a = ap.parse_args()

    member = {}
    for spec in a.cohort:
        name, path = spec.split("=", 1)
        name = re.sub(r"[^A-Za-z0-9]", "", name.upper())
        for s in read_ids(path):
            if s in member and member[s] != name:
                sys.exit(f"ERROR: {s} is in two cohort lists ({member[s]}, {name})")
            member[s] = name
    groups = {}
    with open(a.groups) as fh:
        for l in fh:
            f = l.split()
            if len(f) >= 2:
                g = f[1].upper()
                groups[f[0]] = {"CSA": "SAS", "MID": "MEN"}.get(g, g)
    sgdp = {}
    if a.sgdp_manifest:
        with open(a.sgdp_manifest) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sgdp[r["sample"]] = r["super_pop"]

    out, sizes = [], Counter()
    for s in read_ids(a.samples):
        co = member.get(s, a.default_cohort)
        if co == "SGDP":
            sp = sgdp.get(s, "NA")
        elif co == "ALMARRI":
            sp = "MEN"
        elif co == "MXB":
            sp = "AMR"
        else:
            sp = groups.get(s, "NA")
        g = f"{co}_{sp}"
        out.append(f"{s}\t{g}")
        sizes[g] += 1
    sys.stdout.write("\n".join(out) + "\n")
    print(f"[miss-groups] {len(out)} samples in {len(sizes)} groups:", file=sys.stderr)
    for g, n in sorted(sizes.items()):
        print(f"[miss-groups]   {g:16s} {n:5d}", file=sys.stderr)


if __name__ == "__main__":
    main()
