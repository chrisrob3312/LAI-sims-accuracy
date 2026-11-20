#!/usr/bin/env python3
"""sgdp-raw-sample-map.py -- which SGDP per-sample VCF is which panel sample.

The vcf_variants release names each file, and the one sample inside it, by
the Illumina ID (SS6004478.annotated.nh2.variants.vcf.gz holds sample
SS6004478). The panel names SGDP samples as reference_ids/sgdp_manifest.tsv
does (`sample`, e.g. B_Australian-3.DG), and everything downstream -- the
cohort groups, 4a's labels, the known-duplicate pairs -- keys on that name.
The public SGDP metadata joins the two through SGDP_ID:
    Illumina_ID SS6004478 -> SGDP_ID B_Australian-3 -> manifest sample B_Australian-3.DG

Only manifest samples are mapped (the manifest IS the SGDP candidate list);
each must resolve to exactly one Illumina ID with a VCF and an index on disk.
Anything that does not is reported, never silently dropped.

Output (--out): illumina_id <TAB> panel_sample <TAB> vcf_path, manifest order.
Exit 1 if nothing maps, 2 on bad input.
"""
import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--suffix", default=".annotated.nh2.variants.vcf.gz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.metadata, encoding="utf-8", errors="replace", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        sys.exit(2)
    head = [h.lstrip("#").strip().strip('"') for h in rows[0]]
    try:
        i_ill, i_sg = head.index("Illumina_ID"), head.index("SGDP_ID")
    except ValueError:
        print(f"ERROR: metadata header lacks Illumina_ID / SGDP_ID: {head[:8]}", file=sys.stderr)
        sys.exit(2)
    by_sgdp = {}
    for r in rows[1:]:
        if len(r) <= max(i_ill, i_sg) or not r[i_sg].strip():
            continue
        by_sgdp.setdefault(r[i_sg].strip(), []).append(r[i_ill].strip())

    with open(a.manifest, newline="") as fh:
        man = list(csv.DictReader(fh, delimiter="\t"))
    if not man or "sample" not in man[0] or "sgdp_id" not in man[0]:
        print("ERROR: manifest needs `sample` and `sgdp_id` columns", file=sys.stderr)
        sys.exit(2)

    out, no_meta, multi, no_vcf = [], [], [], []
    for m in man:
        ids = by_sgdp.get(m["sgdp_id"].strip(), [])
        if not ids:
            no_meta.append(m["sgdp_id"])
            continue
        if len(ids) > 1:
            multi.append(f"{m['sgdp_id']}: {','.join(ids)}")
        ill = ids[0]
        vcf = os.path.join(a.vcf_dir, ill + a.suffix)
        if not (os.path.isfile(vcf) and (os.path.isfile(vcf + ".tbi") or os.path.isfile(vcf + ".csi"))):
            no_vcf.append(f"{m['sgdp_id']} ({ill})")
            continue
        out.append((ill, m["sample"].strip(), vcf))

    seen = {}
    for ill, s, _ in out:
        if ill in seen:
            print(f"ERROR: {ill} maps to two panel samples ({seen[ill]}, {s})", file=sys.stderr)
            sys.exit(2)
        seen[ill] = s
    with open(a.out, "w") as g:
        for row in out:
            g.write("\t".join(row) + "\n")
    print(f"[sgdp-map] manifest {len(man)} samples: {len(out)} mapped to a per-sample VCF; "
          f"{len(no_meta)} not in the metadata, {len(no_vcf)} without a VCF + index, "
          f"{len(multi)} with several Illumina IDs (first used)", file=sys.stderr)
    for label, items in (("not in the metadata", no_meta), ("no VCF + index on disk", no_vcf),
                         ("several Illumina IDs", multi)):
        for x in items[:20]:
            print(f"[sgdp-map]   {label}: {x}", file=sys.stderr)
    # Files in the release that no manifest sample uses: public SGDP genomes the
    # panel does not take (not in the candidate manifest), named with their
    # SGDP_ID so they can be looked up.
    used = {ill for ill, _, _ in out}
    by_ill = {ill: sg for sg, ills in by_sgdp.items() for ill in ills}
    unused = sorted(f[: -len(a.suffix)] for f in os.listdir(a.vcf_dir)
                    if f.endswith(a.suffix) and f[: -len(a.suffix)] not in used)
    print(f"[sgdp-map] {len(unused)} per-sample VCF(s) on disk not used by any manifest sample", file=sys.stderr)
    for ill in unused[:20]:
        print(f"[sgdp-map]   unused: {ill} ({by_ill.get(ill, 'not in the metadata')})", file=sys.stderr)
    if not out:
        sys.exit(1)


if __name__ == "__main__":
    main()
