#!/usr/bin/env python3
"""
LASI-DAD (Longitudinal Aging Study in India -- Harmonized Diagnostic Assessment
of Dementia) ingestion helper. Produces SAS anchor sample lists once you have
local access to the LASI-DAD manifest.

Paper: Li et al., HGG Adv 2026 (PMID 41656741, PMC12945573).
       "A reference panel for linkage disequilibrium and genotype imputation
        using whole-genome sequencing data from 2,680 participants across India"
n_samples: 2,680 hi-cov WGS Indian participants -- currently the largest
           public India-focused reference panel. Nationally representative
           (many caste and linguistic groups).

Access reality (not a simple wget like SGDP):
  * Individual-level genotypes: dbGaP phs003013 (LASI-DAD study).
      -> requires IRB approval + Data Access Request (DAR) + DUA. Turnaround
         4-8 weeks. Apply as PI-of-record at your institution.
  * Imputation reference: available on Michigan / TOPMed Imputation Server
      -> lets you impute YOUR samples against LASI-DAD, but does NOT give you
         individual haplotypes -- not usable to expand our LAI reference panel.
  * Aggregate LD blocks / summary stats: paper supplement + LASI-DAD website.

Workflow:
  1. (One-time) Apply for dbGaP phs003013 DAR through your institution.
  2. Once approved, download the LASI-DAD sample manifest and VCFs. Save the
     manifest TSV somewhere and pass its path via --manifest-tsv.
  3. Run this script to emit:
       reference_ids/lasi_dad_sas_anchors.txt   sample IDs
       reference_ids/lasi_dad_manifest.tsv       filtered provenance rows
  4. Wire the LASI-DAD VCFs into 1b_phase_prep_clm.sh (same softunion flow
     used for MXB) -- add another bcftools merge input and rebuild the panel.

Column names expected in the manifest (customizable via --column-map JSON):
    sample_id           default: "Sample_ID"
    sex                 default: "Sex"
    age                 default: "Age"
    caste               default: "Caste"     (optional -- for stratification)
    state               default: "State"     (optional -- for geographic bins)
    coverage_x          default: "Coverage_X" (drop <20X if desired via
                                               --min-coverage-x)

Usage:
    scripts/fetch-lasi-dad-anchors.py \
        --manifest-tsv /path/to/lasi_dad_sample_manifest.tsv \
        [--outdir reference_ids] \
        [--min-coverage-x 20] \
        [--column-map columns.json] \
        [--per-state-cap N]

If --column-map is provided, it should be a JSON object mapping our internal
field names to actual columns in your manifest, e.g. {"sample_id": "IID"}.
"""

import argparse
import csv
import json
import os
import sys
from collections import Counter

DEFAULT_COLS = {
    "sample_id":  "Sample_ID",
    "sex":        "Sex",
    "age":        "Age",
    "caste":      "Caste",
    "state":      "State",
    "coverage_x": "Coverage_X",
}


def parse_manifest(path, cols):
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest-tsv", required=True,
                    help="Local path to LASI-DAD sample manifest TSV "
                         "(obtained after dbGaP phs003013 DAR approval).")
    ap.add_argument("--outdir", default="reference_ids")
    ap.add_argument("--min-coverage-x", type=float, default=None,
                    help="Drop samples below this coverage threshold (e.g. 20).")
    ap.add_argument("--per-state-cap", type=int, default=None,
                    help="Optional cap per Indian state -- useful for balanced "
                         "geographic anchoring so no single region dominates.")
    ap.add_argument("--column-map", default=None,
                    help="JSON file with column-name overrides.")
    args = ap.parse_args()

    cols = dict(DEFAULT_COLS)
    if args.column_map:
        with open(args.column_map) as f:
            cols.update(json.load(f))

    os.makedirs(args.outdir, exist_ok=True)
    rows_kept = []
    n_seen = 0
    n_dropped_cov = 0
    per_state = Counter()

    for r in parse_manifest(args.manifest_tsv, cols):
        n_seen += 1
        sample = r.get(cols["sample_id"], "").strip()
        if not sample:
            continue
        if args.min_coverage_x is not None:
            try:
                cov = float(r.get(cols["coverage_x"], "0"))
            except ValueError:
                cov = 0.0
            if cov < args.min_coverage_x:
                n_dropped_cov += 1
                continue
        state = r.get(cols["state"], "").strip() or "unknown"
        if args.per_state_cap and per_state[state] >= args.per_state_cap:
            continue
        per_state[state] += 1
        rows_kept.append({
            "sample": sample,
            "sex": r.get(cols["sex"], "").strip(),
            "age": r.get(cols["age"], "").strip(),
            "caste": r.get(cols["caste"], "").strip(),
            "state": state,
            "coverage_x": r.get(cols["coverage_x"], "").strip(),
        })

    anchors_path = os.path.join(args.outdir, "lasi_dad_sas_anchors.txt")
    with open(anchors_path, "w") as g:
        g.write("\n".join(sorted({r["sample"] for r in rows_kept}))
                + ("\n" if rows_kept else ""))
    manifest_path = os.path.join(args.outdir, "lasi_dad_manifest.tsv")
    fields = ["sample", "sex", "age", "caste", "state", "coverage_x"]
    with open(manifest_path, "w") as g:
        g.write("\t".join(fields) + "\n")
        for r in sorted(rows_kept, key=lambda x: (x["state"], x["sample"])):
            g.write("\t".join(r.get(k, "") for k in fields) + "\n")

    print(f"[lasi-dad] parsed {n_seen} manifest rows", file=sys.stderr)
    if args.min_coverage_x is not None:
        print(f"[lasi-dad] dropped {n_dropped_cov} below "
              f"{args.min_coverage_x}x coverage", file=sys.stderr)
    print(f"[lasi-dad] wrote {len(rows_kept)} SAS anchors -> {anchors_path}",
          file=sys.stderr)
    top_states = per_state.most_common(6)
    print(f"[lasi-dad] top states: "
          f"{', '.join(f'{s}={n}' for s, n in top_states)}", file=sys.stderr)
    print(f"[lasi-dad] manifest -> {manifest_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
