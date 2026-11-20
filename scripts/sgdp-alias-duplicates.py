#!/usr/bin/env python3
"""
sgdp-alias-duplicates.py

Name the individuals who appear TWICE in the merged panel, from metadata,
before KING ever runs -- and close a hole in the no-1000G-donor rule.

WHAT THE ALIAS COLUMN SAYS. SGDP's manifest carries an `alias` for each
sample: the identifier that individual holds in whatever cohort sequenced them
first. On our 277-sample candidate list, 132 aliases are HGDP ids and 41 are
1000 Genomes ids (HG#####/NA#####). SGDP re-sequenced those people; they are
not new individuals.

TWO SEPARATE PROBLEMS FALL OUT, AND THEY NEED DIFFERENT FIXES.

1. DUPLICATION. Where the aliased individual is ALSO in our HGDP+1000G panel,
   the merged panel holds one person twice, under two ids. ADMIXTURE counts
   them twice; one copy can be a donor while the other is a reference or a
   test sample. 1e (KING) finds these from genotypes. This script finds the
   same thing from metadata, which makes it a CROSS-CHECK rather than a
   replacement: KING should find at least these pairs, and a KING run that
   misses a pair named here is a KING run to distrust.

2. SERVER-PANEL LEAKAGE, which KING does NOT solve. A 1000 Genomes individual
   re-sequenced by SGDP arrives as `S_Name-1.DG`, is recorded as cohort SGDP,
   and is therefore DONOR-ELIGIBLE under a rule that excludes `cohort == 1KG`.
   Their haplotypes are nevertheless inside TOPMed's and Michigan's reference
   panels under their NA/HG identifier. Using them as a simulation donor means
   imputing against their own haplotypes -- exactly the circularity the
   no-1000G-donor rule exists to prevent, entering through the SGDP door.

   De-duplication does not fix this. Even when the 1000 Genomes copy is NOT in
   our panel, so there is no duplicate pair to resolve, the SGDP copy is still
   that person, and that person is still in the server's panel.

OUTPUTS
  known_duplicate_pairs.tsv   every aliased sample, its alias, the cohort the
                              alias belongs to, and whether that alias is in
                              our panel (i.e. whether it is a real duplicate)
  sgdp_server_panel_ids.txt   SGDP ids that are 1000 Genomes people. Treat as
                              1000 Genomes for DONOR ELIGIBILITY regardless of
                              whether a duplicate exists.

Usage:
  python3 -I scripts/sgdp-alias-duplicates.py \\
      --manifest reference_ids/sgdp_manifest.tsv \\
      --panel    reference_ids/lai_ref_panel_samples.tsv \\
      --out      reference_ids
"""
import argparse
import csv
import os
import re
import sys

HGDP_RE = re.compile(r"^HGDP\d+$")
KG_RE = re.compile(r"^(HG|NA)\d{5}$")


def norm_hgdp(s):
    """HGDP00540 / HGDP_00540 / hgdp540 -> '540'. Formats differ by source."""
    d = re.sub(r"\D", "", s)
    return d.lstrip("0") or "0"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="sgdp_manifest.tsv")
    ap.add_argument("--panel", required=True,
                    help="a panel table with sample+cohort columns, or any "
                         "one-column list of the ids already in the panel")
    ap.add_argument("--out", default="reference_ids")
    args = ap.parse_args()

    # --- the ids already in the merged panel ------------------------------
    panel_ids, panel_hgdp_norm = set(), set()
    with open(args.panel, newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        if "\t" in first and "sample" in first.split("\t")[0].lower():
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                s = (r.get("sample") or "").strip()
                if not s:
                    continue
                panel_ids.add(s)
                if HGDP_RE.match(s):
                    panel_hgdp_norm.add(norm_hgdp(s))
        else:
            for line in fh:
                s = line.split("\t")[0].strip()
                if not s:
                    continue
                panel_ids.add(s)
                if HGDP_RE.match(s):
                    panel_hgdp_norm.add(norm_hgdp(s))
    print(f"[alias] panel ids: {len(panel_ids)} ({len(panel_hgdp_norm)} HGDP)")

    # --- the SGDP manifest -------------------------------------------------
    rows = []
    with open(args.manifest, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"sample", "alias"}
        missing = need - set(rd.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: {args.manifest} lacks {sorted(missing)}; "
                     f"found {rd.fieldnames}")
        for r in rd:
            rows.append(r)
    print(f"[alias] SGDP manifest: {len(rows)} samples")

    dup_rows, server_ids = [], []
    n_hgdp = n_kg = n_dup_hgdp = n_dup_kg = 0
    for r in rows:
        sgdp_id = (r.get("sample") or "").strip()
        alias = (r.get("alias") or "").strip()
        sp = (r.get("super_pop") or "").strip()
        pop = (r.get("population") or "").strip()
        if not sgdp_id or not alias:
            continue

        if HGDP_RE.match(alias):
            n_hgdp += 1
            in_panel = norm_hgdp(alias) in panel_hgdp_norm
            n_dup_hgdp += in_panel
            dup_rows.append((sgdp_id, alias, "HGDP", sp, pop,
                             "1" if in_panel else "0", ""))
        elif KG_RE.match(alias):
            n_kg += 1
            in_panel = alias in panel_ids
            n_dup_kg += in_panel
            # Donor-ineligible WHETHER OR NOT the 1000 Genomes copy is in our
            # panel: the server's panel contains that person either way.
            server_ids.append(sgdp_id)
            dup_rows.append((sgdp_id, alias, "1KG", sp, pop,
                             "1" if in_panel else "0", "server_panel"))

    os.makedirs(args.out, exist_ok=True)
    dup_f = os.path.join(args.out, "known_duplicate_pairs.tsv")
    with open(dup_f, "w") as g:
        g.write("sgdp_id\talias\talias_cohort\tsuper_pop\tpopulation\t"
                "alias_in_panel\tnote\n")
        for t in sorted(dup_rows):
            g.write("\t".join(t) + "\n")

    srv_f = os.path.join(args.out, "sgdp_server_panel_ids.txt")
    with open(srv_f, "w") as g:
        g.write("# SGDP ids that are 1000 Genomes individuals re-sequenced.\n")
        g.write("# Their haplotypes are inside the TOPMed and Michigan\n")
        g.write("# reference panels under their NA/HG identifier, so they are\n")
        g.write("# NOT donor-eligible even though their cohort reads SGDP.\n")
        for s in sorted(set(server_ids)):
            g.write(s + "\n")

    print(f"\n[alias] aliases that are HGDP ids       : {n_hgdp}")
    print(f"[alias]   of those, ALSO in our panel   : {n_dup_hgdp}  <- real duplicates")
    print(f"[alias] aliases that are 1000G ids      : {n_kg}")
    print(f"[alias]   of those, ALSO in our panel   : {n_dup_kg}  <- real duplicates")
    print(f"[alias] total named duplicate pairs     : {n_dup_hgdp + n_dup_kg}")
    print(f"[alias] SGDP ids inside a server panel  : {len(set(server_ids))}"
          f"  <- donor-ineligible regardless")
    print(f"\n[alias] {dup_f}")
    print(f"[alias] {srv_f}")
    print(f"""
[alias] These are found from METADATA, so they are a CROSS-CHECK on 1e, not a
[alias] replacement: KING should find at least the {n_dup_hgdp + n_dup_kg} pairs above
[alias] from genotypes, and a KING run that misses one named here is a KING
[alias] run to distrust. Conversely KING finds duplicates no alias records --
[alias] the same person across HGDP and 1000 Genomes, say -- which this
[alias] cannot see.""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
