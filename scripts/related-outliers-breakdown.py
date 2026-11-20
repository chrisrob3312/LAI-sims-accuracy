#!/usr/bin/env python3
"""
related-outliers-breakdown.py -- READ-ONLY. What the lab's related_outliers.txt
holds, read against the gnomAD HGDP+1KG metadata, and which of the curated
reference / donor samples are on it.

Why: the July and November panels dropped every ID on that list in 1c (693
HGDP+1KG genomes). The rebuild keeps them (OUTLIERS=none). The methods need to
say what was dropped before -- how many were relatives, how many PCA outliers,
from which populations -- and the curated lists need to know which of their
samples come back, and who they are related to (2b refuses a donor/reference
pair at KING >= 0.177).

Prints, all raw:
  1. the inputs (paths, line counts, first lines of the list)
  2. metadata columns whose names look like relatedness / outlier flags
  3. flag combinations, on the list vs off it
  4. populations losing the most samples to the list (n on list / n in metadata)
  5. every curated reference/donor sample on the list: lists, population,
     flags, and recorded relationships if the metadata has such a column

Run inside the shapeit5 env, from the repo:
  python3 -I scripts/related-outliers-breakdown.py
  python3 -I scripts/related-outliers-breakdown.py --meta <tsv> --list <txt>
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

SR = "/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall"
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--meta", default=f"{SR}/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv")
ap.add_argument("--list", default=f"{SR}/processed_data/sample_map_files/related_outliers.txt")
ap.add_argument("--refs", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reference_ids"))
ap.add_argument("--top", type=int, default=20)
args = ap.parse_args()

csv.field_size_limit(sys.maxsize)

print("== 1. inputs")
for p in (args.meta, args.list):
    if not os.path.exists(p):
        sys.exit(f"ERROR: {p} does not exist")
    with open(p) as fh:
        n = sum(1 for _ in fh)
    print(f"  {p}: {n} lines")
with open(args.list) as fh:
    head = [next(fh, "").rstrip("\n") for _ in range(3)]
print("  first lines of the list (raw):")
for h in head:
    print(f"    {h!r}")
on_list = {}
with open(args.list) as fh:
    for line in fh:
        w = line.split()
        if len(w) >= 2:
            on_list[w[1]] = w[0]
        elif len(w) == 1:
            on_list[w[0]] = "?"
print(f"  IDs on the list: {len(on_list)}")

with open(args.meta, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    rows = list(rd)
col = {name: i for i, name in enumerate(header)}

def pick(cands):
    for c in cands:
        if c in col:
            return c
    return None

sid_c = pick(["project_meta.sample_id", "s", "sample_id", "sample"])
pop_c = pick(["hgdp_tgp_meta.Population", "population", "pop"])
reg_c = pick(["hgdp_tgp_meta.Genetic.region", "hgdp_tgp_meta.Genetic_region", "super_pop", "region"])
if sid_c is None:
    sys.exit(f"ERROR: no sample-ID column in {args.meta}; header starts {header[:8]}")

flag_re = re.compile(r"(related|outlier)$", re.I)
rel_re = re.compile(r"relationship", re.I)
flag_cols = [h for h in header if flag_re.search(h)]
rel_cols = [h for h in header if rel_re.search(h)]
print("\n== 2. metadata columns used (raw names)")
print(f"  sample: {sid_c}   population: {pop_c}   region: {reg_c}")
print(f"  flag columns: {flag_cols if flag_cols else 'NONE FOUND -- names ending in related/outlier'}")
print(f"  relationship columns: {rel_cols if rel_cols else 'none'}")

meta = {}
for r in rows:
    if len(r) < len(header):
        r = r + [""] * (len(header) - len(r))
    meta[r[col[sid_c]]] = r

print("\n== 3. flag combinations (count  on/off-list  " + "  ".join(flag_cols) + ")")
combos = collections.Counter()
for s, r in meta.items():
    key = ("ON " if s in on_list else "off",) + tuple(r[col[c]] for c in flag_cols)
    combos[key] += 1
for key, n in sorted(combos.items(), key=lambda kv: (kv[0][0], -kv[1])):
    print(f"  {n:5d}  " + "  ".join(key))
missing = sorted(s for s in on_list if s not in meta)
print(f"  on the list but not in the metadata: {len(missing)}" + (f" (first: {missing[:5]})" if missing else ""))

print(f"\n== 4. populations losing the most to the list (top {args.top}): on list / in metadata")
tot = collections.Counter(); lost = collections.Counter()
for s, r in meta.items():
    p = r[col[pop_c]] if pop_c else "?"
    rg = r[col[reg_c]] if reg_c else "?"
    tot[(rg, p)] += 1
    if s in on_list:
        lost[(rg, p)] += 1
for (rg, p), n in lost.most_common(args.top):
    print(f"  {rg:6s} {p:22s} {n:4d} / {tot[(rg, p)]:4d}")
by_region = collections.Counter(); by_region_tot = collections.Counter()
for (rg, p), n in tot.items():
    by_region_tot[rg] += n; by_region[rg] += lost[(rg, p)]
print("  by region: " + "  ".join(f"{rg} {by_region[rg]}/{by_region_tot[rg]}" for rg in sorted(by_region_tot)))

print("\n== 5. curated reference / donor samples on the list")
listed = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(args.refs, "*_rfmix.txt")) + glob.glob(os.path.join(args.refs, "*_simulation.txt"))):
    name = os.path.basename(f)[:-len(".txt")]
    with open(f) as fh:
        for line in fh:
            w = line.split()
            if w:
                listed[w[0]].append(name)
hits = sorted(s for s in listed if s in on_list)
print(f"  {len(hits)} of {len(listed)} listed samples are on the list")
print("  sample  lists  population  " + "  ".join(flag_cols) + ("  " + "  ".join(rel_cols) if rel_cols else ""))
for s in hits:
    r = meta.get(s)
    pop = r[col[pop_c]] if (r and pop_c) else "?"
    flags = "  ".join(r[col[c]] for c in flag_cols) if r else "(not in metadata)"
    rel = ("  " + "  ".join(r[col[c]] for c in rel_cols)) if (r and rel_cols) else ""
    print(f"  {s}  {','.join(listed[s])}  {pop}  {flags}{rel}")

# 6. Recorded relationships among the curated lists, by side. A donor list is
# *_simulation, a reference list *_rfmix. Once relatives are kept, a pair with
# one member on each side is what 2b refuses (at KING >= 0.177) -- this shows
# it from the metadata, before 1e has computed KING on the panel.
print("\n== 6. recorded relationships between curated samples, by side")
if not rel_cols:
    print("  the metadata has no relationship column: 1e's KING table is the first look"
          " (verify-rebuild.sh section 7)")
    sys.exit(0)
typed = re.compile(r'"?([A-Za-z][A-Za-z /_-]*?)"?\s*[:=]\s*[\[{(]([^\]})]*)[\]})]')
def partners(text):
    out = []
    for kind, body in typed.findall(text):
        for t in re.findall(r"[A-Za-z0-9_.-]+", body):
            if t in meta:
                out.append((t, kind.strip()))
    if not out:   # untyped: any known sample ID in the text
        out = [(t, "?") for t in re.findall(r"[A-Za-z0-9_.-]+", text) if t in meta]
    return out
def side(lists):
    d = any(x.endswith("_simulation") for x in lists)
    r = any(x.endswith("_rfmix") for x in lists)
    return "both" if (d and r) else "donor" if d else "reference" if r else "none"
# The column may name only the TYPE of relationship ("[siblings]") and not the
# partner, in which case no pair can be formed from it -- say so rather than
# report "0 pairs", which would read as "no relatives across the split".
EMPTY = {"", "NA", "[]", "{}", "None", "nan"}
n_typed = sum(1 for s in listed if s in meta for c in rel_cols if meta[s][col[c]].strip() not in EMPTY)
n_named = sum(1 for s in listed if s in meta for c in rel_cols
              if any(t != s for t, _ in partners(meta[s][col[c]])))
print(f"  curated samples with a recorded relationship: {n_typed}; naming the partner: {n_named}")
if n_typed and not n_named:
    print("  The column records the TYPE of relationship only (e.g. '[siblings]'), not who with,")
    print("  so pairs -- and which side each member is on -- cannot be read from the metadata.")
    print("  First look: 1e's KING table, verify-rebuild.sh section 7, as soon as 1e has run.")
    sys.exit(0)
seen = set(); n_split = 0; n_pairs = 0
for s in sorted(listed):
    r = meta.get(s)
    if not r:
        continue
    for c in rel_cols:
        for t, kind in partners(r[col[c]]):
            if t == s or t not in listed:
                continue
            key = tuple(sorted((s, t)))
            if key in seen:
                continue
            seen.add(key); n_pairs += 1
            ss, ts = side(listed[s]), side(listed[t])
            split = {ss, ts} == {"donor", "reference"}
            n_split += split
            print(f"  {'SPLIT' if split else 'same '}  {s} ({','.join(listed[s])})  --{kind}--  "
                  f"{t} ({','.join(listed[t])})")
print(f"  {n_pairs} recorded pair(s) among curated samples; {n_split} split across donor/reference"
      " (2b refuses those at KING >= 0.177 -- first degree)")
