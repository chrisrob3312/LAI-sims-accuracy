#!/usr/bin/env python3
"""
split-homog-donors-reference.py

Split the homogeneous panel into a REFERENCE half and a DONOR half, per
super-pop, by a size-dependent rule rather than one fixed ratio.

WHY NOT A FIXED RATIO. "2/3 donors, 1/3 reference" is right for a large pool
and wrong for a small one: at a pool of 40 it leaves a 13-individual reference
panel, which is 26 haplotypes. The two sides do not fail the same way:

  Too few DONORS    simulated individuals reuse haplotypes and become cryptic
                    relatives. The intake's sim-kinship check catches this and
                    the fix is to simulate fewer individuals.
                    -> detectable, and recoverable.

  Too few REFERENCE local ancestry genuinely degrades, and that degradation is
                    confounded with the ancestry effect being measured. The
                    result reads as "local ancestry is poor for this group"
                    when the truth is "this group's reference panel had 13
                    people".
                    -> invisible in the output, and not recoverable.

So the reference side is allocated FIRST, to a target, and the donors take the
remainder. For a large pool that still yields 2/3 or more donors; for a small
one it protects the panel that cannot be repaired afterwards.

    reference = clamp(round(REF_FRACTION * N), REF_FLOOR, REF_TARGET)
    donors    = N - reference

With the defaults (0.375, floor 15, target 50):

    N = 400  ->  50 ref / 350 donors   (88% donors)
    N = 200  ->  50 ref / 150 donors   (75% donors)
    N = 100  ->  38 ref /  62 donors   (62% donors)
    N =  40  ->  15 ref /  25 donors   (62% donors)
    N =  30  ->  15 ref /  15 donors   (50% donors)
    N <  30  ->  refuses to split; see below

REF_TARGET = 50 matches the curated EUR/AFR/EAS/SAS panels, which are 50 each.

A pool below REF_FLOOR + DONOR_FLOOR cannot supply both roles. The script does
not quietly halve it: it reports the group as UNDERSUPPLIED and writes nothing
for it, because the useful options at that point are decisions, not defaults --
use the group as a reference/control only and simulate no model containing it,
or accept a thin panel and say so in the methods.

THE SPLIT IS STRATIFIED BY POPULATION, not random over the super-pop. An
unstratified split lets one population land mostly on the donor side and
another mostly on the reference side, so the benchmark measures
cross-population transfer while reporting it as within-super-pop accuracy.
Both sides therefore get a proportional share of every population, following
split_mxb_samples.R, and the allocation inside a stratum is deterministic
(lexicographic sample ID) so it reproduces without a seeded RNG.

SUPER-POPS WITH A SOURCE-BASED SPLIT ARE NOT SPLIT HERE. For EUR, EAS and AFR
the reference is the 1000 Genomes slice and donors come from everything else,
so disjointness is guaranteed by source and no sample has to be withheld.
Pass them with --source-split to get those two lists instead.

Usage:
  python3 -I scripts/split-homog-donors-reference.py \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k7.tsv \\
      --out reference_ids --suffix b2 \\
      --source-split EUR,EAS,AFR

  # see the allocation without writing anything
  python3 -I scripts/split-homog-donors-reference.py --panel ... --dry-run
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

REF_TARGET_DEFAULT = 65      # cap on the reference panel
REF_FLOOR_DEFAULT = 15
DONOR_FLOOR_DEFAULT = 15
REF_FRACTION_DEFAULT = 0.5   # even split until the cap bites


def plan_sizes(n, ref_fraction, ref_floor, ref_target, donor_floor):
    """-> (n_reference, n_donors, status). status '' when the split is usable.

    Even split, capped. With the defaults (0.5, floor 15, cap 65):

        N = 400 -> 65 / 335   N = 100 -> 50 / 50
        N = 200 -> 65 / 135   N =  40 -> 20 / 20
        N = 130 -> 65 /  65   N =  30 -> 15 / 15

    so a small pool is halved and a large one gives the reference its cap and
    everything else to donors -- which lands above 2/3 donors for any pool
    over ~195 without that being a separate rule.
    """
    if n < ref_floor + donor_floor:
        return 0, 0, "UNDERSUPPLIED"
    ref = int(round(ref_fraction * n))
    ref = max(ref_floor, min(ref_target, ref))
    don = n - ref
    if don < donor_floor:
        don = donor_floor
        ref = n - don
    return ref, don, ""


def allocate_clustered(by_pop, n_ref, cluster_of):
    """Stratified allocation that never splits a relative cluster.

    Members are collapsed to their cluster representative, the clusters are
    allocated by the same largest-remainder rule, and the result is expanded.
    A cluster therefore lands entirely on one side, which is what the intake's
    cross-role rule requires -- and the counts are approximate rather than
    exact, because a cluster of 3 cannot be split 1/2. That is the trade: an
    off-by-two reference panel against a cohort the intake refuses.
    """
    rep_pop, rep_members = {}, defaultdict(list)
    for pop, members in by_pop.items():
        for m in members:
            r = cluster_of(m)
            rep_pop[r] = pop
            rep_members[r].append(m)
    by_pop_rep = defaultdict(list)
    for r, pop in rep_pop.items():
        by_pop_rep[pop].append(r)
    # Weight each representative by its cluster size so the largest-remainder
    # pass is working in units of SAMPLES, not clusters.
    target_reps = {}
    for pop, reps in by_pop_rep.items():
        reps.sort()
    total = sum(len(v) for v in by_pop.values())
    if total <= 0:
        return [], []
    want_frac = n_ref / float(total)
    ref_ids, don_ids = [], []
    for pop in sorted(by_pop_rep):
        reps = sorted(by_pop_rep[pop], key=lambda r: (-len(rep_members[r]), r))
        quota = want_frac * sum(len(rep_members[r]) for r in reps)
        taken = 0
        for r in reps:
            if taken + len(rep_members[r]) <= quota + 1e-9:
                ref_ids.extend(sorted(rep_members[r])); taken += len(rep_members[r])
            else:
                don_ids.extend(sorted(rep_members[r]))
    return sorted(ref_ids), sorted(don_ids)


def allocate_stratified(by_pop, n_ref):
    """Largest-remainder allocation of n_ref reference slots across populations.

    Proportional by population size, so neither side over-represents a
    population. Deterministic throughout: populations are ordered by name and
    samples within a population by ID, so the same panel always yields the
    same split with no seeded RNG.
    """
    total = sum(len(v) for v in by_pop.values())
    pops = sorted(by_pop)
    exact = {p: n_ref * len(by_pop[p]) / total for p in pops}
    base = {p: int(exact[p]) for p in pops}
    # A population must not be wiped off the reference side by flooring alone.
    for p in pops:
        if base[p] == 0 and len(by_pop[p]) > 0:
            base[p] = 1
    over = sum(base.values()) - n_ref
    if over > 0:
        # Took too many: trim from the populations with the largest allocation
        # first, never below 1, so small populations keep representation.
        for p in sorted(pops, key=lambda q: (-base[q], q)):
            while over > 0 and base[p] > 1:
                base[p] -= 1
                over -= 1
            if over == 0:
                break
    elif over < 0:
        # Distribute what's left by largest fractional remainder.
        rema = sorted(pops, key=lambda q: (-(exact[q] - int(exact[q])), q))
        i = 0
        while over < 0:
            p = rema[i % len(rema)]
            if base[p] < len(by_pop[p]):
                base[p] += 1
                over += 1
            i += 1
            if i > 100000:
                break
    ref_ids, don_ids = [], []
    for p in pops:
        members = sorted(by_pop[p])
        k = min(base[p], len(members))
        ref_ids.extend(members[:k])
        don_ids.extend(members[k:])
    return sorted(ref_ids), sorted(don_ids)


HOLDOUT_KEEP = 3          # a population keeps at least this many donors
HOLDOUT_SEED = "20261118"  # fixed order within a population, independent of ID prefixes


def allocate_holdout(cands, k, cap):
    """k picks across populations, proportional to what each population can
    give (its candidates, capped at cap[g]) by largest remainder, with NO
    population forced in (unlike a reference, a 10-genome holdout must not
    take one from every population). Within a population the order is a
    seeded hash of the ID, so a holdout is reproducible but not biased towards
    one cohort's ID prefix (HGDP before LP6005 before S_)."""
    import hashlib
    size = {g: min(len(v), cap.get(g, 0)) for g, v in cands.items()}
    size = {g: n for g, n in size.items() if n > 0}
    tot = sum(size.values())
    k = min(k, tot)
    if k <= 0:
        return []
    exact = {g: k * size[g] / tot for g in size}
    base = {g: int(exact[g]) for g in size}
    left = k - sum(base.values())
    for g in sorted(size, key=lambda q: (-(exact[q] - base[q]), -size[q], q)):
        if left == 0:
            break
        if base[g] < size[g]:
            base[g] += 1
            left -= 1
    def order(s):
        return hashlib.md5((HOLDOUT_SEED + s).encode()).hexdigest()
    picks = []
    for g in sorted(base):
        picks += sorted(cands[g], key=order)[:base[g]]
    return sorted(picks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True,
                    help="lai_ref_panel_samples.k<K>.tsv from 4a")
    ap.add_argument("--out", default="reference_ids")
    ap.add_argument("--suffix", default="b2",
                    help="written as <sp>_reference_<suffix>.txt / <sp>_simulation_<suffix>.txt")
    ap.add_argument("--mode", choices=("benchmark", "lai"), default="benchmark",
                    help="benchmark: 1000 Genomes is never a donor (its "
                         "haplotypes sit inside the imputation servers' "
                         "reference panels). lai: an even stratified split for "
                         "every super-pop with no source rule, because the LAI "
                         "and ADMIXTURE panels are not imputed against those "
                         "servers and so carry no such circularity.")
    ap.add_argument("--source-split", default="",
                    help="deprecated and ignored; --mode now selects the policy")
    ap.add_argument("--kg-cohort", default="1KG",
                    help="value of the cohort column marking 1000 Genomes")
    ap.add_argument("--server-panel-ids",
                    default="reference_ids/sgdp_server_panel_ids.txt",
                    help="extra ids that are inside an imputation server's "
                         "reference panel under ANOTHER identifier, and so "
                         "are donor-ineligible even though their cohort is "
                         "not 1000 Genomes. Written by "
                         "scripts/sgdp-alias-duplicates.py. Missing file = "
                         "empty, with a warning.")
    ap.add_argument("--ref-fraction", type=float, default=REF_FRACTION_DEFAULT)
    ap.add_argument("--ref-floor", type=int, default=REF_FLOOR_DEFAULT)
    ap.add_argument("--ref-target", type=int, default=REF_TARGET_DEFAULT)
    ap.add_argument("--donor-floor", type=int, default=DONOR_FLOOR_DEFAULT)
    ap.add_argument("--ref-share", default="",
                    help="how much of each super-pop's NON-1000G pool goes to "
                         "the reference, as comma-separated terms: "
                         "'default=0,NAT@MXB=0.5,MEN=0.5' (the contract's "
                         "value). default=0 means the reference is 1000 "
                         "Genomes ONLY and every homogeneous non-1KG sample "
                         "is a donor -- which is the point, since 1KG can "
                         "never BE a donor, so spending non-1KG samples on a "
                         "role 1KG already fills wastes the only "
                         "donor-eligible material there is. Without this "
                         "flag the old even split applies. A term may be "
                         "keyed on the super-pop (MEN=0.5) or on a cohort "
                         "within it (NAT@MXB=0.5 -> AMR's MXB samples).")
    ap.add_argument("--ref-min", type=int, default=0,
                    help="smallest REFERENCE panel for a split super-pop, when "
                         "its pool allows it (the rule of plan_simset.R "
                         "--ref-min on the Nextflow side): raised to this many "
                         "where share x pool and the 1000 Genomes top-up fall "
                         "short, but never past half the non-1KG pool when "
                         "that pool is under 2 x ref-min, and never so far "
                         "that donors drop below --donor-floor. MEN, with no "
                         "1000 Genomes population, is the case it exists for. "
                         "Not applied to --no-donor-superpops (OCE), whose "
                         "remainder is the test pool. 0 = off.")
    ap.add_argument("--pop-aliases",
                    default="reference_ids/population_aliases.tsv",
                    help="labels that name the same ancestry under different "
                         "conventions. Gap detection compares population "
                         "labels, and the sources name the same people "
                         "differently -- 1000 Genomes YRI and HGDP Yoruba are "
                         "one population sampled twice. Without this, Yoruba "
                         "reads as uncovered while YRI sits in the reference, "
                         "and --ref-min-per-pop pulls donor-eligible samples "
                         "to fill a gap that does not exist.")
    ap.add_argument("--alias-confidence", choices=("same", "close"),
                    default="same",
                    help="'same' (default) applies only labels that are the "
                         "same population sampled twice. 'close' also treats "
                         "regionally equivalent populations as covering each "
                         "other -- a judgement about the analysis, so it is "
                         "opt-in.")
    ap.add_argument("--ref-min-per-pop", default="0",
                    help="guarantee at least this many REFERENCE samples from "
                         "every population that has donors. Exists because "
                         "--ref-share default=0 leaves sub-populations "
                         "present only in the non-1KG sources with no "
                         "reference representation at all: 1000 Genomes AFR "
                         "has no Biaka, Mbuti, San or Mandenka, EUR no "
                         "Basque, Sardinian or Orcadian, EAS none of the HGDP "
                         "minority populations. Simulating a donor haplotype "
                         "whose ancestry no reference covers measures the gap, "
                         "not the method. 2 or 3 is usually enough to anchor "
                         "a lineage; 0 keeps the strict contract behaviour. "
                         "Either one number for every super-pop, or terms like "
                         "--ref-share's: 'default=3,EUR=0'. A super-pop NAMED "
                         "in the terms is a decision already taken: its "
                         "uncovered populations are reported, not listed "
                         "under NEEDS A DECISION (EUR=0 since 18 Nov: its "
                         "donors are too few to spend on coverage).")
    ap.add_argument("--relatedness", default="",
                    help="KING table (1e's king.kin0). Relatives at or above "
                         "--role-kinship are kept TOGETHER in one role. Since "
                         "2b59278 ingest_simset.sh REFUSES any pair >= 0.177 "
                         "sitting in two different roles (donor / reference / "
                         "real test), and that check runs after the whole "
                         "build. Disjointness by identity is not enough: a "
                         "donor whose sibling is a test genome passes every "
                         "earlier check.")
    ap.add_argument("--role-kinship", type=float, default=0.177,
                    help="kinship at or above which two samples must share a "
                         "role (default 0.177, first degree)")
    ap.add_argument("--anchors-to-reference", action="store_true",
                    help="OPTIONAL. Pin supervised ADMIXTURE anchors (the "
                         "`assigned` column) to the REFERENCE side. Anchors "
                         "are NOT barred from any role -- the supervised run "
                         "only identifies the pools and its Q is never a "
                         "benchmark outcome -- so this is a preference, not a "
                         "correctness fix: it concentrates the purest samples "
                         "in the reference panel at the cost of a smaller "
                         "donor pool, worst where it is already thin. Leave "
                         "it off to let anchors fall wherever the stratified "
                         "split puts them.")
    ap.add_argument("--no-donor-superpops", default="",
                    help="comma-separated super-pops that get a reference but "
                         "NO donors, because no simulation model draws on "
                         "them (OCE is the case this exists for: Table 1 has "
                         "no OCE model, so an OCE donor is withheld for "
                         "nothing). Their non-reference samples are left "
                         "unallocated and so stay test-eligible.")
    ap.add_argument("--test-holdout", default="",
                    help="homogeneous genomes to take back from the DONORS as "
                         "real_homogeneous TEST genomes, per super-pop: "
                         "'AFR=15,EAS=15,EUR=10' (NAT accepted for AMR). With "
                         "--ref-share default=0 every homogeneous genome "
                         "outside 1000 Genomes is a donor, which leaves almost "
                         "no real homogeneous test genome for those groups. "
                         "Taken stratified over populations, never leaving a "
                         "population fewer than 3 donors or the super-pop under "
                         "--donor-floor, and only from donors with no relative "
                         "anywhere in the pool, so no relative can sit in "
                         "another role. Empty = none.")
    ap.add_argument("--pool-column", default="included",
                    help="the panel column (1 = in the pool) to split. 4a writes "
                         "`homogeneous` (the shared pool: duplicates and named "
                         "exclusions out, no project's roles) and `included` "
                         "(the LAI panel: homogeneous minus the LAI simulation "
                         "donors and their relatives). The benchmark split reads "
                         "`homogeneous`; the default `included` keeps the LAI "
                         "mode, and old tables, as they were.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    _ = args.source_split  # deprecated; --mode selects the policy

    # --ref-min-per-pop: one number, or per-super-pop terms ('default=3,EUR=0').
    # NAT is accepted for AMR, as in --ref-share.
    _sp_alias = {"NAT": "AMR"}
    rmpp, rmpp_named = {"__default__": 0}, set()
    for term in str(args.ref_min_per_pop).split(","):
        term = term.strip()
        if not term:
            continue
        key, val = term.split("=", 1) if "=" in term else ("default", term)
        key = key.strip()
        try:
            n_min = int(val.strip())
        except ValueError:
            sys.exit(f"ERROR: --ref-min-per-pop term '{term}': '{val.strip()}' "
                     f"is not a whole number")
        if n_min < 0:
            sys.exit(f"ERROR: --ref-min-per-pop term '{term}' is negative")
        if key.lower() == "default":
            rmpp["__default__"] = n_min
        else:
            key = _sp_alias.get(key.upper(), key.upper())
            rmpp[key] = n_min
            rmpp_named.add(key)

    def ref_min_per_pop_for(sp):
        return rmpp.get(sp, rmpp["__default__"])

    if str(args.ref_min_per_pop).strip() not in ("", "0"):
        print(f"[split] ref-min-per-pop: default={rmpp['__default__']}"
              + "".join(f", {k}={rmpp[k]}" for k in sorted(rmpp_named)))

    # --test-holdout: real homogeneous test genomes taken back from the donors.
    holdout = {}
    for term in args.test_holdout.split(","):
        term = term.strip()
        if not term:
            continue
        if "=" not in term:
            sys.exit(f"ERROR: --test-holdout term '{term}' has no '=' (e.g. AFR=15)")
        key, val = term.split("=", 1)
        try:
            n_h = int(val.strip())
        except ValueError:
            sys.exit(f"ERROR: --test-holdout term '{term}': '{val.strip()}' is not a whole number")
        if n_h < 0:
            sys.exit(f"ERROR: --test-holdout term '{term}' is negative")
        holdout[_sp_alias.get(key.strip().upper(), key.strip().upper())] = n_h
    if holdout:
        if args.mode != "benchmark":
            sys.exit("ERROR: --test-holdout applies to --mode benchmark only")
        print("[split] test holdout (real homogeneous test genomes taken from the donors): "
              + ", ".join(f"{k}={v}" for k, v in sorted(holdout.items())))
        if not args.relatedness:
            print("[split] WARNING: --test-holdout without --relatedness cannot keep a "
                  "held-out genome's relatives out of the other roles", file=sys.stderr)

    # Ids that are inside a server panel under another identifier. Only used
    # in benchmark mode; in lai mode no source rule applies at all.
    server_ids = set()
    if args.mode == "benchmark":
        try:
            with open(args.server_panel_ids) as fh:
                for line in fh:
                    s = line.strip()
                    if s and not s.startswith("#"):
                        server_ids.add(s)
            print(f"[split] server-panel ids (donor-ineligible): {len(server_ids)}")
        except FileNotFoundError:
            print(f"[split] WARNING: {args.server_panel_ids} not found.\n"
                  f"        SGDP re-sequenced 1000 Genomes individuals; without\n"
                  f"        this list they read as cohort SGDP and become\n"
                  f"        DONOR-ELIGIBLE while their haplotypes sit in the\n"
                  f"        imputation servers' panels under an NA/HG id.\n"
                  f"        Generate it: python3 -I scripts/sgdp-alias-duplicates.py",
                  file=sys.stderr)

    # --- read the panel ---------------------------------------------------
    rows = []
    with open(args.panel, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"sample", "cohort", "pop", "super_pop"}
        missing = need - set(rd.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: {args.panel} lacks column(s): {sorted(missing)}\n"
                     f"       found: {rd.fieldnames}")
        # WHICH POOL. 4a writes two flags (18 Nov): `homogeneous`, the pool
        # both projects choose from, and `included`, the LAI project's panel
        # (homogeneous minus ITS simulation donors and their relatives). The
        # benchmark reads `homogeneous`: the LAI project's donors are no reason
        # to withhold a genome from the benchmark, whose own roles are decided
        # here. A table older than the column cannot say who the LAI hold-out
        # removed, so asking for it there is an error, not a silent fallback.
        cols = rd.fieldnames or []
        pool_col = args.pool_column
        if pool_col not in cols:
            if pool_col == "included":
                pool_col = None          # a plain list with no flag: every row
            else:
                sys.exit(f"ERROR: --pool-column {pool_col}: {args.panel} has no such "
                         f"column (it predates 4a's 18 Nov `homogeneous` flag). Rerun "
                         f"4a, or pass --pool-column included to reproduce the old "
                         f"pool, which leaves out the LAI simulation donors and "
                         f"their relatives.")
        n_seen = 0
        for r in rd:
            n_seen += 1
            if pool_col and str(r.get(pool_col, "")).strip() not in ("1", "true", "True"):
                continue
            rows.append(r)
        print(f"[split] pool: {len(rows)} of {n_seen} rows with "
              + (f"{pool_col}=1" if pool_col else "no flag column (every row)"))

    # Supervised anchors, from the `assigned` column (ADMIXTURE's .pop prior).
    # Anchors are ordinary samples for role purposes: the supervised run only
    # identifies the pools, and its Q is never a benchmark outcome, so an
    # anchor is eligible for donor, reference or real test sample under the
    # same rules as everyone else. Tracked only so --anchors-to-reference can
    # express a preference and so the roles can be reported.
    anchor_of = {}
    if "assigned" in (rd.fieldnames or []):
        for r in rows:
            v = (r.get("assigned") or "").strip()
            if v and v not in ("-", "NA", "", "."):
                anchor_of[r["sample"]] = v
        print(f"[split] supervised anchors in the panel: {len(anchor_of)}"
              + (f" (pinning to reference)" if args.anchors_to_reference
                 else " (not pinned; pass --anchors-to-reference to pin)"))
    elif args.anchors_to_reference:
        sys.exit(f"ERROR: --anchors-to-reference given but {args.panel} has no\n"
                 f"       `assigned` column, so anchors cannot be identified.\n"
                 f"       That column is written by 4a from ADMIXTURE's .pop\n"
                 f"       prior; without it, pinning would silently do nothing.")

    if not rows:
        sys.exit(f"ERROR: no samples in the pool of {args.panel}: every "
                 f"{pool_col or 'row'} was 0 -- check 4a's stage C.")
    print(f"[split] homogeneous samples: {len(rows)}")

    by_sp = defaultdict(list)
    for r in rows:
        by_sp[r["super_pop"].strip().upper()].append(r)

    os.makedirs(args.out, exist_ok=True)
    print(f"\n{'super_pop':>9} {'N':>6} {'mode':>12} {'ref':>6} {'donors':>7} "
          f"{'%don':>6}  note")
    print("-" * 72)

    # RELATIVE CLUSTERS. Union-find over pairs at or above --role-kinship, so
    # a chain (A-B, B-C) lands in one cluster even where A and C are not
    # themselves related: the intake's rule is pairwise, and splitting a chain
    # anywhere leaves a cross-role pair behind.
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx
    n_rel_pairs = 0
    if args.relatedness:
        if not os.path.exists(args.relatedness):
            sys.exit(f"ERROR: --relatedness {args.relatedness} does not exist.")
        with open(args.relatedness) as fh:
            head = None
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if head is None:
                    head = [c.lstrip("#") for c in f]
                    def _f(c):
                        for n in c:
                            if n in head:
                                return head.index(n)
                        return -1
                    i1 = _f(("ID1", "IID1", "FID1"))
                    i2 = _f(("ID2", "IID2", "FID2"))
                    ik = _f(("KINSHIP", "Kinship", "KIN", "PropIBD"))
                    if min(i1, i2, ik) < 0:
                        sys.exit(f"ERROR: {args.relatedness} header {head} has "
                                 f"no ID1/ID2/KINSHIP columns")
                    continue
                if len(f) <= max(i1, i2, ik):
                    continue
                try:
                    k = float(f[ik])
                except ValueError:
                    continue
                if k >= args.role_kinship:
                    union(f[i1], f[i2]); n_rel_pairs += 1
        print(f"[split] relatedness: {n_rel_pairs} pair(s) >= "
              f"{args.role_kinship} will be kept within one role")

    # For --test-holdout: a donor turned test genome must have no relative
    # anywhere in the pool, since every pool genome ends up a donor, a
    # reference or (outside 1KG) a test genome itself -- so no relative of a
    # held-out genome can sit in another role. Relatives OUTSIDE the pool
    # (admixed genomes) are test genomes too, which the intake allows.
    pool_ids = {r["sample"] for r in rows}
    _members = defaultdict(set)
    for _s in list(parent):
        _members[find(_s)].add(_s)
    def has_pool_relative(x):
        return x in parent and any(m != x and m in pool_ids for m in _members[find(x)])
    held_report = []

    # label -> ancestry group. Labels with no row are their own group, so an
    # unlisted population is never silently merged with anything.
    alias_group = {}
    if args.pop_aliases and os.path.exists(args.pop_aliases):
        n_alias = 0
        with open(args.pop_aliases) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                f = line.split("\t")
                if len(f) < 3 or f[0] == "label":
                    continue
                if f[2].strip() == "same" or args.alias_confidence == "close":
                    alias_group[f[0].strip()] = f[1].strip()
                    n_alias += 1
        print(f"[split] population aliases: {n_alias} label(s) grouped "
              f"(confidence '{args.alias_confidence}')")
    elif any(v > 0 for v in rmpp.values()):
        print(f"[split] WARNING: no alias file at {args.pop_aliases}. Gap "
              f"detection will compare raw labels, so the same population "
              f"under two names (YRI / Yoruba) reads as a gap and pulls "
              f"donors to fill it.", file=sys.stderr)
    def anc_group(pop):
        return alias_group.get(pop, pop)

    COMPONENT_SP = {"NAT": "AMR", "EUR": "EUR", "AFR": "AFR", "MEN": "MEN",
                    "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}
    ref_share, ref_share_cohort = {}, {}
    if args.ref_share:
        for term in args.ref_share.split(","):
            term = term.strip()
            if not term:
                continue
            if "=" not in term:
                sys.exit(f"ERROR: --ref-share term '{term}' has no '='")
            key, val = term.split("=", 1)
            try:
                frac = float(val)
            except ValueError:
                sys.exit(f"ERROR: --ref-share '{term}' value is not a number")
            if not 0.0 <= frac <= 1.0:
                sys.exit(f"ERROR: --ref-share '{term}' must be in [0,1]")
            key = key.strip()
            if key == "default":
                ref_share["__default__"] = frac
            elif "@" in key:
                sp_key, coh = key.split("@", 1)
                sp_key = COMPONENT_SP.get(sp_key.upper(), sp_key.upper())
                ref_share_cohort[(sp_key, coh.strip())] = frac
            else:
                sp_key = COMPONENT_SP.get(key.upper(), key.upper())
                ref_share[sp_key] = frac
        print(f"[split] ref-share: default="
              f"{ref_share.get('__default__', 'even split')}"
              + ("".join(f", {k}={v}" for k, v in sorted(ref_share.items())
                         if k != "__default__"))
              + ("".join(f", {k[0]}@{k[1]}={v}"
                         for k, v in sorted(ref_share_cohort.items()))))

    written, problems, comp_lines, decided = [], [], [], []
    no_donor_sps = {s.strip().upper() for s in args.no_donor_superpops.split(",")
                    if s.strip()}
    if no_donor_sps:
        print(f"[split] reference-only super-pops (no donors): "
              f"{sorted(no_donor_sps)}")
    for sp in sorted(by_sp):
        members = by_sp[sp]
        n = len(members)

        # 1000 GENOMES IS NEVER A DONOR.
        #
        # The benchmark's imputation reference panels (TOPMed, Michigan, All of
        # Us) all contain 1000 Genomes. A 1000 Genomes sample used as a
        # simulation donor would therefore have its own haplotypes inside the
        # server's panel, and the imputation would be scored on data the server
        # has memorised. That is circular in a way no split can repair, so the
        # exclusion is by SOURCE and applies to every super-pop, not only the
        # three that have enough 1000 Genomes to build a reference from.
        # In lai mode there is no source rule: every sample is eligible for
        # either role, so the "non-1000G" pool is simply everyone and the
        # 1000G top-up pool is empty. The rest of the allocation is identical,
        # which is the point -- one stratified, deterministic code path, with
        # the policy difference expressed only in who is eligible.
        if args.mode == "lai":
            kg, nonkg = [], list(members)
        else:
            # Donor-ineligible = cohort is 1000 Genomes, OR the sample is a
            # 1000 Genomes individual re-sequenced under another id. SGDP
            # re-sequenced 41 of the people on our candidate list; they read
            # as cohort SGDP and would otherwise be donor-eligible while their
            # haplotypes sit in TOPMed's panel under an NA/HG identifier.
            # De-duplication does not catch this: the 1000 Genomes copy need
            # not be in OUR panel for the server to hold that person.
            def _is_kg(r):
                return (r["cohort"].strip().upper() == args.kg_cohort.upper()
                        or r["sample"].strip() in server_ids)
            kg = sorted(r["sample"] for r in members if _is_kg(r))
            nonkg = [r for r in members if not _is_kg(r)]

        # THE 1000 GENOMES TOP-UP IS STRATIFIED BY POPULATION.
        #
        # Until 18 Nov it was kg[:k] -- the first k IDs in alphabetical order --
        # so the reference took whichever 1KG populations have the lowest IDs:
        # EUR got GBR and FIN only (HG00xxx; no CEU, IBS or TSI), AFR got ACB
        # and GWD (no YRI, LWK, ESN, MSL), EAS got CHS. That is not only a
        # skewed panel; it also cost DONORS, because the coverage rule below
        # then found HGDP Yoruba uncovered (no YRI in the reference) and pulled
        # Yoruba donors into the reference to cover them. Now the k slots are
        # spread over the 1KG populations by the same largest-remainder rule
        # as everything else, and a 1KG sample related (>= --role-kinship) to a
        # donor is skipped: the pair would sit in two roles, which the intake
        # refuses.
        pop_of_kg = {r["sample"]: (r["pop"].strip() or "NA") for r in members}
        def kg_topup(k, donors):
            if k <= 0 or not kg:
                return []
            dcl = {find(d) for d in donors} if args.relatedness else set()
            by_pop_kg = defaultdict(list)
            for x in kg:
                if args.relatedness and find(x) in dcl:
                    continue
                by_pop_kg[pop_of_kg.get(x, "NA")].append(x)
            avail = sum(len(v) for v in by_pop_kg.values())
            k = min(k, avail)
            if k <= 0:
                return []
            if k < len(by_pop_kg):
                # Fewer slots than populations: one each from the largest.
                big = sorted(by_pop_kg, key=lambda q: (-len(by_pop_kg[q]), q))[:k]
                return sorted(sorted(by_pop_kg[q])[0] for q in big)
            return allocate_stratified(by_pop_kg, k)[0]

        # A super-pop that no simulation model draws on needs no donors. OCE
        # is the case this exists for: Table 1 has no OCE model, so every OCE
        # sample withheld as a donor is withheld for nothing -- and in a group
        # that thin, those samples are worth far more as test genomes, since
        # OCE and OTH are precisely the ancestries that cannot be pooled into
        # a stratified analysis and so have to be tested directly.
        if sp in no_donor_sps:
            # The reference takes up to its target and the REMAINDER IS LEFT
            # UNALLOCATED, which is the whole point: an unallocated sample is
            # neither donor nor reference, so the cohort assembler picks it up
            # as a test genome. Making all of them reference instead would
            # leave this super-pop with no test samples at all, which for OCE
            # would mean no OCE reading anywhere in the benchmark.
            # NOT min(n, ref_target). For a group of ~20 that takes every
            # sample as reference and leaves NO test genomes -- which for OCE
            # would mean the benchmark has no OCE imputation reading at all,
            # the exact gap that putting OCE in the cohort was meant to close.
            # Use the same sizing policy as a split group, with no donor floor
            # since there are no donors, and leave the remainder unallocated
            # so the assembler picks it up as test samples.
            n_ref_only, _nd, _st = plan_sizes(
                n, args.ref_fraction, args.ref_floor, args.ref_target, 0)
            n_ref_only = min(n_ref_only, n)
            # 1000 GENOMES FIRST HERE TOO (18 Nov, for SAS). A 1KG sample can
            # never be a test genome (it is inside the servers' panels), so it
            # costs nothing as reference; a non-1KG genome left unallocated is
            # a test genome. SAS has ~10 homogeneous non-1KG genomes, too few to
            # simulate from, and the point of making it reference-only is that
            # they are read as real_homogeneous test genomes instead -- so the
            # reference comes from 1KG, spread over its populations, skipping a
            # 1KG sample related to one of them (the pair would sit in two
            # roles). Non-1KG genomes go to the reference only for what 1KG
            # cannot cover: all of it for OCE, which has no 1KG population.
            nonkg_ids = [r["sample"] for r in nonkg]
            ref_kg = kg_topup(n_ref_only, nonkg_ids)
            by_pop_nk = defaultdict(list)
            for r in nonkg:
                by_pop_nk[r["pop"].strip() or "NA"].append(r["sample"])
            ref_nk, rest = (allocate_stratified(by_pop_nk, n_ref_only - len(ref_kg))
                            if n_ref_only > len(ref_kg) and nonkg_ids
                            else ([], sorted(nonkg_ids)))
            ref = sorted(ref_kg + ref_nk)
            print(f"{sp:>9} {n:>6} {'ref-only':>12} {len(ref):>6} {0:>7} "
                  f"{0:>5.0f}%  no model draws on {sp}: 0 donors; ref = "
                  f"{len(ref_nk)} non-1000G + {len(ref_kg)} 1000G; "
                  f"{len(rest)} non-1000G left for the test pool")
            written.append((sp, ref, []))
            continue

        n_ref, n_don, status = plan_sizes(
            n, args.ref_fraction, args.ref_floor, args.ref_target,
            args.donor_floor)
        if status == "UNDERSUPPLIED":
            note = (f"*** {n} < {args.ref_floor}+{args.donor_floor}: cannot "
                    f"supply both roles -- DECIDE, do not default")
            print(f"{sp:>9} {n:>6} {'--':>12} {'-':>6} {'-':>7} {'-':>6}  {note}")
            problems.append((sp, note))
            continue

        # MATCHED COMPOSITION.
        #
        # The non-1000G samples are split EVENLY and stratified by population,
        # so the reference and donor sides carry the same sub-population mix.
        # The reference is then topped up from 1000 Genomes to its target.
        #
        # Splitting the non-1000G half evenly is what makes the two sides
        # comparable: if the reference were taken whole from one source and the
        # donors from another, a difference between them would be a difference
        # between datasets, and the benchmark would read it as an ancestry
        # effect. The 1000 Genomes top-up unavoidably adds populations the
        # donor side cannot have -- that is the price of the no-donor rule
        # above, and the per-population table below reports it rather than
        # hiding it.
        by_pop_nonkg = defaultdict(list)
        for r in nonkg:
            by_pop_nonkg[r["pop"].strip() or "NA"].append(r["sample"])

        n_nonkg = len(nonkg)
        # --- how much of the NON-1KG pool the reference takes -------------
        # Without --ref-share this stays None and the even split below runs.
        share = None
        if args.ref_share:
            share = ref_share.get(sp, ref_share.get("__default__"))
        nonkg_by_cohort = defaultdict(list)
        for r in nonkg:
            nonkg_by_cohort[r["cohort"].strip()].append(r)
        # How many non-1000G go to the reference: half of them, but never more
        # than the reference target, and never so many that donors starve.
        if share is None:
            ref_from_nonkg = min(n_nonkg // 2, n_ref)
            if n_nonkg - ref_from_nonkg < args.donor_floor:
                ref_from_nonkg = max(0, n_nonkg - args.donor_floor)
        else:
            # A cohort-keyed term (NAT@MXB=0.5) applies to THAT cohort only;
            # everything else in the super-pop takes the super-pop or default
            # share. Summed, not averaged: "half the MXB" is a count of MXB
            # samples, and averaging it over a pool that also holds HGDP and
            # SGDP would silently take half of those too.
            want = 0.0
            for coh, members_c in nonkg_by_cohort.items():
                f = ref_share_cohort.get((sp, coh), share)
                want += f * len(members_c)
            ref_from_nonkg = int(round(want))
            # The cap still applies: the reference never exceeds its target.
            ref_from_nonkg = min(ref_from_nonkg, n_ref, n_nonkg)
            # And the donor floor still applies, but in the other direction --
            # with share=0 donors are maximal, so this can only bite when a
            # share is high.
            if n_nonkg - ref_from_nonkg < args.donor_floor:
                ref_from_nonkg = max(0, n_nonkg - args.donor_floor)

        # REFERENCE FLOOR (--ref-min). A share of a small pool is a small
        # panel: MEN at 0.33 of 59 was 18 genomes. Raise the non-1KG part until
        # the whole reference (non-1KG + what 1KG can supply) reaches ref-min,
        # within the limits plan_simset.R uses -- at most half the non-1KG pool
        # when that pool is under 2 x ref-min -- and never below the donor floor.
        cohort_keyed = any(k[0] == sp for k in ref_share_cohort)
        floor_note = ""
        if args.ref_min > 0 and not cohort_keyed:
            kg_room = min(len(kg), max(0, n_ref - ref_from_nonkg))
            if ref_from_nonkg + kg_room < args.ref_min:
                need = args.ref_min - min(len(kg), args.ref_min)
                lim = n_nonkg // 2 if n_nonkg < 2 * args.ref_min else n_nonkg
                lim = min(lim, max(0, n_nonkg - args.donor_floor))
                raised = min(max(ref_from_nonkg, need), lim)
                if raised > ref_from_nonkg:
                    floor_note = f"; reference raised {ref_from_nonkg} -> {raised} non-1000G (--ref-min {args.ref_min})"
                    ref_from_nonkg = raised
                    n_ref = max(n_ref, ref_from_nonkg)

        # ALLOCATE PER COHORT WHEN A COHORT-KEYED SHARE APPLIES.
        #
        # NAT@MXB=0.5 means "half the MXB", not "a number equal to half the
        # MXB, taken from the whole non-1KG pool". Pooling first and
        # allocating proportionally hits the right TOTAL with the wrong
        # COMPOSITION -- on a 30 MXB / 10 HGDP pool it took 11 MXB and 4
        # HGDP, when the contract wants 15 MXB and no HGDP, and those 4 HGDP
        # samples are donor-eligible material spent on a role 1KG already
        # fills. So each cohort is allocated against its own share and the
        # results are concatenated.
        if ref_from_nonkg <= 0:
            ref_nonkg, don = [], sorted(r["sample"] for r in nonkg)
        elif share is not None and cohort_keyed:
            ref_nonkg, don = [], []
            for coh in sorted(nonkg_by_cohort):
                members_c = nonkg_by_cohort[coh]
                f = ref_share_cohort.get((sp, coh), share)
                k_c = int(round(f * len(members_c)))
                by_pop_c = defaultdict(list)
                for r in members_c:
                    by_pop_c[r["pop"].strip() or "NA"].append(r["sample"])
                if k_c <= 0:
                    don.extend(sorted(r["sample"] for r in members_c))
                    continue
                if args.relatedness:
                    rc, dc = allocate_clustered(by_pop_c, k_c, find)
                else:
                    rc, dc = allocate_stratified(by_pop_c, k_c)
                ref_nonkg.extend(rc); don.extend(dc)
            ref_nonkg, don = sorted(ref_nonkg), sorted(don)
        elif args.relatedness:
            ref_nonkg, don = allocate_clustered(by_pop_nonkg, ref_from_nonkg, find)
        else:
            ref_nonkg, don = allocate_stratified(by_pop_nonkg, ref_from_nonkg)

        # ANCHORS PINNED TO THE REFERENCE SIDE.
        #
        # Done AFTER the stratified allocation and as a MOVE, not an addition:
        # adding them would grow the reference past its target, and dropping
        # them from the donor side without compensating would shrink the
        # cohort. Each anchor moved out of the donor half is replaced by a
        # donor-side sample from the same POPULATION where one exists, so the
        # proportional sub-population matching the stratification exists to
        # produce is preserved.
        n_moved = 0
        if args.anchors_to_reference and anchor_of:
            r_set, d_set = set(ref_nonkg), set(don)
            pop_of_s = {r["sample"]: (r["pop"].strip() or "NA") for r in nonkg}
            for aid in sorted(d_set & set(anchor_of)):
                swap_pop = pop_of_s.get(aid, "NA")
                cand = sorted(x for x in r_set
                              if pop_of_s.get(x, "NA") == swap_pop
                              and x not in anchor_of)
                if not cand:
                    # No same-population swap available: move the anchor and
                    # accept a one-sample imbalance rather than break the
                    # population match by swapping across populations.
                    d_set.discard(aid); r_set.add(aid)
                else:
                    victim = cand[0]
                    r_set.discard(victim); d_set.add(victim)
                    d_set.discard(aid); r_set.add(aid)
                n_moved += 1
            ref_nonkg, don = sorted(r_set), sorted(d_set)

        # --- sub-population coverage -------------------------------------
        # With --ref-share default=0 the reference is 1000 Genomes only, and
        # 1KG does not carry the deeply divergent lineages the non-1KG sources
        # do: no Biaka, Mbuti, San or Mandenka in AFR, no Basque, Sardinian or
        # Orcadian in EUR, none of the HGDP minority populations in EAS.
        # Simulating a donor haplotype whose ancestry NO reference covers
        # measures that gap rather than the method, so --ref-min-per-pop pulls
        # a few back per uncovered population.
        n_pulled = 0
        min_per_pop = ref_min_per_pop_for(sp)
        if min_per_pop > 0:
            r_set, d_set = set(ref_nonkg), set(don)
            pop_of_nk = {r["sample"]: (r["pop"].strip() or "NA") for r in nonkg}
            # COVERAGE IS COUNTED ON THE 1000 GENOMES SIDE TOO.
            #
            # The reference is filled from 1KG first, so by the time this runs
            # most of it is 1KG -- and those samples cover ancestry groups the
            # donor-side labels also belong to. Counting only the non-1KG part
            # of the reference would report every HGDP population as uncovered
            # and pull donors to fill gaps 1KG already fills.
            have = defaultdict(int)
            for x in r_set:
                have[anc_group(pop_of_nk.get(x, "NA"))] += 1
            for x in kg_topup(max(0, n_ref - len(ref_nonkg)), d_set):
                have[anc_group(pop_of_kg.get(x, "NA"))] += 1
            by_pop_donor = defaultdict(list)
            for x in d_set:
                by_pop_donor[anc_group(pop_of_nk.get(x, "NA"))].append(x)
            for pop in sorted(by_pop_donor):
                need = min_per_pop - have[pop]
                if need <= 0:
                    continue
                # Never take a population's LAST donors: a population that
                # exists only as a handful of samples is more use as donors
                # for its own component than as an under-powered reference.
                avail = sorted(by_pop_donor[pop])
                take = min(need, max(0, len(avail) - min_per_pop))
                # Coverage never takes the super-pop's donors below the donor
                # floor: where donors are that scarce, a donor is worth more
                # than a lineage's reference coverage (18 Nov).
                take = min(take, max(0, len(d_set) - args.donor_floor))
                # A sample with a relative among the donors stays a donor:
                # moving it alone would put the pair in two roles.
                if args.relatedness:
                    avail = [x for x in avail
                             if not any(y != x and find(y) == find(x) for y in d_set)]
                for x in avail[:take]:
                    d_set.discard(x); r_set.add(x); n_pulled += 1
            if n_pulled:
                ref_nonkg, don = sorted(r_set), sorted(d_set)

        # --- real homogeneous test genomes, taken back from the donors -------
        # (18 Nov) With every homogeneous genome outside 1KG a donor, the real
        # test half had almost no homogeneous genome for AFR, EAS or EUR. A
        # few are taken back here, stratified over populations; they land in
        # neither list, so the assembler reads them as real_homogeneous.
        held = []
        k_hold = holdout.get(sp, 0)
        if k_hold > 0:
            pop_of_d = {r["sample"]: (r["pop"].strip() or "NA") for r in nonkg}
            ndon_g, cands = defaultdict(int), defaultdict(list)
            for x in don:
                g = anc_group(pop_of_d.get(x, "NA"))
                ndon_g[g] += 1
                if not has_pool_relative(x):
                    cands[g].append(x)
            cap = {g: max(0, ndon_g[g] - HOLDOUT_KEEP) for g in ndon_g}
            room = max(0, len(don) - args.donor_floor)
            held = allocate_holdout(cands, min(k_hold, room), cap)
            if held:
                don = sorted(set(don) - set(held))
            held_report.append((sp, k_hold, held, room,
                                sum(min(len(cands[g]), cap[g]) for g in cands),
                                {x: pop_of_d.get(x, "NA") for x in held}))

        # FILL ORDER: 1000 GENOMES FIRST, NON-1KG ONLY FOR WHAT IT CANNOT
        # COVER.
        #
        # 1KG can ONLY ever be reference -- it is barred from the donor side
        # because the imputation servers contain it. A non-1KG sample can be
        # either. So the single-use resource is spent first and the dual-use
        # one only on the remainder, which is what --ref-share default=0
        # expresses: for AFR, EUR, EAS and SAS the 1KG slice fills the
        # reference and every homogeneous non-1KG sample stays a donor.
        #
        # The rule degenerates correctly where 1KG cannot cover the super-pop,
        # and that is exactly where the contract gives a non-zero share:
        #   MEN  1000 Genomes has no Middle Eastern population  -> MEN=0.5
        #        (0.33 until 18 Nov: 18 reference genomes of 59), with --ref-min
        #   NAT  its American populations are ADMIXED, not Indigenous, and
        #        mostly will not clear the homogeneity bar       -> NAT@MXB=0.5
        #   OCE  no Oceanian population at all          -> all reference is non-1KG
        top_up = max(0, n_ref - len(ref_nonkg))
        ref = sorted(ref_nonkg + kg_topup(top_up, don))

        if len(don) < args.donor_floor:
            note = (f"*** only {len(don)} donors available "
                    f"({n_nonkg} non-1000G in total) -- below the donor floor "
                    f"of {args.donor_floor}")
            problems.append((sp, note))
        else:
            note = (f"{len(by_pop_nonkg)} non-1000G populations stratified; "
                    f"ref = {len(ref_nonkg)} non-1000G + {len(ref)-len(ref_nonkg)} 1000G")
            if n_moved:
                note += f"; {n_moved} anchor(s) pinned to reference"
            if share is not None:
                note += f"; ref-share {share:g}"
            if n_pulled:
                note += f"; {n_pulled} pulled for sub-pop coverage"
            if held:
                note += f"; {len(held)} held out as real homogeneous test genomes"
            note += floor_note
            if args.ref_min > 0 and len(ref) < args.ref_min:
                short = ("the pool cannot supply it with the donor floor kept"
                         if ref_from_nonkg + min(len(kg), max(0, n_ref - ref_from_nonkg)) < args.ref_min
                         else "per-population rounding, relatives kept together")
                note += f"; reference {len(ref)} < --ref-min {args.ref_min} ({short})"
        pct = 100.0 * len(don) / (len(ref) + len(don)) if (ref or don) else 0
        print(f"{sp:>9} {n:>6} {'split':>12} {len(ref):>6} {len(don):>7} "
              f"{pct:>5.0f}%  {note}")

        # Composition report. The matched split is the whole point of
        # stratifying, so the evidence that it worked is printed rather than
        # assumed. A population present on one side only is named: that is
        # either a 1000 Genomes population (reference-only by the no-donor
        # rule) or a stratum too small to divide, and both are facts the
        # methods have to state.
        pop_of = {r["sample"]: (r["pop"].strip() or "NA") for r in members}
        # POPULATIONS WITH DONORS BUT NO REFERENCE. check_roles_vs_king.R
        # flags exactly this, so flag it here where it can still be fixed: a
        # simulated haplotype drawn from a population no reference sample
        # covers is scored against the nearest thing in the panel, which
        # measures the panel's gap rather than the caller.
        _ref_pops = {anc_group(pop_of.get(x, "NA")) for x in ref}
        _uncovered = sorted({anc_group(pop_of.get(x, "NA")) for x in don}
                            - _ref_pops)
        if _uncovered:
            # Each gap with the population's donor count: the methods have
            # to say which donor lineages the reference does not carry.
            _ndon = defaultdict(int)
            for x in don:
                _ndon[anc_group(pop_of.get(x, "NA"))] += 1
            def _fmt(gs):
                return (", ".join(f"{g} ({_ndon[g]})" for g in gs[:8])
                        + (" ..." if len(gs) > 8 else ""))
            if sp in rmpp_named and min_per_pop == 0:
                # Named with 0: the coverage pull was declined on purpose
                # (EUR, 18 Nov -- its donors are too few to spend).
                decided.append((sp, f"{len(_uncovered)} population(s) have donors "
                                    f"but no reference sample: {_fmt(_uncovered)}  "
                                    f"(kept as donors: --ref-min-per-pop {sp}=0)"))
            elif min_per_pop == 0:
                problems.append((sp, f"{len(_uncovered)} population(s) have donors "
                                     f"but no reference sample: {_fmt(_uncovered)}  "
                                     f"(--ref-min-per-pop pulls a few back)"))
            else:
                # The coverage rule ran and these are its own limits, not open
                # questions (18 Nov: they had been listed under NEEDS A
                # DECISION): it never takes a population's last donors, never
                # takes the super-pop under the donor floor, and never splits
                # relatives across roles. Reported, not a failure.
                small = [g for g in _uncovered if _ndon[g] <= min_per_pop]
                rest = [g for g in _uncovered if _ndon[g] > min_per_pop]
                if small:
                    decided.append((sp, f"{len(small)} population(s) with "
                                        f"{min_per_pop} or fewer donors keep them "
                                        f"all, so no reference sample: {_fmt(small)}  "
                                        f"(the coverage rule never takes a "
                                        f"population's last donors)"))
                if rest:
                    _why = (f"donors at the floor of {args.donor_floor}"
                            if len(don) <= args.donor_floor
                            else "relatives kept together in one role")
                    decided.append((sp, f"{len(rest)} population(s) not pulled "
                                        f"({_why}): {_fmt(rest)}"))
        cr, cd = defaultdict(int), defaultdict(int)
        for s in ref:
            cr[pop_of[s]] += 1
        for s in don:
            cd[pop_of[s]] += 1
        worst, worst_pop = 0.0, ""
        oneside = []
        for p in sorted(set(cr) | set(cd)):
            pr = 100.0 * cr[p] / len(ref) if ref else 0.0
            pd = 100.0 * cd[p] / len(don) if don else 0.0
            if abs(pr - pd) > worst:
                worst, worst_pop = abs(pr - pd), p
            if (cr[p] == 0) != (cd[p] == 0):
                oneside.append(f"{p}({'ref' if cd[p] == 0 else 'don'} only)")
        comp_lines.append(
            f"  {sp:<5} {len(set(cr) | set(cd)):>3} populations, "
            f"largest share gap {worst:>5.1f}pp ({worst_pop})"
            + (f"; one side only: {', '.join(oneside)}" if oneside else ""))

        overlap = set(ref) & set(don)
        if overlap:
            sys.exit(f"ERROR: {sp} has {len(overlap)} samples on BOTH sides. "
                     f"That is the circularity this script exists to prevent.")
        written.append((sp, ref, don))

    if comp_lines:
        print("\n[split] sub-population composition, reference vs donors:")
        for line in comp_lines:
            print(line)

    if args.dry_run:
        print("\n[split] --dry-run: nothing written")
    else:
        for sp, ref, don in written:
            lc = sp.lower()
            for name, ids in ((f"{lc}_reference_{args.suffix}.txt", ref),
                              (f"{lc}_simulation_{args.suffix}.txt", don)):
                with open(os.path.join(args.out, name), "w") as fh:
                    fh.write("\n".join(ids) + ("\n" if ids else ""))
        print(f"\n[split] wrote {2*len(written)} lists to {args.out}/")

    # A sample must not be a donor for one super-pop and a reference for
    # another: it is one person either way, and the method would be scored on
    # a haplotype it was given.
    ref_all, don_all = set(), set()
    for _sp, ref, don in written:
        ref_all |= set(ref); don_all |= set(don)
    cross = ref_all & don_all
    if cross:
        sys.exit(f"\nERROR: {len(cross)} sample(s) are a reference for one "
                 f"super-pop and a donor for another, e.g. "
                 f"{sorted(cross)[:5]}. Same person, same haplotypes.")
    print(f"[split] disjoint across all super-pops: "
          f"{len(ref_all)} reference, {len(don_all)} donor, 0 shared")

    if held_report:
        print("\n[split] real homogeneous test genomes held out of the donors (--test-holdout):")
        for sp, k_hold, held, room, n_cand, grp_of in held_report:
            byg = defaultdict(int)
            for x in held:
                byg[grp_of[x]] += 1
            why = ""
            if len(held) < k_hold:
                why = (f"  -- short of {k_hold}: " +
                       (f"only {room} above the donor floor of {args.donor_floor}"
                        if room < k_hold and room <= n_cand else
                        f"only {n_cand} donors without a relative in the pool, in "
                        f"populations keeping {HOLDOUT_KEEP}+ donors"))
            print(f"  {sp}: {len(held)} of {k_hold} asked"
                  + (f" ({', '.join(f'{g} {n}' for g, n in sorted(byg.items(), key=lambda kv: (-kv[1], kv[0])))})"
                     if held else "") + why)
    if decided:
        print("\n[split] decided already (reported, not a failure):")
        for sp, note in decided:
            print(f"  {sp}: {note}")
    if problems:
        print("\n[split] NEEDS A DECISION:")
        for sp, note in problems:
            print(f"  {sp}: {note}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
