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
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k8.tsv \\
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
                         "'default=0,NAT@MXB=0.5,MEN=0.33' (the contract's "
                         "value). default=0 means the reference is 1000 "
                         "Genomes ONLY and every homogeneous non-1KG sample "
                         "is a donor -- which is the point, since 1KG can "
                         "never BE a donor, so spending non-1KG samples on a "
                         "role 1KG already fills wastes the only "
                         "donor-eligible material there is. Without this "
                         "flag the old even split applies. A term may be "
                         "keyed on the super-pop (MEN=0.33) or on a cohort "
                         "within it (NAT@MXB=0.5 -> AMR's MXB samples).")
    ap.add_argument("--ref-min-per-pop", type=int, default=0,
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
                         "a lineage; 0 keeps the strict contract behaviour.")
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
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    _ = args.source_split  # deprecated; --mode selects the policy

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
        has_incl = "included" in (rd.fieldnames or [])
        for r in rd:
            if has_incl and str(r.get("included", "")).strip() not in ("1", "true", "True"):
                continue
            rows.append(r)

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
        sys.exit(f"ERROR: no included samples in {args.panel}. If the file has "
                 f"an 'included' column, every row was 0 -- check stage C.")
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

    written, problems, comp_lines = [], [], []
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
            by_pop_all = defaultdict(list)
            for r in members:
                by_pop_all[r["pop"].strip() or "NA"].append(r["sample"])
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
            ref, rest = allocate_stratified(by_pop_all, n_ref_only)
            print(f"{sp:>9} {n:>6} {'ref-only':>12} {len(ref):>6} {0:>7} "
                  f"{0:>5.0f}%  no model draws on {sp}: 0 donors, "
                  f"{len(rest)} left for the test pool")
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
        cohort_keyed = any(k[0] == sp for k in ref_share_cohort)
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
        if args.ref_min_per_pop > 0:
            r_set, d_set = set(ref_nonkg), set(don)
            pop_of_nk = {r["sample"]: (r["pop"].strip() or "NA") for r in nonkg}
            have = defaultdict(int)
            for x in r_set:
                have[pop_of_nk.get(x, "NA")] += 1
            by_pop_donor = defaultdict(list)
            for x in d_set:
                by_pop_donor[pop_of_nk.get(x, "NA")].append(x)
            for pop in sorted(by_pop_donor):
                need = args.ref_min_per_pop - have[pop]
                if need <= 0:
                    continue
                # Never take a population's LAST donors: a population that
                # exists only as a handful of samples is more use as donors
                # for its own component than as an under-powered reference.
                avail = sorted(by_pop_donor[pop])
                take = min(need, max(0, len(avail) - args.ref_min_per_pop))
                for x in avail[:take]:
                    d_set.discard(x); r_set.add(x); n_pulled += 1
            if n_pulled:
                ref_nonkg, don = sorted(r_set), sorted(d_set)

        top_up = max(0, n_ref - len(ref_nonkg))
        ref = sorted(ref_nonkg + kg[:top_up])

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
        _ref_pops = {pop_of.get(x, "NA") for x in ref}
        _uncovered = sorted({pop_of.get(x, "NA") for x in don} - _ref_pops)
        if _uncovered:
            problems.append((sp, f"{len(_uncovered)} population(s) have donors "
                                 f"but no reference sample: "
                                 f"{', '.join(_uncovered[:8])}"
                                 + (" ..." if len(_uncovered) > 8 else "")
                                 + "  (--ref-min-per-pop pulls a few back)"))
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

    if problems:
        print("\n[split] NEEDS A DECISION:")
        for sp, note in problems:
            print(f"  {sp}: {note}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
