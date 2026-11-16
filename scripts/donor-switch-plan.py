#!/usr/bin/env python3
"""
donor-switch-plan.py -- choose simulation donors so that fewer of their
relatives have to be kept out of the reference panel. READ-ONLY on
reference_ids/: proposals are written to --out-dir for review.

THE COST. A donor's first-degree relatives (KING >= --thr, 0.177) cannot be
reference haplotypes -- 4a keeps them out of the homogeneous panel and 2b
refuses a donor whose relative sits in a reference list. In 1KG, where most
samples belong to parent-child trios, the choice of donor decides how many
samples that costs:
    a relative-free donor           costs 0
    both parents of one trio        cost 1 (the child) for 2 donors
    one parent                      costs 1 (the child)
    the child                       costs 2 (both parents)
For each donor list and each population in it, this finds the set of the same
size, from the same population, with the fewest relatives kept out, keeping as
many current donors as that allows. Donors are never relatives of each other
(a parent and child would put the same haplotypes in the donor pool twice).

A candidate must be: in the panel table, in the same population as the donors
it replaces, not a 1e duplicate drop or named exclusion, not in a curated
reference list (reference_ids/*_rfmix.txt other than the *_homog_rfmix.txt
lists 4a writes), not a relative of anyone in those lists, and not a donor in
another list. 4a rebuilds the homogeneous panel without the new donors and
their relatives, so the homogeneous lists are not a constraint here -- but 2b
checks them, so 4a must run after a donor change and before 2b.

Families are the connected components of the relatedness graph; the cost adds
up across them, so each population is solved exactly: every non-related
subset of each family, then a knapsack over families for the target count.

  python3 -I scripts/donor-switch-plan.py \\
      --table <04_homogeneity_panel>/lai_ref_panel_samples.k8.tsv \\
      --kin <04_homogeneity_panel>/king/king.kin0 --refs reference_ids --out-dir results/donor_plan

Writes <out-dir>/<list>.proposed.txt (current donors kept first, in their
order, then the new ones), donor_plan_summary.tsv and donor_plan_changes.tsv.
"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

ENUM_MAX = 16          # candidates per family enumerated exactly (2^16 subsets)


def ids(path):
    try:
        with open(path) as fh:
            return [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        return []


def family_options(cands, rel, countable, current):
    """Every donor subset of one family that holds no related pair:
    (n donors, cost, current kept, subset). cost = relatives kept out."""
    cands = sorted(cands)
    opts = {}
    if len(cands) <= ENUM_MAX:
        subsets = (s for k in range(len(cands) + 1) for s in combinations(cands, k))
    else:                                   # large family: greedy, cheapest first
        order = sorted(cands, key=lambda c: (len(rel[c] & countable), c not in current, c))
        picked, subsets = [], [()]
        for c in order:
            if not (rel[c] & set(picked)):
                picked.append(c)
                subsets.append(tuple(picked))
    for s in subsets:
        ss = set(s)
        if any(rel[a] & ss for a in s):
            continue
        out = set().union(*(rel[a] for a in s)) - ss if s else set()
        val = (len(out & countable), -len(ss & current), tuple(sorted(s)))
        if len(s) not in opts or val < opts[len(s)][0]:
            opts[len(s)] = (val, s)
    return [(k, v[0][0], -v[0][1], v[1]) for k, v in opts.items()]


def solve(families, target):
    """Knapsack over families: exactly `target` donors, fewest relatives kept
    out, then most current donors kept. -> (cost, kept, donors) or None."""
    best = {0: (0, 0, ())}
    for opts in families:
        nxt = {}
        for n, (c, kept, chosen) in best.items():
            for k, cost, kp, s in opts:
                m = n + k
                if m > target:
                    continue
                cand = (c + cost, kept + kp, chosen + tuple(s))
                if m not in nxt or (cand[0], -cand[1], cand[2]) < (nxt[m][0], -nxt[m][1], nxt[m][2]):
                    nxt[m] = cand
        best = nxt
    return best.get(target)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="lai_ref_panel_samples.k<K>.tsv (sample, pop columns)")
    ap.add_argument("--kin", required=True, help="1e's king.kin0 (IID1, IID2, KINSHIP)")
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--thr", type=float, default=0.177)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.table), delimiter="\t"))
    pop = {r["sample"]: r["pop"] for r in rows}
    panel = set(pop)
    drop = set(ids(os.path.join(a.refs, "panel_duplicate_drop_ids.txt"))) \
        | set(ids(os.path.join(a.refs, "panel_exclude_ids.txt")))
    countable = panel - drop                   # a relative only costs if it could be a reference

    curated = {}
    for f in sorted(glob.glob(os.path.join(a.refs, "*_rfmix.txt"))):
        if f.endswith("_homog_rfmix.txt"):
            continue
        for s in ids(f):
            curated.setdefault(s, os.path.basename(f))
    lists = {os.path.basename(f)[:-len(".txt")]: ids(f)
             for f in sorted(glob.glob(os.path.join(a.refs, "*_simulation.txt")))}
    # A list that is exactly the union of others (amr_hgdpmxb = amr + mxb) is
    # not planned on its own: its proposal is the union of theirs.
    derived = {}
    for name, ds in lists.items():
        parts = [m for m, ms in lists.items() if m != name and ms and set(ms) < set(ds)]
        if len(parts) >= 2 and set().union(*(set(lists[m]) for m in parts)) == set(ds):
            derived[name] = parts
    base = {n: ds for n, ds in lists.items() if n not in derived}
    donor_of = {}
    for name, ds in base.items():
        for s in ds:
            donor_of.setdefault(s, name)

    rel = defaultdict(set)
    n_pairs = 0
    with open(a.kin) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        if not all(c in ix for c in ("IID1", "IID2", "KINSHIP")):
            sys.exit(f"ERROR: {a.kin} header {head} lacks IID1/IID2/KINSHIP")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            x, y, k = f[ix["IID1"]], f[ix["IID2"]], float(f[ix["KINSHIP"]])
            if k >= a.thr and x != y:
                rel[x].add(y); rel[y].add(x); n_pairs += 1
    near_curated = {s for c in curated for s in rel[c]}
    print(f"{a.table}: {len(panel)} samples; {n_pairs} pairs at KING >= {a.thr}; "
          f"{len(curated)} curated reference ids; {len(drop)} duplicate drops / named exclusions")

    os.makedirs(a.out_dir, exist_ok=True)
    summ = open(os.path.join(a.out_dir, "donor_plan_summary.tsv"), "w")
    chg = open(os.path.join(a.out_dir, "donor_plan_changes.tsv"), "w")
    summ.write("list\tpop\tdonors\tcandidates\trelative_free_candidates\tcurrent_cost\tcurrent_with_relative\t"
               "current_related_pairs\tbest_cost\tkept\treplaced\tnote\n")
    chg.write("list\tpop\taction\tsample\trelatives_kept_out\n")
    print(f"\ncost = relatives kept out of the reference (duplicate drops and named exclusions not counted)")
    print(f"{'list':22s} {'pop':16s} {'donors':>6s} {'cands':>5s} {'relfree':>7s} {'cost now':>8s} "
          f"{'best':>4s} {'kept':>4s} {'swap':>4s}  note")
    proposals = {}

    def cost_of(donors):
        dset = set(donors)
        return len((set().union(*(rel[s] for s in dset)) - dset) & countable) if dset else 0

    for name, ds in base.items():
        proposal = []
        missing = [s for s in ds if s not in pop]
        by_pop = defaultdict(list)
        for s in ds:
            if s in pop:
                by_pop[pop[s]].append(s)
        for p, cur in sorted(by_pop.items()):
            cur_set = set(cur)
            cands = {s for s in panel if pop[s] == p and s not in drop and s not in curated
                     and s not in near_curated and donor_of.get(s, name) == name}
            outside = [s for s in cur if s not in cands]
            # A donor kept as is still rules out its relatives as fellow donors.
            cands -= set().union(*(rel[s] for s in outside)) if outside else set()
            # Families: candidates and their relatives, joined through any shared relative.
            seen, families = set(), []
            for s in sorted(cands):
                if s in seen:
                    continue
                comp, stack = set(), [s]
                while stack:
                    x = stack.pop()
                    if x in comp:
                        continue
                    comp.add(x)
                    stack.extend(rel[x] - comp)
                seen |= comp
                families.append(family_options(comp & cands, rel, countable, cur_set))
            free = sum(1 for s in cands if not (rel[s] & countable))
            cur_out = set().union(*(rel[s] for s in cur)) - cur_set
            cur_cost = len(cur_out & countable)
            cur_with = sum(1 for s in cur if rel[s] & countable)
            cur_pairs = sum(1 for x, y in combinations(cur, 2) if y in rel[x])
            note = []
            if outside:
                note.append(f"{len(outside)} current donor(s) not eligible as candidates "
                            f"(e.g. {outside[0]}: in/related to a curated list, or a drop) -- kept as is")
            target = len(cur) - len(outside)
            sol = solve(families, target)
            if sol is None:
                note.append(f"no {target}-donor set without related donors exists; current kept")
                chosen = list(cur)
                best_cost, kept = cur_cost, len(cur)
            else:
                best_cost, kept, chosen = sol
                chosen = [s for s in cur if s in set(chosen)] + sorted(set(chosen) - cur_set) + outside
            if sol is not None and cost_of(chosen) >= cur_cost and not cur_pairs:
                chosen = list(cur)
                note.append("current set already as good")
            chosen_set = set(chosen)
            out_ids = [s for s in cur if s not in chosen_set]
            in_ids = [s for s in chosen if s not in cur_set]
            new_cost = cost_of(chosen)
            proposal.extend(chosen)
            for s in out_ids:
                chg.write(f"{name}\t{p}\tout\t{s}\t{','.join(sorted(rel[s] & countable)) or '-'}\n")
            for s in in_ids:
                chg.write(f"{name}\t{p}\tin\t{s}\t{','.join(sorted(rel[s] & countable)) or '-'}\n")
            summ.write(f"{name}\t{p}\t{len(cur)}\t{len(cands)}\t{free}\t{cur_cost}\t{cur_with}\t{cur_pairs}\t"
                       f"{new_cost}\t{len(cur) - len(out_ids)}\t{len(in_ids)}\t{'; '.join(note)}\n")
            print(f"{name:22s} {p:16s} {len(cur):6d} {len(cands):5d} {free:7d} {cur_cost:8d} "
                  f"{new_cost:4d} {len(cur) - len(out_ids):4d} {len(in_ids):4d}  {'; '.join(note)}"
                  + (f"; {cur_pairs} related donor pair(s) now" if cur_pairs else ""))
        if missing:
            print(f"{name:22s} {'(not in table)':16s} {len(missing):6d}  kept as is: {', '.join(missing[:5])}")
            proposal.extend(missing)
        proposals[name] = proposal
    for name, parts in derived.items():
        # Keep the union list's own order for donors that stay; new ones after.
        merged = set().union(*(set(proposals[m]) for m in parts))
        proposals[name] = [s for s in lists[name] if s in merged] \
            + [s for m in parts for s in proposals[m] if s not in set(lists[name])]
        print(f"{name:22s} {'(union)':16s} {len(lists[name]):6d}  = {' + '.join(parts)}: "
              f"{len(proposals[name])} donors, {len(set(proposals[name]) - set(lists[name]))} new")
    for name, proposal in proposals.items():
        with open(os.path.join(a.out_dir, f"{name}.proposed.txt"), "w") as g:
            g.write("\n".join(proposal) + "\n")
    summ.close(); chg.close()
    print(f"\nproposals -> {a.out_dir}/<list>.proposed.txt; table -> donor_plan_summary.tsv; "
          f"every swap and the relatives each donor would keep out -> donor_plan_changes.tsv")


if __name__ == "__main__":
    main()
