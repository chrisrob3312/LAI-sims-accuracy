#!/usr/bin/env python3
"""
cohort-power-summary.py -- how much sample size the benchmark cohort actually
gives a power or GWAS simulation, once relatedness is accounted for.

WHY THIS EXISTS. ~500 real genomes is too few for a power calculation, so the
simulated 2,160 have to be usable for it too -- about 2,660 together. But
simulated individuals drawn from a finite donor pool are cryptic relatives of
each other: two individuals built from the same donor haplotypes share long
stretches by construction, and the thinner the donor pool the worse it gets.
SAS is the binding case, ~24 donors for 200 SASEUR individuals.

That relatedness does NOT have to be a bias. The GWAS uses GENESIS with a GRM,
which models relatedness rather than pruning it, and power simulations run
through the same machinery inherit that. What relatedness does cost is
EFFECTIVE sample size, and that is a number rather than a caveat:

    N_eff  =  N^2 / sum_ij R_ij        R = the GRM (2*kinship off-diagonal)

For unrelated individuals sum_ij R_ij = N and N_eff = N. Every related pair
adds to the denominator, so N_eff falls below N by exactly the amount the
duplication costs. A power calculation should use N_eff, not the roster count,
and this prints both per group so a thin-donor model cannot quietly contribute
200 nominal individuals that behave like 40.

IT DOES NOT COMPUTE KINSHIP. Kinship needs genotypes, so it reads a KING or
plink2 king-table for the combined cohort (6e_cohort_kinship_clm.sh produces
one). Pairs absent from the table are taken as unrelated, which is how both
tools report: they emit pairs above a threshold and omit the rest.

Usage:
  python3 -I scripts/cohort-power-summary.py \\
      --kinship results/cohort_king.kin0 \\
      --manifest results/benchmark_cohort_manifest.tsv \\
      --individuals results/benchmark_cohort_simset_individuals.tsv \\
      --out results/cohort_power_summary.txt
"""
import argparse
import os
import sys
from collections import defaultdict

# KING-robust thresholds, as used by 1e_king_dedup_clm.sh.
KING_BANDS = [
    (0.354, "duplicate / MZ"),
    (0.177, "first degree"),
    (0.0884, "second degree"),
    (0.0442, "third degree"),
]


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None, []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        rows = []
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            p = ln.split("\t")
            if len(p) < len(head):
                p += [""] * (len(head) - len(p))
            rows.append(p)
    return ix, rows


def read_kinship(path):
    """Pairs from a KING .kin0 or plink2 --make-king-table output.

    Both are whitespace-separated with a header. Column names differ between
    tools and versions, so they are resolved BY NAME from a set of known
    spellings rather than by position -- a positional guess would silently
    read the wrong column and report everyone unrelated, which is the
    optimistic direction.
    """
    if not os.path.exists(path):
        sys.exit(f"FATAL: missing {path}\n"
                 f"Run 6e_cohort_kinship_clm.sh to produce it. Without a "
                 f"kinship table this would assume the cohort is unrelated, "
                 f"which is the optimistic direction and exactly the "
                 f"assumption the simulated half breaks.")
    id1_names = ("ID1", "IID1", "#IID1", "FID1")
    id2_names = ("ID2", "IID2", "IID2", "FID2")
    kin_names = ("KINSHIP", "Kinship", "KIN", "PropIBD", "PI_HAT")
    pairs = []
    with open(path) as f:
        head = None
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            p = ln.split()
            if head is None:
                head = p
                def find(cands):
                    for c in cands:
                        if c in head:
                            return head.index(c)
                    return -1
                i1, i2 = find(id1_names), find(id2_names)
                ik = find(kin_names)
                if min(i1, i2, ik) < 0:
                    sys.exit(f"FATAL: {path} header {head} has no recognised "
                             f"id/kinship columns. Expected one of "
                             f"{id1_names} / {id2_names} / {kin_names}.")
                continue
            if len(p) <= max(i1, i2, ik):
                continue
            try:
                k = float(p[ik])
            except ValueError:
                continue
            pairs.append((p[i1], p[i2], k))
    return pairs


def n_eff(ids, kin_by_pair):
    """N^2 / sum_ij R_ij, with R = 2*kinship off-diagonal and 1 on the diagonal.

    Only pairs present in the kinship table contribute; the rest are taken as
    unrelated, which is what both tools' thresholded output means.
    """
    n = len(ids)
    if n == 0:
        return 0.0, 0.0
    s = float(n)                      # the diagonal
    idset = set(ids)
    rel = 0.0
    for (a, b), k in kin_by_pair.items():
        if a in idset and b in idset:
            rel += 2.0 * (2.0 * k)    # both (i,j) and (j,i), R_ij = 2*phi
    s += rel
    return (n * n / s) if s > 0 else 0.0, rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kinship", required=True,
                    help="KING .kin0 or plink2 --make-king-table output for "
                         "the COMBINED cohort")
    ap.add_argument("--manifest", default="",
                    help="benchmark_cohort_manifest.tsv (the real half)")
    ap.add_argument("--individuals", default="",
                    help="benchmark_cohort_simset_individuals.tsv")
    ap.add_argument("--min-neff-ratio", type=float, default=0.5,
                    help="flag a group whose N_eff falls below this fraction "
                         "of its nominal N (default 0.5)")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    pairs = read_kinship(a.kinship)
    kin_by_pair = {}
    for x, y, k in pairs:
        key = (x, y) if x <= y else (y, x)
        kin_by_pair[key] = max(kin_by_pair.get(key, -1.0), k)

    groups = defaultdict(list)   # label -> ids
    half = {}                    # id -> "real" | "simulated"
    mix, mrows = read_tsv(a.manifest)
    if mix is not None:
        for p in mrows:
            s = p[mix["sample"]]
            lab = f"real:{p[mix['sample_type']]}"
            groups[lab].append(s)
            half[s] = "real"
    six, srows = read_tsv(a.individuals)
    if six is not None:
        for p in srows:
            s = p[six["sim_id"]]
            groups[f"sim:{p[six['model']]}"].append(s)
            half[s] = "simulated"
    if not groups:
        sys.exit("FATAL: no cohort members found. Pass --manifest and/or "
                 "--individuals.")

    L = []
    w = L.append
    w("=" * 78)
    w(" COHORT SAMPLE SIZE FOR POWER / GWAS SIMULATIONS")
    w("=" * 78)
    w("")
    w(" Both halves are used. ~500 real genomes alone is too few for a power")
    w(" calculation; the simulated 2,160 make the total workable. The cost is")
    w(" that simulated individuals built from a shared donor pool are cryptic")
    w(" relatives, so the nominal count overstates the information. Use N_eff.")
    w("")
    w(" N_eff = N^2 / sum_ij R_ij, R = GRM (2*kinship off-diagonal, 1 on the")
    w(" diagonal). Unrelated gives N_eff = N; every related pair lowers it.")
    w("")
    w(" The GWAS uses GENESIS with a GRM, which MODELS relatedness rather than")
    w(" pruning it, so this relatedness is handled rather than biasing the")
    w(" estimate -- provided the power simulations use the same GRM. What it")
    w(" still costs is the information N_eff measures.")
    w("")

    w(" KINSHIP PAIRS IN THE COMBINED COHORT")
    w("-" * 78)
    cohort_ids = set(half)
    inside = {k: v for k, v in kin_by_pair.items()
              if k[0] in cohort_ids and k[1] in cohort_ids}
    w(f" pairs in the kinship table: {len(kin_by_pair)}")
    w(f" pairs with both members in the cohort: {len(inside)}")
    w("")
    w(f" {'band':<22}{'pairs':>8}  {'of which':>10}")
    prev = 10.0
    for thr, name in KING_BANDS:
        sel = [k for k, v in inside.items() if thr <= v < prev]
        n_rr = sum(1 for k in sel if half.get(k[0]) == "real"
                   and half.get(k[1]) == "real")
        n_ss = sum(1 for k in sel if half.get(k[0]) == "simulated"
                   and half.get(k[1]) == "simulated")
        n_rs = len(sel) - n_rr - n_ss
        w(f" {name + f' (>={thr})':<22}{len(sel):>8}  "
          f"real-real {n_rr}, sim-sim {n_ss}, cross {n_rs}")
        prev = thr
    below = [k for k, v in inside.items() if v < KING_BANDS[-1][0]]
    w(f" {'below third degree':<22}{len(below):>8}")
    w("")
    # A real-simulated pair above second degree is a LEAK, not just
    # relatedness: it means a donor's own genome is in the test set.
    leak = [k for k, v in inside.items() if v >= 0.0884
            and {half.get(k[0]), half.get(k[1])} == {"real", "simulated"}]
    if leak:
        w(f" *** {len(leak)} real-simulated pair(s) at second degree or closer.")
        w(" That is not ordinary relatedness: a real test genome related to a")
        w(" simulated one means a donor, or a donor's close relative, is in")
        w(" the test set. Check the donor/test exclusions before scoring.")
        for k in leak[:5]:
            w(f"      {k[0]}  {k[1]}  kinship {inside[k]:.4f}")
        w("")

    w(" EFFECTIVE SAMPLE SIZE BY GROUP")
    w("-" * 78)
    w(f" {'group':<22}{'N':>7}{'N_eff':>9}{'ratio':>8}{'rel pairs':>11}")
    flagged = []
    for lab in sorted(groups):
        ids = groups[lab]
        ne, _rel = n_eff(ids, inside)
        idset = set(ids)
        npair = sum(1 for k in inside if k[0] in idset and k[1] in idset)
        ratio = ne / len(ids) if ids else 0.0
        mark = ""
        if ratio < a.min_neff_ratio:
            mark = "  <<<"
            flagged.append((lab, len(ids), ne, ratio))
        w(f" {lab:<22}{len(ids):>7}{ne:>9.1f}{ratio:>7.0%}{npair:>11}{mark}")
    all_ids = [s for ids in groups.values() for s in ids]
    ne_all, _ = n_eff(all_ids, inside)
    w("-" * 78)
    w(f" {'COMBINED':<22}{len(all_ids):>7}{ne_all:>9.1f}"
      f"{(ne_all / len(all_ids) if all_ids else 0):>7.0%}")
    real_ids = [s for s, h in half.items() if h == "real"]
    sim_ids = [s for s, h in half.items() if h == "simulated"]
    ne_r, _ = n_eff(real_ids, inside)
    ne_s, _ = n_eff(sim_ids, inside)
    w(f" {'  real only':<22}{len(real_ids):>7}{ne_r:>9.1f}")
    w(f" {'  simulated only':<22}{len(sim_ids):>7}{ne_s:>9.1f}")
    w("")
    w(" Use COMBINED N_eff in the power calculation, not the roster count.")

    if flagged:
        w("")
        w(" GROUPS WHOSE EFFECTIVE SIZE IS WELL BELOW NOMINAL")
        w("-" * 78)
        for lab, n, ne, ratio in flagged:
            w(f"   {lab}: {n} individuals behave like {ne:.0f} "
              f"({ratio:.0%}) -- donor reuse")
        w("")
        w(" This is the SAS decision in BENCHMARK_COHORT_PLAN.md section 5")
        w(" showing up as a number. The options there were: reduce the model's")
        w(" n, accept the reuse and report this distribution, or drop the")
        w(" model to control-only. Reporting it is now possible, so accepting")
        w(" the reuse is defensible -- as long as the power calculation uses")
        w(" N_eff and the paper states it.")

    text = "\n".join(L) + "\n"
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as g:
            g.write(text)
    print(text, end="")
    # Low N_eff is a finding to report, not an error: it is the expected
    # consequence of a thin donor pool and the decision it informs is
    # Christina's. A real-simulated leak IS an error.
    sys.exit(2 if leak else 0)


if __name__ == "__main__":
    main()
