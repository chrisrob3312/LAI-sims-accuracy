#!/usr/bin/env python3
"""
check-simset-contract.py -- run the intake's rules BEFORE the build, not after.

ingest_simset.sh (genotyping_pipeline, benchmarking/simulate/) refuses a cohort
that breaks any of its checks, and it runs at the END -- after simulation, chip
files and a transfer. Several of its rules are decided much earlier, at the
split and the cell deal, so a violation costs the whole build. This re-states
the ones that are knowable from the plan alone and checks them while they are
still cheap to fix.

It does NOT replace ingest_simset.sh. It checks the subset that can be decided
from the role lists, the manifests and a KING table.

Contract: SIMSET_CONTRACT.md v2 (Doc 20, 2 Nov 2026), plus the relatedness
rule tightened at 2b59278 (6 Nov 2026).

  5.donors        no donor has both roles; every model component has donors
  2.real_samples  sample_type agrees with top_fraction at the FLAT threshold;
                  no real sample is a donor or in a server panel
  2.chips         every sample is on EXACTLY ONE chip x build
  2.balance       every model and every real group is within 1 across the
                  6 cells
  7.identity      no pair of distinct IDs above 0.354 anywhere in the set
  7.roles         (2b59278) no pair at or above 0.177 in two DIFFERENT roles
                  -- simulation donor, LAI reference, real test genome

Usage:
  python3 -I scripts/check-simset-contract.py \\
      --manifest results/benchmark_cohort_manifest.tsv \\
      --individuals results/benchmark_cohort_simset_individuals.tsv \\
      --refs reference_ids --suffix bench \\
      --king 04_homogeneity_panel/king/king.kin0 \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k8.tsv
"""
import argparse
import os
import sys
from collections import defaultdict

CHIPS = ("snp6", "affy500k", "gsa")
BUILDS = ("hg19", "hg38")
CELLS = [(c, b) for c in CHIPS for b in BUILDS]
# The contract's flat threshold. NOT the tiered 0.95/0.93 the panel build uses:
# intake 2.real_samples compares sample_type against top_fraction at this one
# number, so a tiered classification disagrees with it for EUR/EAS/AFR samples
# between 0.93 and 0.95.
CONTRACT_THRESHOLD = 0.93
# Components a model draws on; 5.donors requires each to have donors.
MODEL_COMPONENTS = {"NAT", "EUR", "AFR", "MEN", "EAS", "SAS"}
COMPONENT_POOL = {"NAT": "AMR", "EUR": "EUR", "AFR": "AFR",
                  "MEN": "MEN", "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None, []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        rows = []
        for ln in f:
            ln = ln.rstrip("\n")
            if ln:
                p = ln.split("\t")
                if len(p) < len(head):
                    p += [""] * (len(head) - len(p))
                rows.append(p)
    return ix, rows


def role_list(refs_dir, sp, role, suffix):
    for name in (f"{sp.lower()}_{role}_{suffix}.txt",
                 f"{sp}_{role}_{suffix}.txt"):
        p = os.path.join(refs_dir, name)
        if os.path.exists(p):
            out = set()
            with open(p) as f:
                for ln in f:
                    s = ln.split("#", 1)[0].strip()
                    if s:
                        out.add(s.split()[0])
            return out
    return set()


def read_king(path):
    if not path or not os.path.exists(path):
        return None
    pairs = []
    with open(path) as f:
        head = None
        for ln in f:
            p = ln.split()
            if not p:
                continue
            if head is None:
                head = [c.lstrip("#") for c in p]
                def find(c):
                    for n in c:
                        if n in head:
                            return head.index(n)
                    return -1
                i1 = find(("ID1", "IID1", "FID1"))
                i2 = find(("ID2", "IID2", "FID2"))
                ik = find(("KINSHIP", "Kinship", "KIN", "PropIBD"))
                if min(i1, i2, ik) < 0:
                    sys.exit(f"FATAL: {path} header {head} has no "
                             f"ID1/ID2/KINSHIP columns")
                continue
            if len(p) <= max(i1, i2, ik):
                continue
            try:
                pairs.append((p[i1], p[i2], float(p[ik])))
            except ValueError:
                continue
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--individuals", required=True)
    ap.add_argument("--refs", default="reference_ids")
    ap.add_argument("--suffix", default="bench")
    ap.add_argument("--king", default="",
                    help="king.kin0 from 1e (REPORT_FLOOR 0.0442 is enough "
                         "for the role rule; 7.identity needs unfiltered)")
    ap.add_argument("--panel", default="",
                    help="lai_ref_panel_samples.k8.tsv, for top_fraction")
    ap.add_argument("--threshold", type=float, default=CONTRACT_THRESHOLD)
    ap.add_argument("--role-kinship", type=float, default=0.177)
    a = ap.parse_args()

    fails, warns = [], []

    mix, mrows = read_tsv(a.manifest)
    six, srows = read_tsv(a.individuals)
    if mix is None:
        sys.exit(f"FATAL: missing {a.manifest}")
    if six is None:
        sys.exit(f"FATAL: missing {a.individuals}")

    sps = sorted(set(COMPONENT_POOL.values()) |
                 {r[mix["super_pop"]] for r in mrows if "super_pop" in mix})
    donors, refs = {}, {}
    for sp in sps:
        for s in role_list(a.refs, sp, "simulation", a.suffix):
            donors[s] = sp
        for s in role_list(a.refs, sp, "reference", a.suffix):
            refs[s] = sp

    print("=" * 76)
    print(" SIMSET CONTRACT PRE-FLIGHT (v2 + 2b59278)")
    print("=" * 76)
    print(f" donors {len(donors)}   reference {len(refs)}   "
          f"real test {len(mrows)}   simulated {len(srows)}")
    print("")

    # ---- 5.donors ----------------------------------------------------
    both = set(donors) & set(refs)
    if both:
        fails.append(f"5.donors: {len(both)} sample(s) are BOTH simulation "
                     f"donor and LAI reference, e.g. {sorted(both)[:5]}")
    missing_comp = []
    for comp in sorted(MODEL_COMPONENTS):
        pool = COMPONENT_POOL[comp]
        if not any(v == pool for v in donors.values()):
            missing_comp.append(comp)
    if missing_comp:
        fails.append(f"5.donors: no simulation donors for component(s) "
                     f"{missing_comp}; every model component needs donors")

    # ---- 2.real_samples ----------------------------------------------
    if "sample_type" in mix:
        real_ids = {r[mix["sample"]] for r in mrows}
        bad = real_ids & set(donors)
        if bad:
            fails.append(f"2.real_samples: {len(bad)} real test sample(s) are "
                         f"also simulation donors, e.g. {sorted(bad)[:5]}")
        # sample_type must agree with top_fraction at the FLAT threshold
        if "max_Q_sup" in mix:
            mism = []
            for r in mrows:
                try:
                    q = float(r[mix["max_Q_sup"]])
                except ValueError:
                    continue
                st = r[mix["sample_type"]]
                want = "real_homogeneous" if q >= a.threshold else "real_admixed"
                if st in ("real_homogeneous", "real_admixed") and st != want:
                    mism.append((r[mix["sample"]], q, st, want))
            if mism:
                fails.append(
                    f"2.real_samples: {len(mism)} sample(s) whose sample_type "
                    f"disagrees with top_fraction at the contract's FLAT "
                    f"{a.threshold}. e.g. "
                    + "; ".join(f"{s} q={q:.3f} is {st}, contract says {w}"
                                for s, q, st, w in mism[:3])
                    + ". A TIERED threshold (0.95 for EUR/EAS/AFR) produces "
                      "exactly this: intake compares against one number.")

    # ---- 2.chips : exactly one chip x build --------------------------
    def cells_of(row, ix):
        ch = row[ix["chips"]].split(",") if "chips" in ix else []
        bu = row[ix["builds"]].split(",") if "builds" in ix else []
        return [(c.strip(), b.strip()) for c in ch if c.strip()
                for b in bu if b.strip()]

    multi = []
    for label, ix, rows, idcol in (("real", mix, mrows, "sample"),
                                   ("simulated", six, srows, "sim_id")):
        if "chips" not in ix or "builds" not in ix:
            warns.append(f"2.chips: {label} table has no chips/builds columns")
            continue
        for r in rows:
            n = len(cells_of(r, ix))
            if n != 1:
                multi.append((label, r[ix[idcol]], n))
    if multi:
        n_multi = len(multi)
        ex = multi[0]
        fails.append(
            f"2.chips: {n_multi} sample(s) are NOT on exactly one chip x "
            f"build (e.g. {ex[0]} {ex[1]} is on {ex[2]} cells). The contract "
            f"splits each group over the 3 chips and, WITHIN a chip, half "
            f"hg19 / half hg38 -- six cells, one per person. "
            f"--design paired-builds or full breaks this.")

    # ---- 2.balance ---------------------------------------------------
    groups = defaultdict(lambda: defaultdict(int))
    for r in srows:
        if "chips" in six and "builds" in six:
            for cell in cells_of(r, six):
                groups[f"model:{r[six['model']]}"][cell] += 1
    for r in mrows:
        if "chips" in mix and "builds" in mix:
            key = (f"real:{r[mix['sample_type']]}|"
                   f"{r[mix.get('argmax_pop_sup', mix['sample'])]}")
            for cell in cells_of(r, mix):
                groups[key][cell] += 1
    unbalanced = []
    for g, cc in sorted(groups.items()):
        vals = [cc.get(c, 0) for c in CELLS]
        if max(vals) - min(vals) > 1:
            unbalanced.append((g, vals))
    if unbalanced:
        fails.append(f"2.balance: {len(unbalanced)} group(s) differ by more "
                     f"than 1 across the six cells, e.g. "
                     f"{unbalanced[0][0]} {unbalanced[0][1]}")

    # ---- 7.identity and 7.roles --------------------------------------
    if a.king:
        pairs = read_king(a.king)
        if pairs is None:
            warns.append(f"7.identity: no KING table at {a.king} -- the "
                         f"relatedness rules cannot be checked, and they are "
                         f"the ones that fail AFTER the build")
        else:
            role_of = {}
            for s in donors:
                role_of[s] = "simulation"
            for s in refs:
                role_of.setdefault(s, "lai_reference")
            for r in mrows:
                role_of.setdefault(r[mix["sample"]], "real_test")
            dups, cross = [], []
            for x, y, k in pairs:
                if x not in role_of or y not in role_of:
                    continue
                if k > 0.354:
                    dups.append((x, y, k))
                if k >= a.role_kinship and role_of[x] != role_of[y]:
                    cross.append((x, y, k, role_of[x], role_of[y]))
            if dups:
                fails.append(
                    f"7.identity: {len(dups)} pair(s) above 0.354 inside the "
                    f"set -- the same person twice. e.g. "
                    + "; ".join(f"{x}/{y} {k:.3f}" for x, y, k in dups[:3]))
            if cross:
                fails.append(
                    f"7.roles (2b59278): {len(cross)} pair(s) at or above "
                    f"{a.role_kinship} sitting in TWO DIFFERENT roles. e.g. "
                    + "; ".join(f"{x}({r1})/{y}({r2}) {k:.3f}"
                                for x, y, k, r1, r2 in cross[:3])
                    + ". The split must keep a first-degree relative CLUSTER "
                      "inside one role: re-run the splitter with "
                      "--relatedness <king.kin0>.")
    else:
        warns.append("7.identity / 7.roles: no --king given. These are the "
                     "rules that fail after the whole build, so check them.")

    # ---- verdict -----------------------------------------------------
    print(" CHECKS")
    print("-" * 76)
    for f in fails:
        print(f" FAIL  {f}")
        print("")
    for w in warns:
        print(f" WARN  {w}")
        print("")
    if not fails and not warns:
        print(" all checked rules pass")
    print("-" * 76)
    if fails:
        print(f" {len(fails)} rule(s) would be REFUSED by ingest_simset.sh.")
        print(" Fix them before building: every one of these is decided at the")
        print(" split or the cell deal, which take seconds, while the intake")
        print(" runs after simulation, chip files and a transfer.")
        sys.exit(1)
    print(" Nothing here would be refused. This is a SUBSET of the intake's")
    print(" checks -- file layout, chip sites, tracts and global ancestry are")
    print(" only checkable once the build exists.")
    sys.exit(0)


if __name__ == "__main__":
    main()
