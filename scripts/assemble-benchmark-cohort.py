#!/usr/bin/env python3
"""
assemble-benchmark-cohort.py -- build steps 7 and 8 of BENCHMARK_COHORT_PLAN.md

Step 7: decide which real genomes join the cohort as TEST samples, and label
        each one real_homogeneous or real_admixed.
Step 8: deal every individual -- real and simulated -- into the six
        chip x build cells, balanced within 1.

This needs only the K=8 supervised Q table from 4a and the donor/reference
split. It does NOT need the simulator, so it runs the moment stage 3 lands,
which is the point: BENCHMARK_COHORT_PLAN.md section 5 says the SAS and OCE
decisions are "decided when stage 3 prints the counts", and this is what
prints them.

WHO IS ELIGIBLE AS A TEST SAMPLE, AND WHY EACH EXCLUSION IS THERE.
A real genome may be a test sample only if no method under evaluation has
already seen its haplotypes. Four separate ways that can happen, so four
exclusions:

  donor        its haplotypes were used to simulate admixed individuals
  reference    it is in the LAI / imputation reference panel
  1000 Genomes the imputation servers' own panels contain 1KG
  server panel SGDP re-sequenced individuals who ARE in 1KG under other ids,
               so a cohort column reading "SGDP" does not establish absence
               from a server panel. That is the donor-leak class that
               sgdp_server_panel_ids.txt exists to close, and it applies to
               test samples for exactly the same reason it applies to donors.

Dropping any one of the four silently inflates accuracy.

THE HOMOGENEITY THRESHOLD IS TIERED BY DEFAULT.
BENCHMARK_COHORT_PLAN.md section 3 writes a flat 0.93. The homogeneity panel
that 4a builds uses 0.95 for EUR/EAS/AFR and 0.93 for AMR/OCE/MEN/SAS, because
the deeply sampled super-pops can meet a stricter bar and the others cannot.
Classifying test genomes on a flat bar while the panel uses a tiered one would
put samples in real_homogeneous that the panel itself rejects as insufficiently
homogeneous. So the tiered thresholds are the default here; --flat-threshold
reproduces the plan document's literal reading.

Usage:
  python3 -I scripts/assemble-benchmark-cohort.py \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k8.tsv \\
      --refs  reference_ids --suffix b2 \\
      --out   results/benchmark_cohort
"""
import argparse
import os
import sys
from collections import defaultdict, OrderedDict

# Document 19, Table 1 (contract v2). Kept identical to the roster in
# cohort-composition-report.py: the two scripts must not disagree about what
# is being built.
MODELS = [
    # name,    components,            proportions,        gen, n,   priority
    ("LAT3",   ("NAT", "EUR", "AFR"), (0.50, 0.45, 0.05),  16, 600, "primary"),
    ("CAR3",   ("NAT", "EUR", "AFR"), (0.10, 0.50, 0.40),  16, 200, "secondary"),
    ("AFAM2",  ("AFR", "EUR"),        (0.80, 0.20),         7, 200, "secondary"),
    ("MENA3",  ("EUR", "MEN", "AFR"), (0.30, 0.50, 0.20),  30, 200, "secondary"),
    ("MENEUR", ("MEN", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    ("EASEUR", ("EAS", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    ("SASEUR", ("SAS", "EUR"),        (0.50, 0.50),         2, 200, "secondary"),
    ("NAT1",   ("NAT",),              (1.0,),               0,  60, "control"),
    ("EUR1",   ("EUR",),              (1.0,),               0,  60, "control"),
    ("AFR1",   ("AFR",),              (1.0,),               0,  60, "control"),
    ("EAS1",   ("EAS",),              (1.0,),               0,  60, "control"),
    ("SAS1",   ("SAS",),              (1.0,),               0,  60, "control"),
    ("MEN1",   ("MEN",),              (1.0,),               0,  60, "control"),
]

COMPONENT_POOL = {"NAT": "AMR", "EUR": "EUR", "AFR": "AFR",
                  "MEN": "MEN", "EAS": "EAS", "SAS": "SAS", "OCE": "OCE"}

CHIPS = ("snp6", "affy500k", "gsa")
BUILDS = ("hg19", "hg38")
CELLS = [(c, b) for c in CHIPS for b in BUILDS]   # six, fixed order

THRESH_STRICT = 0.95
STRICT_SUPERPOPS = ("EUR", "EAS", "AFR")
THRESH_RELAXED = 0.93


def read_tsv(path):
    """Header-keyed rows. Columns are resolved by NAME, never by offset: the
    panel table's width changes with K, so an offset that is right at K=5 is
    wrong at K=8 and fails silently by reading a neighbouring column."""
    if not os.path.exists(path):
        sys.exit(f"FATAL: missing {path}")
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(head)}
        rows = []
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            p = ln.split("\t")
            if len(p) < len(head):
                p += [""] * (len(head) - len(p))
            rows.append(p)
    return idx, rows


def read_list(path):
    if not path or not os.path.exists(path):
        return set()
    out = set()
    with open(path) as f:
        for ln in f:
            s = ln.split("#", 1)[0].strip()
            if s:
                out.add(s.split()[0])
    return out


def threshold_for(sp, flat):
    if flat is not None:
        return flat
    return THRESH_STRICT if sp in STRICT_SUPERPOPS else THRESH_RELAXED


def deal(groups):
    """Deal each group round-robin across the six cells, rotating the starting
    cell per group.

    Rotating matters. A fixed start would send every group's remainder to the
    same cell, so with 13 models and many real groups cell 0 would end up
    systematically larger -- and a cell-size imbalance is indistinguishable
    from a chip effect in the results. Rotating spreads the remainders, which
    keeps each group balanced within 1 AND the six totals balanced within 1.
    """
    out = {}
    for gi, (gname, members) in enumerate(sorted(groups.items())):
        off = gi % len(CELLS)
        for i, m in enumerate(sorted(members)):
            out[m] = CELLS[(i + off) % len(CELLS)]
    return out


def deal_counts(n, offset):
    """Split a simulated model's n across the six cells, balanced within 1."""
    base, rem = divmod(n, len(CELLS))
    counts = [base] * len(CELLS)
    for r in range(rem):
        counts[(r + offset) % len(CELLS)] += 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True,
                    help="lai_ref_panel_samples.k<K>.tsv from 4a")
    ap.add_argument("--refs", default="reference_ids",
                    help="directory holding <sp>_reference_<suffix>.txt and "
                         "<sp>_simulation_<suffix>.txt from the splitter")
    ap.add_argument("--suffix", default="b2")
    ap.add_argument("--server-panel-ids",
                    default="reference_ids/sgdp_server_panel_ids.txt",
                    help="ids inside an imputation server's reference panel; "
                         "these cannot be test samples either")
    ap.add_argument("--kg-cohort", default="1KG")
    ap.add_argument("--flat-threshold", type=float, default=None,
                    help="use one homogeneity threshold for every super-pop "
                         "(the plan document's literal 0.93) instead of the "
                         "tiered thresholds 4a applies")
    ap.add_argument("--out", default="results/benchmark_cohort",
                    help="output prefix; three files are written")
    a = ap.parse_args()

    idx, rows = read_tsv(a.panel)
    need = ["sample", "cohort", "pop", "super_pop",
            "max_Q_sup", "argmax_pop_sup", "included"]
    missing = [c for c in need if c not in idx]
    if missing:
        sys.exit(f"FATAL: {a.panel} lacks column(s) {missing}. Columns "
                 f"present: {sorted(idx)}")

    # ---- the four exclusions -------------------------------------------
    donors, refs = set(), set()
    sps = set(r[idx["super_pop"]] for r in rows) | set(COMPONENT_POOL.values())
    for sp in sorted(x for x in sps if x and x != "NA"):
        donors |= read_list(os.path.join(a.refs, f"{sp}_simulation_{a.suffix}.txt"))
        refs   |= read_list(os.path.join(a.refs, f"{sp}_reference_{a.suffix}.txt"))
    server = read_list(a.server_panel_ids)

    if not donors and not refs:
        sys.exit(f"FATAL: no donor or reference lists found under {a.refs} "
                 f"with suffix '{a.suffix}'. Run split-homog-donors-reference.py "
                 f"first -- without it every panel sample looks test-eligible, "
                 f"which would put donors in the test set.")

    overlap = donors & refs
    if overlap:
        sys.exit(f"FATAL: {len(overlap)} sample(s) are BOTH donor and "
                 f"reference, e.g. {sorted(overlap)[:5]}. The split is broken; "
                 f"refusing to build a cohort on it.")

    # ---- step 7: classify the real test genomes ------------------------
    test, excl = [], defaultdict(list)
    for r in rows:
        s = r[idx["sample"]]
        if r[idx["included"]].strip() not in ("1", "TRUE", "True", "true"):
            excl["not in panel"].append(s); continue
        if s in donors:
            excl["donor"].append(s); continue
        if s in refs:
            excl["reference panel"].append(s); continue
        if r[idx["cohort"]].strip() == a.kg_cohort:
            excl["1000 Genomes (in server panels)"].append(s); continue
        if s in server:
            excl["inside a server reference panel"].append(s); continue
        try:
            q = float(r[idx["max_Q_sup"]])
        except ValueError:
            excl["unparseable max_Q_sup"].append(s); continue
        sp = r[idx["super_pop"]]
        thr = threshold_for(sp, a.flat_threshold)
        test.append({
            "sample": s,
            "cohort": r[idx["cohort"]],
            "pop": r[idx["pop"]],
            "super_pop": sp,
            "argmax_pop_sup": r[idx["argmax_pop_sup"]],
            "max_Q_sup": q,
            "threshold": thr,
            "sample_type": "real_homogeneous" if q >= thr else "real_admixed",
        })

    # ---- step 8: deal into the six chip x build cells ------------------
    groups = defaultdict(list)
    for t in test:
        groups[f"{t['sample_type']}|{t['argmax_pop_sup']}"].append(t["sample"])
    assign = deal(groups)
    for t in test:
        t["chip"], t["build"] = assign[t["sample"]]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    man = f"{a.out}_manifest.tsv"
    with open(man, "w") as g:
        cols = ["sample", "cohort", "pop", "super_pop", "argmax_pop_sup",
                "max_Q_sup", "threshold", "sample_type", "chip", "build"]
        g.write("\t".join(cols) + "\n")
        for t in sorted(test, key=lambda x: x["sample"]):
            g.write("\t".join(str(t[c]) for c in cols) + "\n")

    ros = f"{a.out}_simset_roster.tsv"
    with open(ros, "w") as g:
        g.write("model\tcomponents\tproportions\tgenerations\tpriority"
                "\tchip\tbuild\tn\n")
        for mi, (name, comp, prop, gen, n, pri) in enumerate(MODELS):
            counts = deal_counts(n, mi)
            for (chip, build), k in zip(CELLS, counts):
                g.write(f"{name}\t{','.join(comp)}"
                        f"\t{','.join(f'{p:g}' for p in prop)}"
                        f"\t{gen}\t{pri}\t{chip}\t{build}\t{k}\n")

    # Per-individual rows, not just counts. Step 6 simulates a genome per
    # model and step 10 masks it to a chip x build cell, so the two steps must
    # agree about which simulated individual is in which cell. Counts alone
    # cannot express that, and letting each step deal independently would put
    # the same id in different cells in the two steps.
    ind = f"{a.out}_simset_individuals.tsv"
    with open(ind, "w") as g:
        g.write("sim_id\tmodel\tcomponents\tproportions\tgenerations"
                "\tpriority\tchip\tbuild\n")
        for mi, (name, comp, prop, gen, n, pri) in enumerate(MODELS):
            k = 0
            for (chip, build), c in zip(CELLS, deal_counts(n, mi)):
                for _ in range(c):
                    k += 1
                    g.write(f"{name}_{k:04d}\t{name}\t{','.join(comp)}"
                            f"\t{','.join(f'{p:g}' for p in prop)}"
                            f"\t{gen}\t{pri}\t{chip}\t{build}\n")

    # ---- the counts, and the two decisions section 5 defers ------------
    rep = f"{a.out}_counts.txt"
    L = []
    w = L.append
    w("=" * 72)
    w(" BENCHMARK COHORT -- composition after ADMIXTURE K=8")
    w("=" * 72)
    w(f" panel       {a.panel}")
    w(f" split       {a.refs} (suffix {a.suffix})")
    w(f" threshold   " + ("flat %.2f" % a.flat_threshold if a.flat_threshold
                          else "tiered: %.2f for %s, %.2f otherwise"
                               % (THRESH_STRICT, "/".join(STRICT_SUPERPOPS),
                                  THRESH_RELAXED)))
    w("")
    w(" REAL TEST GENOMES")
    w("-" * 72)
    by = defaultdict(lambda: defaultdict(int))
    for t in test:
        by[t["argmax_pop_sup"]][t["sample_type"]] += 1
    w(f" {'component':<12}{'homogeneous':>13}{'admixed':>10}{'total':>8}")
    for c in sorted(by):
        h = by[c]["real_homogeneous"]; m = by[c]["real_admixed"]
        w(f" {c:<12}{h:>13}{m:>10}{h+m:>8}")
    th = sum(1 for t in test if t["sample_type"] == "real_homogeneous")
    tm = len(test) - th
    w(f" {'TOTAL':<12}{th:>13}{tm:>10}{len(test):>8}")
    w("")
    w(" Real genomes are scored for imputation accuracy, calibration and")
    w(" retention ONLY. They have no recorded tracts, so they are never")
    w(" scored for local-ancestry accuracy.")
    w("")
    w(" EXCLUDED FROM THE TEST SET")
    w("-" * 72)
    for k in sorted(excl):
        w(f" {k:<38}{len(excl[k]):>6}")
    w("")
    w(" CHIP x BUILD CELLS")
    w("-" * 72)
    cc = defaultdict(int)
    for t in test:
        cc[(t["chip"], t["build"])] += 1
    sim_cc = defaultdict(int)
    for mi, (name, _c, _p, _g, n, _pr) in enumerate(MODELS):
        for cell, k in zip(CELLS, deal_counts(n, mi)):
            sim_cc[cell] += k
    w(f" {'cell':<18}{'simulated':>11}{'real':>7}{'total':>8}")
    for cell in CELLS:
        w(f" {cell[0] + '/' + cell[1]:<18}{sim_cc[cell]:>11}"
          f"{cc[cell]:>7}{sim_cc[cell] + cc[cell]:>8}")
    w(f" {'TOTAL':<18}{sum(sim_cc.values()):>11}"
      f"{sum(cc.values()):>7}{sum(sim_cc.values()) + sum(cc.values()):>8}")
    spread = max(cc.values()) - min(cc.values()) if cc else 0
    w("")
    w(f" real-sample spread across cells: {spread} "
      f"({'balanced' if spread <= 1 else 'NOT BALANCED -- investigate'})")
    w("")
    w(" THE TWO DECISIONS SECTION 5 DEFERS TO THESE COUNTS")
    w("-" * 72)
    dn = defaultdict(int)
    for sp in sorted(set(COMPONENT_POOL.values())):
        dn[sp] = len(read_list(os.path.join(a.refs,
                                            f"{sp}_simulation_{a.suffix}.txt")))
    need_by_comp = defaultdict(int)
    for name, comp, _p, _g, n, _pr in MODELS:
        for c in comp:
            need_by_comp[c] += n
    w(f" {'component':<10}{'donors':>8}{'haps':>7}{'indiv needing it':>18}"
      f"{'haps/indiv':>12}")
    for c in sorted(need_by_comp):
        pool = COMPONENT_POOL.get(c, c)
        d = dn.get(pool, 0)
        ratio = (2.0 * d / need_by_comp[c]) if need_by_comp[c] else 0.0
        w(f" {c:<10}{d:>8}{2 * d:>7}{need_by_comp[c]:>18}{ratio:>12.3f}")
    w("")
    sas_d = dn.get("SAS", 0)
    w(f" SAS: {sas_d} donors -> {2 * sas_d} haplotypes for 200 SASEUR + 60 SAS1.")
    w("      Options, in the plan's preference order: reduce SASEUR's n;")
    w("      accept the reuse and report the kinship distribution; or drop")
    w("      SASEUR to control-only and rely on real admixed SAS for the")
    w("      imputation readings. NOT DEFAULTED -- decide and record it.")
    w("")
    oce_h = by.get("OCE", {}).get("real_homogeneous", 0)
    oce_m = by.get("OCE", {}).get("real_admixed", 0)
    w(f" OCE: {dn.get('OCE', 0)} donors; no model in Table 1 uses OCE.")
    w(f"      {oce_h} real_homogeneous + {oce_m} real_admixed OCE test samples")
    w("      carry the OCE reading instead. Confirm OCE stays unsimulated.")
    w("")
    w("=" * 72)
    w(" NOT BUILT BY THIS SCRIPT: step 6 (simulate the 2,160), step 10")
    w(" (write chip files per chip x build x chr), step 11 (tables,")
    w(" checksums, provenance). Steps 6 and 10 run in")
    w(" genotyping_pipeline/benchmarking/simulate/.")
    w("=" * 72)

    with open(rep, "w") as g:
        g.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwrote {man}\n      {ros}\n      {ind}\n      {rep}")


if __name__ == "__main__":
    main()
