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


def role_list(refs_dir, sp, role, suffix):
    """Read <sp>_<role>_<suffix>.txt, trying both letter cases.

    split-homog-donors-reference.py lowercases the super-pop in the filename
    (eur_reference_bench.txt). Looking only for the upper-case spelling finds
    nothing, every panel sample then looks test-eligible, and the cohort would
    be built with donors in the test set. Try both rather than depend on which
    convention the splitter happens to use.
    """
    for name in (f"{sp.lower()}_{role}_{suffix}.txt",
                 f"{sp}_{role}_{suffix}.txt",
                 f"{sp.upper()}_{role}_{suffix}.txt"):
        p = os.path.join(refs_dir, name)
        if os.path.exists(p):
            return read_list(p)
    return set()


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
    ap.add_argument("--panel-exclude", default="reference_ids/panel_exclude_ids.txt",
                    help="named exclusions (e.g. the CHMI hydatidiform-mole "
                         "cell line); not test-eligible")
    ap.add_argument("--dup-drop", default="reference_ids/panel_duplicate_drop_ids.txt",
                    help="KING duplicate drops from 1e; not test-eligible")
    ap.add_argument("--thin-reference", type=int, default=25,
                    help="a super-pop whose benchmark reference list is "
                         "smaller than this is marked reference_support=thin")
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
        donors |= role_list(a.refs, sp, "simulation", a.suffix)
        refs   |= role_list(a.refs, sp, "reference",  a.suffix)
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
    #
    # WHY `included` IS NOT THE ELIGIBILITY TEST.
    # 4a sets included=1 only when max_Q >= the homogeneity threshold AND the
    # supervised argmax matches the assigned super-pop. Gating on it therefore
    # admits ONLY homogeneous samples, which makes real_admixed empty by
    # construction -- and real_admixed is the category that contributes real
    # admixture with the ancestry/genotype correlation structure a simulation
    # idealises away. It also drops every OTH sample, since 4a excludes the
    # leftover cluster from the homogeneous panel on purpose.
    #
    # Those are exactly the people a harmonization pipeline has to handle: the
    # multiracial individuals who fit no single reference ancestry, and the
    # thinly represented groups that cannot be pooled into a stratified
    # analysis. Excluding them would test the pipeline only on the cases it
    # finds easy.
    #
    # But included=0 is overloaded -- it also marks donors, super_pop=NA, and
    # (via 4a's exclusion set) named exclusions and KING duplicate drops, and
    # the TSV does not record which branch fired. So eligibility is decided
    # from the role lists directly, and `included` is carried as a reported
    # column rather than used as a filter.
    panel_excl = read_list(a.panel_exclude)
    dup_drop = read_list(a.dup_drop)
    ref_sizes = {}
    for sp in sorted(x for x in sps if x and x != "NA"):
        ref_sizes[sp] = len(role_list(a.refs, sp, "reference", a.suffix))

    def support_for(sp):
        n = ref_sizes.get(sp, 0)
        if n == 0:
            return "none"
        return "thin" if n < a.thin_reference else "full"

    test, excl = [], defaultdict(list)
    for r in rows:
        s = r[idx["sample"]]
        if s in donors:
            excl["donor"].append(s); continue
        if s in refs:
            excl["reference panel"].append(s); continue
        if r[idx["cohort"]].strip() == a.kg_cohort:
            excl["1000 Genomes (in server panels)"].append(s); continue
        if s in server:
            excl["inside a server reference panel"].append(s); continue
        if s in panel_excl:
            excl["named panel exclusion"].append(s); continue
        if s in dup_drop:
            excl["KING duplicate drop"].append(s); continue
        sp = r[idx["super_pop"]]
        if sp in ("", "NA"):
            excl["super_pop unresolved"].append(s); continue
        try:
            q = float(r[idx["max_Q_sup"]])
        except ValueError:
            excl["unparseable max_Q_sup"].append(s); continue
        ap = r[idx["argmax_pop_sup"]]
        thr = threshold_for(sp, a.flat_threshold)
        # OTH is the empirical leftover cluster at K=8. It is not a reference
        # ancestry, so an OTH sample is neither homogeneous (there is no
        # ancestry for it to be homogeneous FOR) nor admixed between named
        # components. It gets its own label.
        if sp == "OTH" or ap == "OTH":
            stype = "real_unassigned"
            # Support is "none" by definition, not whatever the sample's
            # metadata label happens to have. An OTH sample's assigned
            # super-pop is a label the supervised Q contradicts, so reading
            # support off it would report a reference panel that does not
            # cover this person.
            support = "none"
        elif q >= thr and ap == sp:
            stype = "real_homogeneous"
            support = support_for(sp)
        else:
            stype = "real_admixed"
            support = support_for(sp)
        test.append({
            "sample": s,
            "cohort": r[idx["cohort"]],
            "pop": r[idx["pop"]],
            "super_pop": sp,
            "argmax_pop_sup": ap,
            "max_Q_sup": q,
            "threshold": thr,
            "sample_type": stype,
            # The person's ancestry and whether a matching reference panel
            # exists are different facts, and conflating them would make it
            # impossible to ask "how does imputation do when the reference is
            # thin?" separately from "how does it do on people no reference
            # fits?". So support is its own column.
            "reference_support": support,
            "panel_included": r[idx["included"]].strip(),
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
                "max_Q_sup", "threshold", "sample_type", "reference_support",
                "panel_included", "chip", "build"]
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
    sup_of = {}
    for t in test:
        sup_of[t["argmax_pop_sup"]] = t["reference_support"]
    w(f" {'component':<12}{'homog':>8}{'admixed':>9}{'unassigned':>12}"
      f"{'total':>8}{'ref':>7}")
    for c in sorted(by):
        h = by[c]["real_homogeneous"]
        m = by[c]["real_admixed"]
        u = by[c]["real_unassigned"]
        w(f" {c:<12}{h:>8}{m:>9}{u:>12}{h+m+u:>8}{sup_of.get(c, '?'):>7}")
    th = sum(1 for t in test if t["sample_type"] == "real_homogeneous")
    tm = sum(1 for t in test if t["sample_type"] == "real_admixed")
    tu = sum(1 for t in test if t["sample_type"] == "real_unassigned")
    w(f" {'TOTAL':<12}{th:>8}{tm:>9}{tu:>12}{len(test):>8}")
    w("")
    w(" real_homogeneous  top Q >= threshold and argmax == assigned super-pop")
    w(" real_admixed      naturally admixed, or argmax disagrees with the label")
    w(" real_unassigned   OTH, the empirical leftover cluster at K=8 -- no")
    w("                   reference ancestry fits, so it is neither homogeneous")
    w("                   nor admixed between named components. These are the")
    w("                   multiracial-like individuals a real diverse cohort")
    w("                   carries, and every one not in a reference panel is")
    w("                   included here on purpose.")
    w("")
    w(" ref column: full / thin (< %d reference samples) / none. A thin or"
      % a.thin_reference)
    w(" absent reference panel is a property of the PANEL, not of the person,")
    w(" which is why it is separate from sample_type: it lets 'imputation with")
    w(" a thin reference' be asked apart from 'imputation on someone no")
    w(" reference fits'. OCE and OTH are the groups that cannot be pooled into")
    w(" a stratified analysis, and they are in the cohort for that reason.")
    w("")
    w(" Real genomes have no recorded tracts, so they are NEVER scored for")
    w(" local-ancestry accuracy. They are used for:")
    w("   1. imputation accuracy, calibration and retention")
    w("   2. power and GWAS simulations -- phenotypes on real genotypes")
    w(" The second is why the real half cannot be replaced by more simulated")
    w(" individuals: power depends on the LD structure, allele-frequency")
    w(" spectrum and relatedness the genotypes actually carry, and simulated")
    w(" genomes from a finite donor pool have compressed haplotype diversity")
    w(" plus donor-reuse cryptic relatedness, which makes power estimated on")
    w(" them optimistic and uncorrectable after the fact.")
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
        dn[sp] = len(role_list(a.refs, sp, "simulation", a.suffix))
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
    # Donors allocated for a component no model draws on are doing nothing,
    # and in a thin super-pop those are exactly the samples that would be
    # worth more as test genomes.
    unused = []
    for sp in sorted(dn):
        comp = [c for c, p in COMPONENT_POOL.items() if p == sp]
        if dn[sp] > 0 and not any(need_by_comp.get(c, 0) > 0 for c in comp):
            unused.append((sp, dn[sp]))
    if unused:
        w(" DONORS ALLOCATED BUT UNUSED BY ANY MODEL")
        for sp, n_u in unused:
            w(f"      {sp}: {n_u} donors withheld, no model draws on them.")
        w("      Nothing simulates from these, so they are withheld for no")
        w("      benefit while being exactly the samples a thin group needs as")
        w("      test genomes. Re-run the split with")
        w("        --no-donor-superpops " + ",".join(sp for sp, _ in unused))
        w("      to return them to the test pool.")
        w("")
    w("=" * 72)
    w(" NEXT: step 6  6b_simulate_benchmark_cohort_clm.sh  (the 2,160)")
    w("       step 7  6c_extract_real_test_genomes_clm.sh (the real half)")
    w("       transfer 6a_transfer_to_mab_clm.sh benchmark")
    w(" NOT BUILT: step 10 (chip files per chip x build x chr),")
    w("            step 11 (tables, checksums, provenance).")
    w("=" * 72)

    with open(rep, "w") as g:
        g.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwrote {man}\n      {ros}\n      {ind}\n      {rep}")


if __name__ == "__main__":
    main()
