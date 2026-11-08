#!/usr/bin/env python3
"""
cohort-ancestry-summary.py -- mean ancestry proportions across every testing
population in the benchmark cohort, simulated and real, in one table.

WHY THIS IS A SEPARATE CHECK FROM THE TRACT VALIDATOR.
validate-simset-tracts.py counts ancestry SWITCHES along each haplotype. This
counts ancestry FRACTIONS per individual. They are different statistics over
the same truth, so a bug in one does not reproduce in the other, and the
second has its own independent read on generations:

  mean  individual ancestry fraction  ->  should equal the model's proportion,
                                          at any number of generations
  SD    between individuals           ->  shrinks as generations rise

An individual's genome is a mosaic of roughly 2*(g-1)*L segments (two
haplotypes, L Morgans), so the between-individual spread in ancestry fraction
is about

    SD  ~  sqrt( p(1-p) / (2 (g-1) L) )

Per chromosome (L ~ 1.4 M) that is ~0.30 at g=2 and ~0.077 at g=16: a factor
of four. So a cohort where the .dat terminal line was ignored and every model
actually ran one generation shows the RIGHT means -- the pulse proportions are
set on line 1 and survive -- with SDs several times too large, uniformly.
That is exactly the failure the means alone cannot see, which is what makes
checking the means alone a false reassurance.

THE REAL HALF HAS NO TRUTH TRACTS, so its ancestry comes from the supervised
ADMIXTURE Q in the panel table. There is nothing to validate it against -- it
IS the measurement -- so it is reported for inspection rather than tested:
whether the real admixed AMR samples look Latino, whether a group is
unexpectedly homogeneous, whether the cohort spans the ancestry space the
models were chosen to cover.

Usage:
  python3 -I scripts/cohort-ancestry-summary.py \\
      --panel 04_homogeneity_panel/lai_ref_panel_samples.k8.tsv \\
      --manifest results/benchmark_cohort_manifest.tsv \\
      --simset 06_benchmark_simset \\
      --individuals results/benchmark_cohort_simset_individuals.tsv \\
      --chr chr2 --out results/cohort_ancestry_summary.txt
"""
import argparse
import glob
import math
import os
import sys
from collections import defaultdict


def read_tsv(path):
    if not os.path.exists(path):
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


def mean_sd(xs):
    n = len(xs)
    if n == 0:
        return 0.0, 0.0
    m = sum(xs) / n
    if n < 2:
        return m, 0.0
    v = sum((x - m) ** 2 for x in xs) / (n - 1)
    return m, math.sqrt(v)


def snp_span_morgans(path):
    pos = []
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if len(p) >= 3:
                try:
                    pos.append(float(p[2]))
                except ValueError:
                    pass
    if len(pos) < 2:
        return None
    span = pos[-1] - pos[0]
    if span <= 0:
        return None
    # Morgans if the span is small, centiMorgans otherwise. Guessing wrong
    # rescales the predicted SD by 10 and would make this a false-alarm
    # generator; no human chromosome is 300+ Morgans or under 3 cM.
    return span if span < 10 else span / 100.0


def per_individual_fractions(hanc_path, n_anc):
    """Ancestry fraction per INDIVIDUAL, averaging its two haplotypes.

    .hanc is one line per site, one character per haplotype, haplotypes 2j-1
    and 2j belonging to individual j. Counted streaming: a 600-individual
    model over ~100k sites would be 120M characters held at once otherwise.
    """
    counts = None
    n_sites = 0
    with open(hanc_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if counts is None:
                nh = len(ln)
                if nh % 2:
                    return None, f"odd haplotype count {nh} -- not diploid"
                counts = [[0] * n_anc for _ in range(nh // 2)]
            elif len(ln) != len(counts) * 2:
                return None, (f"ragged .hanc at site {n_sites+1}: "
                              f"{len(ln)} vs {len(counts)*2}")
            for i, ch in enumerate(ln):
                if ch.isdigit():
                    k = int(ch)
                    if k < n_anc:
                        counts[i // 2][k] += 1
            n_sites += 1
    if counts is None or n_sites == 0:
        return None, "empty .hanc"
    denom = 2.0 * n_sites
    return [[c / denom for c in ind] for ind in counts], None


def predicted_sd(p, gen, span_M):
    if gen is None or gen <= 1 or span_M is None or span_M <= 0:
        return None
    seg = 2.0 * (gen - 1) * span_M
    if seg <= 0:
        return None
    return math.sqrt(max(p * (1.0 - p), 0.0) / seg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simset", default="", help="06_benchmark_simset")
    ap.add_argument("--individuals", default="",
                    help="benchmark_cohort_simset_individuals.tsv")
    ap.add_argument("--panel", default="",
                    help="lai_ref_panel_samples.k<K>.tsv, for the real half")
    ap.add_argument("--manifest", default="",
                    help="benchmark_cohort_manifest.tsv")
    ap.add_argument("--chr", default="",
                    help="one chromosome for the simulated side; default: all")
    ap.add_argument("--sd-tolerance", type=float, default=2.5,
                    help="flag a model whose observed between-individual SD "
                         "is off by more than this factor (default 2.5)")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    L, flagged, notes = [], [], []
    w = L.append
    w("=" * 78)
    w(" COHORT ANCESTRY PROPORTIONS -- every testing population")
    w("=" * 78)

    # ---------------- simulated half ----------------
    if a.simset and a.individuals:
        ix, rows = read_tsv(a.individuals)
        if ix is None:
            sys.exit(f"FATAL: missing {a.individuals}")
        for c in ("model", "components", "proportions", "generations"):
            if c not in ix:
                sys.exit(f"FATAL: {a.individuals} lacks column '{c}'")
        models = {}
        for p in rows:
            m = p[ix["model"]]
            if m not in models:
                models[m] = (p[ix["components"]].split(","),
                             [float(x) for x in p[ix["proportions"]].split(",")],
                             int(p[ix["generations"]]))
        chrs = ([a.chr] if a.chr else
                sorted(os.path.basename(d) for d in
                       glob.glob(os.path.join(a.simset, "chr*"))
                       if os.path.isdir(d)))
        w("")
        w(" SIMULATED -- truth tracts (.hanc), per individual")
        w("-" * 78)
        w(f" {'model':<9}{'gen':>4}{'comp':>6}{'want':>7}{'mean':>7}"
          f"{'SD obs':>8}{'SD exp':>8}{'ratio':>7}{'n':>6}")
        any_sim = False
        for ch in chrs:
            snp = os.path.join(a.simset, ch, f"{ch}.snp")
            span_M = snp_span_morgans(snp) if os.path.exists(snp) else None
            for m in sorted(models):
                comps, props, gen = models[m]
                if gen == 0:
                    continue        # constant ancestry by construction
                hanc = os.path.join(a.simset, ch, f"{m}.{ch}.hanc")
                if not os.path.exists(hanc):
                    continue
                fr, err = per_individual_fractions(hanc, len(comps))
                if err:
                    notes.append(f"{ch}/{m}: {err}")
                    continue
                any_sim = True
                for k, comp in enumerate(comps):
                    xs = [ind[k] for ind in fr]
                    mo, so = mean_sd(xs)
                    se = predicted_sd(props[k], gen, span_M)
                    ratio = (so / se) if (se and se > 0) else None
                    mark = ""
                    if ratio is not None and (ratio > a.sd_tolerance
                                              or ratio < 1.0 / a.sd_tolerance):
                        mark = " <<<"
                        flagged.append((m, gen, comp, ratio))
                    w(f" {m:<9}{gen:>4}{comp:>6}{props[k]:>7.2f}{mo:>7.3f}"
                      f"{so:>8.3f}"
                      f"{(f'{se:.3f}' if se else '    -'):>8}"
                      f"{(f'{ratio:.2f}' if ratio else '    -'):>7}"
                      f"{len(xs):>6}{mark}")
                w("")
        if not any_sim:
            w(" no simulated output found yet.")
        w("-" * 78)
        w(" mean should match 'want' at ANY number of generations: the pulse")
        w(" proportions are set on the .dat's first line and survive even if")
        w(" the terminal generation line is ignored. So a correct mean is NOT")
        w(" evidence the demography is right -- the SD is what carries that.")
        w(" SD exp = sqrt(p(1-p) / (2(g-1)L)); ratio = obs/exp, want ~1.")

    # ---------------- real half ----------------
    if a.panel and a.manifest:
        pix, prows = read_tsv(a.panel)
        mix, mrows = read_tsv(a.manifest)
        if pix is None:
            sys.exit(f"FATAL: missing {a.panel}")
        if mix is None:
            sys.exit(f"FATAL: missing {a.manifest}")
        qcols = [(c[len("Q_sup_"):], i) for c, i in pix.items()
                 if c.startswith("Q_sup_")]
        qcols.sort(key=lambda t: t[0])
        if not qcols:
            sys.exit(f"FATAL: {a.panel} has no Q_sup_* columns")
        qby = {}
        for p in prows:
            vals = {}
            for name, i in qcols:
                try:
                    vals[name] = float(p[i])
                except (ValueError, IndexError):
                    vals[name] = 0.0
            qby[p[pix["sample"]]] = vals
        groups = defaultdict(list)
        for p in mrows:
            s = p[mix["sample"]]
            if s not in qby:
                notes.append(f"manifest sample {s} absent from the panel table")
                continue
            key = (p[mix["sample_type"]], p[mix["argmax_pop_sup"]])
            groups[key].append(qby[s])
        w("")
        w(" REAL -- supervised ADMIXTURE Q (no truth tracts; for inspection)")
        w("-" * 78)
        if not groups:
            w(" no real test genomes in the manifest.")
        for (stype, ap_) in sorted(groups):
            members = groups[(stype, ap_)]
            w(f" {stype} / {ap_}   n={len(members)}")
            parts = []
            for name, _i in qcols:
                xs = [mm.get(name, 0.0) for mm in members]
                mo, so = mean_sd(xs)
                if mo >= 0.01:
                    parts.append(f"{name} {mo:.3f}+/-{so:.3f}")
            w("    " + ("  ".join(parts) if parts else "(all components <0.01)"))
        w("-" * 78)
        w(" These are a MEASUREMENT, not a truth, so nothing here is pass/fail.")
        w(" Read them for whether each group looks like the ancestry it claims")
        w(" and whether the cohort spans the space the models were chosen to")
        w(" cover. Real genomes are never scored for local ancestry.")

    if notes:
        w("")
        w(" NOTES")
        for n in notes[:20]:
            w(f"   {n}")
        if len(notes) > 20:
            w(f"   ... and {len(notes) - 20} more")

    w("")
    if flagged:
        w(" *** BETWEEN-INDIVIDUAL SPREAD DISAGREES WITH THE INTENDED"
          " GENERATIONS ***")
        hi = [f for f in flagged if f[3] > 1.0]
        if len(hi) == len(flagged) and len({f[1] for f in flagged}) > 1:
            w("")
            w(" Every flagged component is MORE variable between individuals")
            w(" than its generations predict, across models with different")
            w(" intended generations. That is the signature of the .dat")
            w(" terminal line being ignored: each model ran as a single pulse,")
            w(" which leaves the means correct and the spread far too wide.")
            w(" Cross-check with validate-simset-tracts.py -- it should show")
            w(" the same models with near-identical switch rates.")
        else:
            w(" Deviations are not uniform, so the generation lines are being")
            w(" read but at least one model's history is wrong. Check those")
            w(" models' .dat in the provenance directory.")
        w("")

    text = "\n".join(L) + "\n"
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as g:
            g.write(text)
    print(text, end="")
    sys.exit(1 if flagged else 0)


if __name__ == "__main__":
    main()
