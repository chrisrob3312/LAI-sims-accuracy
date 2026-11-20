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
from itertools import permutations
from collections import defaultdict


def hanc_path(d, model, ch):
    """Prefer the site-per-line .hanc2; fall back to .hanc and SAY SO.

    bp2anc.pl writes ONE LINE PER HAPLOTYPE. Every reader here wants one line
    per site, so 6b (like 2b) transposes once into .hanc2. Reading a raw .hanc
    as if it were site-per-line does not fail -- it silently transposes the
    interpretation, turning N_SITES x 2n into 2n x N_SITES, and every tract
    and ancestry number computed from it is meaningless. So the fallback is
    reported, not silent.
    """
    h2 = os.path.join(d, f"{model}.{ch}.hanc2")
    if os.path.exists(h2):
        return h2, None
    h1 = os.path.join(d, f"{model}.{ch}.hanc")
    if os.path.exists(h1):
        return h1, (f"{ch}/{model}: no .hanc2, reading raw .hanc -- if that "
                    f"file is haplotype-per-line (bp2anc's layout) every "
                    f"number from it is transposed. Re-run 6b, which writes "
                    f".hanc2.")
    return None, None


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


def read_dat(path):
    """-> (n, model, [components], [pulse proportions], terminal_gen) or None.

    The .dat is what write_dat ACTUALLY wrote, as opposed to what the roster
    intended. Reading it is what lets a disagreement be localised: roster vs
    .dat blames write_dat, .dat vs .hanc blames mixer.
    """
    try:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip()]
    except OSError:
        return None
    if len(lines) < 2:
        return None
    h = lines[0].split()
    if len(h) < 3:
        return None
    try:
        n = int(h[0])
    except ValueError:
        return None
    model, comps = h[1], h[2:]
    pulse = None
    term_gen = 1
    for ln in lines[1:]:
        f = ln.split()
        if len(f) != len(comps) + 1:
            continue
        try:
            g = int(float(f[0]))
            vals = [float(x) for x in f[1:]]
        except ValueError:
            continue
        if pulse is None and any(v > 0 for v in vals):
            pulse = vals
        term_gen = max(term_gen, g)
    if pulse is None:
        return None
    return n, model, comps, pulse, term_gen


def is_permutation(a, b, tol=1e-6):
    """True when b is a reordering of a but not equal to it."""
    if len(a) != len(b):
        return False
    if all(abs(x - y) < tol for x, y in zip(a, b)):
        return False
    return all(abs(x - y) < tol for x, y in zip(sorted(a), sorted(b)))


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
    ap.add_argument("--provenance", default="",
                    help="directory holding the archived .dat files "
                         "(default <simset>/provenance/<chr>). Comparing the "
                         ".dat against both the roster and the realized "
                         "ancestry is what localises a disagreement: roster "
                         "vs .dat blames write_dat, .dat vs .hanc blames "
                         "mixer.")
    ap.add_argument("--sd-tolerance", type=float, default=2.5,
                    help="flag a model whose observed between-individual SD "
                         "is off by more than this factor (default 2.5)")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    L, flagged, notes = [], [], []
    # problems: anything worth naming. link1 / perm / mean_bad localise
    # WHICH of roster -> .dat -> .hanc broke, which is the whole point
    # of reading the archived .dat rather than only the outcome.
    problems, dat_notes, link1, perm, mean_bad = [], [], [], [], []
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
                hanc, warn = hanc_path(os.path.join(a.simset, ch), m, ch)
                if warn:
                    notes.append(warn)
                if hanc is None:
                    continue
                fr, err = per_individual_fractions(hanc, len(comps))
                if err:
                    notes.append(f"{ch}/{m}: {err}")
                    continue
                any_sim = True

                # --- link 1: roster -> .dat (did write_dat write the plan?) --
                provdir = a.provenance or os.path.join(a.simset,
                                                       "provenance", ch)
                dat = read_dat(os.path.join(provdir, f"{m}.{ch}.dat"))
                if dat is None:
                    dat_notes.append(f"{ch}/{m}: no archived .dat under "
                                     f"{provdir} -- cannot confirm the file "
                                     f"matched the plan, only the outcome")
                    dat_props, dat_comps, dat_gen = None, None, None
                else:
                    _n, _mn, dat_comps, dat_props, dat_gen = dat
                    if dat_comps != comps:
                        problems.append(
                            f"{m}: .dat components {dat_comps} != roster "
                            f"{comps}. The -POP argument order follows the "
                            f"roster, so a different .dat header order would "
                            f"PERMUTE the ancestry indices.")
                        link1.append(m)
                    elif not all(abs(x - y) < 1e-6
                                 for x, y in zip(dat_props, props)):
                        problems.append(
                            f"{m}: .dat pulse {dat_props} != roster {props} "
                            f"-- write_dat did not write the planned "
                            f"proportions")
                        link1.append(m)
                    if dat_gen is not None and dat_gen != gen:
                        problems.append(
                            f"{m}: .dat terminal generation {dat_gen} != "
                            f"roster {gen}")
                        link1.append(m)

                # --- link 2: .dat -> .hanc (did mixer simulate the file?) ----
                # THE MEANS ARE TESTED, not just printed. A variance-aware
                # bound, because realized ancestry is a mean over independent
                # segments whose number scales with generations: a gen-2 model
                # scatters widely by construction, so a flat tolerance would
                # flag every low-generation model while letting a transposed
                # gen-30 model through.
                want = dat_props if dat_props else props
                obs_means, obs_sds = [], []
                n_ind = len(fr)
                for k, comp in enumerate(comps):
                    xs = [ind[k] for ind in fr]
                    mo, so = mean_sd(xs)
                    obs_means.append(mo)
                    se = predicted_sd(props[k], gen, span_M)
                    ratio = (so / se) if (se and se > 0) else None
                    obs_sds.append(so)
                    # THE MEAN'S BOUND COMES FROM THE OBSERVED SPREAD, not
                    # from the intended generations. Deriving it from `gen`
                    # would be circular: whether the simulation ran at `gen`
                    # is exactly what is in question, and data generated at
                    # the wrong generation count is more variable than `gen`
                    # predicts, so a gen-derived bound is too tight and
                    # reports a mean failure for a demography failure. The
                    # standard error of a mean over n individuals is
                    # SD/sqrt(n) whatever produced the SD.
                    m_se = (so / (n_ind ** 0.5)) if n_ind > 1 else 1.0
                    m_bound = max(0.03, 4.0 * m_se)
                    m_off = abs(mo - want[k])
                    mark = ""
                    if m_off > m_bound:
                        mark = " <<<M"
                        mean_bad.append((m, comp, want[k], mo, m_off, m_bound))
                    if ratio is not None and (ratio > a.sd_tolerance
                                              or ratio < 1.0 / a.sd_tolerance):
                        mark += " <<<SD"
                        flagged.append((m, gen, comp, ratio))
                    w(f" {m:<9}{gen:>4}{comp:>6}{want[k]:>7.2f}{mo:>7.3f}"
                      f"{so:>8.3f}"
                      f"{(f'{se:.3f}' if se else '    -'):>8}"
                      f"{(f'{ratio:.2f}' if ratio else '    -'):>7}"
                      f"{len(xs):>6}{mark}")
                # A permutation is a component-ORDER bug, not a demography
                # one, and for LAT3 (0.50/0.45/0.05) swapping the first two
                # is a change small enough to hide in noise unless named.
                # A permutation claim is only meaningful when the identity
                # assignment actually FAILS and some reordering fits clearly
                # better. Asserting it on noise alone mislabels a demography
                # failure as a component-order bug -- the proportions of LAT3
                # (0.50/0.45/0.05) sit closer together than the noise at
                # modest n, so sorted-value agreement proves nothing by
                # itself.
                if any(b[0] == m for b in mean_bad):
                    mb = max(0.03, 4.0 * (max(obs_sds) / (n_ind ** 0.5))
                             if n_ind > 1 else 1.0)
                    err_id = max(abs(o - wv)
                                 for o, wv in zip(obs_means, want))
                    best = None
                    for pm in permutations(range(len(want))):
                        if list(pm) == list(range(len(want))):
                            continue
                        e = max(abs(obs_means[i] - want[pm[i]])
                                for i in range(len(want)))
                        best = e if best is None else min(best, e)
                    if best is not None and err_id > mb and best < mb / 2.0:
                        problems.append(
                            f"{m}: realized ancestry is a PERMUTATION of the "
                            f"planned proportions ({want} -> "
                            f"{[round(x, 3) for x in obs_means]}); identity "
                            f"misfit {err_id:.3f} vs best reordering "
                            f"{best:.3f}. That is a component-order bug "
                            f"between the .dat header and the -POP arguments, "
                            f"not a wrong demography.")
                        perm.append(m)
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

    if dat_notes:
        w("")
        w(" .dat PROVENANCE")
        for n in dat_notes[:10]:
            w(f"   {n}")

    w("")
    w(" DID IT SIMULATE WHAT WE PLANNED?")
    w("-" * 78)
    w(" Three things must agree, and naming the failing link is the point:")
    w("   roster  (the plan)  ->  .dat  (what write_dat wrote)  ->  .hanc")
    w("")
    if link1:
        w(f" *** LINK 1 BROKEN for {sorted(set(link1))}: the archived .dat does")
        w(" NOT match the roster. write_dat is wrong -- mixer simulated the")
        w(" file it was given faithfully. Fix write_dat in 6b; re-simulating")
        w(" without fixing it reproduces the same error.")
    elif perm:
        w(f" *** PERMUTED for {sorted(set(perm))}: the .dat matched the plan but")
        w(" the realized ancestry is a reordering of it. The .dat header order")
        w(" and the -POP argument order disagree, so the ancestry INDICES in")
        w(" .hanc point at the wrong components. Not a demography error.")
    elif mean_bad:
        w(" *** LINK 2: the .dat matched the plan but the realized means do")
        w(" not match the .dat:")
        for mm, comp, wv, ov, off, bnd in mean_bad[:12]:
            w(f"      {mm} {comp}: planned {wv:.3f}, got {ov:.3f} "
              f"(off {off:.3f} > {bnd:.3f})")
        w(" mixer did not produce the proportions in its own input file.")
    elif flagged:
        w(" Proportions are correct. The SPREAD is not -- see below. That")
        w(" combination is the signature of the demography being wrong while")
        w(" the pulse proportions are right, because the pulse comes off the")
        w(" .dat's first line and survives a mis-parsed terminal line.")
    else:
        w(" ok -- the .dat matches the roster, and the realized ancestry")
        w(" matches the .dat in both mean and spread.")
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
    sys.exit(1 if (flagged or link1 or perm or mean_bad) else 0)


if __name__ == "__main__":
    main()
