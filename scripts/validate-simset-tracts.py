#!/usr/bin/env python3
"""
validate-simset-tracts.py -- did the simulated cohort get the demography it
was asked for?

WHY THIS EXISTS, AND WHY IT BEATS INSPECTING THE .dat FILE.
The .dat body is the admixture history: a founding pulse, then a terminal line
whose only job is to set how many generations of random mating follow. The
`mixer` binary parses it and simulates whatever history it reads. It has no
notion of what was intended, so a misshapen history does not fail -- it
produces well-formed .bp and .hanc files for the WRONG demography.

Generations set tract length, and tract length is the main determinant of how
hard local-ancestry inference is. A gen-16 model that actually ran one
generation yields near-whole-chromosome tracts, ancestry proportions still
close to target, and RFMix and FLARE both scoring far better than they should
-- uniformly, with nothing in the output to say so. The benchmark would be
systematically too easy and still look right.

Reading the .dat cannot catch that, because the .dat is the input whose
interpretation is in question. The TRUTH TRACTS can, because they have a sharp
prediction. For a single pulse g generations ago with ancestry proportions
p_i, recombination events ancestral to a haplotype arrive at (g-1) per Morgan,
and each switches ancestry only when the two sides come from different
sources, with probability (1 - sum p_i^2):

    switches per Morgan  ~  (g - 1) * (1 - sum_i p_i^2)
    mean tract (cM)      ~  100 / that

gen 2 predicts ~200 cM and gen 30 predicts ~6 cM -- an order of magnitude
apart. So if the generation lines are not taking effect, every model reads as
gen 1-2 whatever its .dat says, and that is unmissable in one table.

This also catches the quieter failure: proportions silently transposed,
because the realized ancestry fractions are reported against the intended
ones.

Usage:
  python3 -I scripts/validate-simset-tracts.py \\
      --simset 06_benchmark_simset \\
      --individuals results/benchmark_cohort_simset_individuals.tsv \\
      --chr chr2
"""
import argparse
import glob
import os
import sys
from collections import defaultdict


def read_snp_positions(path):
    """Genetic positions in cM from an admix-simu .snp (EIGENSTRAT snp).

    Column 3 is the genetic position. Some writers use Morgans and some use
    centiMorgans, and guessing wrong rescales every tract length by 100, which
    would turn this check into a generator of false alarms. Decide from the
    span: no human chromosome is 300+ Morgans, and none is under 3 cM.
    """
    pos = []
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 3:
                continue
            try:
                pos.append(float(p[2]))
            except ValueError:
                continue
    if len(pos) < 2:
        return None, None
    span = pos[-1] - pos[0]
    if span <= 0:
        return None, None
    unit = "M" if span < 10 else "cM"
    if unit == "M":
        pos = [p * 100.0 for p in pos]
    return pos, unit


def tract_stats(hanc_path, pos_cm):
    """Switch count per haplotype and realized ancestry fractions.

    .hanc is one line per site, one character per haplotype, the character
    being the ancestry index. Read it streaming: a 600-individual model over
    ~100k sites is 120M characters, and holding it all would be pointless.
    """
    n_sites = len(pos_cm)
    prev = None
    switches = None
    counts = None
    seen = 0
    with open(hanc_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if switches is None:
                switches = [0] * len(ln)
                counts = [defaultdict(int) for _ in range(len(ln))]
            elif len(ln) != len(switches):
                return None, None, (f"ragged .hanc: site {seen+1} has "
                                    f"{len(ln)} haplotypes, expected "
                                    f"{len(switches)}"), 0, 0.0
            for i, ch in enumerate(ln):
                counts[i][ch] += 1
                if prev is not None and prev[i] != ch:
                    switches[i] += 1
            prev = ln
            seen += 1
    if switches is None:
        return None, None, "empty .hanc", 0, 0.0
    if seen != n_sites:
        return None, None, (f".hanc has {seen} sites but .snp has {n_sites}; "
                            f"they must describe the same marker set"), 0, 0.0
    span_cm = pos_cm[-1] - pos_cm[0]
    per_morgan = [s / (span_cm / 100.0) for s in switches]
    n_hap = len(switches)
    tot = defaultdict(int)
    for c in counts:
        for k, v in c.items():
            tot[k] += v
    grand = sum(tot.values()) or 1
    frac = {k: v / grand for k, v in sorted(tot.items())}
    return per_morgan, frac, None, n_hap, span_cm


def predicted_switches(gen, props):
    if gen is None or gen <= 1:
        return 0.0
    return (gen - 1) * (1.0 - sum(p * p for p in props))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simset", required=True,
                    help="06_benchmark_simset (holds chr*/ subdirs)")
    ap.add_argument("--individuals", required=True,
                    help="benchmark_cohort_simset_individuals.tsv")
    ap.add_argument("--chr", default="",
                    help="one chromosome, e.g. chr2. Default: every chr* "
                         "directory present.")
    ap.add_argument("--tolerance", type=float, default=0.35,
                    help="relative deviation in switches/Morgan tolerated "
                         "before a model is flagged (default 0.35)")
    a = ap.parse_args()

    models = {}
    with open(a.individuals) as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for c in ("model", "components", "proportions", "generations"):
            if c not in ix:
                sys.exit(f"FATAL: {a.individuals} lacks column '{c}'")
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(head):
                continue
            m = p[ix["model"]]
            if m in models:
                continue
            models[m] = {
                "components": p[ix["components"]].split(","),
                "props": [float(x) for x in p[ix["proportions"]].split(",")],
                "gen": int(p[ix["generations"]]),
            }
    if not models:
        sys.exit(f"FATAL: no models in {a.individuals}")

    chrs = ([a.chr] if a.chr
            else sorted(os.path.basename(d) for d in
                        glob.glob(os.path.join(a.simset, "chr*"))
                        if os.path.isdir(d)))
    if not chrs:
        sys.exit(f"FATAL: no chr* directories under {a.simset}")

    rows, problems, skipped = [], [], []
    for ch in chrs:
        d = os.path.join(a.simset, ch)
        snp = os.path.join(d, f"{ch}.snp")
        if not os.path.exists(snp):
            skipped.append(f"{ch}: no {ch}.snp"); continue
        pos, unit = read_snp_positions(snp)
        if pos is None:
            skipped.append(f"{ch}: unusable genetic positions in {ch}.snp")
            continue
        for m, info in sorted(models.items()):
            if info["gen"] == 0:
                continue   # controls have constant ancestry by construction
            hanc = os.path.join(d, f"{m}.{ch}.hanc")
            if not os.path.exists(hanc):
                skipped.append(f"{ch}/{m}: no .hanc"); continue
            per_m, frac, err, n_hap, span_cm = tract_stats(hanc, pos)
            if err:
                problems.append(f"{ch}/{m}: {err}"); continue
            obs = sum(per_m) / len(per_m)
            exp = predicted_switches(info["gen"], info["props"])
            obs_cm = (100.0 / obs) if obs > 0 else float("inf")
            exp_cm = (100.0 / exp) if exp > 0 else float("inf")
            rel = abs(obs - exp) / exp if exp > 0 else 0.0
            rows.append((ch, m, info["gen"], exp, obs, exp_cm, obs_cm, rel,
                         frac, info, n_hap, span_cm))

    print("=" * 78)
    print(" SIMSET TRACT VALIDATION -- did the .dat generations take effect?")
    print("=" * 78)
    if skipped:
        print(f" skipped {len(skipped)} unit(s) with no output yet "
              f"(e.g. {skipped[0]})")
        print("")
    print(f" {'chr':<7}{'model':<9}{'gen':>4}{'sw/M exp':>10}{'sw/M obs':>10}"
          f"{'cM exp':>9}{'cM obs':>9}{'dev':>7}")
    print("-" * 78)
    flagged = []
    for (ch, m, g, exp, obs, ecm, ocm, rel, frac, info, n_hap, span_cm) in rows:
        mark = ""
        if rel > a.tolerance:
            mark = "  <<<"
            flagged.append((ch, m, g, exp, obs))
        print(f" {ch:<7}{m:<9}{g:>4}{exp:>10.2f}{obs:>10.2f}"
              f"{ecm:>9.1f}{ocm:>9.1f}{rel:>6.0%}{mark}")
    print("-" * 78)

    # Realized vs intended ancestry fractions. Catches transposed proportions,
    # which tract length alone cannot see.
    print("")
    print(" REALIZED ANCESTRY FRACTIONS (intended -> observed)")
    print("-" * 78)
    seen_models = set()
    ancestry_bad = []
    for (ch, m, g, exp, obs, ecm, ocm, rel, frac, info, n_hap, span_cm) in rows:
        if m in seen_models:
            continue
        seen_models.add(m)
        # VARIANCE-AWARE BOUND, not a flat tolerance.
        #
        # Realized ancestry is a mean over independent ancestry segments, and
        # the number of those scales with generations: a gen-2 model has under
        # one switch per haplotype, so its realized fractions scatter widely
        # around the target by construction. A flat +/-0.05 would flag every
        # low-generation model as broken while letting a genuinely transposed
        # gen-30 model through -- wrong in both directions.
        #
        # Effective independent draws ~ haplotypes * segments per haplotype,
        # and segments ~ max(1, (g-1) * span in Morgans).
        span_M = span_cm / 100.0
        n_eff = max(1.0, n_hap * max(1.0, (g - 1) * span_M))
        parts, bad_comp = [], []
        for i, comp in enumerate(info["components"]):
            want = info["props"][i]
            got = frac.get(str(i), 0.0)
            se = (want * (1.0 - want) / n_eff) ** 0.5
            bound = max(0.03, 4.0 * se)
            off = abs(got - want)
            flag = "!" if off > bound else ""
            if flag:
                bad_comp.append(f"{comp} {want:.2f}->{got:.2f} "
                                f"(off {off:.3f} > {bound:.3f})")
            parts.append(f"{comp} {want:.2f}->{got:.2f}{flag}")
        bad = "  <<<" if bad_comp else ""
        print(f" {m:<9}({ch})  " + "  ".join(parts) + bad)
        if bad_comp:
            problems.append(f"{m}: realized ancestry outside 4 SE -- "
                            + "; ".join(bad_comp))
            ancestry_bad.append(m)
    print("-" * 78)
    print("")

    if problems:
        print(" PROBLEMS")
        for p in problems[:20]:
            print(f"   {p}")
        print("")

    if flagged:
        print(" *** TRACT LENGTHS DISAGREE WITH THE INTENDED GENERATIONS ***")
        spread = {g for (_c, _m, g, _e, _o) in flagged}
        obs_vals = [o for (_c, _m, _g, _e, o) in flagged]
        print("")
        if len(spread) > 1 and max(obs_vals) - min(obs_vals) < 0.5:
            print(" Every flagged model shows nearly the SAME switch rate")
            print(" despite different intended generations. That is the")
            print(" signature of the generation lines not taking effect at")
            print(" all -- the terminal line of each .dat is being ignored or")
            print(" misparsed, so every model ran as a single pulse.")
            print("")
            print(" Do NOT score this cohort. Fix write_dat in")
            print(" 6b_simulate_benchmark_cohort_clm.sh against the admix-simu")
            print(" .dat format, re-simulate one chromosome, and re-run this.")
        else:
            print(" Deviations differ by model, so the generation lines are")
            print(" being read but at least one model's history is wrong.")
            print(" Check those models' .dat in the provenance directory.")
        print("")
        sys.exit(1)

    if not rows:
        print(" nothing to validate yet -- no simulated output found.")
        sys.exit(0)
    if problems:
        # Printing PROBLEMS and then "ok" with exit 0 is worse than either
        # alone: it trains you to ignore the section.
        print(" NOT OK -- see PROBLEMS above. Tract lengths track the intended")
        print(" generations, but realized ancestry does not match the models'")
        print(" proportions, which points at the .dat PULSE line rather than")
        print(" its terminal generation line.")
        sys.exit(1)
    print(" ok -- observed tract lengths track the intended generations, and")
    print(" realized ancestry matches the models' proportions.")
    sys.exit(0)


if __name__ == "__main__":
    main()
