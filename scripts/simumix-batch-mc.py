#!/usr/bin/env python3
"""Monte Carlo of 6b's simulation units, built from the two programs' code.

mixer (admix-simu main.cc): a population of max(10000, 20*N) haplotypes per
generation, so the N output haplotypes are effectively independent. Gen 1 is
a pulse: each haplotype is ONE source pop (chosen by the proportions), with
same-pop recombination breakpoints at 1 per Morgan. Each later generation a
haplotype recombines (1 per Morgan) between two random haplotypes of the
previous generation, keeping their segment labels and breakpoints.

simu-mix (scripts/simu-mix-seeded.pl): at every breakpoint, copy a new donor
haplotype of that segment's pop, drawn uniformly among those NOT being copied
by another haplotype of the same run (batch) at that marker, and not the one
just left if the pop is the same. If none is free: die "used all haplotypes".

Part A: how often a (model, chromosome) task dies at 6b's batch cap, and at
alternatives. Part B: pairwise kinship between simulated individuals for the
two gen-2 models and the gen-0 controls.

The numbers in 6b's FEAS_Z comment, WORKFLOW.md and BENCHMARK_COHORT_PLAN.md
(20 Nov) come from `a 200` and `b`. Donor counts and models are the 18 Nov
K=7 split and reference_ids/models.tsv; edit DONORS / MODELS to re-run for
another split. Needs numpy (not in the shapeit5 env; any Python 3 with numpy).

Usage:
  python3 -I scripts/simumix-batch-mc.py a 200   # Part A, 200 replicates
  python3 -I scripts/simumix-batch-mc.py b       # Part B
"""
import math
import sys

import numpy as np

CHR_M = [2.86, 2.68, 2.23, 2.14, 2.04, 1.92, 1.87, 1.68, 1.66, 1.81, 1.58,
         1.75, 1.26, 1.19, 1.41, 1.34, 1.29, 1.17, 1.07, 1.08, 0.62, 0.73]
STEP = 0.005  # Morgans (0.5 cM grid)
DONORS = {"NAT": 76, "EUR": 32, "AFR": 82, "MEN": 39, "EAS": 113}
MODELS = [  # name, comps, props, gen, n   (reference_ids/models.tsv, 20 Nov)
    ("LAT3", ("NAT", "EUR", "AFR"), (0.50, 0.45, 0.05), 16, 300),
    ("CAR3", ("NAT", "EUR", "AFR"), (0.10, 0.50, 0.40), 16, 200),
    ("AFAM2", ("AFR", "EUR"), (0.80, 0.20), 7, 200),
    ("MENA3", ("EUR", "MEN", "AFR"), (0.30, 0.50, 0.20), 30, 200),
    ("MENEUR", ("MEN", "EUR"), (0.50, 0.50), 2, 200),
    ("EASEUR", ("EAS", "EUR"), (0.50, 0.50), 2, 200),
]


def cap_margin(comps, props, m=0.8):
    """6b feasible_batch today: n such that 2*n*p <= m * donor haplotypes."""
    best = 0
    for c, p in zip(comps, props):
        if p <= 0:
            continue
        h = 2 * DONORS[c]
        v = max(1, int(m * h / (2 * p)))
        best = v if best == 0 or v < best else best
    return best


def cap_z(comps, props, z=5.0):
    """Largest n with mean + z*sd of the per-marker count <= haplotypes - 1."""
    best = 0
    for c, p in zip(comps, props):
        if p <= 0:
            continue
        H = 2 * DONORS[c]
        if p >= 1:
            v = (H - 1) // 2
        else:
            b = z * math.sqrt(2 * p * (1 - p))
            x = (-b + math.sqrt(b * b + 8 * p * (H - 1))) / (4 * p)
            v = max(1, int(x * x))
        best = v if best == 0 or v < best else best
    return best


def cap_n1(comps, props):
    """admix-simu README guarantee: N haplotypes need N+1 donors per pop."""
    return min(DONORS[c] - 1 for c, p in zip(comps, props) if p > 0)


def batches(n, cap):
    out = [cap] * (n // cap)
    if n % cap:
        out.append(n % cap)
    return out


def sim_haps(rng, nhap, length_m, props, gen):
    """Segment labels on a grid for nhap independent haplotypes of one
    chromosome. Returns (pop index per marker x hap, breakpoint flags)."""
    nm = int(length_m / STEP) + 1
    k = len(props)
    pops = np.empty((nm, nhap), dtype=np.int8)
    brk = np.zeros((nm, nhap), dtype=bool)
    p = np.asarray(props) / sum(props)

    def gen1(count):
        # pure haplotypes, same-pop breakpoints at 1/M
        lab = rng.choice(k, size=count, p=p)
        P = np.repeat(lab[None, :], nm, axis=0)
        B = rng.random((nm, count)) < (1 - math.exp(-STEP))
        B[0, :] = True
        return P, B

    if gen <= 1:
        P, B = gen1(nhap)
        return P, B
    # Exact recursion is 2^(g-1) ancestors per haplotype; do it exactly for
    # small g, and for larger g use the standard approximation: ancestry
    # junctions at (g-1)/M with iid labels, plus same-pop breakpoints so the
    # total breakpoint density is g/M.
    if gen <= 4:
        P, B = gen1(nhap * 2 ** (gen - 1))
        width = P.shape[1]
        while width > nhap:
            half = width // 2
            a_P, b_P = P[:, :half], P[:, half:]
            a_B, b_B = B[:, :half], B[:, half:]
            sw = rng.random((nm, half)) < (1 - math.exp(-STEP))
            sw[0, :] = False
            start = rng.random(half) < 0.5
            use_b = (np.cumsum(sw, axis=0) % 2 == 1) ^ start[None, :]
            P = np.where(use_b, b_P, a_P)
            B = np.where(use_b, b_B, a_B) | sw
            width = half
        B[0, :] = True
        return P, B
    lam_j = (gen - 1) * STEP
    junction = rng.random((nm, nhap)) < (1 - math.exp(-lam_j))
    junction[0, :] = True
    seg = np.cumsum(junction, axis=0) - 1
    nseg = seg.max(axis=0) + 1
    labels = rng.choice(k, size=(int(nseg.max()), nhap), p=p)
    P = labels[seg, np.arange(nhap)[None, :]].astype(np.int8)
    extra = rng.random((nm, nhap)) < (1 - math.exp(-STEP))
    B = junction | extra
    return P, B


def run_simumix(rng, P, B, comps, want_codes=False):
    """Draw donors marker by marker; return (died, codes or None)."""
    nm, nh = P.shape
    ndon = [2 * DONORS[c] for c in comps]
    busy = [np.zeros(n, dtype=bool) for n in ndon]
    cur = np.full(nh, -1, dtype=np.int64)
    curpop = np.full(nh, -1, dtype=np.int64)
    codes = np.empty((nm, nh), dtype=np.int32) if want_codes else None
    for m in range(nm):
        need = np.nonzero(B[m] | (P[m] != curpop))[0]
        # release first (as simu-mix does per haplotype, in order; releasing
        # all first is the more permissive reading, so deaths here are a floor)
        for h in need:
            if cur[h] >= 0:
                busy[curpop[h]][cur[h]] = False
        for h in need:
            pp = int(P[m, h])
            free = np.nonzero(~busy[pp])[0]
            if curpop[h] == pp and cur[h] >= 0:
                free = free[free != cur[h]]
            if free.size == 0:
                return True, None
            d = int(free[rng.integers(free.size)])
            busy[pp][d] = True
            cur[h] = d
            curpop[h] = pp
        if want_codes:
            codes[m] = curpop * 1000 + cur
    return False, codes


def die_check(rng, comps, props, gen, size, length_m):
    """Fast check without donor identities: a run dies when the number of
    haplotypes copying a pop at a marker exceeds its donor haplotypes, or
    equals it where one of them breaks within the pop (it may not re-take
    its own donor)."""
    P, B = sim_haps(rng, 2 * size, length_m, props, gen)
    prev = np.vstack([np.full((1, P.shape[1]), -1, dtype=P.dtype), P[:-1]])
    for i, c in enumerate(comps):
        H = 2 * DONORS[c]
        cnt = (P == i).sum(axis=1)
        if (cnt > H).any():
            return True
        # a breakpoint that stays in pop i may not re-take its own donor
        same_brk = ((P == i) & (prev == i) & B).any(axis=1)
        if ((cnt == H) & same_brk).any():
            return True
    return False


def part_a(reps, seed=1):
    rng = np.random.default_rng(seed)
    print("PART A -- does a 6b chromosome task die? (reps per rule: %d)" % reps)
    print(f"{'model':<8}{'rule':<10}{'cap':>5}{'batches':>8}"
          f"{'P(task dies) per chr, mean':>28}{'P(any of 22 dies)':>19}")
    for name, comps, props, gen, n in MODELS:
        for rule, cap in (("now 0.8", cap_margin(comps, props)),
                          ("z=5", cap_z(comps, props)),
                          ("N+1", cap_n1(comps, props))):
            bs = batches(n, cap)
            any_die = 0
            per_chr = np.zeros(22)
            for _ in range(reps):
                died_run = False
                for ci, L in enumerate(CHR_M):
                    if any(die_check(rng, comps, props, gen, s, L) for s in bs):
                        per_chr[ci] += 1
                        died_run = True
                any_die += died_run
            print(f"{name:<8}{rule:<10}{cap:>5}{len(bs):>8}"
                  f"{per_chr.mean() / reps:>28.4f}{any_die / reps:>19.3f}")
        sys.stdout.flush()


def kinship_model(rng, comps, props, gen, n, cap):
    bs = batches(n, cap)
    codes_all, width_all = [], []
    for L in CHR_M:
        cols = []
        for s in bs:
            for attempt in range(50):
                P, B = sim_haps(rng, 2 * s, L, props, gen)
                died, codes = run_simumix(rng, P, B, comps, want_codes=True)
                if not died:
                    break
            else:
                raise SystemExit("could not simulate a batch in 50 tries")
            cols.append(codes)
        codes_all.append(np.concatenate(cols, axis=1))
    C = np.concatenate(codes_all, axis=0)
    total = C.shape[0]
    nh = C.shape[1]
    f = np.zeros((nh, nh))
    for h in range(nh - 1):
        eq = (C[:, h][:, None] == C[:, h + 1:])
        f[h, h + 1:] = eq.sum(axis=0) / total
    f = f + f.T
    ind = np.arange(nh) // 2
    ni = nh // 2
    phi = np.zeros((ni, ni))
    for a in range(2):
        for b in range(2):
            phi += f[a::2, :][:, b::2]
    phi /= 4.0
    np.fill_diagonal(phi, np.nan)
    batch_of = np.concatenate([[i] * s for i, s in enumerate(bs)])
    return phi, batch_of


def kinship_controls(rng, comp, n):
    H = 2 * DONORS[comp]
    tot = sum(CHR_M)
    phi = np.zeros((n, n))
    for L in CHR_M:
        a = rng.integers(H, size=n)
        b = rng.integers(H, size=n)
        while (a == b).any():
            same = a == b
            b[same] = rng.integers(H, size=same.sum())
        haps = np.stack([a, b], axis=1)
        m = np.zeros((n, n))
        for i in range(2):
            for j in range(2):
                m += haps[:, i][:, None] == haps[:, j][None, :]
        phi += m * (L / tot)
    phi /= 4.0
    np.fill_diagonal(phi, np.nan)
    return phi


def summarise(label, phi, batch_of=None):
    iu = np.triu_indices_from(phi, k=1)
    v = phi[iu]
    line = (f"{label:<28} pairs {v.size:>6}  mean {np.nanmean(v):.4f}  "
            f"p99 {np.nanpercentile(v, 99):.4f}  max {np.nanmax(v):.4f}  "
            f">=0.0221 {int((v >= 0.0221).sum()):>4}  >=0.0442 "
            f"{int((v >= 0.0442).sum()):>3}  >=0.0884 {int((v >= 0.0884).sum()):>2}  "
            f">=0.177 {int((v >= 0.177).sum()):>2}")
    print(line)
    if batch_of is not None:
        same = batch_of[iu[0]] == batch_of[iu[1]]
        print(f"{'':<28} within-batch mean {np.nanmean(v[same]):.4f}, "
              f"between-batch mean {np.nanmean(v[~same]):.4f}")


def part_b(seed=2):
    rng = np.random.default_rng(seed)
    print("\nPART B -- kinship between simulated individuals (one full run,"
          " all 22 chromosomes)")
    print("KING bins: >=0.0442 3rd degree, >=0.0884 2nd, >=0.177 1st;"
          " 0.0221 is the 4th-degree edge")
    for name, comps, props, gen, n in MODELS:
        if name not in ("MENEUR", "EASEUR"):
            continue
        for rule, cap in (("now 0.8", cap_margin(comps, props)),
                          ("z=5", cap_z(comps, props))):
            phi, bo = kinship_model(rng, comps, props, gen, n, cap)
            summarise(f"{name} gen{gen} {rule} cap {cap}", phi, bo)
            sys.stdout.flush()
    for comp in ("EUR", "MEN", "AFR", "EAS", "NAT"):
        phi = kinship_controls(rng, comp, 60)
        summarise(f"{comp}1 gen0 controls n=60", phi)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "ab"
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    if "a" in what:
        part_a(reps)
    if "b" in what:
        part_b()
