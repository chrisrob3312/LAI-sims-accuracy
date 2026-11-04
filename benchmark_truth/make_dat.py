#!/usr/bin/env python3
"""Write an admix-simu .dat file for one cohort.

admix-simu (williamslab/admix-simu) .dat layout, as used by simu-mix.pl:
    <n_admixed_samples> ADMIX <POP1> <POP2> ... <POPk>
    <generation> <prop1> <prop2> ... <propk>
    (one row per generation block; a single row = one pulse of admixture
    <generation> generations ago)

K=1 (a "homogeneous" cohort) is the same format with one population and
proportion 1.0: admix-simu then recombines haplotypes from that single donor
pool for <generation> generations, i.e. within-population splicing.

Verify the header line against an existing .dat in $PROJECT_ROOT/02_simulations
(e.g. Brasa.dat) the first time this is run; the population labels must match
the -LABEL flags given to simu-mix.pl (5b does this automatically).

Usage:
    make_dat.py --n 100 --pops EUR AMR AFR --props 0.5 0.45 0.05 --gen 12 --out ADM_MEX.dat
"""
import argparse

p = argparse.ArgumentParser()
p.add_argument("--n", type=int, required=True)
p.add_argument("--pops", nargs="+", required=True)
p.add_argument("--props", nargs="+", type=float, required=True)
p.add_argument("--gen", type=int, required=True)
p.add_argument("--out", required=True)
a = p.parse_args()

if len(a.pops) != len(a.props):
    raise SystemExit("pops and props must have the same length")
tot = sum(a.props)
if abs(tot - 1.0) > 1e-6:
    raise SystemExit(f"proportions must sum to 1 (got {tot})")

with open(a.out, "w") as f:
    f.write(f"{a.n} ADMIX " + " ".join(a.pops) + "\n")
    f.write(f"{a.gen} " + " ".join(f"{x:g}" for x in a.props) + "\n")
print(f"[make_dat] wrote {a.out}: n={a.n} pops={a.pops} props={a.props} gen={a.gen}")
