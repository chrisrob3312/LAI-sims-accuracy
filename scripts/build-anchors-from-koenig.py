#!/usr/bin/env python3
"""Build supervised-ADMIXTURE anchor lists from Koenig/gnomAD flat metadata.

Filters:
  - high_quality == true
  - population_inference.prob_<pop> >= thresh
  - cap N per super-pop

--amr-mode koenig-plus-mxb : gnomAD-AMR top N ∪ MXB samples (previous default)
--amr-mode mxb-only        : ONLY MXB samples (cleaner Amerindigenous centroid)

Maps: afr->AFR  amr->AMR  eas->EAS  nfe->EUR  sas->SAS
"""
import sys, os, argparse
from collections import defaultdict

POP2SUPER = {"afr": "AFR", "amr": "AMR", "eas": "EAS", "nfe": "EUR", "sas": "SAS"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta_tsv")
    ap.add_argument("--out-dir", default="reference_ids")
    ap.add_argument("--thresh", type=float, default=0.95)
    ap.add_argument("--cap", type=int, default=50)
    ap.add_argument("--mxb", default="reference_ids/MXB50genomes_popinfo.tsv")
    ap.add_argument("--amr-mode", choices=["koenig-plus-mxb", "mxb-only"],
                    default="mxb-only")
    args = ap.parse_args()

    by_pop = defaultdict(list)
    with open(args.meta_tsv) as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        c_s = idx["s"]
        c_hq = idx["high_quality"]
        c_pop = idx["population_inference.pop"]
        c_prob = {p: idx[f"population_inference.prob_{p}"] for p in POP2SUPER}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[c_hq].lower() != "true": continue
            pop = f[c_pop].lower()
            if pop not in POP2SUPER: continue
            try: prob = float(f[c_prob[pop]])
            except (ValueError, IndexError): continue
            if prob < args.thresh: continue
            by_pop[POP2SUPER[pop]].append((f[c_s], prob))

    os.makedirs(args.out_dir, exist_ok=True)
    mxb_ids = []
    if os.path.isfile(args.mxb):
        with open(args.mxb) as fh:
            _ = next(fh)
            mxb_ids = [ln.split("\t")[0].strip() for ln in fh if ln.strip()]
        print(f"MXB: {len(mxb_ids)} samples available", file=sys.stderr)

    print(f"{'super-pop':10s} {'qualified':>10s} {'written':>10s}", file=sys.stderr)
    for sp in ("EUR","AFR","AMR","EAS","SAS"):
        cands = sorted(by_pop[sp], key=lambda x: -x[1])
        if sp == "AMR":
            if args.amr_mode == "mxb-only":
                chosen = sorted(mxb_ids)
            else:
                chosen = sorted(set([s for s,_ in cands[:args.cap]]) | set(mxb_ids))
        else:
            chosen = [s for s,_ in cands[:args.cap]]
        outp = os.path.join(args.out_dir, f"{sp.lower()}_rfmix.txt")
        with open(outp, "w") as fh:
            for s in chosen: fh.write(s + "\n")
        print(f"{sp:10s} {len(cands):>10d} {len(chosen):>10d}", file=sys.stderr)

if __name__ == "__main__": main()
