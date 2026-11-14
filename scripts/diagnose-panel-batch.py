#!/usr/bin/env python3
"""diagnose-panel-batch.py -- is a cohort, not an ancestry, driving the K=8 fit?

Read-only. Every mode reads text that bcftools / plink2 / ADMIXTURE already
wrote, so the raw inputs can be inspected next to the summary. Nothing here
writes into the pipeline's directories.

WHY THIS EXISTS. 4a job 3747721 found zero anchors. Its discovery log shows
SGDP-only label pools (the seed_for bug) AND one cluster carrying a large share
of SGDP samples from every continent: AFR meanQ 0.551 + 0.402 on cluster 7, AMR
0.721 + 0.266 on cluster 7, OCE 0.747 on cluster 7. SGDP is the only cohort
that arrives through its own pipeline (Reich phased_data2021, hg19, Picard
liftover, +fixref). A cross-continental cluster that follows a cohort is the
signature of a processing batch, not of ancestry. These modes test that
directly, most decisive first:

  q        Q profile by super-pop x cohort, who loads on each cluster, the
           curated lists' Q, a preview of what FIXED 4a discovery would find
           from the EXISTING unsupervised Q, and -- the sharpest test -- the
           same person present twice (two pipelines) compared side by side.
  psc      per-sample missingness / heterozygosity from `bcftools stats -s -`
           on the PRE-phasing soft-union BCF (raw calls, before SHAPEIT5
           imputes anything), by super-pop x cohort.
  het      per-sample F from `plink2 --het` on the ADMIXTURE input (what
           ADMIXTURE actually saw, after imputation), by super-pop x cohort.
  pair-ids IDs of both copies of every duplicate pair (feed to bcftools -S).
  pairs-gt genotype concordance between the two copies, from
           `bcftools query -H -f '%POS[\\t%GT]\\n'` on stdin.

Cohort = explicit membership lists first (--cohort SGDP=ids.txt ...), then ID
form: ^(HG|NA)<digit> 1KG, ^HGDP HGDP, ^(LP6005|SS600) HGDP_LPSS (HGDP
individuals sequenced by SGDP and called by gnomAD -- they come in with the
HGDP+1KG source), ^MXB MXB, *.DG SGDP, longranger/APPG ALMARRI, ^CHMI CHMI,
else OTHER (listed, never silently bucketed).

Super-pop = SGDP manifest for SGDP samples, sample_groups.tsv for the rest
(CSA->SAS, MID->MEN), MEN for ALMARRI, AMR for MXB, else NA (listed).
"""
import argparse
import csv
import re
import statistics
import sys
from collections import Counter, defaultdict

BASE_COHORTS = ("1KG", "HGDP", "HGDP_LPSS", "MXB")
OTHER_PIPELINES = ("SGDP", "ALMARRI")
SP_ORDER = ["AFR", "AMR", "EAS", "EUR", "MEN", "OCE", "SAS", "NA"]
ALMARRI_RE = re.compile(r"(^longranger\d*_wgs_|(^|_)APPG\d+(_|$))")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def read_ids(path):
    try:
        with open(path) as fh:
            return [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        die(f"not found: {path}")


# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------
class Labels:
    def __init__(self, a):
        self.explicit = []  # (name, set) in the order given
        for spec in a.cohort or []:
            if "=" not in spec:
                die(f"--cohort wants NAME=FILE, got {spec!r}")
            name, path = spec.split("=", 1)
            ids = set(read_ids(path))
            self.explicit.append((name.upper(), ids))
            print(f"[labels] cohort {name.upper()}: {len(ids)} ids from {path}", file=sys.stderr)
        self.sgdp_sp = {}
        if a.sgdp_manifest:
            with open(a.sgdp_manifest) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    self.sgdp_sp[r["sample"]] = r["super_pop"]
            print(f"[labels] SGDP manifest: {len(self.sgdp_sp)} samples", file=sys.stderr)
        self.groups = {}
        if a.groups:
            with open(a.groups) as fh:
                for line in fh:
                    ff = line.split()
                    if len(ff) >= 2:
                        g = ff[1].upper()
                        self.groups[ff[0]] = {"CSA": "SAS", "MID": "MEN"}.get(g, g)
            print(f"[labels] sample groups: {len(self.groups)} samples "
                  f"{dict(Counter(self.groups.values()))}", file=sys.stderr)

    def cohort(self, s):
        for name, ids in self.explicit:
            if s in ids:
                return name
        if re.match(r"^(HG|NA)\d", s): return "1KG"
        if s.startswith("HGDP"): return "HGDP"
        if s.startswith(("LP6005", "SS600")): return "HGDP_LPSS"
        if s.startswith("MXB"): return "MXB"
        if s.endswith(".DG"): return "SGDP"
        if ALMARRI_RE.search(s): return "ALMARRI"
        if s.startswith("CHMI"): return "CHMI"
        return "OTHER"

    def superpop(self, s, cohort=None):
        cohort = cohort or self.cohort(s)
        if s in self.sgdp_sp: return self.sgdp_sp[s]
        if s in self.groups: return self.groups[s]
        if cohort == "ALMARRI": return "MEN"
        if cohort == "MXB": return "AMR"
        return "NA"


def sp_key(sp):
    return (SP_ORDER.index(sp) if sp in SP_ORDER else len(SP_ORDER), sp)


def label_table(samples, L):
    co = {s: L.cohort(s) for s in samples}
    sp = {s: L.superpop(s, co[s]) for s in samples}
    tab = Counter((co[s], sp[s]) for s in samples)
    cohorts = sorted({c for c, _ in tab})
    sps = sorted({p for _, p in tab}, key=sp_key)
    print("\n== 1. LABELS: cohort x super-pop (every sample in the panel) ==")
    print(f"{'cohort':10s}" + "".join(f"{p:>6s}" for p in sps) + f"{'total':>7s}")
    for c in cohorts:
        row = [tab[(c, p)] for p in sps]
        print(f"{c:10s}" + "".join(f"{n:6d}" for n in row) + f"{sum(row):7d}")
    na = [s for s in samples if sp[s] == "NA"]
    if na:
        print(f"  super-pop NA: {len(na)} sample(s) -- unlabelled, so excluded downstream:")
        byc = defaultdict(list)
        for s in na: byc[co[s]].append(s)
        for c, ss in sorted(byc.items()):
            print(f"    {c}: {len(ss)}  e.g. {', '.join(ss[:4])}")
    other = [s for s in samples if co[s] == "OTHER"]
    if other:
        print(f"  cohort OTHER (matched no rule): {len(other)}  e.g. {', '.join(other[:6])}")
    return co, sp


# --------------------------------------------------------------------------
# q mode
# --------------------------------------------------------------------------
def fmt_q(v):
    return " ".join(f"{x:4.2f}" for x in v)


def mean_vec(rows, K):
    if not rows: return [0.0] * K
    return [sum(r[c] for r in rows) / len(rows) for c in range(K)]


def greedy_map(label_pools, samples, Q, K):
    """4a's cluster naming, verbatim in logic: exemplars in exactly one label
    pool, mean-Q per (super-pop, cluster), one-to-one greedy by max cell."""
    owners = defaultdict(set)
    for sp, ids in label_pools.items():
        for s in ids: owners[s].add(sp)
    ex = {s: next(iter(v)) for s, v in owners.items() if len(v) == 1}
    idx = {s: i for i, s in enumerate(samples)}
    sums = {sp: [0.0] * K for sp in label_pools}; cnt = Counter()
    for s, sp in ex.items():
        if s in idx:
            cnt[sp] += 1
            for c in range(K): sums[sp][c] += Q[idx[s]][c]
    meanQ = {sp: [sums[sp][c] / cnt[sp] if cnt[sp] else 0.0 for c in range(K)] for sp in label_pools}
    cells = sorted(((meanQ[sp][c], sp, c) for sp in label_pools for c in range(K)), reverse=True)
    m, used_c, used_s = {}, set(), set()
    for v, sp, c in cells:
        if sp in used_s or c in used_c: continue
        m[sp] = c; used_s.add(sp); used_c.add(c)
    return m, meanQ, cnt


def load_pairs(a, present):
    """Unordered pairs (id1, id2) -> info, from KING's table and the metadata
    alias table. Only pairs with both copies in `present` are kept."""
    pairs = {}
    if a.dups:
        with open(a.dups) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                x, y = r["id1"], r["id2"]
                k = tuple(sorted((x, y)))
                pairs.setdefault(k, {"kinship": float(r["kinship"]), "src": set()})["src"].add("KING")
    if a.known_dups:
        with open(a.known_dups) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                x, y = r["sgdp_id"], r["alias"]
                k = tuple(sorted((x, y)))
                pairs.setdefault(k, {"kinship": None, "src": set()})["src"].add("metadata")
    kin = {}
    if a.kin0:
        want = set(pairs)
        with open(a.kin0) as fh:
            h = fh.readline().lstrip("#").rstrip("\n").split("\t")
            ix = {n: i for i, n in enumerate(h)}
            i1 = ix.get("IID1", ix.get("ID1")); i2 = ix.get("IID2", ix.get("ID2")); ik = ix.get("KINSHIP")
            if None in (i1, i2, ik): die(f"{a.kin0}: no IID1/IID2/KINSHIP in header {h}")
            for line in fh:
                ff = line.rstrip("\n").split("\t")
                k = tuple(sorted((ff[i1], ff[i2])))
                if k in want: kin[k] = float(ff[ik])
    out = {}
    for k, v in pairs.items():
        if k[0] in present and k[1] in present:
            if v["kinship"] is None and k in kin: v["kinship"] = kin[k]
            out[k] = v
    return out


def orient(x, y, co):
    """(base copy, other copy): other = the SGDP/ALMARRI pipeline copy."""
    cx, cy = co.get(x, "?"), co.get(y, "?")
    if cx in OTHER_PIPELINES and cy not in OTHER_PIPELINES: return y, x
    if cy in OTHER_PIPELINES and cx not in OTHER_PIPELINES: return x, y
    return (x, y) if (cx, x) <= (cy, y) else (y, x)


def mode_q(a):
    L = Labels(a)
    fam = [l.split()[1] for l in open(a.fam) if l.strip()]
    Q = [list(map(float, l.split())) for l in open(a.q) if l.strip()]
    if len(fam) != len(Q): die(f"{a.fam} has {len(fam)} rows but {a.q} has {len(Q)}")
    K = len(Q[0])
    th = [float(t) for t in a.thresholds.split(",")]
    print(f"[q] {len(fam)} samples, K={K}, thresholds {th}")
    print(f"[q] raw first rows of {a.q}:")
    for s, q in list(zip(fam, Q))[:3]:
        print(f"      {s:40s} {fmt_q(q)}")
    co, sp = label_table(fam, L)
    idx = {s: i for i, s in enumerate(fam)}
    excl = set(read_ids(a.exclude)) if a.exclude else set()

    # 2. profile
    print(f"\n== 2. MEAN Q PROFILE by super-pop x cohort (clusters c0..c{K-1}) ==")
    print("   A cohort row that loads on a cluster its same-super-pop rows from")
    print("   other cohorts do NOT load on is the batch signature.")
    hdr = " ".join(f"  c{c}" for c in range(K))
    print(f"{'sp':4s} {'cohort':10s} {'n':>5s}  {hdr}  " + " ".join(f">={t:<4}" for t in th) + "  top")
    groups = defaultdict(list)
    for s in fam: groups[(sp[s], co[s])].append(Q[idx[s]])
    for (p, c), rows in sorted(groups.items(), key=lambda kv: (sp_key(kv[0][0]), kv[0][1])):
        mv = mean_vec(rows, K)
        top = max(range(K), key=lambda j: mv[j])
        npass = [sum(1 for r in rows if max(r) >= t) for t in th]
        print(f"{p:4s} {c:10s} {len(rows):5d}  {fmt_q(mv)}  " + " ".join(f"{n:6d}" for n in npass) + f"  c{top}")

    # 3. cluster membership
    print("\n== 3. WHO LOADS ON EACH CLUSTER (samples with Q_c >= 0.5) ==")
    for c in range(K):
        hit = [s for s in fam if Q[idx[s]][c] >= 0.5]
        cc = Counter(co[s] for s in hit); ss = Counter(sp[s] for s in hit)
        mean_by_c = {k: statistics.mean(Q[idx[s]][c] for s in fam if co[s] == k)
                     for k in sorted(set(co.values()))}
        print(f"  c{c}: {len(hit):5d} samples  cohorts {dict(cc.most_common())}")
        print(f"        super-pops {dict(ss.most_common(5))}")
        print(f"        mean Q_c by cohort: " + ", ".join(f"{k} {v:.2f}" for k, v in mean_by_c.items()))

    # 4. curated lists
    curated = {}
    for spec in a.curated or []:
        name, path = spec.split("=", 1)
        curated[name.upper()] = read_ids(path)
    if curated:
        print("\n== 4. CURATED *_rfmix.txt LISTS (the samples job 3747721 could not see) ==")
        print(f"{'list':5s} {'n':>3s} {'inPanel':>7s}  {hdr}  best  " + " ".join(f">={t:<4}" for t in th))
        for name, ids in curated.items():
            got = [Q[idx[s]] for s in ids if s in idx]
            mv = mean_vec(got, K)
            b = max(range(K), key=lambda j: mv[j]) if got else 0
            npass = [sum(1 for r in got if r[b] >= t) for t in th]
            print(f"{name:5s} {len(ids):3d} {len(got):7d}  {fmt_q(mv)}  c{b}    " + " ".join(f"{n:6d}" for n in npass))

    # 5. preview of fixed discovery
    if curated:
        print("\n== 5. PREVIEW: what FIXED 4a discovery would find from THIS unsupervised Q ==")
        print(f"   label pools: curated lists (+MXB for AMR); OCE/MEN: every panel member")
        print(f"   draw pools : curated + every panel member with that metadata super-pop")
        print(f"   anchors    : top {a.cap} at Q >= {th[0]} on the mapped cluster; floor {a.floor}")
        named = ["EUR", "AFR", "AMR", "EAS", "SAS", "OCE", "MEN"]
        label = {}
        for p in named:
            if p in curated:
                ids = set(curated[p])
                if p == "AMR": ids |= {s for s in fam if co[s] == "MXB"}
            else:
                ids = {s for s in fam if sp[s] == p}
            label[p] = ids - excl
        m, meanQ, cnt = greedy_map(label, fam, Q, K)
        print(f"   cluster map: { {p: f'c{c}' for p, c in sorted(m.items(), key=lambda kv: kv[1])} }")
        for variant, drop_sgdp in (("all cohorts", False), ("SGDP removed from draw pool", True)):
            print(f"\n   -- draw pool: {variant} --")
            print(f"   {'sp':4s} {'clus':4s} {'meanQ':>6s} {'runner-up':>12s} {'pool':>5s} "
                  f"{'>='+str(th[0]):>7s} {'>='+str(th[-1]):>7s} {'kept':>5s}  kept-by-cohort")
            for p in named:
                if p not in m:
                    print(f"   {p:4s}  --  claimed no cluster"); continue
                c = m[p]
                ru = sorted(((meanQ[p][j], j) for j in range(K) if j != c), reverse=True)[0]
                pool = (set(curated.get(p, [])) | {s for s in fam if sp[s] == p}) - excl
                if drop_sgdp: pool = {s for s in pool if co.get(s) != "SGDP"}
                pool = [s for s in pool if s in idx]
                q = sorted(((Q[idx[s]][c], s) for s in pool), reverse=True)
                p0 = [x for x in q if x[0] >= th[0]]; p1 = [x for x in q if x[0] >= th[-1]]
                kept = p0[:a.cap]
                kc = Counter(co[s] for _, s in kept)
                flag = "" if len(kept) >= a.floor else "  <-- below floor"
                print(f"   {p:4s} c{c:<3d} {meanQ[p][c]:6.3f} {'c'+str(ru[1])+'@'+format(ru[0], '.3f'):>12s} "
                      f"{len(pool):5d} {len(p0):7d} {len(p1):7d} {len(kept):5d}  {dict(kc.most_common())}{flag}")
            left = [j for j in range(K) if j not in m.values()]
            for j in left:
                q = sorted(((Q[idx[s]][j], s) for s in fam if s not in excl
                            and not (drop_sgdp and co[s] == "SGDP")), reverse=True)
                p0 = [x for x in q if x[0] >= th[0]]
                print(f"   OTH  c{j:<3d} leftover: {len(p0)} samples >= {th[0]}; "
                      f"by cohort {dict(Counter(co[s] for _, s in p0[:a.cap]).most_common())}; "
                      f"by super-pop {dict(Counter(sp[s] for _, s in p0[:a.cap]).most_common(4))}")

    # 6. same individual, two pipelines
    if a.dups or a.known_dups:
        pairs = load_pairs(a, set(idx))
        print(f"\n== 6. SAME PERSON TWICE: Q of the two copies ({len(pairs)} pairs with both in the panel) ==")
        print("   TVD = half the summed |Q difference| (0 = identical, 1 = disjoint).")
        print("   Same genome through two pipelines should give TVD ~0.0x. A large TVD")
        print("   with one cluster gained by the SGDP/Almarri copy = processing batch.")
        bytype = defaultdict(list)
        rows = []
        for (x, y), v in pairs.items():
            b, o = orient(x, y, co)
            qb, qo = Q[idx[b]], Q[idx[o]]
            d = [qo[c] - qb[c] for c in range(K)]
            tvd = 0.5 * sum(abs(t) for t in d)
            j = max(range(K), key=lambda c: abs(d[c]))
            t = f"{co[o]}~{co[b]}"
            bytype[t].append(d)
            rows.append((tvd, b, o, v, d, j, t))
        for tvd, b, o, v, d, j, t in sorted(rows, reverse=True)[: a.show_pairs]:
            k = "n/a" if v["kinship"] is None else f"{v['kinship']:.3f}"
            print(f"   {tvd:5.2f}  {o} ({co[o]}) vs {b} ({co[b]})  kin {k}  [{'+'.join(sorted(v['src']))}]"
                  f"  biggest: c{j} {d[j]:+.2f}")
            print(f"          {co[o]:>9s}: {fmt_q(Q[idx[o]])}")
            print(f"          {co[b]:>9s}: {fmt_q(Q[idx[b]])}")
        if len(rows) > a.show_pairs:
            print(f"   ... {len(rows) - a.show_pairs} more pairs (raise --show-pairs to see them)")
        print("\n   summary by pair type (mean signed Q change, other copy minus base copy):")
        print(f"   {'type':22s} {'n':>3s} {'meanTVD':>8s} {'medTVD':>7s}  " + " ".join(f"   c{c}" for c in range(K)))
        for t, ds in sorted(bytype.items()):
            tv = [0.5 * sum(abs(x) for x in d) for d in ds]
            mv = [sum(d[c] for d in ds) / len(ds) for c in range(K)]
            print(f"   {t:22s} {len(ds):3d} {statistics.mean(tv):8.3f} {statistics.median(tv):7.3f}  "
                  + " ".join(f"{x:+5.2f}" for x in mv))
        meta_only = [(k, v) for k, v in pairs.items() if v["src"] == {"metadata"}]
        if meta_only:
            print(f"\n   metadata says SAME person but KING did NOT call a duplicate ({len(meta_only)}):")
            print("   (kinship near 0.5 = KING table was filtered; well below 0.354 = the")
            print("    two copies' genotypes disagree, or the alias is wrong)")
            for (x, y), v in sorted(meta_only, key=lambda kv: (kv[1]["kinship"] is None, kv[1]["kinship"] or 0)):
                k = "not in kin0" if v["kinship"] is None else f"{v['kinship']:.3f}"
                print(f"     {x} ({co[x]}) ~ {y} ({co[y]})  kinship {k}")

    # 7. marker set vs July
    if a.old_bim and a.new_bim:
        def bimset(p):
            s, sp2 = set(), set()
            with open(p) as fh:
                for l in fh:
                    f = l.split()
                    ch = f[0][3:] if f[0].startswith("chr") else f[0]
                    s.add((ch, f[3], tuple(sorted((f[4], f[5])))))
                    sp2.add((ch, f[3]))
            return s, sp2
        o, op = bimset(a.old_bim); n, np_ = bimset(a.new_bim)
        print("\n== 7. ADMIXTURE MARKER SET: now vs the July K=5 run ==")
        print(f"   July {len(o):,} SNPs | now {len(n):,} SNPs | shared (pos+alleles) {len(o & n):,}"
              f" | shared positions {len(op & np_):,}")
        print(f"   {100 * len(o & n) / max(len(n), 1):.1f}% of today's markers were in July's set")


# --------------------------------------------------------------------------
# psc / het
# --------------------------------------------------------------------------
def summarize(per, L, metric_cols, title, extremes):
    co = {s: L.cohort(s) for s in per}
    sp = {s: L.superpop(s, co[s]) for s in per}
    g = defaultdict(list)
    for s, v in per.items(): g[(sp[s], co[s])].append(v)
    print(f"\n== {title} ==")
    print(f"{'sp':4s} {'cohort':10s} {'n':>5s}  " + "  ".join(f"{name:>14s}" for name, _ in metric_cols))
    for (p, c), vs in sorted(g.items(), key=lambda kv: (sp_key(kv[0][0]), kv[0][1])):
        cells = []
        for name, key in metric_cols:
            xs = [v[key] for v in vs if v[key] is not None]
            cells.append(f"{statistics.mean(xs):8.4f}+-{(statistics.pstdev(xs) if len(xs) > 1 else 0):.3f}"
                         if xs else f"{'n/a':>14s}")
        print(f"{p:4s} {c:10s} {len(vs):5d}  " + "  ".join(f"{x:>14s}" for x in cells))
    for name, key, rev in extremes:
        xs = sorted(((v[key], s) for s, v in per.items() if v[key] is not None), reverse=rev)[:8]
        print(f"  {'highest' if rev else 'lowest'} {name}: " +
              ", ".join(f"{s} ({co[s]}/{sp[s]}) {x:.4f}" for x, s in xs))


def mode_psc(a):
    L = Labels(a)
    n_records = None; cols = None; per = {}
    with open(a.stats) as fh:
        for line in fh:
            if line.startswith("SN\t") and "number of records:" in line:
                n_records = int(line.rstrip("\n").split("\t")[-1])
            elif line.startswith("# PSC\t"):
                cols = [re.sub(r"^\[\d+\]", "", x) for x in line.rstrip("\n").split("\t")]
                print(f"[psc] PSC header: {cols}")
            elif line.startswith("PSC\t"):
                if cols is None: die("PSC rows before the '# PSC' header -- not bcftools stats -s output?")
                f = dict(zip(cols, line.rstrip("\n").split("\t")))
                rr, ra, hh = int(f["nRefHom"]), int(f["nNonRefHom"]), int(f["nHets"])
                called = rr + ra + hh
                if "nMissing" in f:
                    miss = int(f["nMissing"])
                elif n_records is not None:
                    miss = n_records - called
                else:
                    miss = None
                tot = called + (miss or 0)
                per[f["sample"]] = {
                    "miss": (miss / tot) if (miss is not None and tot) else None,
                    "het": hh / called if called else None,
                    "homalt": ra / called if called else None,
                }
    if not per: die(f"no PSC rows in {a.stats}")
    print(f"[psc] {len(per)} samples; SN number of records: {n_records}")
    first = next(iter(per.items()))
    print(f"[psc] first sample parsed: {first[0]} {first[1]}")
    summarize(per, L, [("missing", "miss"), ("het/called", "het"), ("homalt/called", "homalt")],
              "PRE-PHASING per-sample rates by super-pop x cohort (mean +- sd)",
              [("missing", "miss", True), ("het", "het", True), ("het", "het", False)])


def mode_het(a):
    L = Labels(a)
    per = {}
    with open(a.het) as fh:
        h = fh.readline().lstrip("#").split()
        print(f"[het] header: {h}")
        ix = {n: i for i, n in enumerate(h)}
        for need in ("IID", "O(HOM)", "OBS_CT", "F"):
            if need not in ix: die(f"{a.het}: no {need} column in {h}")
        for line in fh:
            f = line.split()
            if not f: continue
            ohom, obs = float(f[ix["O(HOM)"]]), float(f[ix["OBS_CT"]])
            per[f[ix["IID"]]] = {"F": float(f[ix["F"]]), "het": 1 - ohom / obs if obs else None}
    print(f"[het] {len(per)} samples")
    summarize(per, L, [("F", "F"), ("het rate", "het")],
              "ADMIXTURE INPUT (post-phasing) F and heterozygosity by super-pop x cohort (mean +- sd)",
              [("F", "F", True), ("F", "F", False)])


# --------------------------------------------------------------------------
# pairs
# --------------------------------------------------------------------------
def mode_pair_ids(a):
    fam = {l.split()[1] for l in open(a.fam) if l.strip()}
    pairs = load_pairs(a, fam)
    ids = sorted({x for k in pairs for x in k})
    print(f"[pair-ids] {len(pairs)} pairs, {len(ids)} ids (both copies in the panel)", file=sys.stderr)
    sys.stdout.write("\n".join(ids) + "\n")


GT_CODE = {}


def gt_code(g):
    c = GT_CODE.get(g)
    if c is None:
        al = re.split(r"[/|]", g)
        c = None if any(x in (".", "") for x in al) else sum(1 for x in al if x != "0")
        GT_CODE[g] = c
    return c


def mode_pairs_gt(a):
    L = Labels(a)
    hdr = sys.stdin.readline()
    if not hdr.startswith("#"):
        die("first line is not a bcftools query -H header; run bcftools query with -H")
    cols = []
    for f in hdr.rstrip("\n").split("\t"):
        f = re.sub(r"^#\s*", "", f).strip()
        f = re.sub(r"^\[\d+\]", "", f)
        cols.append(f[:-3] if f.endswith(":GT") else f)
    names = cols[1:]
    print(f"[pairs-gt] header: {cols[:4]} ... ({len(names)} samples)")
    pos = {s: i for i, s in enumerate(names)}
    pairs = load_pairs(a, set(pos))
    co = {s: L.cohort(s) for s in names}
    P = []
    for (x, y) in pairs:
        b, o = orient(x, y, co)
        P.append((b, o, pos[b], pos[o]))
    if not P: die("no pair has both copies in the bcftools output")
    T = [Counter() for _ in P]
    n = 0
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")[1:]
        n += 1
        for t, (_, _, ib, io) in zip(T, P):
            t[(gt_code(f[ib]), gt_code(f[io]))] += 1
    print(f"[pairs-gt] {n:,} sites, {len(P)} pairs")
    print("   base = 1KG/HGDP/MXB copy, other = SGDP/Almarri copy. disc = discordant among")
    print("   sites both called. swap = 0/0 vs 1/1 (allele flip). Replicates of one genome")
    print("   through one pipeline sit well under 1% disc and ~0 swaps.")
    print("   bHet/oHom = base het, other hom (other lost a het);  bHom/oHet = the reverse.")
    print(f"   {'other copy':26s} {'base copy':18s} {'both':>8s} {'disc%':>6s} {'swap':>6s} "
          f"{'bHet/oHom':>9s} {'bHom/oHet':>9s} {'missO%':>7s} {'missB%':>7s}")
    bytype = defaultdict(list)
    for (b, o, _, _), t in sorted(zip(P, T), key=lambda z: -sum(v for (gb, go), v in z[1].items()
                                                                if gb is not None and go is not None and gb != go)):
        both = sum(v for (gb, go), v in t.items() if gb is not None and go is not None)
        disc = sum(v for (gb, go), v in t.items() if gb is not None and go is not None and gb != go)
        swap = t[(0, 2)] + t[(2, 0)]
        # base het, other hom  /  base hom, other het
        h2h = t[(1, 0)] + t[(1, 2)]
        hom2het = t[(0, 1)] + t[(2, 1)]
        mo = sum(v for (gb, go), v in t.items() if go is None)
        mb = sum(v for (gb, go), v in t.items() if gb is None)
        dp = 100 * disc / both if both else float("nan")
        bytype[f"{co[o]}~{co[b]}"].append((dp, 100 * swap / both if both else 0, 100 * mo / n, 100 * mb / n))
        print(f"   {o[:26]:26s} {b[:18]:18s} {both:8d} {dp:6.2f} {swap:6d} {h2h:9d} {hom2het:9d} "
              f"{100 * mo / n:7.2f} {100 * mb / n:7.2f}")
    print("\n   by pair type: mean disc% | mean swap% | mean missing% other copy | base copy")
    for k, v in sorted(bytype.items()):
        print(f"   {k:22s} n={len(v):3d}  disc {statistics.mean(x[0] for x in v):6.2f}%  "
              f"swap {statistics.mean(x[1] for x in v):6.3f}%  miss {statistics.mean(x[2] for x in v):6.2f}% "
              f"| {statistics.mean(x[3] for x in v):6.2f}%")

    # WHAT DOES "ABSENT" MEAN? Where one copy has no call (its source file had
    # no record there), the other copy of the SAME person says what the truth
    # was. Mostly 0/0 -> the source lists only variant sites, so absent means
    # hom-ref and filling it with 0/0 is right. Plenty of 0/1 and 1/1 -> the
    # source dropped sites the person does carry, so absent is truly missing.
    def absent_view(selector, label):
        tot = Counter()
        for (b, o, _, _), t in zip(P, T):
            for (gb, go), v in t.items():
                key = selector(gb, go)
                if key is not None: tot[key] += v
        n_all = sum(tot.values())
        if not n_all:
            print(f"   {label}: none"); return
        parts = ", ".join(f"{g} {100 * tot[c] / n_all:5.1f}%" for c, g in ((0, "0/0"), (1, "0/1"), (2, "1/1")))
        print(f"   {label}: {n_all:,} genotypes -> {parts}")
    print("\n   WHAT 'ABSENT' MEANS (pooled over pairs):")
    absent_view(lambda gb, go: gb if (go is None and gb is not None) else None,
                "other copy missing -> base copy's genotype")
    absent_view(lambda gb, go: go if (gb is None and go is not None) else None,
                "base copy missing  -> other copy's genotype")
    print("   mostly 0/0 = absent means hom-ref (fill with 0/0); "
          "real 0/1 + 1/1 share = sites were filtered (keep missing)")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    def lab(p):
        p.add_argument("--cohort", action="append", help="NAME=FILE of ids (repeatable)")
        p.add_argument("--groups", help="sample_groups.tsv (sample<TAB>superpop)")
        p.add_argument("--sgdp-manifest", help="reference_ids/sgdp_manifest.tsv")

    def dup(p):
        p.add_argument("--dups", help="04_homogeneity_panel/king/duplicate_pairs.tsv")
        p.add_argument("--known-dups", help="reference_ids/known_duplicate_pairs.tsv")
        p.add_argument("--kin0", help="04_homogeneity_panel/king/king.kin0 (kinship for metadata pairs)")

    q = sub.add_parser("q"); lab(q); dup(q)
    q.add_argument("--fam", required=True); q.add_argument("--q", required=True)
    q.add_argument("--curated", action="append", help="SP=FILE, e.g. EUR=reference_ids/eur_rfmix.txt")
    q.add_argument("--exclude", help="ids never eligible (panel_exclude_ids.txt)")
    q.add_argument("--thresholds", default="0.93,0.95")
    q.add_argument("--cap", type=int, default=50); q.add_argument("--floor", type=int, default=10)
    q.add_argument("--show-pairs", type=int, default=12)
    q.add_argument("--old-bim"); q.add_argument("--new-bim")
    ps = sub.add_parser("psc"); lab(ps); ps.add_argument("--stats", required=True)
    he = sub.add_parser("het"); lab(he); he.add_argument("--het", required=True)
    pi = sub.add_parser("pair-ids"); dup(pi); pi.add_argument("--fam", required=True)
    pg = sub.add_parser("pairs-gt"); lab(pg); dup(pg)
    a = ap.parse_args()
    {"q": mode_q, "psc": mode_psc, "het": mode_het,
     "pair-ids": mode_pair_ids, "pairs-gt": mode_pairs_gt}[a.mode](a)


if __name__ == "__main__":
    main()
