#!/usr/bin/env python3
"""
variant-accounting.py -- variants kept and lost at each step, per chromosome
and genome-wide, for the methods. Run by scripts/merge-sites-collect.sh
after it has gathered 1c's tables into --out-dir. READ-ONLY on the pipeline.

Every number has a source column: a file's record count (bcftools index -s,
from the wrapper), a table 1c wrote (provenance/, gathered into --out-dir),
or a line a stage printed to its own log (the newest log of that stage and
chromosome newer than the rebuild epoch; the line itself is kept in
variant_accounting_loglines.tsv).

Steps, in run order:
  HGDP+1KG  release records (SNVs + indels) -> SNP records 1c merged
  MXB       records after liftover to GRCh38 (rejects: mxb_lifted.rejected)
  SGDP      raw union -> sites masked -> liftover rejects -> prepped
  Almarri   merged -> sites masked -> prepped
  panel     union of the four -> after plink2 site QC (biallelic ACGT SNPs,
            duplicate IDs out) -> positions with >1 record dropped ->
            fill-contradicted positions dropped -> soft-union (MAF >= 0.005
            in >= 1 super-pop group) -> ADMIXTURE site list; phased; RFMix
  4a        genome-wide: ADMIXTURE candidates -> LD-pruned SNPs
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import OrderedDict, defaultdict

COHORTS = ["HGDP1KG", "MXB", "SGDP", "ALMARRI"]
STEPS = OrderedDict([
    # key: (order, label)
    ("hgdp_release", (10, "HGDP+1KG release records (SNVs + indels, all samples)")),
    ("mxb_lifted", (20, "MXB records after liftover to GRCh38")),
    ("mxb_rejects", (21, "MXB liftover rejects (hg19 contig)")),
    ("sgdp_raw", (30, "SGDP raw records (union of per-sample calls)")),
    ("sgdp_low", (31, "SGDP records with a call below the filter level")),
    ("sgdp_masked", (32, "SGDP sites set missing (more failing carriers than passing)")),
    ("sgdp_rejects", (33, "SGDP liftover rejects")),
    ("sgdp_prepped", (34, "SGDP records after prep (GRCh38)")),
    ("almarri_merged", (40, "Almarri records after merging per-sample calls")),
    ("almarri_masked", (41, "Almarri sites set missing (more missing than ALT)")),
    ("almarri_prepped", (42, "Almarri records after prep")),
    ("cohort_file", (50, "records each cohort's file brings to the merge")),
    ("union", (60, "panel: union of the four cohort files")),
    ("after_plink_qc", (61, "panel: after plink2 site QC (biallelic ACGT SNPs, duplicate IDs out)")),
    ("multi_positions", (62, "dropped: positions holding more than one record")),
    ("multi_records", (63, "dropped: records at those positions")),
    ("fill_positions", (64, "dropped: positions where a 0/0 fill was contradicted")),
    ("after_drops", (65, "panel: after those drops")),
    ("softunion", (66, "panel: soft-union (MAF >= 0.005 in at least one super-pop group)")),
    ("admix_sites", (67, "panel: ADMIXTURE site list (missingness rule)")),
    ("phased", (70, "panel: phased (SHAPEIT5, ligated)")),
    ("rechr", (71, "panel: RFMix input (.rechr)")),
    ("admix_all", (80, "4a: ADMIXTURE candidate sites, all chromosomes")),
    ("pruned", (81, "4a: LD-pruned SNPs fitted (MAF >= 0.01, r2 < 0.1 in 500-SNP windows)")),
])

NUM = r"([\d,]+)"
RX = {
    "sgdp_raw": re.compile(r"\] " + NUM + r" records; " + NUM + r" have at least one call below FL \S+; " + NUM + " have more failing"),
    "sgdp_done": re.compile(r"DONE\. \d+ samples, " + NUM + r" records on chr\w+, " + NUM + r" liftover rejects"),
    "almarri_merged": re.compile(r"merged: \d+ samples, " + NUM + r" records, missing GT fraction"),
    "almarri_mask": re.compile(r"site mask: all of Almarri set missing at " + NUM + r" of " + NUM + " sites"),
    "almarri_nomask": re.compile(r"ALMARRI_SITE_MASK=0: no site mask \(" + NUM + " of " + NUM),
    "prov_panel": re.compile(r"\[provenance \S+ (merged|qc|softunion)\] panel: " + NUM + r" records"),
    "dropped": re.compile(r"dropped " + NUM + r" positions: " + NUM + r" holding more than one record \(" + NUM
                          + r" records: split multi-allelic sites\); " + NUM + " where a 0/0 fill was contradicted"),
    "softunion": re.compile(r"\[softunion[^\]]*\] " + NUM + r" records in; " + NUM + r" kept"),
    "admix": re.compile(r"\[admix-sites\] " + NUM + r" records at " + NUM + r" positions; kept " + NUM + r" positions"),
    "plink_remaining": re.compile(NUM + r" variants? remaining after main filters"),
}


def n(s):
    return int(s.replace(",", ""))


def newest_log(logdir, stage, chrom, epoch):
    """The newest <stage>_<chrom>_<job>.out of this rebuild, and its .err."""
    hits = [p for p in glob.glob(os.path.join(logdir, f"{stage}_{chrom}_*.out"))
            if epoch is None or os.path.getmtime(p) > epoch]
    if not hits:
        return None, None
    out = max(hits, key=os.path.getmtime)
    err = out[:-4] + ".err"
    return out, (err if os.path.exists(err) else None)


def read_lines(*paths):
    for p in paths:
        if p:
            with open(p, errors="replace") as fh:
                for line in fh:
                    yield p, line.rstrip("\n")


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    if not rows:
        return []
    h = rows[0]
    return [dict(zip(h, r)) for r in rows[1:]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, help="results/merge_sites, already holding the collector's tables")
    ap.add_argument("--file-counts", required=True, help="chrom step value source, from the wrapper")
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--epoch", default="")
    ap.add_argument("--chroms", default=" ".join(str(i) for i in range(1, 23)))
    a = ap.parse_args()

    epoch = os.path.getmtime(a.epoch) if a.epoch and os.path.exists(a.epoch) else None
    chroms = [f"chr{c.replace('chr', '')}" for c in a.chroms.split()]
    rows = []          # chrom, step, cohort, value, source
    loglines = []      # chrom, stage, log, line

    def add(chrom, step, value, source, cohort="-"):
        if value is None or value == "" or value == "NA":
            return
        rows.append((chrom, step, cohort, int(value), source))

    for r in read_tsv(a.file_counts):
        add(r["chrom"], r["step"], r["value"], r["source"])

    # 1c's own tables, as the collector gathered them
    pat = read_tsv(os.path.join(a.out_dir, "patterns.tsv"))
    stage_sum = defaultdict(int)
    for r in pat:
        stage_sum[(r["chrom"], r["stage"])] += int(r["n_records"])
    for (ch, st), v in stage_sum.items():
        step = {"merged": "union", "qc": "after_drops", "softunion": "softunion"}.get(st)
        if step:
            add(ch, step, v, f"patterns.tsv ({st} stage, summed over cohort combinations)")
    for r in read_tsv(os.path.join(a.out_dir, "cohorts.tsv")):
        if r["stage"] == "merged":
            add(r["chrom"], "cohort_file", r["file_alleles"], "cohorts.tsv (merged stage, file_alleles)", r["cohort"])
    for r in read_tsv(os.path.join(a.out_dir, "fill_check_totals.tsv")):
        add(r["chrom"], "after_plink_qc", r["records"], "fill_check_totals.tsv (records the fill check judged)")
        add(r["chrom"], "fill_positions", r["positions_dropped"], "fill_check_totals.tsv (positions_dropped)")
    for r in read_tsv(os.path.join(a.out_dir, "chrom_totals.tsv")):
        add(r["chrom"], "admix_sites", r["admix_positions"], "chrom_totals.tsv (admix_sites.bed lines)")

    # the stages' own logs
    have = {(c, s) for c, s, *_ in rows}
    for ch in chroms:
        for stage in ("1b_sgdp_raw", "1d_almarri", "1c_phase"):
            out, err = newest_log(a.logdir, stage, ch, epoch)
            if not out:
                loglines.append((ch, stage, "-", f"no {stage}_{ch}_*.out newer than the epoch in {a.logdir}"))
                continue
            for p, line in read_lines(out, err):
                hit = None
                for key, rx in RX.items():
                    m = rx.search(line)
                    if m:
                        hit = (key, m)
                        break
                if not hit:
                    continue
                key, m = hit
                loglines.append((ch, stage, os.path.basename(p), line.strip()))
                src = f"{os.path.basename(p)}: {line.strip()[:160]}"
                if key == "sgdp_raw":
                    add(ch, "sgdp_raw", n(m.group(1)), src)
                    add(ch, "sgdp_low", n(m.group(2)), src)
                    add(ch, "sgdp_masked", n(m.group(3)), src)
                elif key == "sgdp_done":
                    if (ch, "sgdp_prepped") not in have:
                        add(ch, "sgdp_prepped", n(m.group(1)), src)
                    add(ch, "sgdp_rejects", n(m.group(2)), src)
                elif key == "almarri_merged":
                    add(ch, "almarri_merged", n(m.group(1)), src)
                elif key in ("almarri_mask", "almarri_nomask"):
                    add(ch, "almarri_masked", n(m.group(1)) if key == "almarri_mask" else 0, src)
                elif key == "prov_panel" and m.group(1) == "merged" and (ch, "union") not in have:
                    add(ch, "union", n(m.group(2)), src)
                elif key == "dropped":
                    add(ch, "multi_positions", n(m.group(2)), src)
                    add(ch, "multi_records", n(m.group(3)), src)
                elif key == "softunion" and (ch, "softunion") not in have:
                    add(ch, "softunion", n(m.group(2)), src)
                elif key == "admix" and (ch, "admix_sites") not in have:
                    add(ch, "admix_sites", n(m.group(3)), src)

    # one value per (chrom, step, cohort): the first source wins (files, then tables, then logs)
    seen, final = set(), []
    for r in rows:
        k = r[:3]
        if k in seen:
            continue
        seen.add(k)
        final.append(r)
    unknown = sorted({r[1] for r in final if r[1] not in STEPS})
    if unknown:
        sys.exit(f"ERROR: steps with no label: {unknown}")
    final.sort(key=lambda r: (STEPS[r[1]][0], r[2], (int(r[0][3:]) if r[0][3:].isdigit() else 99)))

    out_path = os.path.join(a.out_dir, "variant_accounting.tsv")
    with open(out_path, "w") as g:
        g.write("chrom\torder\tstep\tlabel\tcohort\tvalue\tsource\n")
        for ch, step, co, v, src in final:
            g.write(f"{ch}\t{STEPS[step][0]}\t{step}\t{STEPS[step][1]}\t{co}\t{v}\t{src}\n")
    with open(os.path.join(a.out_dir, "variant_accounting_loglines.tsv"), "w") as g:
        g.write("chrom\tstage\tlog\tline\n")
        for r in loglines:
            g.write("\t".join(r) + "\n")

    # genome-wide: sum over chromosomes; a step missing on some shows how many it has
    tot = OrderedDict()
    for ch, step, co, v, src in final:
        k = (step, co)
        t = tot.setdefault(k, [0, set()])
        t[0] += v
        t[1].add(ch)
    gpath = os.path.join(a.out_dir, "variant_accounting_genome.tsv")
    with open(gpath, "w") as g:
        g.write("order\tstep\tlabel\tcohort\ttotal\tchromosomes\n")
        print(f"\n[accounting] genome-wide (sum over chromosomes; 'genome' rows are one figure)")
        print(f"[accounting] {'step':58s} {'cohort':8s} {'total':>14s}  chromosomes")
        for (step, co), (v, chs) in sorted(tot.items(), key=lambda kv: (STEPS[kv[0][0]][0], kv[0][1])):
            nch = "genome" if chs == {"genome"} else str(len(chs))
            g.write(f"{STEPS[step][0]}\t{step}\t{STEPS[step][1]}\t{co}\t{v}\t{nch}\n")
            print(f"[accounting] {STEPS[step][1][:58]:58s} {co:8s} {v:>14,d}  {nch}")
    missing = [ch for ch in chroms if not any(r[0] == ch and r[1] == "softunion" for r in final)]
    if missing:
        print(f"[accounting] WARNING: no soft-union count for {' '.join(missing)}")
    print(f"[accounting] per chromosome: {out_path}; log lines used: variant_accounting_loglines.tsv")


if __name__ == "__main__":
    main()
