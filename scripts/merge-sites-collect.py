#!/usr/bin/env python3
"""
merge-sites-collect.py -- the merge's site tables, all chromosomes, one place.

1c writes, per chromosome, what happened to every cohort's sites in the merge
(scripts/site-provenance.py, fill-contradictions.py, softunion-sites.py,
admix-site-filter.py), and scripts/cohort-site-report.sh adds per-site
missingness and allele frequency. Those are ~10 small files per chromosome
under ${PANEL_DIR}. This gathers them into a few tidy TSVs (a `chrom` column
first) that can be committed and plotted off the cluster. READ-ONLY on the
panel; writes only --out-dir.

  cohorts.tsv            per stage (merged / qc / softunion) and cohort: file
                         records, records reaching the panel, records unique
                         to the cohort, and for the panel records the cohort
                         lacks, why (absent, swapped, flipped, other allele)
  patterns.tsv           per stage: which COMBINATION of cohort files carries
                         each panel record ("HGDP1KG+SGDP", ...), with counts
  fill_check.tsv         the 0/0-fill guard, long: cohort x super-pop, how
                         many of the cohort's absent records another cohort of
                         the same super-pop contradicted (position dropped)
  fill_check_totals.tsv  per chromosome: records judged, positions dropped
  softunion_groups.tsv   per super-pop group: records with MAF >= the cut in it
  admix_groups.tsv / admix_totals.tsv   the ADMIXTURE site list's missingness rule
  chrom_totals.tsv       per chromosome: samples, soft-union records, ADMIXTURE
                         positions, fill-dropped positions, the 1c stamp
and, when chrN.site_report.tsv.gz exists (scripts/cohort-site-report.sh):
  site_report_src.tsv    per cohort: soft-union records by where they stood in
                         the cohort's own file (P present, A absent -> filled
                         0/0, S/F/FS/X/R/U allele mismatches), all and ADMIXTURE
  site_report_miss.tsv   per cohort: records by fraction of genotypes missing
  site_report_af2d.tsv   per cohort pair: 10 x 10 histogram of the two ALT
                         frequencies at records both files carry (P and P)
  site_report_bins.tsv   per 1-Mb bin and cohort, long: absent, other alleles,
                         mean missingness
  MANIFEST.tsv           every output: rows, chromosomes, and whether each
                         source file is newer than the rebuild epoch
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import OrderedDict, defaultdict
from datetime import datetime

STAGES = ["merged", "qc", "softunion"]
MISS_BINS = [(0.0, 0.0, "0"), (0.0, 0.01, "(0,0.01]"), (0.01, 0.05, "(0.01,0.05]"), (0.05, 0.10, "(0.05,0.10]"),
             (0.10, 0.15, "(0.10,0.15]"), (0.15, 0.30, "(0.15,0.30]"), (0.30, 0.50, "(0.30,0.50]"),
             (0.50, 1.00, "(0.50,1]")]


def miss_bin(m):
    if m <= 0.0:
        return "0"
    for lo, hi, lab in MISS_BINS[1:]:
        if lo < m <= hi:
            return lab
    return MISS_BINS[-1][2]


def af_bin(x):
    return min(int(x * 10), 9)


class Out:
    """One output TSV: header written once, rows counted, chromosomes noted."""

    def __init__(self, d, name, header):
        self.path = os.path.join(d, name)
        self.name, self.header = name, header
        self.fh = open(self.path + ".partial", "w", newline="")
        self.w = csv.writer(self.fh, delimiter="\t", lineterminator="\n")
        self.w.writerow(header)
        self.rows, self.chroms = 0, set()

    def row(self, chrom, values):
        self.w.writerow([chrom] + list(values))
        self.rows += 1
        self.chroms.add(chrom)

    def close(self):
        self.fh.close()
        os.replace(self.path + ".partial", self.path)


def read_tsv(path):
    """(comment lines, header, rows) of a '#'-commented TSV."""
    comments, header, rows = [], None, []
    with open(path, newline="") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line.lstrip("# "))
            elif header is None:
                header = line.split("\t")
            else:
                rows.append(line.split("\t"))
    return comments, header, rows


def kv(comment_lines):
    d = {}
    for c in comment_lines:
        for m in re.finditer(r"(\w+)=(\S+)", c):
            d[m.group(1)] = m.group(2)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--totals", required=True, help="chrom samples softunion_records admix_positions softunion_new recipe")
    ap.add_argument("--epoch", default="", help="rebuild_epoch file: sources older than it are flagged")
    ap.add_argument("--chroms", default=" ".join(str(i) for i in range(1, 23)))
    ap.add_argument("--no-site-report", action="store_true", help="skip the per-site tables (slow: ~27M rows)")
    a = ap.parse_args()

    prov = os.path.join(a.panel_dir, "provenance")
    os.makedirs(a.out_dir, exist_ok=True)
    epoch = os.path.getmtime(a.epoch) if a.epoch and os.path.exists(a.epoch) else None
    sources = defaultdict(lambda: [0, 0, 0])   # output -> [files read, newer than epoch, older]

    def used(out_name, path):
        s = sources[out_name]
        s[0] += 1
        if epoch is not None:
            s[1 if os.path.getmtime(path) > epoch else 2] += 1

    chroms = [f"chr{c.replace('chr', '')}" for c in a.chroms.split()]
    o_co = o_pat = None
    o_fc = Out(a.out_dir, "fill_check.tsv", ["chrom", "cohort", "compared_on", "records_absent_from_file",
                                              "contradicted", "superpop", "contradicted_superpop"])
    o_fct = Out(a.out_dir, "fill_check_totals.tsv", ["chrom", "min_n", "min_ac", "max_p", "records", "positions_dropped"])
    o_sug = Out(a.out_dir, "softunion_groups.tsv", ["chrom", "group", "samples", "judged", "passing_records"])
    o_adg = o_adt = None
    missing = defaultdict(list)
    fill_dropped = {}

    for ch in chroms:
        n = ch[3:]
        for st in STAGES:
            p = os.path.join(prov, f"{ch}.{st}.cohorts.tsv")
            if os.path.exists(p):
                _, h, rows = read_tsv(p)
                if o_co is None:
                    o_co = Out(a.out_dir, "cohorts.tsv", ["chrom"] + h)
                elif h != o_co.header[1:]:
                    sys.exit(f"ERROR: {p}: header differs from the first cohorts table:\n  {h}\n  {o_co.header[1:]}")
                for r in rows:
                    o_co.row(ch, r)
                used("cohorts.tsv", p)
            else:
                missing["cohorts.tsv"].append(f"{ch}.{st}")
            p = os.path.join(prov, f"{ch}.{st}.patterns.tsv")
            if os.path.exists(p):
                _, h, rows = read_tsv(p)
                if o_pat is None:
                    o_pat = Out(a.out_dir, "patterns.tsv", ["chrom"] + h)
                for r in rows:
                    o_pat.row(ch, r)
                used("patterns.tsv", p)
            else:
                missing["patterns.tsv"].append(f"{ch}.{st}")

        p = os.path.join(prov, f"{ch}.fill_check.tsv")
        if os.path.exists(p):
            com, h, rows = read_tsv(p)
            t = kv(com)
            o_fct.row(ch, [t.get(k, "") for k in ("min_n", "min_ac", "max_p", "records", "positions_dropped")])
            fill_dropped[ch] = t.get("positions_dropped", "")
            used("fill_check_totals.tsv", p)
            ix = {c: i for i, c in enumerate(h)}
            for r in rows:
                base = [r[ix["cohort"]], r[ix["compared_on"]], r[ix["records_absent_from_file"]], r[ix["contradicted"]]]
                bysp = r[ix["contradicted_by_superpop"]] if len(r) > ix["contradicted_by_superpop"] else ""
                parts = [x for x in bysp.split(",") if ":" in x]
                if not parts:
                    o_fc.row(ch, base + ["", "0"])
                for x in parts:
                    sp, c = x.split(":", 1)
                    o_fc.row(ch, base + [sp, c])
            used("fill_check.tsv", p)
        else:
            missing["fill_check.tsv"].append(ch)

        p = os.path.join(prov, f"{ch}.softunion_groups.tsv")
        if os.path.exists(p):
            _, h, rows = read_tsv(p)
            for r in rows:
                o_sug.row(ch, (r + [""] * 4)[:4])
            used("softunion_groups.tsv", p)
        else:
            missing["softunion_groups.tsv"].append(ch)

        p = os.path.join(a.panel_dir, f"merged_{ch}.admix_sites.summary.tsv")
        if os.path.exists(p):
            com, h, rows = read_tsv(p)
            t = kv(com)
            tk = ["max_group_miss", "max_total_miss", "min_group_n", "records", "records_kept", "positions",
                  "positions_kept", "records_failing_overall"]
            if o_adt is None:
                o_adt = Out(a.out_dir, "admix_totals.tsv", ["chrom"] + tk)
                o_adg = Out(a.out_dir, "admix_groups.tsv", ["chrom"] + h)
            o_adt.row(ch, [t.get(k, "") for k in tk])
            for r in rows:
                o_adg.row(ch, r)
            used("admix_groups.tsv", p)
            used("admix_totals.tsv", p)
        else:
            missing["admix_groups.tsv"].append(ch)
            missing["admix_totals.tsv"].append(ch)

    # per-chromosome totals (from the published files, via the wrapper's bcftools)
    o_ct = Out(a.out_dir, "chrom_totals.tsv", ["chrom", "samples", "softunion_records", "admix_positions",
                                               "fill_positions_dropped", "softunion_newer_than_epoch", "recipe"])
    _, h, rows = read_tsv(a.totals)
    ix = {c: i for i, c in enumerate(h)}
    for r in rows:
        ch = r[ix["chrom"]]
        o_ct.row(ch, [r[ix["samples"]], r[ix["softunion_records"]], r[ix["admix_positions"]],
                      fill_dropped.get(ch, "NA"), r[ix["softunion_new"]], r[ix["recipe"]] if len(r) > ix["recipe"] else ""])

    # per-site report (scripts/cohort-site-report.sh), when it has run
    outs = [o for o in (o_co, o_pat, o_fc, o_fct, o_sug, o_adg, o_adt, o_ct) if o is not None]
    if not a.no_site_report:
        o_src = o_miss = o_af = o_bins = None
        for ch in chroms:
            p = os.path.join(prov, f"{ch}.site_report.tsv.gz")
            if not os.path.exists(p):
                missing["site_report_*.tsv"].append(ch)
                continue
            src = defaultdict(int); src_ad = defaultdict(int)
            miss = defaultdict(int); miss_ad = defaultdict(int)
            af2 = defaultdict(int)
            with gzip.open(p, "rt") as fh:
                head = fh.readline().rstrip("\n").split("\t")
                labels = [x[4:] for x in head if x.startswith("src_")]
                ci = {c: (head.index(f"src_{c}"), head.index(f"miss_{c}"), head.index(f"af_{c}")) for c in labels}
                iad = head.index("admix")
                pairs = [(x, y) for i, x in enumerate(labels) for y in labels[i + 1:]]
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    ad = f[iad] == "1"
                    for c in labels:
                        s_i, m_i, _ = ci[c]
                        s, mb = f[s_i], miss_bin(float(f[m_i]))
                        src[(c, s)] += 1
                        miss[(c, mb)] += 1
                        if ad:
                            src_ad[(c, s)] += 1
                            miss_ad[(c, mb)] += 1
                    for x, y in pairs:
                        if f[ci[x][0]] == "P" and f[ci[y][0]] == "P" and f[ci[x][2]] != "NA" and f[ci[y][2]] != "NA":
                            af2[(x, y, af_bin(float(f[ci[x][2]])), af_bin(float(f[ci[y][2]])))] += 1
            if o_src is None:
                o_src = Out(a.out_dir, "site_report_src.tsv", ["chrom", "cohort", "src", "records", "admix_records"])
                o_miss = Out(a.out_dir, "site_report_miss.tsv", ["chrom", "cohort", "miss_bin", "records", "admix_records"])
                o_af = Out(a.out_dir, "site_report_af2d.tsv", ["chrom", "cohort_x", "cohort_y", "af_x_bin", "af_y_bin", "records"])
            for (c, s), v in sorted(src.items()):
                o_src.row(ch, [c, s, v, src_ad.get((c, s), 0)])
            order = {lab: i for i, (_, _, lab) in enumerate(MISS_BINS)}
            for (c, mb), v in sorted(miss.items(), key=lambda kv_: (kv_[0][0], order[kv_[0][1]])):
                o_miss.row(ch, [c, mb, v, miss_ad.get((c, mb), 0)])
            for (x, y, bx, by), v in sorted(af2.items()):
                o_af.row(ch, [x, y, f"{bx / 10:.1f}-{(bx + 1) / 10:.1f}", f"{by / 10:.1f}-{(by + 1) / 10:.1f}", v])
            for o in (o_src, o_miss, o_af):
                used(o.name, p)
            pb = os.path.join(prov, f"{ch}.site_report.bins.tsv")
            if os.path.exists(pb):
                _, hb, rows = read_tsv(pb)
                labs = [x[len("absent_"):] for x in hb if x.startswith("absent_")]
                if o_bins is None:
                    o_bins = Out(a.out_dir, "site_report_bins.tsv", ["chrom", "bin_start", "bin_end", "records", "admix_records",
                                                                     "cohort", "absent", "other_alleles", "mean_miss"])
                ib = {c: i for i, c in enumerate(hb)}
                for r in rows:
                    for c in labs:
                        o_bins.row(ch, [r[ib["bin_start"]], r[ib["bin_end"]], r[ib["records"]], r[ib["admix_records"]], c,
                                        r[ib[f"absent_{c}"]], r[ib[f"other_alleles_{c}"]], r[ib[f"mean_miss_{c}"]]])
                used("site_report_bins.tsv", pb)
            print(f"[merge-sites] {ch}: site report summarised ({len(labels)} cohorts)")
        outs += [o for o in (o_src, o_miss, o_af, o_bins) if o is not None]

    for o in outs:
        o.close()
    with open(os.path.join(a.out_dir, "MANIFEST.tsv"), "w") as g:
        g.write(f"# merge-sites-collect {datetime.now():%Y-%m-%d %H:%M}  panel_dir={a.panel_dir}  "
                f"epoch={a.epoch or 'none'}\n")
        g.write("file\trows\tchromosomes\tsource_files\tsources_newer_than_epoch\tsources_older\tnot_found\n")
        for o in outs:
            s = sources.get(o.name, [0, 0, 0])
            if o.name.startswith("site_report"):
                nf = missing.get("site_report_*.tsv", [])
            else:
                nf = missing.get(o.name, [])
            # "chr1.merged,chr1.qc,chr1.softunion" -> "chr1"; one stage -> "chr1(qc)"
            by_ch = OrderedDict()
            for x in nf:
                c, _, st = x.partition(".")
                by_ch.setdefault(c, []).append(st)
            nf = [c if (not sts[0] or len(sts) == len(STAGES)) else f"{c}({'+'.join(sts)})" for c, sts in by_ch.items()]
            if o.name == "fill_check_totals.tsv":
                nf = missing.get("fill_check.tsv", [])
            src_n = s[0] if o.name != "chrom_totals.tsv" else "bcftools"
            g.write(f"{o.name}\t{o.rows}\t{len(o.chroms)}\t{src_n}\t{s[1] if epoch else 'NA'}\t"
                    f"{s[2] if epoch else 'NA'}\t{','.join(nf) if nf else '-'}\n")
        if not a.no_site_report and missing.get("site_report_*.tsv"):
            g.write(f"site_report_*.tsv\t0\t0\t0\tNA\tNA\t{','.join(missing['site_report_*.tsv'])}\n")
    with open(os.path.join(a.out_dir, "MANIFEST.tsv")) as fh:
        print(fh.read(), end="")
    stale = sum(s[2] for s in sources.values())
    if epoch is not None and stale:
        print(f"[merge-sites] WARNING: {stale} source file(s) are older than the rebuild epoch -- from the previous build")


if __name__ == "__main__":
    main()
