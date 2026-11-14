#!/usr/bin/env python3
"""
build-chip-manifest-sites.py

Turn a vendor array manifest pair (GRCh37 + GRCh38, SAME array version) into the
two site lists the benchmark contract asks for, plus a UCSC custom track.

WHY A PAIR AND NOT A LIFTOVER. Document 20, step 9: the hg19 cells of the
benchmark exist to test the pipeline's OWN liftover path. If the hg19 site list
were produced by down-lifting the hg38 list, it would contain exactly the sites
that map back cleanly -- by construction -- and would hide the probes that real
hg19 array data carry and that the pipeline's CrossMap step drops, moves or
strand-flips. So both builds come from the vendor's own manifests for that array
version, joined on the PROBE ID, and a probe is kept only when both builds agree
on the chromosome.

WHAT THIS EMITS (Document 20, step 9 and appendix section 2)

  manifests/<chip>.hg38.sites.tsv   chrom pos ref alt
  manifests/<chip>.hg19.sites.tsv   chrom pos ref alt grch38_site grch38_strand

  <chip>.hg38.annot.tsv             companion: probe, rsid, palindromic, orient
  <chip>.hg19.annot.tsv             (NOT part of the contract; the intake
                                     check reads only the manifests/ files, so
                                     extra columns must not go in them)

  hg38/chr<N>.snps  hg19/chr<N>.snps   chr:pos drop-ins for CHIP_SNPS_TPL
  <chip>.<build>.ucsc.bed              UCSC custom track (0-based, 4 columns)
  build_report.tsv                     every drop reason, counted
  PROVENANCE.md                        source file name, bytes and sha256

THE COORDINATE RULE, STATED ONCE. Illumina `MapInfo` is 1-BASED. UCSC BED is
0-BASED HALF-OPEN. So a BED line for a SNP is
    chromStart = MapInfo - 1      chromEnd = MapInfo
and the 1-based position is chromEnd. This script reads MapInfo as-is from an
annotated manifest and reads chromEnd from a vendor BED; both give the same
1-based number. It never applies the -1 twice. (Verified on GSA-24v1-0_C2:
MapInfo 792461 <-> BED 792460 792461.)

WHY REF/ALT ARE READ FROM THE GENOME, NOT FROM THE VENDOR. Illumina's `Alleles`
column is on the design strand, not the genome plus strand. On GSA v1.0 the
column holds 246,322 [T/C] and 245,940 [A/G] -- the same class, reported both
ways -- so roughly half the probes are on the minus strand of the genome and
the column alone cannot say which. This script therefore reads the reference
base out of each build's FASTA and orients the pair against it:
    direct      one allele equals the reference base
    revcomp     matched only after reverse-complementing
    unresolved  neither matched -- dropped, and counted
That measures the strand convention instead of assuming one.

PALINDROMES. A/T and C/G probes cannot be oriented from alleles at all. REF/ALT
are still well defined for them in a single build (ref A forces alt T either
way), so they are KEPT and the manifests are correct. What is undecidable is
which physical allele a GENOTYPE CALL refers to, which is a merge-time problem,
not a site-list problem. They are flagged in the .annot.tsv companion and
counted in the report. Per Document 20, grch38_strand is taken as + for them.

Usage (GSA v1.0, the only GSA version with both native builds in hand):

  python3 -I scripts/build-chip-manifest-sites.py \\
      --chip gsa \\
      --hg38-annot <dir>/GSA-24v1-0_C2.hg38.annotated.txt \\
      --hg19-bed   <dir>/GSA-24v1-0_C1.bed \\
      --rsids      <dir>/GSA-24v1-0_C1_b150_rsids.txt \\
      --hg38-fasta /storage/.../Homo_sapiens_assembly38.fasta \\
      --hg19-fasta /storage/.../hg19/hg19.fa \\
      --out reference_ids/chip_gsa

Without --hg38-fasta/--hg19-fasta the join, the BED tracks, the .snps lists and
the report are still written; the manifests/ files are not, because they cannot
be written honestly without the reference sequence. The script says so and exits
non-zero so a pipeline cannot mistake a partial run for a complete one.
"""
import argparse
import hashlib
import os
import subprocess
import sys
from collections import Counter

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AUTOSOMES = [str(c) for c in range(1, 23)]
PALINDROMIC_PAIRS = {frozenset("AT"), frozenset("CG")}


def rc(b):
    return COMP.get(b, "N")


def norm_chrom(c):
    """'chr7' / '7' / 'Chr7' -> '7'. Returns '' for anything unusable."""
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def parse_alleles(field):
    """'[A/G]' or 'A/G' -> ('A','G'); returns None when not a biallelic SNV."""
    f = field.strip().strip("[]")
    parts = f.split("/")
    if len(parts) != 2:
        return None
    a1, a2 = parts[0].strip().upper(), parts[1].strip().upper()
    if a1 not in COMP or a2 not in COMP or a1 == a2:
        return None          # I/D, D/I, indels, N, or a degenerate pair
    return a1, a2


def read_annotated(path):
    """Illumina annotated manifest -> {probe: (chrom, pos, (a1,a2) or None)}."""
    out = {}
    with open(path, "r", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = header.index("Name")
            i_chr = header.index("Chr")
            i_pos = header.index("MapInfo")
            i_all = header.index("Alleles")
        except ValueError:
            sys.exit(f"ERROR: {path} is missing Name/Chr/MapInfo/Alleles columns; "
                     f"found: {header[:8]}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_all:
                continue
            try:
                pos = int(f[i_pos])
            except ValueError:
                continue
            if pos <= 0:
                continue                      # unmapped in this build
            out[f[i_name]] = (norm_chrom(f[i_chr]), pos, parse_alleles(f[i_all]))
    return out


def read_bed(path):
    """Vendor BED -> {probe: (chrom, pos_1based)}. pos = chromEnd."""
    out = {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                end = int(f[2])
            except ValueError:
                continue
            if end <= 0:
                continue
            out[f[3]] = (norm_chrom(f[0]), end)
    return out


def read_rsids(path):
    """Illumina b150 rsid file -> {probe_id: rsid}.

    The 'Name' column is the PROBE ID as it appears in the manifest, not a
    position: on GSA-24v1-0_C1 it is 295,713 rs-numbers, 302,310 GSA-* ids and
    610 chr:pos ids, and all 618,540 match the manifest's Name column exactly.
    Keying this on chr:pos instead silently matches only the 610.
    """
    out = {}
    if not path:
        return out
    with open(path, "r", errors="replace") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].startswith("rs"):
                out[f[0]] = f[1]
    return out


def ref_bases(fasta, chrom_prefixed, positions, tag):
    """One samtools faidx call for a whole chromosome's worth of positions.

    Returns a list of single uppercase bases, aligned with `positions`. A
    position samtools cannot serve comes back as 'N' rather than shifting
    every later base by one, which is the failure mode that matters here.
    """
    if not positions:
        return []
    regions = "\n".join(f"{chrom_prefixed}:{p}-{p}" for p in positions)
    try:
        proc = subprocess.run(
            ["samtools", "faidx", fasta, "-r", "-"],
            input=regions, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        sys.exit("ERROR: samtools not found on PATH.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR: samtools faidx failed on {tag} {chrom_prefixed}:\n"
                 f"{e.stderr.strip()[:800]}")
    bases, cur = [], None
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if cur is not None:
                bases.append(cur)
            cur = ""
        elif cur is not None:
            cur += line.strip().upper()
    if cur is not None:
        bases.append(cur)
    if len(bases) != len(positions):
        sys.exit(f"ERROR: {tag} {chrom_prefixed}: asked for {len(positions)} bases, "
                 f"samtools returned {len(bases)}. Refusing to emit misaligned "
                 f"REF alleles. Check that the FASTA contig naming matches "
                 f"'{chrom_prefixed}'.")
    return [b if b in COMP else "N" for b in bases]


def orient(a1, a2, ref):
    """(alt, orientation) against a reference base, or (None, 'unresolved')."""
    if ref not in COMP:
        return None, "no_ref"
    if a1 == ref:
        return a2, "direct"
    if a2 == ref:
        return a1, "direct"
    r1, r2 = rc(a1), rc(a2)
    if r1 == ref:
        return r2, "revcomp"
    if r2 == ref:
        return r1, "revcomp"
    return None, "unresolved"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chip", required=True,
                    help="short chip name used in filenames, e.g. gsa")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--hg38-annot", help="Illumina annotated manifest, GRCh38")
    ap.add_argument("--hg38-bed", help="vendor BED, GRCh38 (alternative/extra)")
    ap.add_argument("--hg19-annot", help="Illumina annotated manifest, GRCh37")
    ap.add_argument("--hg19-bed", help="vendor BED, GRCh37")
    ap.add_argument("--rsids", help="Illumina b150 rsid file keyed on hg19 chr:pos")
    ap.add_argument("--hg38-fasta", help="GRCh38 FASTA (indexed)")
    ap.add_argument("--hg19-fasta", help="GRCh37 FASTA (indexed)")
    ap.add_argument("--hg38-contig", default="chr",
                    help="contig prefix in the hg38 FASTA (default 'chr')")
    ap.add_argument("--hg19-contig", default="chr",
                    help="contig prefix in the hg19 FASTA (default 'chr')")
    ap.add_argument("--track-color", default="0,0,150")
    args = ap.parse_args()

    if not (args.hg38_annot or args.hg38_bed):
        sys.exit("ERROR: need --hg38-annot or --hg38-bed.")
    if not (args.hg19_annot or args.hg19_bed):
        sys.exit("ERROR: need --hg19-annot or --hg19-bed.")

    # ---- load -----------------------------------------------------------
    # Alleles are a property of the probe DESIGN, so the same array version
    # reports the same pair in either build's manifest. Take them from
    # whichever annotated file was given; prefer hg38.
    hg38 = read_annotated(args.hg38_annot) if args.hg38_annot else {}
    hg19 = read_annotated(args.hg19_annot) if args.hg19_annot else {}
    if args.hg38_bed:
        for probe, (c, p) in read_bed(args.hg38_bed).items():
            if probe not in hg38:
                hg38[probe] = (c, p, None)
    if args.hg19_bed:
        for probe, (c, p) in read_bed(args.hg19_bed).items():
            if probe not in hg19:
                hg19[probe] = (c, p, None)

    alleles = {}
    for src in (hg38, hg19):
        for probe, (_, _, al) in src.items():
            if al and probe not in alleles:
                alleles[probe] = al

    rsids = read_rsids(args.rsids)
    print(f"[{args.chip}] hg38 probes : {len(hg38)}")
    print(f"[{args.chip}] hg19 probes : {len(hg19)}")
    print(f"[{args.chip}] with alleles: {len(alleles)}")
    print(f"[{args.chip}] rsids       : {len(rsids)}")

    # ---- join and filter ------------------------------------------------
    drops = Counter()
    kept = []                     # (chrom, pos38, pos19, a1, a2, probe)
    for probe, (c38, p38, _) in hg38.items():
        if probe not in hg19:
            drops["not_in_hg19_manifest"] += 1
            continue
        c19, p19 = hg19[probe][0], hg19[probe][1]
        if c38 not in AUTOSOMES:
            drops["hg38_not_autosomal"] += 1
            continue
        if c19 != c38:
            drops["chromosome_differs_between_builds"] += 1
            continue
        al = alleles.get(probe)
        if al is None:
            drops["not_a_biallelic_snv"] += 1
            continue
        kept.append((c38, p38, p19, al[0], al[1], probe))
    for probe in hg19:
        if probe not in hg38:
            drops["not_in_hg38_manifest"] += 1
    print(f"[{args.chip}] joined, autosomal, biallelic SNV: {len(kept)}")

    # ---- one probe per position, in BOTH builds -------------------------
    # Arrays carry more than one probe at a multi-allelic site: Illumina names
    # the second one `seq-rsNNNN.1` beside `seq-rsNNNN`, and the pair reports
    # different allele pairs at the same coordinate (GSA v1.0: 567 such
    # positions). The contract says every variant is a BIALLELIC SNV, so a
    # position whose probes disagree on the alleles is not one site and the
    # whole group is dropped. A position whose probes agree is one site seen
    # twice, and is collapsed to one row.
    #
    # Both builds are checked, because a pair can collide in one build and not
    # the other, and a duplicate left in either list breaks the chip VCF's
    # "manifest lists exactly the sites" check and crops the position twice.
    for build_idx, build_name in ((1, "hg38"), (2, "hg19")):
        groups = {}
        for r in kept:
            groups.setdefault((r[0], r[build_idx]), []).append(r)
        survivors, seen = [], set()
        for key, grp in groups.items():
            if len(grp) == 1:
                survivors.append(grp[0])
                continue
            pairs = {frozenset((g[3], g[4])) for g in grp}
            if len(pairs) > 1:
                drops[f"duplicate_position_multiallelic_{build_name}"] += len(grp)
                continue
            drops[f"duplicate_position_collapsed_{build_name}"] += len(grp) - 1
            survivors.append(sorted(grp, key=lambda g: g[5])[0])
        kept = survivors
    print(f"[{args.chip}] one probe per position in both builds: {len(kept)}")

    os.makedirs(args.out, exist_ok=True)
    os.makedirs(os.path.join(args.out, "manifests"), exist_ok=True)

    # ---- resolve REF/ALT against each build's FASTA ---------------------
    # This runs BEFORE anything is written, because it drops sites (a probe
    # whose alleles match neither the reference base nor its complement). If
    # the .snps lists and the BED were written from the pre-resolution set they
    # would name sites the manifests do not, and the intake's 3.manifest_sites
    # check ("each chip x build manifest lists exactly the sites of its VCFs")
    # would fail on a set that looked fine.
    have_fasta = bool(args.hg38_fasta and args.hg19_fasta)
    resolved = []        # kept-tuple + (r38, alt38, r19, alt19, strand, pal, o38, o19)
    orient_counts = Counter()
    strand_counts = Counter()
    if have_fasta:
        for c in AUTOSOMES:
            sub = sorted((r for r in kept if r[0] == c), key=lambda r: r[1])
            if not sub:
                continue
            b38 = ref_bases(args.hg38_fasta, f"{args.hg38_contig}{c}",
                            [r[1] for r in sub], "hg38")
            b19 = ref_bases(args.hg19_fasta, f"{args.hg19_contig}{c}",
                            [r[2] for r in sub], "hg19")
            for rec, r38, r19 in zip(sub, b38, b19):
                _, _, _, a1, a2, _ = rec
                alt38, o38 = orient(a1, a2, r38)
                alt19, o19 = orient(a1, a2, r19)
                if alt38 is None:
                    orient_counts[f"hg38_{o38}"] += 1
                    drops[f"hg38_ref_{o38}"] += 1
                    continue
                if alt19 is None:
                    orient_counts[f"hg19_{o19}"] += 1
                    drops[f"hg19_ref_{o19}"] += 1
                    continue
                orient_counts[f"hg38_{o38}"] += 1
                orient_counts[f"hg19_{o19}"] += 1

                is_pal = frozenset((a1, a2)) in PALINDROMIC_PAIRS
                # grch38_strand: does the GRCh37 + strand read as the GRCh38 +
                # strand here? Compare each build's {REF,ALT} set. Equal -> +,
                # reverse complement of each other -> -. A palindrome's set is
                # its own complement and satisfies both, so Document 20 takes
                # those as +.
                s38 = frozenset((r38, alt38))
                s19 = frozenset((r19, alt19))
                if is_pal or s19 == s38:
                    strand = "+"
                elif s19 == frozenset(rc(b) for b in s38):
                    strand = "-"
                else:
                    drops["ref_alt_disagree_between_builds"] += 1
                    continue
                strand_counts[strand] += 1
                resolved.append(rec + (r38, alt38, r19, alt19, strand,
                                       is_pal, o38, o19))
        emit = [r[:6] for r in resolved]
        print(f"[{args.chip}] resolved against both references: {len(emit)}")
    else:
        emit = kept

    # ---- UCSC custom tracks ---------------------------------------------
    # 4-column BED, 0-based start = 1-based position - 1. The vendor's own BED
    # has a TAB between 'track' and 'name=', which UCSC will not parse; the
    # track line here is space-separated, as UCSC requires.
    # The feature name is the probe's dbSNP b150 rsID where the vendor gives
    # one, else the probe ID. (Both files key on the probe ID.)
    for build, idx in (("hg38", 1), ("hg19", 2)):
        path = os.path.join(args.out, f"{args.chip}.{build}.ucsc.bed")
        rows = sorted(((r[0], r[idx], rsids.get(r[5], r[5])) for r in emit),
                      key=lambda t: (int(t[0]), t[1]))
        with open(path, "w") as fh:
            fh.write(f'track name="{args.chip.upper()}_{build}" '
                     f'description="{args.chip.upper()} array sites ({build}, '
                     f'autosomal biallelic SNVs, both builds agree)" '
                     f'visibility=pack color={args.track_color}\n')
            for c, p, name in rows:
                fh.write(f"chr{c}\t{p - 1}\t{p}\t{name}\n")
        print(f"[{args.chip}] UCSC track -> {path}  ({len(rows)} features)")

    # ---- per-chromosome .snps drop-ins (no FASTA needed) ----------------
    for build, idx in (("hg38", 1), ("hg19", 2)):
        d = os.path.join(args.out, build)
        os.makedirs(d, exist_ok=True)
        by_chr = {c: [] for c in AUTOSOMES}
        for r in emit:
            by_chr[r[0]].append(r[idx])
        for c in AUTOSOMES:
            with open(os.path.join(d, f"chr{c}.snps"), "w") as fh:
                for p in sorted(set(by_chr[c])):
                    fh.write(f"{c}:{p}\n")
        print(f"[{args.chip}] {build} .snps lists -> {d}")

    # ---- report ---------------------------------------------------------
    pal = sum(1 for r in emit if frozenset((r[3], r[4])) in PALINDROMIC_PAIRS)
    report = os.path.join(args.out, "build_report.tsv")
    with open(report, "w") as fh:
        fh.write("stage\tcount\n")
        fh.write(f"hg38_manifest_probes\t{len(hg38)}\n")
        fh.write(f"hg19_manifest_probes\t{len(hg19)}\n")
        for k in sorted(drops):
            fh.write(f"dropped_{k}\t{drops[k]}\n")
        fh.write(f"joined_autosomal_biallelic_snv\t{len(kept)}\n")
        for k in sorted(orient_counts):
            fh.write(f"orient_{k}\t{orient_counts[k]}\n")
        for k in sorted(strand_counts):
            fh.write(f"grch38_strand_{k}\t{strand_counts[k]}\n")
        fh.write(f"final_sites\t{len(emit)}\n")
        fh.write(f"palindromic_at_cg\t{pal}\n")
    for k in sorted(drops):
        print(f"[{args.chip}]   dropped {k}: {drops[k]}")
    for k in sorted(orient_counts):
        print(f"[{args.chip}]   {k}: {orient_counts[k]}")
    for k in sorted(strand_counts):
        print(f"[{args.chip}]   grch38_strand {k}: {strand_counts[k]}")
    pct = 100.0 * pal / len(emit) if emit else 0.0
    print(f"[{args.chip}] palindromic (A/T, C/G): {pal} ({pct:.2f}%)")

    # ---- provenance ------------------------------------------------------
    prov = os.path.join(args.out, "PROVENANCE.md")
    with open(prov, "w") as fh:
        fh.write(f"# {args.chip} manifest provenance\n\n")
        fh.write("| role | file | bytes | sha256 |\n|---|---|---|---|\n")
        for role, p in (("hg38 annotated", args.hg38_annot),
                        ("hg38 bed", args.hg38_bed),
                        ("hg19 annotated", args.hg19_annot),
                        ("hg19 bed", args.hg19_bed),
                        ("rsids", args.rsids)):
            if p:
                fh.write(f"| {role} | `{os.path.basename(p)}` | "
                         f"{os.path.getsize(p)} | `{sha256(p)}` |\n")
    print(f"[{args.chip}] provenance -> {prov}")

    # ---- the contract manifests (FASTA required) ------------------------
    if not have_fasta:
        print(f"\n[{args.chip}] NO FASTA GIVEN -- manifests/ not written, and")
        print(f"[{args.chip}] the .snps lists and BED above are PROVISIONAL:")
        print(f"[{args.chip}] they still contain probes whose alleles may match")
        print(f"[{args.chip}] neither reference base. REF/ALT cannot be decided")
        print(f"[{args.chip}] without the reference sequence, because the vendor")
        print(f"[{args.chip}] allele column is on the design strand, not the")
        print(f"[{args.chip}] genome strand. Re-run on the cluster with")
        print(f"[{args.chip}] --hg38-fasta and --hg19-fasta.")
        return 3

    rows38, rows19, annot38, annot19 = [], [], [], []
    for rec in resolved:
        chrom, p38, p19, a1, a2, probe = rec[:6]
        r38, alt38, r19, alt19, strand, is_pal = rec[6:12]
        o38, o19 = rec[12], rec[13]
        g38 = f"chr{chrom}:{p38}:{r38}:{alt38}"
        rows38.append((chrom, p38, r38, alt38))
        rows19.append((chrom, p19, r19, alt19, g38, strand))
        rs = rsids.get(probe, "")
        annot38.append((chrom, p38, probe, rs, int(is_pal), o38))
        annot19.append((chrom, p19, probe, rs, int(is_pal), o19))
    rows38.sort(key=lambda t: (int(t[0]), t[1]))
    rows19.sort(key=lambda t: (int(t[0]), t[1]))
    annot38.sort(key=lambda t: (int(t[0]), t[1]))
    annot19.sort(key=lambda t: (int(t[0]), t[1]))

    m38 = os.path.join(args.out, "manifests", f"{args.chip}.hg38.sites.tsv")
    with open(m38, "w") as fh:
        fh.write("chrom\tpos\tref\talt\n")
        for c, p, r, a in rows38:
            fh.write(f"chr{c}\t{p}\t{r}\t{a}\n")
    m19 = os.path.join(args.out, "manifests", f"{args.chip}.hg19.sites.tsv")
    with open(m19, "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgrch38_site\tgrch38_strand\n")
        for c, p, r, a, g, s in rows19:
            fh.write(f"{c}\t{p}\t{r}\t{a}\t{g}\t{s}\n")
    for name, rows in ((f"{args.chip}.hg38.annot.tsv", annot38),
                       (f"{args.chip}.hg19.annot.tsv", annot19)):
        with open(os.path.join(args.out, name), "w") as fh:
            fh.write("chrom\tpos\tprobe\trsid\tpalindromic\torientation\n")
            for c, p, probe, rs, pl, o in rows:
                fh.write(f"{c}\t{p}\t{probe}\t{rs}\t{pl}\t{o}\n")

    print(f"\n[{args.chip}] manifests written: {len(rows38)} sites")
    print(f"[{args.chip}]   {m38}")
    print(f"[{args.chip}]   {m19}")
    # Document 20 note: hg19 contigs are written bare (1..22) and hg38 with the
    # chr prefix, matching "one contig naming per file" in the contract.
    return 0


if __name__ == "__main__":
    sys.exit(main())
