#!/usr/bin/env python3
"""
inventory-chip-manifests.py -- what chip manifests are on disk, what build each
covers, and which chips can actually be built as a GRCh37+GRCh38 pair.

WHY THIS IS NEEDED. build-chip-manifest-sites.py joins a GRCh37 manifest and a
GRCh38 manifest on probe id, so a chip is buildable only if BOTH native builds
are present for it. Its reader is also Illumina-specific: read_annotated()
requires tab-separated Name / Chr / MapInfo / Alleles and exits otherwise, so
an Affymetrix annot.csv -- comma-separated, with Probe Set ID / Chromosome /
Physical Position / Allele A / Allele B -- cannot be read at all yet.

So "do we have the manifests?" has three different answers per chip, and only
one of them means work can start:

  PAIR READY        both builds, both in a readable format -> build it
  ONE BUILD ONLY    the other build's vendor manifest is missing. Do NOT
                    substitute a liftover of the one you have: the hg19 cells
                    exist to exercise the pipeline's own liftover, and a
                    down-lifted manifest contains only probes that map back
                    cleanly, which hides exactly what is being measured.
  NEEDS A READER    the file is present but in a format the builder does not
                    parse.

READ-ONLY. It opens files to sniff headers and never writes, moves or extracts
anything. Vendor downloads are untrusted input, so run it with `python3 -I`
and pass the directory as an argument rather than running from inside it.

Usage:
  python3 -I scripts/inventory-chip-manifests.py ~/chip_manifests ~/manifests
"""
import gzip
import os
import re
import sys
import zipfile
from collections import defaultdict

CAND_EXT = (".csv", ".txt", ".tsv", ".bpm", ".bed", ".annot",
            ".csv.gz", ".txt.gz", ".tsv.gz", ".zip")

BUILD_PATTERNS = [
    (re.compile(r"(?:^|[^0-9a-z])(?:hg38|grch38|b38|build38)(?:[^0-9a-z]|$)", re.I), "GRCh38"),
    (re.compile(r"(?:^|[^0-9a-z])(?:hg19|grch37|b37|build37)(?:[^0-9a-z]|$)", re.I), "GRCh37"),
    (re.compile(r"na3[6-9]", re.I), "GRCh38?"),   # Affy netaffx release, verify
    (re.compile(r"na3[0-5]", re.I), "GRCh37?"),
]

HEADER_BUILD = [
    (re.compile(r"genome[-_ ]?version[^0-9]*(3[78])", re.I), None),
    (re.compile(r"GenomeBuild[^0-9]*(3[78])", re.I), None),
    (re.compile(r"\bGRCh(3[78])\b", re.I), None),
    (re.compile(r"\bhg(19|38)\b", re.I), None),
]


def open_text(path, inner=None):
    if path.lower().endswith(".zip"):
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist()
                 if n.lower().endswith((".csv", ".txt", ".tsv", ".bpm", ".annot"))]
        if not names:
            return None, None
        pick = inner or names[0]
        return (x.decode("utf-8", "replace")
                for x in z.open(pick).read().split(b"\n")), f"{path}!{pick}"
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace"), path
    return open(path, "r", errors="replace"), path


def sniff(path):
    """-> dict describing the file, or None if it is not a manifest."""
    try:
        fh, shown = open_text(path)
    except Exception as e:
        return {"path": path, "fmt": f"unreadable ({type(e).__name__})",
                "build": "-", "probe_col": "-", "rows": "-", "readable": False}
    if fh is None:
        return None
    pre, header, nrow = [], None, 0
    try:
        for i, raw in enumerate(fh):
            line = raw.rstrip("\r\n")
            if i > 400 and header is None:
                break
            if header is None:
                if not line.strip():
                    continue
                low = line.lower()
                # Illumina manifests put the real header after [Assay];
                # Affy prefixes its comments with #% or #.
                if line.startswith(("#", "!")) or low.startswith("["):
                    pre.append(line)
                    continue
                cells_t = line.split("\t")
                cells_c = next(iter(__import__("csv").reader([line])), [])
                for cells, sep in ((cells_t, "tab"), (cells_c, "comma")):
                    names = [c.strip().strip('"') for c in cells]
                    if len(names) < 3:
                        continue
                    header = (names, sep)
                    break
                if header is None:
                    pre.append(line)
                continue
            if line.strip():
                nrow += 1
            if nrow > 5000:
                nrow = ">5000"
                break
    finally:
        try:
            fh.close()
        except Exception:
            pass
    if header is None:
        return None
    names, sep = header
    lower = [n.lower() for n in names]

    fmt, probe = "unknown", "-"
    readable = False
    if sep == "tab" and {"name", "chr", "mapinfo", "alleles"} <= set(lower):
        fmt, probe, readable = "Illumina annotated (TSV)", "Name", True
    elif {"name", "chr", "mapinfo"} <= set(lower):
        fmt, probe = f"Illumina manifest ({sep})", "Name"
        readable = (sep == "tab")
    elif {"probe set id", "chromosome", "physical position"} <= set(lower):
        fmt, probe = f"Affymetrix annot ({sep})", "Probe Set ID"
    elif {"probeset_id", "chromosome", "position"} <= set(lower):
        fmt, probe = f"Affymetrix annot alt ({sep})", "probeset_id"
    elif lower[:3] == ["chrom", "chromstart", "chromend"] or (
            len(names) >= 4 and lower[0].startswith("chr") and sep == "tab"):
        fmt, probe = "BED", "col4"
    elif {"chrom", "pos", "ref", "alt"} <= set(lower):
        fmt, probe, readable = "contract sites.tsv (already built)", "-", True

    build = "-"
    for pat, b in BUILD_PATTERNS:
        if pat.search(os.path.basename(path)):
            build = b
            break
    if build in ("-", "GRCh38?", "GRCh37?"):
        for line in pre[:60]:
            for pat, _ in HEADER_BUILD:
                m = pat.search(line)
                if m:
                    v = m.group(1)
                    build = ("GRCh38" if v in ("38",) else
                             "GRCh37" if v in ("37", "19") else build)
                    break
            if build in ("GRCh37", "GRCh38"):
                break
    return {"path": shown or path, "fmt": fmt, "build": build,
            "probe_col": probe, "rows": nrow, "readable": readable,
            "cols": names[:6]}


def chip_key(path):
    b = os.path.basename(path).lower()
    b = re.sub(r"\.(csv|txt|tsv|bpm|bed|annot|gz|zip)$", "", b)
    for pat, _ in BUILD_PATTERNS:
        b = pat.sub("_", b)
    b = re.sub(r"[^a-z0-9]+", "_", b).strip("_")
    b = re.sub(r"_(a1|a2|b1|b2|v\d+|rev\w*|manifest|annot\w*|sites?)$", "", b)
    return b or "unknown"


def main():
    dirs = sys.argv[1:]
    if not dirs:
        sys.exit(__doc__.strip())
    found = []
    for d in dirs:
        d = os.path.expanduser(d)
        if not os.path.isdir(d):
            print(f"  (not a directory, skipped: {d})")
            continue
        for root, _sub, files in os.walk(d):
            for fn in sorted(files):
                if not fn.lower().endswith(CAND_EXT):
                    continue
                p = os.path.join(root, fn)
                try:
                    if os.path.getsize(p) < 200:
                        continue
                except OSError:
                    continue
                info = sniff(p)
                if info:
                    info["chip"] = chip_key(fn)
                    found.append(info)

    if not found:
        print("No manifest-like files found in: " + ", ".join(dirs))
        print("Candidate extensions: " + ", ".join(CAND_EXT))
        return

    print("=" * 100)
    print(" CHIP MANIFEST INVENTORY")
    print("=" * 100)
    print(f" {'chip':<22}{'build':<9}{'format':<30}{'probe':<14}{'rows':>8}  file")
    print("-" * 100)
    for i in sorted(found, key=lambda x: (x["chip"], x["build"])):
        print(f" {i['chip']:<22}{i['build']:<9}{i['fmt']:<30}"
              f"{i['probe_col']:<14}{str(i['rows']):>8}  "
              f"{os.path.basename(i['path'])}")
    print("-" * 100)

    by_chip = defaultdict(list)
    for i in found:
        by_chip[i["chip"]].append(i)

    print("")
    print(" VERDICT PER CHIP")
    print("-" * 100)
    for chip in sorted(by_chip):
        items = by_chip[chip]
        builds = {i["build"] for i in items if i["build"] in ("GRCh37", "GRCh38")}
        unreadable = [i for i in items if not i["readable"]]
        has_pair = {"GRCh37", "GRCh38"} <= builds
        if has_pair and not unreadable:
            print(f"  {chip:<22} PAIR READY -- build it")
        elif has_pair:
            print(f"  {chip:<22} both builds present but NEEDS A READER:")
            for i in unreadable:
                print(f"  {'':<22}   {i['fmt']} -- {os.path.basename(i['path'])}")
                print(f"  {'':<22}   columns: {', '.join(i['cols'])}")
        elif builds:
            missing = ({"GRCh37", "GRCh38"} - builds).pop()
            print(f"  {chip:<22} ONE BUILD ONLY ({', '.join(sorted(builds))}) "
                  f"-- {missing} vendor manifest needed")
            print(f"  {'':<22}   do NOT substitute a liftover: the hg19 cells")
            print(f"  {'':<22}   exist to exercise the pipeline's own liftover")
            # Report BOTH blockers. Otherwise the missing build gets chased
            # first and the format wall is only discovered afterwards.
            if unreadable:
                print(f"  {'':<22}   AND needs a reader even once that arrives:")
                for i in unreadable:
                    print(f"  {'':<22}     {i['fmt']} -- "
                          f"{os.path.basename(i['path'])}")
                    print(f"  {'':<22}     columns: {', '.join(i['cols'])}")
        else:
            print(f"  {chip:<22} build not determined from filename or header")
            for i in items:
                print(f"  {'':<22}   {i['build']}  {i['fmt']}  "
                      f"{os.path.basename(i['path'])}")
    print("-" * 100)
    print("")
    print(" build-chip-manifest-sites.py reads Illumina annotated TSV only")
    print(" (Name / Chr / MapInfo / Alleles). Anything marked NEEDS A READER")
    print(" has the data but not in that shape -- send me the column list")
    print(" above and I will add the reader.")


if __name__ == "__main__":
    main()
