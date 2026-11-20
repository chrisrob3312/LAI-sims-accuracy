#!/usr/bin/env python3
"""
lift-cell-to-hg19.py -- rewrite a GRCh38 cell VCF onto the chip's NATIVE hg19
coordinates, flipping GT where the two builds disagree about which allele is
reference.

WHY THIS IS NOT A LIFTOVER. The hg19 cells of the benchmark must use the
chip's own hg19 manifest, never a down-lift of the hg38 sites -- a down-lift
contains only the probes that map back cleanly, by construction, and so hides
every probe the pipeline's own CrossMap drops, moves or strand-flips, which is
part of what the benchmark exists to measure. So the mapping here is the
vendor's: probe -> hg19 position AND probe -> hg38 position, joined on probe
id by build-chip-manifest-sites.py, which emits grch38_site alongside each
hg19 row. This script only transcribes genotypes from the hg38 site to the
hg19 site that the same probe interrogates.

THE PART THAT CAN GO SILENTLY WRONG IS THE ALLELE ORIENTATION. GT is an INDEX
into REF/ALT, so a genotype is only meaningful relative to the REF/ALT pair
written beside it. Between builds a site can:

  same      (R19,A19) == (R38,A38)              -> GT unchanged
  flip      (R19,A19) == complement(R38,A38)    -> GT unchanged. A strand flip
                                                   complements both alleles,
                                                   so the INDICES are
                                                   untouched.
  swap      (R19,A19) == (A38,R38)              -> GT FLIPPED 0<->1. The
                                                   reference allele changed
                                                   between builds, which is
                                                   common and is the case that
                                                   silently inverts dosage if
                                                   missed.
  flip+swap (R19,A19) == complement(A38,R38)    -> GT FLIPPED
  none      no correspondence                   -> site DROPPED

Getting "swap" wrong does not fail: it writes a valid VCF in which the
alternate allele dosage is inverted at those sites, and an imputation
benchmark run on it would report inflated error concentrated at exactly the
sites where the builds disagree -- a plausible-looking result that would be
read as a pipeline problem.

PALINDROMIC SITES ARE IRREDUCIBLY AMBIGUOUS HERE. At an A/T or C/G site,
complement(R,A) == (A,R), so "flip" and "swap" are indistinguishable from the
alleles alone and the choice changes the answer. Those sites are resolved from
the manifest's grch38_strand column instead, and COUNTED SEPARATELY in the
report so the ambiguity is visible rather than absorbed. GSA has 6,422 of them
(1.08% of 596,591).

Usage:
  bcftools view -Ov cell.hg38.bcf \\
    | python3 -I scripts/lift-cell-to-hg19.py \\
        --manifest manifests/gsa.hg19.sites.tsv --report lift.chr1.txt \\
    | bcftools sort -Ob -o cell.hg19.bcf
"""
import argparse
import sys
from collections import Counter

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def comp(a):
    try:
        return "".join(COMP[c] for c in reversed(a.upper()))
    except KeyError:
        return None


def is_palindromic(r, a):
    return len(r) == 1 and len(a) == 1 and comp(r) == a.upper()


def classify(r38, a38, r19, a19, strand):
    """-> (flip_gt, how). flip_gt True means GT indices must be swapped."""
    r38, a38 = r38.upper(), a38.upper()
    r19, a19 = r19.upper(), a19.upper()

    if is_palindromic(r38, a38) and is_palindromic(r19, a19):
        # Alleles cannot settle it: complement == swap at a palindrome. Fall
        # back to the vendor strand, and say so in the report rather than
        # letting it pass as an allele-resolved site.
        if strand in ("-", "minus", "REV", "rev"):
            return False, "palindromic_strand_minus"
        return False, "palindromic_strand_plus"

    if (r19, a19) == (r38, a38):
        return False, "same"
    c38r, c38a = comp(r38), comp(a38)
    if c38r and (r19, a19) == (c38r, c38a):
        return False, "flip"
    if (r19, a19) == (a38, r38):
        return True, "swap"
    if c38r and (r19, a19) == (c38a, c38r):
        return True, "flip_swap"
    return None, "unresolvable"


def flip_gt_field(fmt_keys, sample_field):
    """Swap allele indices 0<->1 in the GT subfield only."""
    try:
        gi = fmt_keys.index("GT")
    except ValueError:
        return sample_field
    parts = sample_field.split(":")
    if gi >= len(parts):
        return sample_field
    gt = parts[gi]
    out = []
    for ch in gt:
        if ch == "0":
            out.append("1")
        elif ch == "1":
            out.append("0")
        else:
            out.append(ch)       # '.', '|', '/' and anything else untouched
    parts[gi] = "".join(out)
    return ":".join(parts)


def load_manifest(path):
    """hg19 site rows keyed by the GRCh38 site they correspond to."""
    by38 = {}
    dup = 0
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for c in ("chrom", "pos", "ref", "alt", "grch38_site"):
            if c not in ix:
                sys.exit(f"FATAL: {path} lacks column '{c}'. This must be the "
                         f"hg19 manifest from build-chip-manifest-sites.py, "
                         f"which emits grch38_site; the hg38 manifest has no "
                         f"such column and cannot drive this rewrite.")
        i_strand = ix.get("grch38_strand", -1)
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(head):
                continue
            key = p[ix["grch38_site"]].strip()
            if not key:
                continue
            if key in by38:
                dup += 1
                continue
            by38[key] = (p[ix["chrom"]], p[ix["pos"]], p[ix["ref"]],
                         p[ix["alt"]],
                         p[i_strand] if 0 <= i_strand < len(p) else "")
    return by38, dup


def norm_key(chrom, pos):
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return f"{c}:{pos}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="<chip>.hg19.sites.tsv (must carry grch38_site)")
    ap.add_argument("--report", default="",
                    help="where to write the per-class counts; stderr if unset")
    ap.add_argument("--hg19-chr-prefix", default="auto",
                    choices=("auto", "yes", "no"),
                    help="emit hg19 CHROM with or without a 'chr' prefix; "
                         "auto follows the manifest")
    a = ap.parse_args()

    by38, dup = load_manifest(a.manifest)
    if not by38:
        sys.exit(f"FATAL: no usable rows in {a.manifest}")

    counts = Counter()
    n_in = n_out = 0
    for line in sys.stdin:
        if line.startswith("##"):
            # Contig lines describe GRCh38 and would be wrong on the output.
            # Drop them and let the caller's `bcftools sort`/reheader supply
            # the hg19 contigs: emitting stale contig lengths would make the
            # file look hg19-valid while carrying hg38 bounds.
            if line.startswith("##contig="):
                counts["dropped_hg38_contig_header"] += 1
                continue
            sys.stdout.write(line)
            continue
        if line.startswith("#CHROM"):
            sys.stdout.write("##liftCellToHg19=native-hg19-manifest;"
                             "source=GRCh38;gt_reoriented=1\n")
            sys.stdout.write(line)
            continue
        n_in += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            counts["malformed"] += 1
            continue
        key = norm_key(f[0], f[1])
        row = by38.get(key)
        if row is None:
            counts["not_on_chip"] += 1
            continue
        c19, p19, r19, a19, strand = row
        if "," in f[4]:
            # Chip sites are biallelic by construction: the manifest builder
            # drops multi-allelic probe pairs. One here means the upstream
            # subset did not use the chip site list.
            counts["multiallelic_dropped"] += 1
            continue
        flip, how = classify(f[3], f[4], r19, a19, strand)
        counts[how] += 1
        if flip is None:
            continue
        if a.hg19_chr_prefix == "yes":
            out_chrom = c19 if c19.lower().startswith("chr") else "chr" + c19
        elif a.hg19_chr_prefix == "no":
            out_chrom = c19[3:] if c19.lower().startswith("chr") else c19
        else:
            out_chrom = c19
        f[0], f[1], f[3], f[4] = out_chrom, p19, r19, a19
        if flip:
            fmt = f[8].split(":")
            f[9:] = [flip_gt_field(fmt, s) for s in f[9:]]
        sys.stdout.write("\t".join(f) + "\n")
        n_out += 1

    lines = ["GRCh38 -> native hg19 rewrite",
             f"  manifest rows (unique grch38_site): {len(by38)}",
             f"  duplicate grch38_site in manifest:  {dup}",
             f"  input records:  {n_in}",
             f"  output records: {n_out}",
             "",
             "  GT unchanged:"]
    for k in ("same", "flip"):
        lines.append(f"    {k:<28}{counts[k]}")
    lines.append("  GT FLIPPED (reference allele differs between builds):")
    for k in ("swap", "flip_swap"):
        lines.append(f"    {k:<28}{counts[k]}")
    lines.append("  palindromic -- resolved by vendor strand, NOT by alleles:")
    for k in ("palindromic_strand_plus", "palindromic_strand_minus"):
        lines.append(f"    {k:<28}{counts[k]}")
    lines.append("  dropped:")
    for k in ("unresolvable", "not_on_chip", "multiallelic_dropped",
              "malformed"):
        lines.append(f"    {k:<28}{counts[k]}")
    pal = (counts["palindromic_strand_plus"]
           + counts["palindromic_strand_minus"])
    if pal:
        lines += ["",
                  f"  {pal} site(s) are palindromic (A/T or C/G): complement "
                  f"and swap are",
                  "  indistinguishable from the alleles there, so orientation "
                  "came from the",
                  "  vendor strand column alone. A wrong strand at those sites "
                  "inverts dosage",
                  "  without any signal in this file."]
    if counts["unresolvable"]:
        lines += ["",
                  f"  WARNING: {counts['unresolvable']} site(s) had no allele "
                  f"correspondence between",
                  "  builds and were dropped. A large count means the manifest "
                  "and the VCF",
                  "  disagree about which reference the hg38 side was called "
                  "against."]
    text = "\n".join(lines) + "\n"
    if a.report:
        with open(a.report, "w") as g:
            g.write(text)
    else:
        sys.stderr.write(text)


if __name__ == "__main__":
    main()
