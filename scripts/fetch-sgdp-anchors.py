#!/usr/bin/env python3
"""
Fetch SGDP (Simons Genome Diversity Project) metadata and emit anchor sample
lists for underrepresented super-populations (SAS, Oceania, Middle East).

Purpose:
  Complementary anchors for a SECONDARY unsupervised-ADMIXTURE run on all
  refs + additions -> global proportions for downstream cohorts that may
  carry SAS/OCE/MEN ancestry. Not the primary LAI panel (that stays at K=5
  HGDP+1KG+MXB).

Sources:
  * SGDP public metadata TSV (Reich lab / Harvard Medical School)
      https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt
  * SGDP phased hg38 VCFs (EBI-EMBL European Nucleotide Archive):
      https://www.ebi.ac.uk/ena/browser/view/PRJEB9586
    Individual FTP path pattern per sample under
      ftp://ftp.sra.ebi.ac.uk/vol1/fastq/... (per-sample; script only writes
      the SAMPLE-ID list; you handle FTP/download separately or point to a
      local reprocessed release).

Output:
  reference_ids/sgdp_sas_anchors.txt        SAS samples (per subpopulation majority)
  reference_ids/sgdp_oceania_anchors.txt    Papuan / Melanesian / Australian
  reference_ids/sgdp_men_anchors.txt        Middle Eastern (Bedouin, Iranian, ...)
  reference_ids/sgdp_amr_anchors.txt        Native American (Karitiana, Surui,
                                            Mayan, Pima, Mixe, Zapotec, Piapoco...)
  reference_ids/sgdp_manifest.tsv           full filtered metadata rows for all
                                            four groups above -- provenance for
                                            the paper + easy re-download

Companion sources (each needs its own ingestion script, not scraped here):
  * Jimenez-Kaufmann et al. (Moreno-Estrada) rare-variant imputation reference
    panel: publicly available InMEGEN / MAIS Amerindigenous sample IDs are in
    the supplementary info. Individual genotypes are dbGaP-restricted; the
    Atkinson lab has a DUA with Moreno-Estrada for the MXB portion, so an
    InMEGEN ingestion script is feasible.
  * LASI-DAD (Li et al., HGG Adv 2026): 2,680 Indian WGS participants,
    publicly available, 69.5M variants. Highest-yield SAS panel currently
    available (SGDP SAS n~30 is a floor; LASI-DAD is ~100x that). Access
    documented at PMC12945573.
  * SG10K_Health (Chang et al., HMG 2026): ~10K Singaporean (Chinese, Malay,
    Indian) WGS. Consortium panel -- access via SG10K_Health MTA. Higher
    ROI for EAS refinement than SAS.

Usage:
  scripts/fetch-sgdp-anchors.py \
      [--metadata-tsv PATH] \
      [--outdir reference_ids] \
      [--per-pop-cap N] \
      [--include-signed-letter] \
      [--include-fan]

  --metadata-tsv     local path to SGDP metadata TSV. If omitted, script
                     downloads the public version to a temp file.
  --per-pop-cap N    optional cap per SGDP sub-population (Papuan, Sindhi, ...).
                     Default: no cap.
  --include-signed-letter  include the 21 samples requiring a signed data-use
                     letter (default: exclude, only public 279)
  --include-fan      include the 44 Fan samples (default: exclude)

Notes:
  * SGDP designates samples in three tiers: 279 fully public, 21 requiring a
    signed letter, 44 "Fan" release. This script defaults to fully public.
  * The metadata TSV has one row per sample. Columns of interest:
      Sample_ID, Population_ID, Region, Country, Latitude, Longitude,
      Sex, Sequencing_source, Contributor
  * "Region" is the super-pop key we use to bucket into SAS / OCE / MEN.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict, Counter

DEFAULT_METADATA_URL = (
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
    "SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"
)

# Region tag -> our super-population bucket.
# SGDP "Region" values seen in the public metadata:
#   Africa, America, CentralAsiaSiberia, EastAsia, WestEurasia, SouthAsia, Oceania
# "WestEurasia" is Europe + Middle East -> refine by Population_ID.
# "America" is Native American -> whitelist homogeneous NAT-like populations.
SGDP_TO_SUPERPOP = {
    "SouthAsia":         "SAS",
    "Oceania":           "OCE",
}

# SGDP "America" populations that are homogeneous NAT (per Reich lab curation
# in Mallick et al. 2016 Table S1). Excludes Mexican-American / admixed cohorts.
AMR_POPULATIONS = {
    "Karitiana", "Surui", "Mayan", "Pima", "Mixe",
    "Zapotec", "Piapoco", "Quechua", "Chane",
}

# WestEurasia populations that we bucket as Middle Eastern.
# From the SGDP paper (Mallick et al., Nature 2016) supplementary Table S1.
MEN_POPULATIONS = {
    "Bedouin",       "Palestinian", "Druze", "Mozabite", "Saharawi",
    "Iranian",       "Yemenite",    "Jordanian",
    "Assyrian",      "Georgian",    "Armenian",
    "Turkish",       "Kurdish",     "Iraqi_Jew",
    "Samaritan",     "Egyptian",    "Lebanese",
    "Chuvash",       # borderline Eurasian; drop if you want strict ME.
}

def download_metadata(dst_path: str) -> None:
    """Fetch via urllib; fall back to curl/wget when the cluster's Python has no
    CA bundle (common HPC issue: SSLCertVerificationError). curl/wget on the
    cluster inherits the system trust store and Just Works."""
    print(f"[fetch] downloading SGDP metadata -> {dst_path}", file=sys.stderr)
    try:
        urllib.request.urlretrieve(DEFAULT_METADATA_URL, dst_path)
        return
    except Exception as e:
        print(f"[fetch] urllib failed ({type(e).__name__}: {e}); "
              f"trying curl/wget", file=sys.stderr)
    for tool, argv in (("curl", ["curl", "-sSL", "-o", dst_path, DEFAULT_METADATA_URL]),
                       ("wget", ["wget", "-q", "-O", dst_path, DEFAULT_METADATA_URL])):
        if shutil.which(tool) is None:
            continue
        r = subprocess.run(argv)
        if r.returncode == 0 and os.path.getsize(dst_path) > 0:
            print(f"[fetch] downloaded via {tool}", file=sys.stderr)
            return
    raise RuntimeError(
        "All download attempts failed. Fetch the file manually with:\n"
        f"  curl -sSL -o sgdp_metadata.tsv '{DEFAULT_METADATA_URL}'\n"
        "and pass it with --metadata-tsv sgdp_metadata.tsv")


def parse_metadata(path: str, include_signed_letter: bool, include_fan: bool):
    """Yield sample dict rows from SGDP metadata TSV, respecting tier filters."""
    with open(path, encoding="utf-8", errors="replace") as f:
        # File is whitespace/tab-delimited with a header comment block.
        # Find the header row starting with 'Sample_ID' or similar.
        reader = csv.DictReader(
            (line for line in f if not line.startswith("#") and line.strip()),
            delimiter="\t")
        for row in reader:
            # SGDP uses one of these column names; normalize.
            sample = (row.get("Sample_ID") or row.get("SGDP_ID")
                      or row.get("SampleID") or row.get("Illumina_ID"))
            pop    = row.get("Population_ID") or row.get("Population")
            region = row.get("Region")
            source = (row.get("Sequencing_source") or row.get("Panel") or "").strip()
            if not sample or not pop or not region:
                continue
            src_lower = source.lower()
            if "letter" in src_lower and not include_signed_letter:
                continue
            if "fan" in src_lower and not include_fan:
                continue
            yield {
                "sample": sample.strip(),
                "population": pop.strip(),
                "region": region.strip(),
                "sequencing_source": source,
                "country": (row.get("Country") or "").strip(),
                "sex": (row.get("Sex") or "").strip(),
                "latitude": (row.get("Latitude") or "").strip(),
                "longitude": (row.get("Longitude") or "").strip(),
            }


def bucket(rows, cap_per_pop=None):
    """Assign each row to SAS/OCE/MEN/AMR or drop. Optional per-population cap."""
    per_pop_seen = Counter()
    buckets = defaultdict(list)
    for r in rows:
        region = r["region"]; pop = r["population"]
        super_pop = SGDP_TO_SUPERPOP.get(region)
        if super_pop is None:
            if region == "WestEurasia" and pop in MEN_POPULATIONS:
                super_pop = "MEN"
            elif region == "America" and pop in AMR_POPULATIONS:
                super_pop = "AMR"
            else:
                continue
        if cap_per_pop and per_pop_seen[(super_pop, pop)] >= cap_per_pop:
            continue
        per_pop_seen[(super_pop, pop)] += 1
        r_out = {"super_pop": super_pop, **r}
        buckets[super_pop].append(r_out)
    return buckets


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata-tsv", default=None,
                    help="Local path to SGDP metadata TSV. If omitted, download it.")
    ap.add_argument("--outdir", default="reference_ids",
                    help="Where to write anchor lists (default: reference_ids/).")
    ap.add_argument("--per-pop-cap", type=int, default=None,
                    help="Cap per SGDP sub-population (default: no cap).")
    ap.add_argument("--include-signed-letter", action="store_true",
                    help="Include the 21-sample signed-letter tier.")
    ap.add_argument("--include-fan", action="store_true",
                    help="Include the 44 Fan release samples.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tsv = args.metadata_tsv
    if tsv is None:
        tsv = os.path.join(args.outdir, ".sgdp_metadata.tsv")
        if not os.path.exists(tsv):
            download_metadata(tsv)

    rows = list(parse_metadata(tsv,
                               include_signed_letter=args.include_signed_letter,
                               include_fan=args.include_fan))
    print(f"[fetch] parsed {len(rows)} eligible SGDP rows", file=sys.stderr)

    buckets = bucket(rows, cap_per_pop=args.per_pop_cap)

    # Anchor lists (one sample per line, sorted for deterministic output).
    for sp_key, out_name in [("SAS", "sgdp_sas_anchors.txt"),
                             ("OCE", "sgdp_oceania_anchors.txt"),
                             ("MEN", "sgdp_men_anchors.txt"),
                             ("AMR", "sgdp_amr_anchors.txt")]:
        out_path = os.path.join(args.outdir, out_name)
        samples = sorted({r["sample"] for r in buckets.get(sp_key, [])})
        with open(out_path, "w") as g:
            g.write("\n".join(samples) + ("\n" if samples else ""))
        pops = Counter(r["population"] for r in buckets.get(sp_key, []))
        pop_summary = ", ".join(f"{p}={n}" for p, n in pops.most_common())
        print(f"[fetch] {sp_key:>3}: wrote {len(samples):3d} anchors "
              f"-> {out_path}  ({pop_summary or 'no matches'})", file=sys.stderr)

    # Combined manifest for provenance.
    manifest_path = os.path.join(args.outdir, "sgdp_manifest.tsv")
    fields = ["super_pop", "sample", "population", "region", "country",
              "sex", "latitude", "longitude", "sequencing_source"]
    with open(manifest_path, "w") as g:
        g.write("\t".join(fields) + "\n")
        for sp_key in ("SAS", "OCE", "MEN", "AMR"):
            for r in sorted(buckets.get(sp_key, []),
                            key=lambda x: (x["population"], x["sample"])):
                g.write("\t".join(str(r.get(k, "")) for k in fields) + "\n")
    print(f"[fetch] wrote manifest -> {manifest_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
