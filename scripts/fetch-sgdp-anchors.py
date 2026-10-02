#!/usr/bin/env python3
"""
Fetch SGDP (Simons Genome Diversity Project) metadata and emit CANDIDATE
sample lists for underrepresented super-populations across all 7 regions.

IMPORTANT -- these are CANDIDATES, not verified-homogeneous samples:
  The lists this script writes are selected purely by SGDP population-label
  majority from the public metadata (e.g. Region=SouthAsia -> SAS candidate,
  Population_ID in AMR_POPULATIONS -> AMR candidate). That is a metadata
  pick, NOT a genetic homogeneity test. Whether a given candidate actually
  clears the Q>=0.95 homogeneity threshold is decided DOWNSTREAM, after these
  are merged into the panel (1b) and run through supervised ADMIXTURE in
  4a_build_homogeneity_panel_clm.sh (K=5 now; K=7 once OCE/MEN are added).
  The final homogeneous reference set is whatever 4a marks included=1 -- it
  will be a SUBSET of these candidates. Do not describe the lists below as
  "the homogeneous SGDP set" in the methods; they are the pre-filter pool.

  The name "anchors" refers only to their ROLE as known-population priors
  that seed supervised ADMIXTURE -- it is not a homogeneity claim.

Purpose:
  Candidate priors for ADMIXTURE on all refs + additions -> global
  proportions for downstream cohorts that may carry SAS/OCE/MEN ancestry.
  Not the primary LAI panel (that stays at K=5 HGDP+1KG+MXB until 4a
  re-decides inclusion over the merged panel).

Sources:
  * SGDP public metadata TSV (Reich lab / Harvard Medical School)
      https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt
  * SGDP phased hg38 VCFs (EBI-EMBL European Nucleotide Archive):
      https://www.ebi.ac.uk/ena/browser/view/PRJEB9586
    Individual FTP path pattern per sample under
      ftp://ftp.sra.ebi.ac.uk/vol1/fastq/... (per-sample; script only writes
      the SAMPLE-ID list; you handle FTP/download separately or point to a
      local reprocessed release).

Output (per-super-pop CANDIDATE lists; one sample ID per line; each is a
pre-homogeneity-filter pool, NOT the final homogeneous set):
  reference_ids/sgdp_sas_anchors.txt        SAS (South Asia)
  reference_ids/sgdp_oceania_anchors.txt    OCE (Papuan / Melanesian / Australian)
  reference_ids/sgdp_men_anchors.txt        MEN (Middle Eastern: Bedouin, Iranian, ...)
  reference_ids/sgdp_amr_anchors.txt        AMR (Karitiana, Surui, Mayan, Pima, Mixe,
                                            Zapotec, Piapoco, Quechua, Chane)
  reference_ids/sgdp_afr_anchors.txt        AFR (sub-Saharan African)
  reference_ids/sgdp_eur_anchors.txt        EUR (WestEurasia minus MEN populations)
  reference_ids/sgdp_eas_anchors.txt        EAS (EastAsia + CentralAsiaSiberia)
  reference_ids/sgdp_combined_anchors.txt   union of all 7 candidates — fed into
                                            1c_prep_sgdp_clm.sh as the download
                                            keep-list (who to pull from SGDP);
                                            final inclusion decided later by 4a
  reference_ids/sgdp_manifest.tsv           full filtered metadata rows for all
                                            seven groups above -- provenance for
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
# "CentralAsiaSiberia" is bucketed into EAS (closer to EAS than SAS/EUR in PCA).
SGDP_TO_SUPERPOP = {
    "SouthAsia":           "SAS",
    "Oceania":             "OCE",
    "Africa":              "AFR",
    "EastAsia":            "EAS",
    "CentralAsiaSiberia":  "EAS",   # Siberian / Central Asian cluster w/ EAS
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
    """Yield sample dict rows from SGDP metadata TSV, respecting tier filters.

    SGDP header line starts with '#Sequencing_Panel' -- treat it as the header,
    not a comment. Tier is the 'Embargo' column with values FullyPublic,
    SignedLetterNoDelay, SignedLetterDelay (or DO_NOT_USE)."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        if not lines:
            return
        # Strip leading '#' from the header row so DictReader keys match.
        if lines[0].startswith("#"):
            lines[0] = lines[0][1:]
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            # Real SGDP columns: Sample_ID, Population_ID, Region, Gender, Embargo, ...
            sample = (row.get("Sample_ID") or row.get("SGDP_ID")
                      or row.get("SampleID") or row.get("Illumina_ID"))
            pop    = row.get("Population_ID") or row.get("Population")
            region = row.get("Region")
            embargo = (row.get("Embargo") or "").strip()
            if not sample or not pop or not region:
                continue
            if embargo == "DO_NOT_USE":
                continue
            if embargo == "SignedLetterNoDelay" and not include_signed_letter:
                continue
            if embargo == "SignedLetterDelay" and not (include_signed_letter and include_fan):
                continue
            yield {
                "sample": sample.strip(),
                "population": pop.strip(),
                "region": region.strip(),
                "sequencing_source": embargo,
                "country": (row.get("Country") or "").strip(),
                "sex": (row.get("Gender") or row.get("Sex") or "").strip(),
                "latitude": (row.get("Latitude") or "").strip(),
                "longitude": (row.get("Longitude") or "").strip(),
            }


def bucket(rows, cap_per_pop=None):
    """Assign each row to SAS/OCE/MEN/AMR/AFR/EUR/EAS or drop.

    WestEurasia samples are split:
      - populations in MEN_POPULATIONS -> MEN
      - all other WestEurasia -> EUR (European)
    America samples are kept only if population is in AMR_POPULATIONS
    (homogeneous NAT per Mallick 2016; drops admixed Mex-Am etc.).
    """
    per_pop_seen = Counter()
    buckets = defaultdict(list)
    for r in rows:
        region = r["region"]; pop = r["population"]
        super_pop = SGDP_TO_SUPERPOP.get(region)
        if super_pop is None:
            if region == "WestEurasia":
                super_pop = "MEN" if pop in MEN_POPULATIONS else "EUR"
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
                             ("AMR", "sgdp_amr_anchors.txt"),
                             ("AFR", "sgdp_afr_anchors.txt"),
                             ("EUR", "sgdp_eur_anchors.txt"),
                             ("EAS", "sgdp_eas_anchors.txt")]:
        out_path = os.path.join(args.outdir, out_name)
        samples = sorted({r["sample"] for r in buckets.get(sp_key, [])})
        with open(out_path, "w") as g:
            g.write("\n".join(samples) + ("\n" if samples else ""))
        pops = Counter(r["population"] for r in buckets.get(sp_key, []))
        pop_summary = ", ".join(f"{p}={n}" for p, n in pops.most_common())
        print(f"[fetch] {sp_key:>3}: wrote {len(samples):3d} candidates "
              f"-> {out_path}  ({pop_summary or 'no matches'})", file=sys.stderr)

    # Combined manifest for provenance.
    manifest_path = os.path.join(args.outdir, "sgdp_manifest.tsv")
    fields = ["super_pop", "sample", "population", "region", "country",
              "sex", "latitude", "longitude", "sequencing_source"]
    all_superpops = ("SAS", "OCE", "MEN", "AMR", "AFR", "EUR", "EAS")
    with open(manifest_path, "w") as g:
        g.write("\t".join(fields) + "\n")
        for sp_key in all_superpops:
            for r in sorted(buckets.get(sp_key, []),
                            key=lambda x: (x["population"], x["sample"])):
                g.write("\t".join(str(r.get(k, "")) for k in fields) + "\n")
    print(f"[fetch] wrote manifest -> {manifest_path}", file=sys.stderr)

    # Combined candidate list — union across all 7 super-pops, deduped + sorted.
    # 1c_prep_sgdp_clm.sh reads this as the keep-list for sample-filtering the
    # downloaded SGDP VCFs. NOT the final homogeneous set: 4a's Q>=0.95 filter
    # decides which of these actually make the published reference panel.
    combined_path = os.path.join(args.outdir, "sgdp_combined_anchors.txt")
    all_samples = sorted({r["sample"]
                          for sp_key in all_superpops
                          for r in buckets.get(sp_key, [])})
    with open(combined_path, "w") as g:
        g.write("\n".join(all_samples) + ("\n" if all_samples else ""))
    print(f"[fetch] wrote combined candidate list ({len(all_samples)} samples) "
          f"-> {combined_path}", file=sys.stderr)
    print(f"[fetch] NOTE: these are population-majority CANDIDATES, not a "
          f"verified homogeneous set — 4a (ADMIXTURE Q>=0.95) decides final "
          f"inclusion downstream.", file=sys.stderr)


if __name__ == "__main__":
    main()
