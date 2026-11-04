#!/usr/bin/env bash
# ============================================================================
# 5a_build_donor_lists.sh  (one-shot, no SLURM needed; seconds)
#
# Build the DONOR table for the simulated benchmarking truth set.
#
# Rule: donors are HGDP and MX Biobank samples ONLY. 1KG samples are never
# donors because they sit inside the imputation reference panels the
# benchmark compares (Michigan 1000G, HRC [which contains 1KG Phase 3],
# TOPMed+1KG). HGDP and MXB are in none of the panels, so genomes spliced
# from them are non-circular for every arm.
#
# Source detection is by sample-ID prefix (robust across metadata versions):
#   HGDP*  -> HGDP        MXB_*  -> MXB        anything else -> 1KG (excluded)
#
# Optional filters:
#   HOMOG_TSV   lai_ref_panel_samples.tsv from 4a (Q_unsup max, included flag).
#               When present, homogeneous-cohort donors are restricted to
#               included==1 (max_Q >= HOMOG_THRESHOLD). Admixed cohorts use the
#               same donors, so the ancestral "sources" are clean.
#   EXCLUDE_RFMIX_REFS=1 (default) drops every sample that appears in any
#               reference_ids/*_rfmix.txt list, so the genotyping pipeline's
#               Module 7 (RFMix with those references) cannot copy its own
#               donors. Set to 0 to widen the pool; the report shows the cost.
#
# Inputs
#   sample_groups.tsv          sample -> superpop (make_sample_groups.sh)
#   gnomad_meta_updated.tsv    for population labels (hgdp_tgp_meta.Population)
#   reference_ids/MXB50genomes_popinfo.tsv  MXB IDs + inferred cluster
#   (optional) $HOMOG_TSV
#
# Outputs ($OUT_DIR, default $PROJECT_ROOT/06_benchmark_truth/donors/)
#   donors.tsv                 sample_id  population  superpop  source  q_max  included  rfmix_ref
#   donors.<POPULATION>.txt    keep-lists (one per population, included==1 only)
#   donor_counts.tsv           per population / superpop counts (read this before 5b)
# ============================================================================
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SAMPLE_GROUPS="${SAMPLE_GROUPS:-${PROJECT_ROOT}/sample_groups.tsv}"
GNOMAD_META="${GNOMAD_META:-/storage/atkinson/shared_resources/reference/ReferencePanels/TGP_HGDP_jointcall/processed_data/phased_haplotypes_v2_filter1/gnomad_meta_updated.tsv}"
MXB_POPINFO="${MXB_POPINFO:-${REPO_DIR}/reference_ids/MXB50genomes_popinfo.tsv}"
HOMOG_TSV="${HOMOG_TSV:-${PROJECT_ROOT}/04_homogeneity_panel/lai_ref_panel_samples.tsv}"
HOMOG_THRESHOLD="${HOMOG_THRESHOLD:-0.93}"
EXCLUDE_RFMIX_REFS="${EXCLUDE_RFMIX_REFS:-1}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/06_benchmark_truth/donors}"

[[ -s "$SAMPLE_GROUPS" ]] || { echo "ERROR: missing $SAMPLE_GROUPS (run init-project-tree.sh)"; exit 1; }
[[ -s "$GNOMAD_META"   ]] || { echo "ERROR: missing $GNOMAD_META"; exit 1; }
[[ -s "$MXB_POPINFO"   ]] || { echo "ERROR: missing $MXB_POPINFO"; exit 1; }
mkdir -p "$OUT_DIR"

# RFMix reference IDs (union of every *_rfmix.txt in reference_ids/)
RFMIX_IDS="${OUT_DIR}/rfmix_reference_ids.txt"
cat "${REFS}"/*_rfmix.txt 2>/dev/null | awk 'NF{print $1}' | sort -u > "$RFMIX_IDS" || : > "$RFMIX_IDS"

python3 - "$SAMPLE_GROUPS" "$GNOMAD_META" "$MXB_POPINFO" "$HOMOG_TSV" "$HOMOG_THRESHOLD" \
          "$EXCLUDE_RFMIX_REFS" "$RFMIX_IDS" "$OUT_DIR" <<'PYEOF'
import csv, os, sys
from collections import Counter, defaultdict
sg, meta, mxb, homog, thr, excl_ref, rfmix_f, out = sys.argv[1:9]
thr = float(thr); excl_ref = excl_ref == "1"

# sample -> superpop
superpop = {}
for line in open(sg):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2: superpop[p[0]] = p[1]

# sample -> population (HGDP/1KG from gnomAD meta; MXB from popinfo cluster)
population = {}
with open(meta) as f:
    r = csv.DictReader(f, delimiter="\t")
    pop_col = "hgdp_tgp_meta.Population" if "hgdp_tgp_meta.Population" in r.fieldnames else None
    for row in r:
        sid = row.get("project_meta.sample_id") or row.get("s")
        if sid: population[sid] = (row.get(pop_col) if pop_col else "NA") or "NA"
with open(mxb) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        population[p[0]] = "MXB_" + (p[3] if len(p) > 3 and p[3] else "unknown")
        superpop.setdefault(p[0], "AMR")

# homogeneity Q (optional)
qmax, included = {}, {}
if os.path.exists(homog):
    with open(homog) as f:
        r = csv.DictReader(f, delimiter="\t")
        qcols = [c for c in r.fieldnames if c.startswith("Q_unsup_")]
        for row in r:
            sid = row.get("sample_id") or row.get("sample") or row.get("IID")
            qs = [float(row[c]) for c in qcols if row.get(c) not in (None, "", "NA")]
            qmax[sid] = max(qs) if qs else float("nan")
            inc = row.get("included")
            included[sid] = (inc == "1") if inc is not None else (qmax[sid] >= thr)
rfmix = set(l.strip() for l in open(rfmix_f) if l.strip())

def source(sid):
    if sid.startswith("HGDP"): return "HGDP"
    if sid.startswith("MXB"):  return "MXB"
    return "1KG"

rows = []
for sid, sp in superpop.items():
    src = source(sid)
    if src == "1KG":  # never a donor
        continue
    is_ref = sid in rfmix
    inc = included.get(sid, True)  # no homogeneity table -> keep
    if excl_ref and is_ref: inc = False
    rows.append((sid, population.get(sid, "NA"), sp, src, qmax.get(sid, float("nan")), int(inc), int(is_ref)))
rows.sort(key=lambda r: (r[2], r[1], r[0]))

with open(os.path.join(out, "donors.tsv"), "w") as g:
    g.write("sample_id\tpopulation\tsuperpop\tsource\tq_max\tincluded\trfmix_ref\n")
    for r in rows:
        g.write("\t".join(str(x) if not isinstance(x, float) else ("NA" if x != x else f"{x:.3f}") for x in r) + "\n")

bypop = defaultdict(list)
for r in rows:
    if r[5] == 1: bypop[(r[2], r[1])].append(r[0])
for (sp, pop), ids in bypop.items():
    safe = pop.replace(" ", "_").replace("/", "_")
    with open(os.path.join(out, f"donors.{sp}.{safe}.txt"), "w") as g:
        g.write("\n".join(ids) + "\n")

cnt_all = Counter((r[2], r[1]) for r in rows)
cnt_inc = Counter((r[2], r[1]) for r in rows if r[5] == 1)
cnt_ref = Counter((r[2], r[1]) for r in rows if r[6] == 1)
with open(os.path.join(out, "donor_counts.tsv"), "w") as g:
    g.write("superpop\tpopulation\tn_all\tn_rfmix_ref\tn_included\n")
    for k in sorted(cnt_all):
        g.write(f"{k[0]}\t{k[1]}\t{cnt_all[k]}\t{cnt_ref[k]}\t{cnt_inc[k]}\n")
    g.write("\n# superpop totals (included)\n")
    tot = Counter()
    for k, v in cnt_inc.items(): tot[k[0]] += v
    for sp in ("AFR", "EUR", "EAS", "CSA", "SAS", "AMR", "OCE", "MEN"):
        g.write(f"# {sp}\t{tot.get(sp, 0)}\n")
print(open(os.path.join(out, "donor_counts.tsv")).read())
PYEOF

echo "[5a] donors written to $OUT_DIR"
echo "[5a] READ donor_counts.tsv before running 5b: superpop totals cap the homogeneous cohort sizes"
echo "     (OCE and AMR are small - cohorts.tsv already caps them at 50)."
