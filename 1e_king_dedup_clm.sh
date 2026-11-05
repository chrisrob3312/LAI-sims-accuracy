#!/usr/bin/env bash
#SBATCH --job-name=1e_king_dedup
#SBATCH --partition=atkinson,mhgcp
#SBATCH --exclude=mhgcp-c02,mhgcp-t01,mhgcp-t02,mhgcp-t03,mhgcp-t04,mhgcp-t05,mhgcp-t06,mhgcp-t07,mhgcp-t08,mhgcp-t09,mhgcp-t10,mhgcp-t11,mhgcp-t12
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1e_king_dedup_%j.out
#SBATCH --error=/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank/logs/1e_king_dedup_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Christina.magyar@bcm.edu
# ============================================================================
# 1e_king_dedup_clm.sh   --   runs ONCE, after 1c, before 4a
#
# Find the same person appearing twice in the merged panel under two different
# IDs, and write the drop list that keeps one copy.
#
# WHY THIS CANNOT BE LEFT TO THE GENOTYPING PIPELINE'S KING. That one removes
# duplicates and MZ twins from a COHORT. Here the two copies are not both
# cohort members: one can be a simulation DONOR and the other a reference
# sample or a real test genome. Nothing in a cohort-level duplicate filter
# looks across those roles, so the pair survives, and the method is then
# scored on a test sample whose haplotypes it was handed as reference. That is
# the circularity every other check in this project exists to prevent, and it
# is invisible in the output: it looks like unusually good accuracy.
#
# WHY IT RUNS BEFORE 4a, NOT AFTER.
#   1. A duplicated individual is counted twice by ADMIXTURE, so its
#      population pulls the cluster centroid twice as hard. The Q matrix
#      everything downstream reads is therefore wrong before the split is even
#      considered.
#   2. The homogeneity cut and the donor/reference split both read that Q. A
#      duplicate removed afterwards would mean redoing both.
#
# WHY IT IS NOT 1c's JOB. 1c deduplicates by SAMPLE ID at the merge -- it
# catches the same identifier in two cohorts, which is what `bcftools merge`
# would abort on. It cannot catch one person carrying two different IDs across
# HGDP, 1KG, MXB, SGDP and Almarri, which is the common case: SGDP re-sequenced
# HGDP individuals, and the same cell line can sit in two releases under two
# names. Only genotypes catch that.
#
# OUTPUT
#   king/king.kin0                     every pair above the reporting floor
#   king/duplicate_pairs.tsv           pairs above DUP_THRESHOLD, annotated
#   king/related_pairs_summary.tsv     counts by degree, for the methods
#   reference_ids/panel_duplicate_drop_ids.txt   the IDs to drop (4a reads it)
#
# Usage:
#   sbatch 1e_king_dedup_clm.sh
#   DUP_THRESHOLD=0.354 sbatch 1e_king_dedup_clm.sh
# ============================================================================
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-shapeit5}"
PROJECT_ROOT="${PROJECT_ROOT:-/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
REFS="${REFS:-${REPO_DIR}/reference_ids}"
PHASED_DIR="${PHASED_DIR:-${PROJECT_ROOT}/01_merged_phased_panel}"
PANEL_TPL="${PANEL_TPL:-${PHASED_DIR}/merged_chr%s.shapeit5_phased.softunion_maf005.bcf}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/04_homogeneity_panel/king}"

# KING-robust kinship coefficients. These are the standard cut-points:
#   > 0.354   duplicate or monozygotic twin   <- what we drop
#   > 0.177   first degree  (parent-child, full sib)
#   > 0.0884  second degree
#   > 0.0442  third degree
# Only the duplicate tier is dropped. Relatives are REPORTED, not removed:
# a panel of unrelated-only individuals is a different panel, and relatedness
# among real population samples is a property of the data the methods section
# should state rather than something to quietly delete.
DUP_THRESHOLD="${DUP_THRESHOLD:-0.354}"
REPORT_FLOOR="${REPORT_FLOOR:-0.0442}"

# LD pruning: KING wants common, roughly independent sites.
MAF="${MAF:-0.05}"
PRUNE_WINDOW="${PRUNE_WINDOW:-200}"
PRUNE_STEP="${PRUNE_STEP:-50}"
PRUNE_R2="${PRUNE_R2:-0.2}"

# Which copy of a duplicate pair to KEEP, highest priority first.
#
# Non-1000G sources are preferred because a 1000 Genomes sample can never be a
# simulation donor (the imputation servers' panels contain 1000 Genomes), so
# keeping the non-1000G copy preserves the sample's eligibility for both
# roles. MXB and Almarri come first because they are the scarcest and the
# hardest to replace. Ties fall back to lexicographic sample ID, so the
# outcome is deterministic without a seed.
SOURCE_PRIORITY="${SOURCE_PRIORITY:-MXB ALMARRI SGDP HGDP 1KG}"

DROP_LIST="${DROP_LIST:-${REFS}/panel_duplicate_drop_ids.txt}"
THREADS=${SLURM_CPUS_PER_TASK:-16}

module load anaconda3/2024.06
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$CONDA_ENV"
set -u
unset PYTHONHOME PYTHONPATH
export PATH="${CONDA_PREFIX}/bin:${PATH}"

mkdir -p "$OUTDIR"
command -v plink2 >/dev/null || { echo "ERROR: plink2 not on PATH"; exit 1; }

for c in $(seq 1 22); do
    f=$(printf "$PANEL_TPL" "$c")
    [[ -s "$f" ]] || { echo "ERROR: missing $f -- run 1c (stage 2) first"; exit 1; }
done

echo "[$(date +%T)] KING de-duplication over the merged phased panel"
echo "[$(date +%T)]   duplicate threshold : ${DUP_THRESHOLD}"
echo "[$(date +%T)]   reporting floor     : ${REPORT_FLOOR}"
echo "[$(date +%T)]   keep priority       : ${SOURCE_PRIORITY}"

# --------------------------------------------------------------------------
# 1. Concatenate and import. KING is a genome-wide statistic: run per
#    chromosome it would estimate kinship from a twenty-second of the genome
#    and the duplicate tier would not separate cleanly from first degree.
# --------------------------------------------------------------------------
CONCAT="${OUTDIR}/all_chr.bcf"
if [[ ! -s "${CONCAT}.csi" ]]; then
    echo "[$(date +%T)] concatenating 22 chromosomes"
    LIST="${OUTDIR}/concat_list.txt"; : > "$LIST"
    for c in $(seq 1 22); do printf "$PANEL_TPL\n" "$c" >> "$LIST"; done
    bcftools concat -f "$LIST" --threads "$THREADS" -Ob -o "${CONCAT}.partial"
    bcftools index --threads "$THREADS" "${CONCAT}.partial"
    mv "${CONCAT}.partial.csi" "${CONCAT}.csi"
    mv "${CONCAT}.partial" "$CONCAT"
fi
N_PANEL=$(bcftools query -l "$CONCAT" | wc -l)
echo "[$(date +%T)] panel: ${N_PANEL} samples"

# --------------------------------------------------------------------------
# 2. LD-prune on common sites.
# --------------------------------------------------------------------------
PRUNE="${OUTDIR}/pruned"
if [[ ! -s "${PRUNE}.prune.in" ]]; then
    echo "[$(date +%T)] plink2 import + LD prune (MAF>=${MAF}, r2<${PRUNE_R2})"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' --rm-dup exclude-all \
           --max-alleles 2 --snps-only just-acgt \
           --maf "$MAF" \
           --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_R2" \
           --out "$PRUNE"
fi
N_PRUNE=$(wc -l < "${PRUNE}.prune.in")
echo "[$(date +%T)] pruned to ${N_PRUNE} variants"
# KING-robust is unreliable on too few markers; the duplicate tier blurs into
# first degree and the drop list becomes a guess.
[[ "$N_PRUNE" -ge 20000 ]] || {
    echo "ERROR: only ${N_PRUNE} pruned variants -- too few for reliable KING." >&2
    echo "       Lower MAF (currently ${MAF}) or relax PRUNE_R2 (${PRUNE_R2})." >&2
    exit 2; }

# --------------------------------------------------------------------------
# 3. KING-robust kinship.
#    --king-table-filter keeps the output small; it is set to the REPORTING
#    floor, not the duplicate threshold, so the degree summary below can be
#    built from the same file.
# --------------------------------------------------------------------------
KING="${OUTDIR}/king"
if [[ ! -s "${KING}.kin0" ]]; then
    echo "[$(date +%T)] plink2 --make-king-table"
    plink2 --bcf "$CONCAT" --threads "$THREADS" \
           --set-missing-var-ids '@:#[b38]' --rm-dup exclude-all \
           --max-alleles 2 --snps-only just-acgt \
           --extract "${PRUNE}.prune.in" \
           --make-king-table --king-table-filter "$REPORT_FLOOR" \
           --out "$KING"
fi
[[ -s "${KING}.kin0" ]] || { echo "ERROR: plink2 produced no ${KING}.kin0"; exit 1; }

# --------------------------------------------------------------------------
# 4. Classify pairs, choose which copy of each duplicate to keep, write the
#    drop list. Columns are read BY NAME: plink2 emits #FID1/IID1/.../KINSHIP
#    or IID1/IID2/... depending on whether the input had family IDs, and
#    taking them positionally silently reads the wrong column.
# --------------------------------------------------------------------------
echo "[$(date +%T)] classifying pairs and choosing copies to keep"
python3 -I - "${KING}.kin0" "$DUP_THRESHOLD" "$REPORT_FLOOR" \
    "${OUTDIR}/duplicate_pairs.tsv" "${OUTDIR}/related_pairs_summary.tsv" \
    "$DROP_LIST" "$SOURCE_PRIORITY" "${REFS}/MXB50genomes_popinfo.tsv" \
    "${REFS}/sgdp_manifest.tsv" <<'PYEOF'
import sys
from collections import defaultdict

(kin_f, dup_t, rep_t, dup_out, summ_out, drop_out, prio_s,
 mxb_f, sgdp_f) = sys.argv[1:10]
dup_t = float(dup_t); rep_t = float(rep_t)
PRIO = prio_s.split()
rank = {s: i for i, s in enumerate(PRIO)}

def load_ids(path, col=0, skip_header=True):
    out = set()
    try:
        with open(path) as fh:
            if skip_header: next(fh, None)
            for line in fh:
                ff = line.rstrip("\n").split("\t")
                if ff and ff[col].strip(): out.add(ff[col].strip())
    except FileNotFoundError:
        pass
    return out

mxb = load_ids(mxb_f)
sgdp = load_ids(sgdp_f, col=1)

def source_of(s):
    if s in mxb: return "MXB"
    if s in sgdp or s.endswith(".DG"): return "SGDP"
    if s.startswith("APPG"): return "ALMARRI"
    if s.startswith("HGDP"): return "HGDP"
    return "1KG"

def prio_of(s):
    return rank.get(source_of(s), len(PRIO))

# --- read kin0 by column NAME ---
with open(kin_f) as fh:
    header = fh.readline().lstrip("#").rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    def col(*names):
        for n in names:
            if n in idx: return idx[n]
        sys.exit(f"ERROR: {kin_f} has none of {names}; header={header}")
    i1, i2 = col("IID1", "ID1"), col("IID2", "ID2")
    ik = col("KINSHIP")
    pairs = []
    for line in fh:
        ff = line.rstrip("\n").split("\t")
        if len(ff) <= max(i1, i2, ik): continue
        try: k = float(ff[ik])
        except ValueError: continue
        pairs.append((ff[i1], ff[i2], k))

print(f"[king] {len(pairs)} pairs at or above kinship {rep_t}", file=sys.stderr)

TIERS = [(0.354, "duplicate_or_MZ"), (0.177, "first_degree"),
         (0.0884, "second_degree"), (0.0442, "third_degree")]
def tier(k):
    for t, name in TIERS:
        if k > t: return name
    return "below_third"

counts = defaultdict(int)
for _a, _b, k in pairs: counts[tier(k)] += 1
with open(summ_out, "w") as g:
    g.write("degree\tn_pairs\n")
    for _t, name in TIERS: g.write(f"{name}\t{counts[name]}\n")
    g.write(f"below_third\t{counts['below_third']}\n")
for _t, name in TIERS:
    print(f"[king]   {name:18s} {counts[name]}", file=sys.stderr)

# --- duplicates: choose a copy to drop, greedily, highest kinship first ---
dups = sorted([p for p in pairs if p[2] > dup_t], key=lambda p: -p[2])
dropped = set()
rows = []
for a, b, k in dups:
    # If either is already dropped the pair is resolved; recording it anyway
    # keeps the table a complete account of what was found.
    if a in dropped or b in dropped:
        rows.append((a, b, k, source_of(a), source_of(b), "already_resolved", ""))
        continue
    # Lower priority number wins. Lexicographic tie-break keeps it
    # deterministic without a seed.
    if (prio_of(a), a) <= (prio_of(b), b): keep, drop = a, b
    else: keep, drop = b, a
    dropped.add(drop)
    rows.append((a, b, k, source_of(a), source_of(b), f"keep:{keep}", drop))

with open(dup_out, "w") as g:
    g.write("id1\tid2\tkinship\tsource1\tsource2\tdecision\tdropped\n")
    for r in rows:
        g.write("\t".join(str(x) for x in r) + "\n")
with open(drop_out, "w") as g:
    for s in sorted(dropped): g.write(s + "\n")

print(f"[king] duplicate/MZ pairs: {len(dups)}", file=sys.stderr)
print(f"[king] samples to drop   : {len(dropped)} -> {drop_out}", file=sys.stderr)
for r in rows[:20]:
    print(f"[king]   {r[0]} ({r[3]}) ~ {r[1]} ({r[4]})  kinship={r[2]:.4f}  {r[5]}",
          file=sys.stderr)
if len(rows) > 20:
    print(f"[king]   ... {len(rows)-20} more in {dup_out}", file=sys.stderr)
PYEOF

N_DROP=$(awk 'NF' "$DROP_LIST" 2>/dev/null | wc -l)
echo
echo "[$(date +%T)] DONE."
echo "  pair table     : ${OUTDIR}/duplicate_pairs.tsv"
echo "  degree summary : ${OUTDIR}/related_pairs_summary.tsv"
echo "  drop list      : ${DROP_LIST} (${N_DROP} samples)"
if [[ "$N_DROP" -gt 0 ]]; then
cat <<EOF

  ${N_DROP} duplicate(s) found. 4a excludes anything in the drop list when it
  builds the ADMIXTURE panel, so stage 3 onwards never sees them. Read
  duplicate_pairs.tsv before stage 3: each row names both IDs, both sources
  and which copy was kept, and the choice is a judgement (non-1000G first,
  then scarcest source, then alphabetical) rather than a fact.
EOF
else
cat <<EOF

  No duplicate pairs above ${DUP_THRESHOLD}. Relatives below that tier are
  reported in related_pairs_summary.tsv and are NOT removed: a panel of
  unrelated-only individuals is a different panel, and relatedness among real
  population samples belongs in the methods, not in a silent filter.
EOF
fi
