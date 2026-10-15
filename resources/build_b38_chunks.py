#!/usr/bin/env python3
# ============================================================================
# build_b38_chunks.py
#
# Emit per-(chr, chunk) phasing regions for SHAPEIT5_phase_common, computed
# from the SHAPEIT4-format b38 genetic maps. Canonical SHAPEIT5 recipe:
#   - 8 cM phasing windows
#   - 2 cM overlap between adjacent windows
# Output (TSV, header'd):
#   chr_num  chr   chunk_idx  phase_region          ligate_region
#   1        chr1  0          chr1:13380-3500000    chr1:13380-3050000
#   1        chr1  1          chr1:3050000-7800000  chr1:3500000-7300000
#   ...
# `phase_region` is the bigger window passed to phase_common's --region.
# `ligate_region` is informational (the non-overlap "core"); SHAPEIT5_ligate
# does the actual stitching from BCF coordinates.
#
# Usage:
#   build_b38_chunks.py <gmap_dir> <out.tsv> [--phase-cm 8] [--overlap-cm 2]
# ============================================================================

import argparse
import bisect
import gzip
import sys
from pathlib import Path


def load_gmap(path):
    pos, cm = [], []
    with gzip.open(path, "rt") as f:
        next(f)  # header
        for line in f:
            parts = line.split()
            # SHAPEIT4 gmap format: pos chr cM
            pos.append(int(parts[0]))
            cm.append(float(parts[-1]))
    return pos, cm


def cm_to_bp(target_cm, pos, cm):
    i = bisect.bisect_left(cm, target_cm)
    if i == 0:
        return pos[0]
    if i >= len(cm):
        return pos[-1]
    frac = (target_cm - cm[i - 1]) / (cm[i] - cm[i - 1] + 1e-9)
    return int(pos[i - 1] + frac * (pos[i] - pos[i - 1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gmap_dir")
    ap.add_argument("out_tsv")
    ap.add_argument("--phase-cm", type=float, default=8.0)
    ap.add_argument("--overlap-cm", type=float, default=2.0)
    args = ap.parse_args()

    gmap_dir = Path(args.gmap_dir)
    with open(args.out_tsv, "w") as out:
        out.write("chr_num\tchr\tchunk_idx\tphase_region\tligate_region\n")
        for chrN in range(1, 23):
            gmap = gmap_dir / f"chr{chrN}.b38.gmap.gz"
            if not gmap.is_file():
                sys.exit(f"missing {gmap}")
            pos, cm = load_gmap(gmap)
            chr_str = f"chr{chrN}"
            cm_min, cm_max = cm[0], cm[-1]
            start_cm = cm_min
            chunk_idx = 0
            while start_cm < cm_max - 1e-3:
                end_cm = min(start_cm + args.phase_cm, cm_max)
                phase_s = cm_to_bp(start_cm, pos, cm)
                phase_e = cm_to_bp(end_cm, pos, cm)
                lig_s_cm = start_cm if chunk_idx == 0 else start_cm + args.overlap_cm / 2
                lig_e_cm = end_cm if end_cm >= cm_max - 1e-3 else end_cm - args.overlap_cm / 2
                lig_s = cm_to_bp(lig_s_cm, pos, cm)
                lig_e = cm_to_bp(lig_e_cm, pos, cm)
                out.write(
                    f"{chrN}\t{chr_str}\t{chunk_idx}\t"
                    f"{chr_str}:{phase_s}-{phase_e}\t"
                    f"{chr_str}:{lig_s}-{lig_e}\n"
                )
                chunk_idx += 1
                if end_cm >= cm_max - 1e-3:
                    break
                start_cm = end_cm - args.overlap_cm
            print(f"chr{chrN}: {chunk_idx} chunks", file=sys.stderr)


if __name__ == "__main__":
    main()
