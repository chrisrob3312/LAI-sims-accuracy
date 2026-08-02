#!/usr/bin/env python3
"""
Render 3 committee-slide figures as PNG:
  fig_wgs_panels.png    grouped bars: NAT recall per panel per sim track (WGS)
  fig_chip_vs_wgs.png   grouped bars: chip GSA vs WGS at the winning panel
  fig_amr_disparity.png bar + within-AMR pie: reference-panel gap

Reads results/accuracy_pilot/accuracy_summary.tsv
and results/accuracy_pilot_chip_gsa/accuracy_summary.tsv.

Numbers are chr20-22 pilot (as committed on branch).
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "slides" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- data readers
def read_summary(tsv):
    """Return {(track, panel, ancestry): weighted_concordance}."""
    d = {}
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = (row["track"], row["panel"], row["ancestry"])
            d[k] = float(row["weighted_concordance"])
    return d

wgs  = read_summary(REPO / "results/accuracy_pilot/accuracy_summary.tsv")
chip = read_summary(REPO / "results/accuracy_pilot_chip_gsa/accuracy_summary.tsv")

# ---------------------------------------------------------------- colors
# Palette matches the Spring TAC pptx theme accent scheme (theme3.xml):
#   accent1 #156082 (teal)   accent2 #E97132 (orange)  accent3 #196B24 (forest)
#   accent4 #0F9ED5 (cyan)   accent5 #A02B93 (plum)    accent6 #4EA72E (lime)
#   dk2     #0E2841 (navy)
NAVY   = "#0E2841"
TEAL   = "#156082"
ORANGE = "#E97132"
FOREST = "#196B24"
CYAN   = "#0F9ED5"
PLUM   = "#A02B93"
LIME   = "#4EA72E"
TEAL_LT= "#91D9CE"

# track colors used in every bar chart
BRASA_COL = ORANGE   # Brazilian-like cohort
MXB_COL   = TEAL     # Mexican-like cohort

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

# ---------------------------------------------------------------- FIG 1: WGS
def fig_wgs_panels():
    # (friendly label, panel key)
    panels = [
        ("HGDP + 1KG",     "NAT_HGDP"),
        ("+ 25 MXB",       "NAT_HGDPMXB"),
        ("+ 50 MXB",       "NAT_HGDPMXB_FULL"),
        ("Homog. Q≥0.95",  "NAT_HOMOG"),
        ("+ PEL (Peru)",   "NAT_PEL"),
        ("+ PEL + EAS",    "NAT_PEL_EAS"),
    ]
    tracks = [("NAT",    "Brazilian-like (Brasa) admixed",  BRASA_COL),
              ("NATMXB", "Mexican-like (MXB) admixed",      MXB_COL)]

    fig, ax = plt.subplots(figsize=(15, 7.2), dpi=180)
    # roomier left margin (y-label breathing room) + roomier bottom (x-labels + legend gap)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.88, bottom=0.36)

    x = np.arange(len(panels))
    w = 0.36
    for i, (track, label, color) in enumerate(tracks):
        vals = [wgs.get((track, key, "NAT"), np.nan) for _, key in panels]
        bars = ax.bar(x + (i - 0.5) * w, vals, w, label=label,
                      color=color, edgecolor="white", linewidth=0.7)
        for xi, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(xi.get_x() + xi.get_width()/2, v + 0.004,
                        f"{v*100:.0f}%", ha="center", va="bottom",
                        fontsize=11, color="#111", fontweight="bold")

    # single-line, concise tick labels; bold
    ax.set_xticks(x)
    ax.set_xticklabels([lab for lab, _ in panels], fontsize=11.5, fontweight="bold")
    for lbl in ax.get_xticklabels():
        lbl.set_multialignment("center")
    ax.tick_params(axis="x", pad=8)
    ax.tick_params(axis="y", pad=6, labelsize=11)

    ax.set_ylim(0.85, 1.00)
    ax.set_ylabel("Weighted AMR recall\n(chr20-22, per-hap concordance)",
                  fontsize=13, fontweight="bold", labelpad=16)

    ax.set_title("Local Ancestry Inference Simulation Accuracy (AMR tracts)",
                 loc="center", fontweight="bold", fontsize=14, pad=14)

    # legend: below plot with clear gap from x-labels
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24),
                    ncol=2, fontsize=13, frameon=False)
    for t in leg.get_texts():
        t.set_fontweight("bold")

    out = OUT / "fig_wgs_panels.png"
    fig.savefig(out)
    plt.close(fig)
    return out

# ---------------------------------------------------------------- FIG 2: chip
def fig_chip_vs_wgs():
    panels = ["NAT_HGDP", "NAT_HGDPMXB", "NAT_HGDPMXB_FULL",
              "NAT_PEL", "NAT_PEL_EAS"]
    track = "NATMXB"      # Mexican-like cohort — the recommendation slide
    wgs_v  = [wgs.get( (track, p, "NAT"), np.nan) for p in panels]
    chip_v = [chip.get((track, p, "NAT"), np.nan) for p in panels]

    fig, ax = plt.subplots(figsize=(11.5, 6.0), dpi=180)
    x = np.arange(len(panels))
    w = 0.36
    b1 = ax.bar(x - w/2, wgs_v, w, label="WGS density",
                color=TEAL, edgecolor="white", linewidth=0.6)
    b2 = ax.bar(x + w/2, chip_v, w, label="Illumina GSA chip (~385k autosomal)",
                color=TEAL_LT, edgecolor="white", linewidth=0.6)
    for bars, vals in ((b1, wgs_v), (b2, chip_v)):
        for xi, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(xi.get_x() + xi.get_width()/2, v + 0.003,
                        f"{v:.3f}", ha="center", va="bottom",
                        fontsize=11, color="#111", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("NAT_", "").replace("_", "\n")
                        for p in panels], fontsize=11, fontweight="bold")
    for lbl in ax.get_xticklabels():
        lbl.set_multialignment("center")
    ax.set_ylim(0.85, 0.96)
    ax.set_ylabel("Weighted AMR recall",
                  fontsize=12, fontweight="bold", labelpad=12)
    ax.tick_params(axis="y", pad=6, labelsize=10.5)
    ax.tick_params(axis="x", pad=8)

    ax.set_title("LAI Accuracy: Genotyping (unimputed) vs. Whole-genome sequencing",
                 loc="center", fontweight="bold", fontsize=13.5, pad=14)

    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
                    ncol=2, fontsize=12, frameon=False)
    for t in leg.get_texts():
        t.set_fontweight("bold")

    fig.subplots_adjust(left=0.11, right=0.98, top=0.86, bottom=0.28)
    out = OUT / "fig_chip_vs_wgs.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out

# ---------------------------------------------------------------- FIG 3: gap
def fig_amr_disparity():
    fig = plt.figure(figsize=(13.5, 6.6), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], wspace=0.32,
                          bottom=0.28, top=0.90, left=0.09, right=0.98)

    # left: bars, samples per continental group in public references
    ax1 = fig.add_subplot(gs[0, 0])
    groups = ["EAS", "AFR", "EUR", "SAS", "AMR\nHGDP+1KG", "AMR\n+MXB"]
    counts = [ 667,   634,   620,   48,   38,             88]
    colors = [TEAL, ORANGE, FOREST, CYAN, PLUM, PLUM]
    alphas = [1,1,1,1,0.55,1.0]
    bars = ax1.bar(groups, counts, color=colors, alpha=None)
    for b, c, a in zip(bars, colors, alphas):
        b.set_alpha(a)
    for b, n in zip(bars, counts):
        ax1.text(b.get_x()+b.get_width()/2, n+8, f"{n}",
                 ha="center", va="bottom", fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Reference haplotypes\navailable in public HGDP+1KG panel",
                   fontweight="bold", fontsize=11.5, labelpad=16)
    ax1.tick_params(axis="y", pad=6, labelsize=10.5)
    ax1.set_ylim(0, 760)
    ax1.set_title("Available Homogeneous Sample Counts in 1KG-HGDP-MXB",
                  loc="center", fontweight="bold", fontsize=13, pad=12)
    ax1.tick_params(axis="x", labelsize=10.5, pad=6)
    for lbl in ax1.get_xticklabels():
        lbl.set_fontweight("bold")
        lbl.set_multialignment("center")

    # right: within-AMR pie showing geographic gaps
    ax2 = fig.add_subplot(gs[0, 1])
    parts = ["Mexico (MXB, +this study)",
             "Mexico (HGDP)",
             "Peru (PEL)",
             "Amazonia + Colombia (HGDP)"]
    sizes = [50, 13, 9, 16]
    cols  = [PLUM, "#c68abd", FOREST, ORANGE]
    wedges, txts, autotxt = ax2.pie(
        sizes, labels=None, colors=cols, autopct="%1.0f%%",
        startangle=90, pctdistance=0.72,
        wedgeprops=dict(edgecolor="white", linewidth=1.5),
        textprops=dict(fontsize=9.5, color="white", fontweight="bold"),
    )
    ax2.set_title("Within AMR: geographic composition after adding MXB",
                  loc="center", fontweight="bold", fontsize=12)
    leg = ax2.legend(wedges, [f"{p} (n={s})" for p, s in zip(parts, sizes)],
                     loc="upper center", bbox_to_anchor=(0.5, -0.02),
                     fontsize=11, frameon=False, ncol=1,
                     handletextpad=0.6, labelspacing=0.4)
    for t in leg.get_texts():
        t.set_fontweight("bold")

    fig.text(0.06, 0.03,
             "Gaps not addressed by MXB: Central America, Southern Cone, "
             "eastern Amazon, Caribbean",
             fontsize=10, color="#555", style="italic")
    out = OUT / "fig_amr_disparity.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out

# ---------------------------------------------------------------- FIG 4: pipeline diagram
def fig_pipeline():
    fig, ax = plt.subplots(figsize=(12, 4.2), dpi=180)
    ax.set_xlim(0, 12); ax.set_ylim(0, 4); ax.axis("off")

    modules = [
        ("Module 1\nPre-imp QC", "harmonize\nfixref\nliftover"),
        ("Module 2\nImputation", "TOPMed\nMichigan\nAoU AnVIL"),
        ("Module 3\nPost-imp QC", "MagicalRsq-X\ncall-rate"),
        ("Module 4\nMerge platforms", "intersection\nstrand-fix"),
        ("Module 5\nRe-imputation", "fill gaps\n(optional)"),
        ("Module 6\nFinal QC", "HWE / MAF\nrelatedness\nPCA"),
        ("Module 7\nAncestry", "GRAF-anc\nADMIXTURE\nRFMix / FLARE"),
        ("Module 8\nBenchmark", "truth vs imp\nsimulation"),
    ]
    # cool tones for QC / imputation modules, warm accent for LAI + benchmark
    _cool_a, _cool_b = "#d6e6ee", "#b6d1de"           # tints of TEAL
    _warm_a, _warm_b = "#f8dccb", "#f0c1a0"           # tints of ORANGE
    colors = [_cool_a, _cool_b, _cool_a, _cool_b,
              _cool_a, _cool_b, _warm_a, _warm_a]

    n = len(modules)
    w = 1.35; h = 2.0
    x0 = 0.15
    gap = (12 - x0 - n*w) / (n - 1)
    for i, ((head, body), c) in enumerate(zip(modules, colors)):
        x = x0 + i * (w + gap)
        y = 1.0
        rect = plt.Rectangle((x, y), w, h, facecolor=c, edgecolor=NAVY,
                             linewidth=0.8, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h - 0.35, head, ha="center", va="top",
                fontsize=9, fontweight="bold", zorder=3)
        ax.text(x + w/2, y + 0.85, body, ha="center", va="top",
                fontsize=7.8, color="#333", zorder=3)
        if i < n - 1:
            ax.annotate("", xy=(x + w + gap - 0.05, y + h/2),
                        xytext=(x + w + 0.02, y + h/2),
                        arrowprops=dict(arrowstyle="->", color="#666", lw=1.2))

    # highlight what this thesis validates
    hi_x = x0 + 6 * (w + gap) - 0.15
    hi_w = w * 2 + gap + 0.30
    ax.add_patch(plt.Rectangle((hi_x, 0.55), hi_w, h + 0.85,
                               facecolor="none", edgecolor=ORANGE,
                               linewidth=1.6, linestyle="--", zorder=1))
    ax.text(hi_x + hi_w/2, 0.30, "LAI accuracy validated in this thesis (Modules 7–8)",
            ha="center", fontsize=9, color=ORANGE, fontweight="bold")

    ax.text(0.05, 3.85, "8-module Nextflow preprocessing pipeline",
            fontsize=13, fontweight="bold")
    ax.text(0.05, 3.55, "containerized · parameterized · parallelized across chromosomes",
            fontsize=9.5, color="#555")

    out = OUT / "fig_pipeline.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out

# ---------------------------------------------------------------- FIG 5: 3-ancestry range (simplified)
def fig_wgs_all_ancestry():
    """
    Small multiples: AMR / EUR / AFR side-by-side, same y-scale.
    Only 3 panels — HGDP baseline, +25 MXB, and PEL — to keep it uncluttered.
    """
    panels = [
        ("HGDP + 1KG",   "NAT_HGDP"),
        ("+ MXB",        "NAT_HGDPMXB"),
        ("+ PEL (Peru)", "NAT_PEL"),
    ]
    tracks = [("NAT",    "Brazilian-like cohort", BRASA_COL),
              ("NATMXB", "Mexican-like cohort",   MXB_COL)]
    # alphabetical: African, Amerindigenous, European
    ancestries = [("AFR", "AFR (African)"),
                  ("NAT", "AMR (Amerindigenous)"),
                  ("EUR", "EUR (European)")]

    fig, axes = plt.subplots(1, 3, figsize=(13, 6.4), dpi=180, sharey=True)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.83, bottom=0.24, wspace=0.20)

    x = np.arange(len(panels))
    w = 0.36

    for ax, (anc, anc_label) in zip(axes, ancestries):
        for i, (track, tlabel, color) in enumerate(tracks):
            vals = [wgs.get((track, key, anc), np.nan) for _, key in panels]
            bars = ax.bar(x + (i - 0.5) * w, vals, w,
                          label=tlabel, color=color,
                          edgecolor="white", linewidth=0.7)
            for xi, v in zip(bars, vals):
                if not np.isnan(v):
                    ax.text(xi.get_x() + xi.get_width()/2, v + 0.005,
                            f"{v*100:.0f}%", ha="center", va="bottom",
                            fontsize=11, color="#111", fontweight="bold")

        ax.set_title(anc_label, fontweight="bold", fontsize=15, pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels([lab for lab, _ in panels], fontsize=11, fontweight="bold")
        for lbl in ax.get_xticklabels():
            lbl.set_multialignment("center")
        ax.tick_params(axis="y", labelsize=10)

    axes[0].set_ylim(0.85, 1.00)
    axes[0].set_ylabel("Weighted per-hap recall  (chr20-22)",
                       fontsize=12.5, fontweight="bold")

    fig.suptitle("Local Ancestry Inference Accuracy",
                 fontweight="bold", fontsize=16, y=0.965)

    leg = axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.24),
                         ncol=2, fontsize=13, frameon=False)
    for t in leg.get_texts():
        t.set_fontweight("bold")

    out = OUT / "fig_wgs_all_ancestry.png"
    fig.savefig(out)
    plt.close(fig)
    return out

if __name__ == "__main__":
    for f in (fig_wgs_panels, fig_wgs_all_ancestry, fig_chip_vs_wgs,
              fig_amr_disparity, fig_pipeline):
        p = f()
        print(f"wrote {p}  ({p.stat().st_size/1024:.0f} kB)")
