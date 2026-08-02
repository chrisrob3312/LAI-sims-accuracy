#!/usr/bin/env python3
"""
Produce publication-ready metrics tables + rendered PNG for the committee deck.

Outputs:
  results/tables/publication_table_full.tsv     — every panel×cohort×density×ancestry
  results/tables/publication_table_amr.tsv      — AMR-only summary
  slides/tables/publication_table_amr.png       — rendered PNG for slide use
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO   = Path(__file__).resolve().parent.parent
OUT_T  = REPO / "results/tables";  OUT_T.mkdir(parents=True, exist_ok=True)
OUT_P  = REPO / "slides/tables";   OUT_P.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- read
def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

wgs  = read_tsv(REPO / "results/accuracy_pilot/accuracy_summary.tsv")
chip = read_tsv(REPO / "results/accuracy_pilot_chip_gsa/accuracy_summary.tsv")

# ---------------------------------------------------------------- panel ordering + friendly labels
PANEL_ORDER = [
    ("NAT_HGDP",         "HGDP + 1KG"),
    ("NAT_HGDPMXB",      "+ 25 MXB"),
    ("NAT_HGDPMXB_FULL", "+ 50 MXB"),
    ("NAT_HOMOG",        "Homog. Q≥0.95"),
    ("NAT_PEL",          "+ PEL (Peru)"),
    ("NAT_PEL_EAS",      "+ PEL + EAS"),
]
TRACK_LABEL = {"NAT": "Brazilian-like", "NATMXB": "Mexican-like"}
ANC_LABEL   = {"AFR": "African", "NAT": "Amerindigenous", "EUR": "European"}

# ---------------------------------------------------------------- full TSV (long tidy)
def write_full():
    fieldnames = ["density", "cohort", "panel_key", "panel_label",
                  "ancestry", "n_haps", "n_sites",
                  "tp", "fp", "fn",
                  "recall", "precision", "f1"]
    with open(OUT_T / "publication_table_full.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for density, rows in [("WGS", wgs), ("GSA_chip", chip)]:
            for r in rows:
                panel_lbl = dict(PANEL_ORDER).get(r["panel"], r["panel"])
                w.writerow(dict(
                    density=density,
                    cohort=TRACK_LABEL.get(r["track"], r["track"]),
                    panel_key=r["panel"],
                    panel_label=panel_lbl,
                    ancestry=ANC_LABEL.get(r["ancestry"], r["ancestry"]),
                    n_haps=r["n_haps_chr"],
                    n_sites=r["total_sites"],
                    tp=r["tp"], fp=r["fp"], fn=r["fn"],
                    recall=f"{float(r['weighted_recall']):.4f}",
                    precision=f"{float(r['weighted_precision']):.4f}",
                    f1=f"{float(r['weighted_f1']):.4f}",
                ))
    print(f"wrote {OUT_T / 'publication_table_full.tsv'}")

# ---------------------------------------------------------------- AMR-focused compact TSV
def build_amr_rows():
    """Return list of dicts: one row per (panel, cohort, density) for AMR only."""
    def by_key(rows, panel_key, track):
        for r in rows:
            if r["panel"] == panel_key and r["track"] == track and r["ancestry"] == "NAT":
                return r
        return None
    out = []
    for panel_key, panel_lbl in PANEL_ORDER:
        for track_key, track_lbl in TRACK_LABEL.items():
            row = dict(panel=panel_lbl, cohort=track_lbl)
            for density, rows in [("WGS", wgs), ("chip", chip)]:
                r = by_key(rows, panel_key, track_key)
                if r is None:
                    row[f"{density}_recall"]    = None
                    row[f"{density}_precision"] = None
                    row[f"{density}_f1"]        = None
                else:
                    row[f"{density}_recall"]    = float(r["weighted_recall"])
                    row[f"{density}_precision"] = float(r["weighted_precision"])
                    row[f"{density}_f1"]        = float(r["weighted_f1"])
            out.append(row)
    return out

def write_amr_tsv(rows):
    fieldnames = ["panel", "cohort",
                  "WGS_recall", "WGS_precision", "WGS_f1",
                  "chip_recall", "chip_precision", "chip_f1"]
    with open(OUT_T / "publication_table_amr.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({
                k: ("" if r[k] is None else
                    (r[k] if k in ("panel", "cohort") else f"{r[k]:.4f}"))
                for k in fieldnames
            })
    print(f"wrote {OUT_T / 'publication_table_amr.tsv'}")

# ---------------------------------------------------------------- rendered PNG
def render_amr_png(rows):
    """
    Publication-style table.

    Layout:
      Columns:   Panel | Brazilian-like (Recall, Prec., F1) | Mexican-like (Recall, Prec., F1)
      Row-band 1: WGS density   — 6 panel rows
      Row-band 2: GSA chip      — 5-6 panel rows (missing panels shown as —)

    Winners bolded per (density band × cohort × metric).
    """
    NAVY   = "#0E2841"
    TEAL   = "#156082"    # Mexican-like
    ORANGE = "#E97132"    # Brazilian-like
    BAND_W = "#F0F6FA"    # WGS band background
    BAND_C = "#FFF3E9"    # chip band background
    HDR_BG = "#D6E6EE"
    BOLD_C = "#B34A17"    # deep orange for winner cells

    # regroup rows -> {(density, panel_label, cohort): (recall, prec, f1)}
    def cell(row, density, metric):
        return row.get(f"{density}_{metric}")

    # winners per (density × cohort × metric)
    winners = {}
    for density in ("WGS", "chip"):
        for cohort in ("Brazilian-like", "Mexican-like"):
            subset = [r for r in rows if r["cohort"] == cohort]
            for metric in ("recall", "precision", "f1"):
                vals = [cell(r, density, metric) for r in subset
                        if cell(r, density, metric) is not None]
                if vals:
                    winners[(density, cohort, metric)] = max(vals)

    # rows per density band: one row per panel, with 6 metric cols (Brasa | Mex)
    def band_rows(density):
        out = []
        for panel_key, panel_lbl in PANEL_ORDER:
            brasa = next((r for r in rows
                          if r["cohort"] == "Brazilian-like"
                          and r["panel"] == panel_lbl), None)
            mex   = next((r for r in rows
                          if r["cohort"] == "Mexican-like"
                          and r["panel"] == panel_lbl), None)
            out.append({
                "panel": panel_lbl,
                "brasa_recall":    (cell(brasa, density, "recall")    if brasa else None),
                "brasa_precision": (cell(brasa, density, "precision") if brasa else None),
                "brasa_f1":        (cell(brasa, density, "f1")        if brasa else None),
                "mex_recall":      (cell(mex,   density, "recall")    if mex   else None),
                "mex_precision":   (cell(mex,   density, "precision") if mex   else None),
                "mex_f1":          (cell(mex,   density, "f1")        if mex   else None),
            })
        return out

    wgs_rows  = band_rows("WGS")
    chip_rows = band_rows("chip")

    # only keep panels that have at least one non-None value in the band
    def prune(band): return [r for r in band if any(
        r[k] is not None for k in
        ("brasa_recall","brasa_precision","brasa_f1","mex_recall","mex_precision","mex_f1"))]
    wgs_rows, chip_rows = prune(wgs_rows), prune(chip_rows)

    col_headers = ["Panel",
                   "Recall", "Prec.", "F1",
                   "Recall", "Prec.", "F1"]
    n_cols = 7
    # column widths (must sum to 1)
    col_w = [0.22,  0.135, 0.125, 0.125,  0.135, 0.125, 0.135]
    assert abs(sum(col_w) - 1.0) < 1e-6, f"col_w sums to {sum(col_w)}"
    col_x = [sum(col_w[:i]) for i in range(n_cols + 1)]

    # layout: header row 1 (cohort banners), header row 2 (col names),
    # WGS band label, WGS rows, chip band label, chip rows
    n_rows = 2 + 1 + len(wgs_rows) + 1 + len(chip_rows)
    row_h  = 1.0 / n_rows
    def row_y(i): return 1.0 - (i + 1) * row_h

    fig, ax = plt.subplots(figsize=(11.5, 0.44 * n_rows + 1.2), dpi=200)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # -------- header row 1: cohort banners
    ax.add_patch(plt.Rectangle((0, row_y(0)), col_x[1], row_h,
                               facecolor=NAVY, edgecolor="none"))
    ax.text(col_x[1]/2, row_y(0) + row_h/2, "AMR-tract calling accuracy",
            ha="center", va="center", color="white",
            fontsize=11.5, fontweight="bold")
    # Brasa banner
    ax.add_patch(plt.Rectangle((col_x[1], row_y(0)), col_x[4] - col_x[1], row_h,
                               facecolor=ORANGE, edgecolor="none"))
    ax.text((col_x[1] + col_x[4])/2, row_y(0) + row_h/2, "Brazilian-like cohort",
            ha="center", va="center", color="white",
            fontsize=12, fontweight="bold")
    # Mexican banner
    ax.add_patch(plt.Rectangle((col_x[4], row_y(0)), col_x[7] - col_x[4], row_h,
                               facecolor=TEAL, edgecolor="none"))
    ax.text((col_x[4] + col_x[7])/2, row_y(0) + row_h/2, "Mexican-like cohort",
            ha="center", va="center", color="white",
            fontsize=12, fontweight="bold")

    # -------- header row 2: metric names
    for c in range(n_cols):
        ax.add_patch(plt.Rectangle((col_x[c], row_y(1)), col_w[c], row_h,
                                   facecolor=HDR_BG, edgecolor="white", linewidth=0.5))
        ax.text(col_x[c] + col_w[c]/2, row_y(1) + row_h/2, col_headers[c],
                ha="center", va="center",
                fontsize=10.5, fontweight="bold", color=NAVY)

    # -------- helper to draw one data-row
    def draw_row(i, band_label, r, band_bg):
        y = row_y(i)
        ax.add_patch(plt.Rectangle((0, y), 1.0, row_h,
                                   facecolor=band_bg, edgecolor="none"))
        text_cells = [
            r["panel"],
            _fmt(r["brasa_recall"]),  _fmt(r["brasa_precision"]),  _fmt(r["brasa_f1"]),
            _fmt(r["mex_recall"]),    _fmt(r["mex_precision"]),    _fmt(r["mex_f1"]),
        ]
        cell_keys = [None,
                     ("brasa_recall","Brazilian-like","recall"),
                     ("brasa_precision","Brazilian-like","precision"),
                     ("brasa_f1","Brazilian-like","f1"),
                     ("mex_recall","Mexican-like","recall"),
                     ("mex_precision","Mexican-like","precision"),
                     ("mex_f1","Mexican-like","f1")]
        for c in range(n_cols):
            val = text_cells[c]
            key = cell_keys[c]
            bold = False; color = "#111"
            if key is not None:
                col_key, cohort, metric = key
                v = r[col_key]
                w = winners.get((band_label, cohort, metric))
                if v is not None and w is not None and v == w:
                    bold = True; color = BOLD_C
            ha = "left" if c == 0 else "center"
            xtext = (col_x[c] + 0.008) if ha == "left" else (col_x[c] + col_w[c]/2)
            ax.text(xtext, y + row_h/2, val,
                    ha=ha, va="center",
                    fontsize=10.5 if bold else 10,
                    fontweight=("bold" if bold else "normal"),
                    color=color)

    # -------- WGS band label row
    y = row_y(2)
    ax.add_patch(plt.Rectangle((0, y), 1.0, row_h,
                               facecolor=NAVY, edgecolor="none"))
    ax.text(0.01, y + row_h/2, "  WGS density",
            ha="left", va="center", color="white",
            fontsize=11, fontweight="bold")

    for i, r in enumerate(wgs_rows):
        draw_row(3 + i, "WGS", r, BAND_W)

    # -------- Chip band label row
    chip_hdr_i = 3 + len(wgs_rows)
    y = row_y(chip_hdr_i)
    ax.add_patch(plt.Rectangle((0, y), 1.0, row_h,
                               facecolor=NAVY, edgecolor="none"))
    ax.text(0.01, y + row_h/2, "  GSA chip density (unimputed)",
            ha="left", va="center", color="white",
            fontsize=11, fontweight="bold")

    for i, r in enumerate(chip_rows):
        draw_row(chip_hdr_i + 1 + i, "chip", r, BAND_C)

    # dividers
    ax.plot([0, 1], [row_y(2), row_y(2)], color=NAVY, lw=0.8)
    ax.plot([0, 1], [row_y(chip_hdr_i), row_y(chip_hdr_i)], color=NAVY, lw=0.8)
    # column vertical rule between cohorts
    ax.plot([col_x[4], col_x[4]], [row_y(1), 0], color=NAVY, lw=0.8)

    # footer
    ax.text(0, -0.02,
            "Weighted per-hap metrics · winner per density×cohort×metric bolded (orange) · "
            "chr1-8+20-22 for full-coverage WGS panels; chip pilot chr20-22.",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9, color="#555", style="italic")

    out = OUT_P / "publication_table_amr.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

def _fmt(v):
    return "—" if v is None else f"{v:.3f}"

if __name__ == "__main__":
    write_full()
    amr = build_amr_rows()
    write_amr_tsv(amr)
    render_amr_png(amr)
