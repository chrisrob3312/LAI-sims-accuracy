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
    Publication-style table. Rows grouped by cohort, then panel.
    Winning cell in each cohort×density×metric column is bolded.
    """
    NAVY   = "#0E2841"
    TEAL   = "#156082"
    ORANGE = "#E97132"
    ROW_A  = "#FFFFFF"
    ROW_B  = "#F4F6F8"
    HDR_BG = "#D6E6EE"

    # sort: cohort first, then panel order
    cohort_order = ["Mexican-like", "Brazilian-like"]
    rows_sorted = []
    for coh in cohort_order:
        for panel_key, panel_lbl in PANEL_ORDER:
            for r in rows:
                if r["cohort"] == coh and r["panel"] == panel_lbl:
                    rows_sorted.append(r); break

    # winners per (cohort × density × metric)
    winners = {}
    for coh in cohort_order:
        subset = [r for r in rows_sorted if r["cohort"] == coh]
        for density in ("WGS", "chip"):
            for metric in ("recall", "precision", "f1"):
                key = f"{density}_{metric}"
                vals = [r[key] for r in subset if r[key] is not None]
                if vals:
                    winners[(coh, key)] = max(vals)

    # header layout: 2 rows of headers
    #   [ Panel | Cohort | WGS(recall|prec|F1) | Chip(recall|prec|F1) ]
    col_headers_row2 = ["Panel", "Cohort",
                        "Recall", "Prec.", "F1",
                        "Recall", "Prec.", "F1"]
    n_cols = len(col_headers_row2)
    n_rows = len(rows_sorted) + 2  # 2 header rows

    fig, ax = plt.subplots(figsize=(12.5, 0.42 * n_rows + 1.4), dpi=200)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # column widths (relative, must sum to 1)
    col_w = [0.22, 0.15,  0.115, 0.105, 0.105,  0.115, 0.095, 0.095]
    assert abs(sum(col_w) - 1.0) < 1e-6, f"col_w sums to {sum(col_w)}"
    col_x = [sum(col_w[:i]) for i in range(n_cols + 1)]

    row_h = 1.0 / n_rows
    # y from top -> row index 0 at top
    def row_y(i): return 1.0 - (i + 1) * row_h

    # header row 1 (density banner across cols 2-4 and 5-7)
    ax.add_patch(plt.Rectangle((0, row_y(0)), col_x[2], row_h,
                               facecolor=NAVY, edgecolor="none"))
    ax.text((0 + col_x[2]) / 2, row_y(0) + row_h/2, "AMR-tract calling accuracy",
            ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")
    ax.add_patch(plt.Rectangle((col_x[2], row_y(0)), col_x[5] - col_x[2], row_h,
                               facecolor=TEAL, edgecolor="none"))
    ax.text((col_x[2] + col_x[5]) / 2, row_y(0) + row_h/2, "WGS density",
            ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")
    ax.add_patch(plt.Rectangle((col_x[5], row_y(0)), col_x[8] - col_x[5], row_h,
                               facecolor=ORANGE, edgecolor="none"))
    ax.text((col_x[5] + col_x[8]) / 2, row_y(0) + row_h/2, "GSA chip density",
            ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")

    # header row 2 (per-column labels)
    for c in range(n_cols):
        ax.add_patch(plt.Rectangle((col_x[c], row_y(1)), col_w[c], row_h,
                                   facecolor=HDR_BG, edgecolor="white", linewidth=0.5))
        ax.text(col_x[c] + col_w[c]/2, row_y(1) + row_h/2, col_headers_row2[c],
                ha="center", va="center", fontsize=10.5, fontweight="bold", color=NAVY)

    # data rows
    for i, r in enumerate(rows_sorted):
        y = row_y(i + 2)
        # zebra
        stripe = ROW_A if i % 2 == 0 else ROW_B
        ax.add_patch(plt.Rectangle((0, y), 1.0, row_h,
                                   facecolor=stripe, edgecolor="none"))
        # cohort banner change
        text_cells = [
            r["panel"],
            r["cohort"],
            _fmt(r["WGS_recall"]),  _fmt(r["WGS_precision"]),  _fmt(r["WGS_f1"]),
            _fmt(r["chip_recall"]), _fmt(r["chip_precision"]), _fmt(r["chip_f1"]),
        ]
        cell_keys = [None, None,
                     "WGS_recall", "WGS_precision", "WGS_f1",
                     "chip_recall", "chip_precision", "chip_f1"]
        for c in range(n_cols):
            val = text_cells[c]
            key = cell_keys[c]
            bold = False
            color = "#111"
            if key is not None:
                wkey = (r["cohort"], key)
                if wkey in winners and r[key] == winners[wkey]:
                    bold = True; color = ORANGE
            ha = "left" if c == 0 else ("center" if c > 1 else "left")
            xpad = 0.008
            xtext = (col_x[c] + xpad) if ha == "left" else (col_x[c] + col_w[c]/2)
            ax.text(xtext, y + row_h/2, val,
                    ha=ha, va="center",
                    fontsize=10 if not bold else 10.5,
                    fontweight=("bold" if bold else "normal"),
                    color=color)

    # thin dividers between cohorts + columns
    ax.plot([0, 1], [row_y(2), row_y(2)], color=NAVY, lw=0.9)  # under header row 2
    # cohort divider
    cutoff_i = sum(1 for r in rows_sorted if r["cohort"] == cohort_order[0])
    y_cut = row_y(cutoff_i + 2)
    ax.plot([0, 1], [y_cut, y_cut], color=NAVY, lw=0.6, linestyle="--")

    # footer note
    ax.text(0, -0.02,
            "Weighted per-hap metrics · chr20–22 pilot · winner per cohort×density×metric bolded (orange)",
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
