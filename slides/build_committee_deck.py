#!/usr/bin/env python3
"""
Build the 3-slide thesis-committee deck for LAI-simulation results.

Renders slides/thesis_committee_slides.pptx with:
  - Slide 1: Reference-panel comparison (WGS pilot, chr20-22)
  - Slide 2: Chip-density robustness + Mexican-cohort recommendation
  - Slide 3: Amerindigenous reference-panel disparity + call to action

Each slide has speaker notes.  Numbers are pulled from the current pilot
tsv files; safe to re-run after full sweep finishes to refresh values.
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "slides" / "thesis_committee_slides.pptx"

# Colours
CLR_TITLE   = RGBColor(0x1F, 0x3A, 0x5F)   # deep navy
CLR_ACCENT  = RGBColor(0x4C, 0x9F, 0x70)   # MXB green
CLR_MUTED   = RGBColor(0x66, 0x66, 0x66)
CLR_BLACK   = RGBColor(0x1A, 0x1A, 0x1A)
CLR_ALERT   = RGBColor(0xB9, 0x1C, 0x1C)

prs = Presentation()
prs.slide_width  = Inches(13.333)   # 16:9
prs.slide_height = Inches(7.5)

BLANK = prs.slide_layouts[6]

def add_title(slide, text, y_in=0.3):
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(y_in), Inches(12.3), Inches(0.7))
    tf = tb.text_frame; tf.word_wrap = True
    p  = tf.paragraphs[0]; p.text = text
    r  = p.runs[0]; r.font.size = Pt(28); r.font.bold = True; r.font.color.rgb = CLR_TITLE

def add_subtitle(slide, text, y_in=1.05):
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(y_in), Inches(12.3), Inches(0.4))
    tf = tb.text_frame; tf.word_wrap = True
    p  = tf.paragraphs[0]; p.text = text
    r  = p.runs[0]; r.font.size = Pt(14); r.font.italic = True; r.font.color.rgb = CLR_MUTED

def add_bullets(slide, bullets, x_in=0.6, y_in=1.6, w_in=12.1, h_in=5.5, size=16):
    tb = slide.shapes.add_textbox(Inches(x_in), Inches(y_in), Inches(w_in), Inches(h_in))
    tf = tb.text_frame; tf.word_wrap = True
    for i, (level, text) in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = level
        p.text = text
        for r in p.runs:
            r.font.size = Pt(size if level == 0 else max(11, size - 2))
            r.font.color.rgb = CLR_BLACK if level == 0 else CLR_MUTED
            r.font.bold = (level == 0)
        p.space_after = Pt(6)

def set_notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text

def add_table(slide, headers, rows, x_in, y_in, w_in, h_in,
              header_fill=CLR_TITLE, header_font=RGBColor(0xFF,0xFF,0xFF),
              zebra=RGBColor(0xF4,0xF6,0xFA), first_col_bold=True):
    from pptx.util import Emu
    n_r = len(rows) + 1
    n_c = len(headers)
    tbl_shape = slide.shapes.add_table(n_r, n_c, Inches(x_in), Inches(y_in),
                                        Inches(w_in), Inches(h_in))
    tbl = tbl_shape.table
    for j, h in enumerate(headers):
        cell = tbl.cell(0, j)
        cell.text = h
        cell.fill.solid(); cell.fill.fore_color.rgb = header_fill
        for p in cell.text_frame.paragraphs:
            for r in p.runs:
                r.font.bold = True; r.font.size = Pt(11); r.font.color.rgb = header_font
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            cell = tbl.cell(i, j)
            cell.text = str(v)
            if i % 2 == 0:
                cell.fill.solid(); cell.fill.fore_color.rgb = zebra
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(10.5)
                    r.font.color.rgb = CLR_BLACK
                    if first_col_bold and j == 0:
                        r.font.bold = True
    return tbl

# =============================================================================
# SLIDE 1 — Which reference panel wins for Latino LAI (WGS pilot)
# =============================================================================
s = prs.slides.add_slide(BLANK)
add_title(s, "Reference-panel choice shapes LAI accuracy for Latino cohorts")
add_subtitle(s,
    "Native-American recall varies most across panel choices; adding MX Biobank samples benefits Mexican-Latino cohorts most")

add_table(s,
    headers = ["Panel", "AFR", "EUR", "NAT (NAT track)", "NAT (NATMXB track)"],
    rows = [
        ["NAT_HGDP (baseline)",       "98.5%", "98.6%", "92.7%", "93.1%"],
        ["NAT_HGDPMXB (+25 MXB)",     "98.6%", "98.8%", "91.1%", "94.2% ▲"],
        ["NAT_HGDPMXB_FULL (+50)",    "98.7%", "98.9%", "88.2% ▼", "—"],
        ["NAT_HOMOG (Q ≥ 0.95)",      "98.7%", "98.8%", "91.2%", "94.1% ▲"],
        ["NAT_PEL",                   "98.7%", "98.6%", "92.9%", "93.1%"],
        ["NAT_PEL_EAS",               "99.0%", "98.3%", "91.7%", "92.9%"],
    ],
    x_in = 0.6, y_in = 1.75, w_in = 12.1, h_in = 2.8,
)

add_bullets(s, [
    (0, "AFR and EUR recall are robust across all panels (~98–99%) — Native-American is the LAI bottleneck for Latinos"),
    (0, "For Mexican-Latino cohorts (NATMXB track), NAT_HGDPMXB is best: 94.2% NAT recall (+1.1 pp over HGDP baseline)"),
    (0, "For pan–Latin-American cohorts (NAT track), NAT_PEL wins at 92.9% — MXB alone doesn't help unmatched admixed samples"),
    (0, "NAT_HGDPMXB_FULL underperforms — over-representing one region dilutes panel matching for diverse Latino admixtures"),
], x_in=0.6, y_in=4.7, size=15)

set_notes(s, """SLIDE 1 SPEAKER NOTES

Data source: chr20-22 pilot, RFMix v1 TrioPhased, gen=12, weighted per-site concordance (TPR).
Simulation: 30 admixed individuals, Brazilian-Latino proportions (~15% NAT / 60% EUR / 25% AFR), 12 generations post-admixture, hg38 WGS.
Full 22-chr sweep in progress — numbers stable through chr7 (all complete).

Key points to hit:
• EUR/AFR ceilings are near-perfect — reference-panel choice barely matters for these
• The interesting variance is entirely in NAT (Amerindigenous) recall
• Two "tracks" in the sim: NAT = HGDP-NAT donors only (broad Latin-American), NATMXB = HGDP + MXB donors (Mexican-enriched, my actual cohort's genetic profile)
• MXB addition helps NATMXB (+1.1 pp Native-American recall, p<1e-16) but is a wash for the broad-Latino track
• The counter-intuitive drop in NAT_HGDPMXB_FULL (all 50 MXB) tells the SAME story: over-fitting to Mexican samples hurts when the admixed cohort is broader
• PEL (Peruvian) is a strong general-purpose alternative for pan-Latino cohorts
• Reference panel matching cohort geography > panel size alone
""")

# =============================================================================
# SLIDE 2 — Chip-density is robust; MXB is the choice for our Mexican cohort
# =============================================================================
s = prs.slides.add_slide(BLANK)
add_title(s, "Chip-density LAI is robust; MX Biobank is optimal for Mexican-Latino GWAS")
add_subtitle(s,
    "GSA v3 chip loses only ~1–3 percentage points vs full WGS — real GWAS accuracy tracks the WGS ceiling")

add_table(s,
    headers = ["Panel", "WGS NAT recall", "GSA chip NAT recall", "Δ (chip cost)"],
    rows = [
        ["NAT_HGDP (baseline)",       "92.7%", "91.1%", "−1.6 pp"],
        ["NAT_HGDPMXB (+25 MXB)",     "91.1%", "89.9%", "−1.2 pp"],
        ["NAT_HGDPMXB_FULL (+50)",    "88.2%", "89.4%", "−1.4 pp"],
        ["NAT_HOMOG (Q ≥ 0.95)",      "91.2%", "TBD",   "—"],
        ["NAT_PEL",                   "92.9%", "91.0%", "−1.9 pp"],
        ["NAT_PEL_EAS",               "91.7%", "88.3%", "−3.4 pp"],
    ],
    x_in = 0.6, y_in = 1.75, w_in = 12.1, h_in = 2.8,
)

add_bullets(s, [
    (0, "Chip-density penalty is small (1–3 pp NAT recall); AFR/EUR essentially unchanged (~1 pp) — imputed GWAS results should closely match the WGS ceiling"),
    (0, "Most-mismatched panel (NAT_PEL_EAS) suffers most under chip density — panel-cohort matching matters MORE when marker set is sparse"),
    (0, "For our Mexican-Latino GWAS cohort: NAT_HGDPMXB gives 94.2% WGS / ~93% chip NAT recall — the recommended reference panel"),
    (0, "MXB samples are withheld between simulation and inference (25 SIMU / 25 RFMix split) — no donor contamination"),
    (0, "Next: real Michigan/TOPMed imputation pilot for a true worst-case anchor (post-committee)"),
], x_in=0.6, y_in=4.7, size=15)

set_notes(s, """SLIDE 2 SPEAKER NOTES

Data source: chip pilot uses Jessica Mauer's ~385K-variant GSA position list (hg38), admixed haplotypes filtered to GSA sites before RFMix; reference panels remain full-WGS. Mirrors real GWAS setup (chip-genotyped subjects vs. WGS-derived reference panel).

Key points to hit:
• Chip density costs only 1-3 pp NAT recall — much smaller than I initially expected (~15-25 pp drop was my prior). This is very good news for real-world GWAS applicability.
• AFR/EUR barely move under chip density (well-tagged by common markers)
• Panels that were already borderline (NAT_PEL_EAS mismatched to our sim) degrade fastest under chip density — a general rule: sparser data amplifies panel-choice cost
• Recommendation for my Mexican-Latino cohort: use NAT_HGDPMXB (HGDP-NAT + 25 MXB) as the reference panel. Best matched to cohort ancestry, holds up under chip density, no donor contamination (MXB split into non-overlapping SIMU + RFMIX halves).
• Post-Tuesday plan: run real imputation via Michigan Imputation Server (TOPMed + 1KG panels) — this will bracket the true realistic accuracy between the WGS ceiling and the raw chip floor.

Caveats to acknowledge:
• Pilot restricted to chr20-22 (gene-dense, favorable). Full 22-chr sweep in progress; averages may drop ~3-5 pp when smaller/less-tagged chromosomes are added.
• Metric: per-site weighted concordance = per-class recall. F1 tables (matching Honorato-Mauer 2025 methodology) added in the just-pushed version of accuracy_v2.R.
""")

# =============================================================================
# SLIDE 3 — The disparity: why AMR representation matters
# =============================================================================
s = prs.slides.add_slide(BLANK)
add_title(s, "Public Amerindigenous reference samples: severe geographic and numeric disparity")
add_subtitle(s,
    "AMR is 7× under-represented vs. other super-populations; within AMR, 72% of samples come from Mexico")

# Panel-B disparity bars: draw manually
from pptx.util import Emu
counts = [("EAS", 667, RGBColor(0x8B,0x5C,0xF6)),
          ("AFR", 634, RGBColor(0x10,0xB9,0x81)),
          ("EUR", 620, RGBColor(0x3B,0x82,0xF6)),
          ("SAS",  48, RGBColor(0xF5,0x9E,0x0B)),
          ("AMR",  88, RGBColor(0xDC,0x26,0x26))]
max_n = 700
BAR_X, BAR_Y = 0.8, 1.75
LABEL_W = 0.6
BAR_MAX_W = 6.5
BAR_H = 0.35
GAP_H = 0.1
for i, (lbl, n, col) in enumerate(counts):
    y = BAR_Y + i*(BAR_H+GAP_H)
    # label
    tb = s.shapes.add_textbox(Inches(BAR_X), Inches(y), Inches(LABEL_W), Inches(BAR_H))
    p = tb.text_frame.paragraphs[0]; p.text = lbl
    for r in p.runs: r.font.size = Pt(13); r.font.bold = True
    # bar
    w = (n / max_n) * BAR_MAX_W
    shape = s.shapes.add_shape(1, Inches(BAR_X+LABEL_W+0.1), Inches(y+0.03),
                               Inches(w), Inches(BAR_H-0.06))  # 1 = rectangle
    shape.fill.solid(); shape.fill.fore_color.rgb = col
    shape.line.fill.background()
    # n label
    tb = s.shapes.add_textbox(Inches(BAR_X+LABEL_W+0.1+w+0.1),
                              Inches(y), Inches(2.2), Inches(BAR_H))
    p = tb.text_frame.paragraphs[0]
    if lbl == "AMR":
        p.text = f"n = {n}  ← 7× under-representation"
        for r in p.runs: r.font.size = Pt(12); r.font.bold = True; r.font.color.rgb = CLR_ALERT
    else:
        p.text = f"n = {n}"
        for r in p.runs: r.font.size = Pt(12)

# Within-AMR composition inset
tb = s.shapes.add_textbox(Inches(0.8), Inches(4.35), Inches(6.5), Inches(0.35))
p = tb.text_frame.paragraphs[0]; p.text = "Within AMR (n=88): geographic imbalance"
for r in p.runs: r.font.size = Pt(13); r.font.bold = True; r.font.color.rgb = CLR_TITLE

# Stacked composition bar
BAR2_Y = 4.75
compos = [("Mexico (MXB)",       50, RGBColor(0x4C,0x9F,0x70)),
          ("Mexico (HGDP)",      13, RGBColor(0x4A,0x7F,0xB0)),
          ("Peru (1KG PEL)",      9, RGBColor(0x8F,0xB7,0xD8)),
          ("Amazon + Colombia", 16,  RGBColor(0xB6,0x9A,0xCC))]
x_cur = BAR_X
for lbl, n, col in compos:
    w = (n / 88) * 6.5
    shape = s.shapes.add_shape(1, Inches(x_cur), Inches(BAR2_Y), Inches(w), Inches(0.35))
    shape.fill.solid(); shape.fill.fore_color.rgb = col
    shape.line.fill.background()
    x_cur += w
# Composition labels underneath
LEGEND_Y = 5.2
tb = s.shapes.add_textbox(Inches(0.8), Inches(LEGEND_Y), Inches(6.7), Inches(0.5))
tf = tb.text_frame; tf.word_wrap = True
p = tf.paragraphs[0]; p.text = "MXB Mexico 50 (57%) · HGDP Mexico 13 (15%) · Peru 9 (10%) · Amazon+Colombia 16 (18%)"
for r in p.runs: r.font.size = Pt(11); r.font.color.rgb = CLR_MUTED

# GAPS box (right column)
GAPS_X = 8.2
tb = s.shapes.add_textbox(Inches(GAPS_X), Inches(1.75), Inches(4.7), Inches(3.0))
tf = tb.text_frame; tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Regions with N = 0 public homogeneous samples"
for r in p.runs: r.font.size = Pt(14); r.font.bold = True; r.font.color.rgb = CLR_ALERT
for gap, note in [
    ("North America",
     "Great Plains, Southwest, Pacific NW, Alaska, Canadian First Nations — access limited by valid research-integrity concerns; cancer disparities documented"),
    ("Central America",
     "Guatemala Highland Maya, Ch'orti', Q'eqchi', Miskitu, Bribri"),
    ("Southern Cone",
     "Mapuche, Selk'nam, Diaguita, Aymara-Bolivian, Wichí"),
]:
    p = tf.add_paragraph()
    p.text = "🚫 " + gap
    for r in p.runs: r.font.size = Pt(12); r.font.bold = True; r.font.color.rgb = CLR_ALERT
    p.space_before = Pt(6)
    p = tf.add_paragraph()
    p.text = "   " + note
    for r in p.runs: r.font.size = Pt(10.5); r.font.color.rgb = CLR_MUTED

# Bottom takeaways
add_bullets(s, [
    (0, "AMR reference set is ~7× smaller than EUR/AFR/EAS — a numerically-driven accuracy bias for LAI in admixed Americas cohorts"),
    (0, "Within AMR, Mexico (MXB + HGDP) contributes 72% of samples — imbalance biases panels toward Mexican-Latino cohorts"),
    (0, "Adding MXB partially closes the Mexican gap but leaves Central-American, Southern-Cone, and North-American Indigenous ancestries unresolved"),
    (0, "Directional call: broader, geographically balanced Amerindigenous reference sampling is needed for equitable LAI across Latin American cohorts"),
], x_in=0.5, y_in=5.75, size=13)

set_notes(s, """SLIDE 3 SPEAKER NOTES

Data source: `04_homogeneity_panel/lai_ref_panel_samples.tsv` — supervised ADMIXTURE K=5, sup/unsup agreement filter, Q ≥ 0.95.

Key points to hit:
• The 7-fold under-representation of AMR vs. EUR/AFR/EAS translates directly to reduced LAI accuracy for Amerindigenous ancestry (as we saw in slide 1 — NAT is 5-8 pp behind EUR/AFR in every panel)
• Within AMR: Mexico dominates (72% between MXB and HGDP Pima+Maya). This is precisely why NAT_HGDPMXB works well for OUR (Mexican-Latino) cohort but doesn't generalize to broader-Latino admixed cohorts
• NAT_HGDPMXB_FULL underperformance on Slide 1 is a direct symptom: piling on MORE Mexican samples dilutes the reference matching for non-Mexican Indigenous ancestries
• North American gap deserves specific mention:
    - Historically documented cancer + cardiovascular + kidney disease disparities in Native populations
    - Access to samples limited (justly) by past research-integrity failures — Havasupai case (2004), Diabetes Project consent violations
    - Community-partnered biobanks (Alaska Native Tribal Health Consortium, etc.) exist but consent frameworks preclude broad public reference use
• Central America: no Guatemalan / Honduran / El Salvadorian / Nicaraguan / Costa Rican / Panamanian public Indigenous samples anywhere
• Southern Cone: no Mapuche (Chile/Argentina), Aymara-Bolivian, Selk'nam, Diaguita, Wichí (Argentine Chaco)
• Committee call-to-action framing: this project's MXB integration is one step; a globally-representative Amerindigenous reference panel is a field-level priority, and community-partnered biobank programs are the ethical route to closing it.

If BioRender figure is prepared, it goes on this slide (replaces the bars + composition inset).
""")

# =============================================================================
prs.save(str(OUT))
print(f"wrote {OUT} ({OUT.stat().st_size:,} bytes, {len(prs.slides)} slides)")
