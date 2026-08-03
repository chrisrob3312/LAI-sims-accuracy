#!/usr/bin/env python3
"""
Build the thesis-committee deck for LAI-simulation + preprocessing pipeline update.

Style: matches Spring TAC deck — 16:9, sentence-headline on top, one dominant
figure, tiny callouts. Text depth lives in speaker notes.

Slides:
  1. Title
  2. Motivation — Amerindigenous reference disparity
  3. Methods — simulation + LAI accuracy pipeline
  4. Result — WGS panel choice depends on cohort composition
  5. Result — chip-density LAI holds within ~1% of WGS
  6. Preprocessing pipeline update — 8-module Nextflow + status
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from datetime import date

REPO = Path(__file__).resolve().parent.parent
FIGDIR = REPO / "slides" / "figures"
OUT    = REPO / "slides" / "thesis_committee_slides.pptx"

# palette matches the Spring TAC theme (theme3.xml)
NAVY   = RGBColor(0x0E, 0x28, 0x41)   # dk2
INK    = RGBColor(0x1A, 0x1A, 0x1A)
MUTED  = RGBColor(0x66, 0x66, 0x66)
TEAL   = RGBColor(0x15, 0x60, 0x82)   # accent1
ORANGE = RGBColor(0xE9, 0x71, 0x32)   # accent2
FOREST = RGBColor(0x19, 0x6B, 0x24)   # accent3
CYAN   = RGBColor(0x0F, 0x9E, 0xD5)   # accent4
PLUM   = RGBColor(0xA0, 0x2B, 0x93)   # accent5
LIME   = RGBColor(0x4E, 0xA7, 0x2E)   # accent6
PILL   = RGBColor(0xE8, 0xE8, 0xE8)   # lt2

# semantic aliases used in slide code below
ACCENT = ORANGE
BLUE   = TEAL
GREEN  = FOREST

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]

# ---------------------------------------------------------------- helpers
def add_headline(slide, text, y=0.32, size=24, color=NAVY):
    tb = slide.shapes.add_textbox(Inches(0.55), Inches(y), Inches(12.2), Inches(0.9))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = text
    r = p.runs[0]; r.font.size = Pt(size); r.font.bold = True; r.font.color.rgb = color
    return tb

def add_text(slide, x, y, w, h, text, size=11, italic=False, bold=False,
             color=INK, align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = ln; p.alignment = align
        for r in p.runs:
            r.font.size = Pt(size); r.font.italic = italic; r.font.bold = bold
            r.font.color.rgb = color
    return tb

def add_picture(slide, path, x, y, w=None, h=None):
    kwargs = dict(left=Inches(x), top=Inches(y))
    if w: kwargs["width"]  = Inches(w)
    if h: kwargs["height"] = Inches(h)
    return slide.shapes.add_picture(str(path), **kwargs)

def add_pill(slide, x, y, w, h, text, fill=PILL, ink=INK, size=10, bold=True):
    sh = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                Inches(x), Inches(y), Inches(w), Inches(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = fill
    tf = sh.text_frame; tf.margin_left = Emu(30000); tf.margin_right = Emu(30000)
    tf.margin_top = Emu(20000); tf.margin_bottom = Emu(20000)
    p = tf.paragraphs[0]; p.text = text; p.alignment = PP_ALIGN.CENTER
    r = p.runs[0]; r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = ink
    return sh

def add_footer(slide, n, total, section=""):
    add_text(slide, 0.55, 7.15, 8, 0.3,
             f"Christina Magyar  ·  Thesis committee, Aug 2026" +
             (f"  ·  {section}" if section else ""),
             size=9, color=MUTED)
    add_text(slide, 12.4, 7.15, 0.6, 0.3, f"{n} / {total}",
             size=9, color=MUTED, align=PP_ALIGN.RIGHT)

def set_notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text

# ---------------------------------------------------------------- SLIDE 1
def slide_title():
    s = prs.slides.add_slide(BLANK)
    # accent bar
    bar = s.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                             Inches(0), Inches(3.35), Inches(13.333), Inches(0.05))
    bar.fill.solid(); bar.fill.fore_color.rgb = ACCENT
    bar.line.fill.background()

    add_text(s, 0.7, 1.7, 12, 1.4,
             "Ancestry-informed analysis of germline genomic contributors\n"
             "to childhood B-ALL clinical outcomes",
             size=30, bold=True, color=NAVY)
    add_text(s, 0.7, 3.7, 12, 0.5,
             "Christina Magyar  ·  GS4 Thesis Committee Update  ·  Aims 1–3",
             size=16, color=INK)
    add_text(s, 0.7, 4.25, 12, 0.4,
             f"{date.today().strftime('%B %Y')}  ·  BCM MSTP · Genetics and Genomics",
             size=12, color=MUTED, italic=True)
    add_text(s, 0.7, 5.2, 12, 0.35,
             "Advisors", size=10, bold=True, color=MUTED)
    add_text(s, 0.7, 5.55, 12, 0.9,
             "Philip J. Lupo, PhD  (Thesis advisor)\n"
             "Elizabeth Atkinson, PhD  (Local advisor)",
             size=13, color=INK)
    set_notes(s,
        "Two-part update today: (1) LAI accuracy benchmarking pilot using an "
        "expanded Mexican Biobank reference panel; (2) status of the 8-module "
        "Nextflow preprocessing pipeline. Numbers shown are the chr20-22 pilot; "
        "full 22-chr sweep is running and will refresh these before the meeting.")
    return s

# ---------------------------------------------------------------- SLIDE 2 (NEW): aims overview
def slide_aims_overview():
    s = prs.slides.add_slide(BLANK)
    add_headline(s, "Three aims — this update focuses on Aim 2.")

    aims = [
        ("Aim 1",
         "Demographic & clinical\npredictors of B-ALL outcomes",
         "REDIAL cohort · MRD-stratified survival · manuscript to Leukemia",
         "Submitted / under revision",
         MUTED, PILL),
        ("Aim 2",
         "Building the preprocessing\n+ LAI pipeline",
         "Nextflow · MXB reference · LAI accuracy benchmarking · manuscript-ready",
         "Focus of today's update",
         ACCENT, RGBColor(0xFF, 0xF1, 0xE4)),
        ("Aim 3",
         "Genetic association study\nof B-ALL outcomes",
         "Trans-ancestry GWAS · Tractor local-ancestry-informed · PolyFun-SuSiE fine-mapping",
         "Next 6-12 months",
         MUTED, PILL),
    ]

    card_w = 4.0; card_h = 4.7; gap = 0.35
    total_w = 3 * card_w + 2 * gap
    x0 = (13.333 - total_w) / 2
    y0 = 1.55

    for i, (aim, title, body, status, border, fill) in enumerate(aims):
        x = x0 + i * (card_w + gap)
        card = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(y0), Inches(card_w), Inches(card_h))
        card.fill.solid(); card.fill.fore_color.rgb = fill
        card.line.color.rgb = border; card.line.width = Pt(2 if aim == "Aim 2" else 0.75)

        add_text(s, x + 0.25, y0 + 0.25, card_w - 0.5, 0.5,
                 aim, size=18, bold=True, color=border)
        add_text(s, x + 0.25, y0 + 0.95, card_w - 0.5, 1.4,
                 title, size=15, bold=True, color=NAVY)
        add_text(s, x + 0.25, y0 + 2.4, card_w - 0.5, 1.7,
                 body, size=11, color=INK)
        add_text(s, x + 0.25, y0 + card_h - 0.7, card_w - 0.5, 0.5,
                 status, size=11, italic=True, bold=True, color=border)

    add_text(s, 0.55, 6.65, 12.2, 0.35,
             "Aim 2 is the technical bridge — the LAI accuracy work "
             "picks the reference panel that Aim 3's local-ancestry GWAS will use.",
             size=12, italic=True, color=INK, align=PP_ALIGN.CENTER)
    add_footer(s, 2, 12, "Aims overview")
    set_notes(s,
        "The three-aim structure hasn't changed since the Spring TAC. Aim 1 is "
        "the REDIAL demographic+clinical outcomes paper (Onwuka/Magyar co-first, "
        "under revision for Leukemia). Aim 2 is what this update focuses on — "
        "the Nextflow preprocessing pipeline plus the LAI accuracy benchmarking "
        "that decides which Amerindigenous reference panel to lock in. Aim 3 is "
        "the downstream GWAS: trans-ancestry SAIGE for the shared variants, "
        "Tractor for local-ancestry-informed effect estimation in the admixed "
        "REDIAL/COG cohorts, PolyFun-SuSiE for fine-mapping. Aim 2 numbers "
        "coming today directly feed Aim 3's panel choice.")
    return s

# ---------------------------------------------------------------- SLIDE 3 (NEW): Aim 1
def slide_aim1_redial():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Aim 1: MRD-negative Latino & NL-Black children still relapse more than NL-White.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "REDIAL cohort · Cox PH adjusted for age, sex, WBC, cytogenetics, CNS, trial",
             size=12, italic=True, color=MUTED)

    # two side-by-side summary boxes
    box_w = 5.9; box_h = 4.3; box_y = 1.75
    for i, (mrd, hr_line1, hr_line2, note) in enumerate([
        ("MRD-negative B-ALL cases",
         "Latino:      HR ~1.4 for disease-free survival vs NL-White",
         "NL-Black:  HR ~1.6 for disease-free survival vs NL-White",
         "Despite the same favorable MRD status, disparity persists."),
        ("MRD-positive B-ALL cases",
         "Latino:      HR consistent with prior COG reports",
         "NL-Black:  HR consistent with prior COG reports",
         "Confirms known MRD-positive disparity in a contemporary cohort."),
    ]):
        x = 0.55 + i * (box_w + 0.4)
        card = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(box_y), Inches(box_w), Inches(box_h))
        card.fill.solid(); card.fill.fore_color.rgb = PILL
        card.line.color.rgb = (ACCENT if i == 0 else TEAL); card.line.width = Pt(1.5)
        add_text(s, x + 0.3, box_y + 0.25, box_w - 0.6, 0.5,
                 mrd, size=15, bold=True, color=(ACCENT if i == 0 else TEAL))
        add_text(s, x + 0.3, box_y + 0.95, box_w - 0.6, 2.2,
                 f"•  {hr_line1}\n•  {hr_line2}",
                 size=12, color=INK)
        add_text(s, x + 0.3, box_y + 3.3, box_w - 0.6, 0.9,
                 note, size=11.5, italic=True, color=MUTED)

    add_text(s, 0.55, 6.3, 12.2, 0.4,
             "Onwuka & Magyar et al. — REDIAL Ethnic and Racial Survival Disparities by MRD "
             "Status — under revision for Leukemia (submission May 2026)",
             size=11, italic=True, color=NAVY, align=PP_ALIGN.CENTER)
    add_footer(s, 3, 12, "Aim 1 · REDIAL")
    set_notes(s,
        "Aim 1 recap. REDIAL analysis found that even among MRD-negative "
        "children — historically the low-risk group — self-identified Latino "
        "and non-Latino Black kids show elevated disease-free-survival hazard "
        "vs non-Latino White. This is important because MRD-negativity is "
        "widely used as a de-escalation threshold in ALL trials, and if the "
        "same MRD-negative status carries different underlying relapse risk "
        "by race/ethnicity, we may be systematically under-treating minority "
        "kids. Motivates Aim 2/3: find the biological substrate. Manuscript "
        "is under revision — HR values shown here are placeholders until the "
        "manuscript numbers are locked; ask me for the current table.\n\n"
        "STATUS: reviewed and updated primary analysis; first draft of "
        "manuscript with figures/tables/supplemental complete; co-author "
        "review; target submission May 2026.")
    return s

# ---------------------------------------------------------------- SLIDE 4 (was 2): motivation
def slide_motivation():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Public Amerindigenous references are ~7× smaller than other continents.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "HGDP + 1000 Genomes reference haplotypes, by continental group",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_amr_disparity.png", x=0.4, y=1.65, w=12.5)
    add_text(s, 0.55, 6.65, 12.2, 0.45,
             "Adding 50 MX Biobank WGS samples brings AMR to 88 — still leaves "
             "Central America, Southern Cone, Amazon, Caribbean underrepresented.",
             size=11, color=INK)
    add_footer(s, 4, 12, "Aim 2 · Motivation")
    set_notes(s,
        "The gap slide. Public HGDP+1KG panels give us ~620 EUR haps, ~634 AFR, "
        "~667 EAS, but only ~88 AMR — a 7× disparity. The MX Biobank donation "
        "of 50 WGS Mexican samples nearly doubles the AMR pool to 138, but the "
        "geographic composition (pie) shows we are still Mexico-heavy and missing "
        "Central America, the Southern Cone, the eastern Amazon, and the "
        "Caribbean. This is why panel-choice matters and why we're benchmarking "
        "before we lock a panel down for the REDIAL/COG cohorts.")
    return s

# ---------------------------------------------------------------- SLIDE 3
def slide_methods():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "We simulated admixed cohorts to measure per-hap LAI recall.")

    # 4-step flow across the middle
    steps = [
        ("Reference haps",
         "HGDP + 1KG\n+ 50 MXB WGS\n(6 panel combos)",
         PILL),
        ("Phasing",
         "SHAPEIT5\npopulation phasing",
         PILL),
        ("Simulate admixed",
         "admix-simu\ngen=12, n=30\nBrasa & MXB tracks",
         PILL),
        ("Local ancestry",
         "RFMix v1 TrioPhased\n22 chromosomes",
         PILL),
        ("Score",
         "per-hap recall\nvs known truth\nweighted concordance / F1",
         PILL),
    ]
    n = len(steps); w = 2.35; gap = 0.15
    x0 = (13.333 - (n*w + (n-1)*gap)) / 2
    y0 = 2.1
    for i, (head, body, fill) in enumerate(steps):
        x = x0 + i*(w+gap)
        card = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(y0), Inches(w), Inches(2.1))
        card.fill.solid(); card.fill.fore_color.rgb = fill
        card.line.color.rgb = MUTED; card.line.width = Pt(0.5)
        tf = card.text_frame; tf.word_wrap = True
        tf.margin_top = Emu(80000); tf.margin_left = Emu(60000); tf.margin_right = Emu(60000)
        p = tf.paragraphs[0]; p.text = head; p.alignment = PP_ALIGN.CENTER
        r = p.runs[0]; r.font.size = Pt(13); r.font.bold = True; r.font.color.rgb = NAVY
        p2 = tf.add_paragraph(); p2.text = body; p2.alignment = PP_ALIGN.CENTER
        for rr in p2.runs:
            rr.font.size = Pt(10); rr.font.color.rgb = INK
        if i < n-1:
            arr = s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW,
                                     Inches(x+w+0.005), Inches(y0+0.92),
                                     Inches(gap-0.01), Inches(0.28))
            arr.fill.solid(); arr.fill.fore_color.rgb = MUTED
            arr.line.fill.background()

    # 6 panels tested
    add_text(s, 0.55, 4.55, 12.2, 0.35,
             "Six reference-panel combinations tested",
             size=13, bold=True, color=NAVY)
    panels = [
        ("NAT_HGDP",         "26 HGDP + 5 LP",           MUTED),
        ("NAT_HGDPMXB",      "+ 25 MXB",                 ACCENT),
        ("NAT_HGDPMXB_FULL", "+ 50 MXB",                 ACCENT),
        ("NAT_HOMOG",        "ADMIXTURE-homog Q≥0.95",   BLUE),
        ("NAT_PEL",          "adds 9 PEL (Peru)",        GREEN),
        ("NAT_PEL_EAS",      "PEL + EAS outgroup",       GREEN),
    ]
    px = 0.55; py = 5.0; pw = 2.05; ph = 1.05
    for i, (name, sub, col) in enumerate(panels):
        r = i // 3; c = i % 3
        x = px + c*(pw + 0.15) + r*6.55*0
        # 2-row grid
        col_i = i % 3; row_i = i // 3
        x = px + col_i*(pw + 0.15)
        y = py + row_i*(ph + 0.05)
        # wider layout: 6 cards in 2 rows of 3
        pill = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(y), Inches(pw), Inches(ph))
        pill.fill.solid(); pill.fill.fore_color.rgb = PILL
        pill.line.color.rgb = col; pill.line.width = Pt(1.4)
        tf = pill.text_frame
        tf.margin_top = Emu(50000); tf.margin_left = Emu(80000); tf.margin_right = Emu(60000)
        p = tf.paragraphs[0]; p.text = name; p.alignment = PP_ALIGN.LEFT
        for rr in p.runs:
            rr.font.size = Pt(11); rr.font.bold = True; rr.font.color.rgb = col
        p2 = tf.add_paragraph(); p2.text = sub; p2.alignment = PP_ALIGN.LEFT
        for rr in p2.runs:
            rr.font.size = Pt(9.5); rr.font.color.rgb = INK

    # method callouts on right
    add_text(s, 7.15, 5.0, 5.7, 0.35,
             "Two simulated admixed cohorts",
             size=13, bold=True, color=NAVY)
    add_text(s, 7.15, 5.4, 5.7, 1.65,
             "•  Brasa track — 15% NAT / 60% EUR / 25% AFR\n"
             "     (Brazilian-like ancestry proportions)\n"
             "•  MXB track — Mexican-like ancestry proportions\n"
             "     (tests whether MXB helps a Mexican cohort)",
             size=11, color=INK)

    add_footer(s, 5, 12, "Aim 2 · Methods")
    set_notes(s,
        "Simulation harness: SHAPEIT5-phased references feed admix-simu to create "
        "30 admixed individuals per cohort at generation 12 (long enough for "
        "chromosome-scale tracts to be resolvable). RFMix v1 TrioPhased mode "
        "matches Honorato-Mauer 2025. Recall = per-haplotype per-ancestry "
        "concordance against the admix-simu-provided truth. We ran 6 candidate "
        "reference panels crossed with 2 simulated cohort tracks (Brazilian-like "
        "and Mexican-like) to see whether the best panel depends on the target "
        "cohort's own ancestry composition.")
    return s

# ---------------------------------------------------------------- SLIDE 4
def slide_ancestry_overview():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Panel choice moves AMR recall by ~3%; EUR and AFR sit near ceiling.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "Weighted per-hap recall by ancestry, WGS density, chr20-22 pilot",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_wgs_all_ancestry.png", x=0.55, y=1.55, w=12.2)
    add_text(s, 0.55, 6.75, 12.2, 0.35,
             "→  AMR is where reference-panel choice actually matters — "
             "so the rest of the deck focuses on AMR recall.",
             size=12, italic=True, color=INK)
    add_footer(s, 6, 12, "Aim 2 · Result overview")
    set_notes(s,
        "Overview slide before we zoom in. Three ancestries side by side, same "
        "panels, same y-scale. EUR and AFR sit near ceiling regardless of which "
        "AMR reference we pick — those calls are easy because the EUR and AFR "
        "reference pools are already large and diverse. AMR is the ancestry that "
        "moves: ~3% swing across panel choice on the Mexican-like cohort. "
        "Motivates zooming in on AMR in the next slide.")
    return s

# ---------------------------------------------------------------- SLIDE 5
def slide_wgs_result():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Best reference panel depends on the target cohort's composition.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "Weighted NAT-ancestry recall, WGS density, chr20-22 pilot",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_wgs_panels.png", x=0.55, y=1.65, w=8.2)

    # takeaway callouts on right
    add_text(s, 9.05, 1.75, 3.9, 0.4,
             "Takeaways", size=13, bold=True, color=NAVY)
    add_text(s, 9.05, 2.15, 3.95, 4.4,
             "For a Mexican-like cohort (blue bars):\n"
             "•  Adding 25 MXB samples\n    lifts NAT recall  0.931 → 0.942\n\n"
             "For a Brazilian-like cohort (red bars):\n"
             "•  PEL (Peru) panel wins on NAT\n    HGDPMXB actually loses ~2%\n"
             "•  Adding non-Mexican ancestry\n    to a Mexican reference hurts\n\n"
             "→  Panel-choice must match the\n     target cohort, not the largest\n     available sample source.",
             size=10.5, color=INK)
    add_footer(s, 7, 12, "Aim 2 · Result · AMR zoom")
    set_notes(s,
        "HEADLINE RESULT. Extended pilot numbers (chr 1-8 + 20-22 for the 4 full-"
        "coverage panels; chr20-22 only for HOMOG and PEL_EAS — flagged with *). "
        "On the Mexican-like cohort (NATMXB track, teal): +25 MXB lifts NAT "
        "recall from 0.954 → 0.958 and F1 from 0.969 → 0.974 (cleanest single "
        "winner). On the Brazilian-like cohort (NAT track, orange): +25 MXB "
        "costs ~1% NAT recall (0.930 → 0.920) but bumps precision (0.981 → "
        "0.990), so F1 is roughly flat. +50 MXB (HGDPMXB_FULL) costs more "
        "recall (0.902) with only marginal F1 benefit — 25 MXB is the sweet "
        "spot. PEL matches or beats HGDP on the Brazilian cohort. NAT_HOMOG "
        "(Q≥0.95 filter, donor-safe) performs comparably to HGDPMXB in the "
        "3-chr subset we have — worth re-running with full coverage. Punchline: "
        "adding MXB is a clear win for a Mexican-cohort GWAS; matched-source "
        "reference (PEL for Brazilian) matters more than panel size.\n\n"
        "──────── METRICS CHEAT SHEET ────────\n"
        "• TP / FP / FN: sites where called ancestry matches truth (TP), calls "
        "were wrong-class (FP), or truth was this class but we missed (FN). "
        "TN is ignored — this is a small-positive-class setting.\n"
        "• Recall = Sensitivity = TPR = TP / (TP + FN). Of all true NAT sites, "
        "how many we recovered. Currently plotted on the hero bars.\n"
        "• Precision = PPV = TP / (TP + FP). Of the sites we CALLED NAT, how "
        "many were actually NAT. Guards against over-calling.\n"
        "• F1 = 2·P·R / (P + R). Harmonic mean of precision and recall — the "
        "class-imbalance-robust summary. NOT the same as balanced accuracy "
        "(which uses specificity); F1 ignores TN by design.\n"
        "• Weighted (in these tables): concordance / recall / precision are "
        "weighted by (n_haps × n_sites) per chromosome so long chromosomes "
        "and larger haplotype counts contribute more.\n"
        "• Concordance in these tables (== weighted_recall in current R impl) "
        "is per-site agreement to truth, per ancestry, one-vs-rest.\n"
        "• SE: standard error of the per-chromosome mean (not a bootstrap CI).\n"
        "Publication table with all three metrics is committed to "
        "results/tables/publication_table_amr.tsv and rendered to "
        "slides/tables/publication_table_amr.png.")
    return s

# ---------------------------------------------------------------- SLIDE 5
def slide_chip_result():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Chip-density LAI holds within ~1% of WGS on a Mexican-like cohort.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "Illumina GSA v3 chip vs full WGS density, same reference panels",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_chip_vs_wgs.png", x=0.55, y=1.65, w=8.4)

    add_text(s, 9.15, 1.75, 3.9, 0.4,
             "Implications for GWAS", size=13, bold=True, color=NAVY)
    add_text(s, 9.15, 2.15, 3.95, 4.4,
             "•  REDIAL / COG cohorts are\n    genotyped on chip, not WGS —\n"
             "    LAI accuracy holds if we\n    use the MXB-augmented panel\n\n"
             "•  Chip → TOPMed imputation is\n    the next validation step\n    (Michigan / 1KG post-meeting)\n\n"
             "•  Panel ranking preserved\n    across densities:\n"
             "     HGDPMXB > HGDP > PEL_EAS\n\n"
             "•  Recommend HGDPMXB for the\n    Mexican-Latino GWAS pipeline",
             size=10.5, color=INK)
    add_footer(s, 8, 12, "Aim 2 · Result · Chip")
    set_notes(s,
        "Complementary result. Same reference panels, but the admixed cohort's "
        "haplotypes were down-sampled to the ~385k autosomal Illumina GSA v3 "
        "positions before calling. NAT recall drops by only ~0.7-1.2% across "
        "panels — the panel ranking is preserved and the direction of every "
        "effect matches WGS. This matters because our REDIAL/COG cohorts are "
        "chip-genotyped, not WGS — so the WGS result generalizes. Michigan/1KG "
        "imputation is the post-committee validation step; today we can already "
        "recommend HGDPMXB (+25 MXB) for the Mexican-Latino LAI pipeline. "
        "NAT_HOMOG is a donor-safe backup that performs equivalently.")
    return s

# ---------------------------------------------------------------- SLIDE 6
def slide_pipeline():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "8-module Nextflow preprocessing pipeline built; benchmarking underway.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "Containerized · parameterized · parallelized across chromosomes and platforms",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_pipeline.png", x=0.4, y=1.55, w=12.5)

    # status pills bottom
    add_text(s, 0.55, 5.55, 12.2, 0.35,
             "Status since Spring TAC", size=13, bold=True, color=NAVY)
    status = [
        ("BUILT",      "9-module Nextflow pipeline",                       GREEN),
        ("BUILT",      "Pilot ran on TOPMed test data (end of Module 8)",  GREEN),
        ("IN PROG.",   "Michigan & AoU AnVIL patched in (Q3 2026)",        BLUE),
        ("IN PROG.",   "Full benchmarking sweep (target Apr 2026 done)",   BLUE),
        ("NEXT",       "Ancestry-stratified GWAS on REDIAL/COG",           ACCENT),
    ]
    x = 0.55; y = 5.95; row_h = 0.32
    for tag, txt, col in status:
        pill = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(y), Inches(0.95), Inches(row_h-0.02))
        pill.fill.solid(); pill.fill.fore_color.rgb = col
        pill.line.fill.background()
        tf = pill.text_frame
        tf.margin_top = Emu(0); tf.margin_bottom = Emu(0)
        p = tf.paragraphs[0]; p.text = tag; p.alignment = PP_ALIGN.CENTER
        for rr in p.runs:
            rr.font.size = Pt(9); rr.font.bold = True
            rr.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        add_text(s, x + 1.05, y + 0.02, 10.5, row_h,
                 txt, size=11, color=INK)
        y += row_h

    add_footer(s, 9, 12, "Aim 2 · Pipeline update")
    set_notes(s,
        "Two-part update: what shipped since the Spring TAC and what's queued. "
        "The 8 modules (colleague-led) now run end-to-end on test data with the "
        "TOPMed imputation server — Module 8 benchmarking plots render. Michigan "
        "and AoU AnVIL back-ends are patched in but not yet run on the full "
        "benchmark set. Full multi-platform benchmarking is expected to complete "
        "before the fall term. Downstream: once Modules 6-7 hand off a QC'd, "
        "ancestry-labeled cohort, the LAI accuracy work presented today directly "
        "feeds panel selection for Module 7's LAI step, and the ancestry-"
        "stratified GWAS begins on REDIAL/COG. The dashed red box on the diagram "
        "marks exactly the modules this thesis has empirically validated.")
    return s

# ---------------------------------------------------------------- SLIDE 10 (NEW): Aim 3 methods
def slide_aim3_methods():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Aim 3: trans-ancestry + local-ancestry-informed GWAS on ~5,400 B-ALL cases.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "Two parallel tracks: ancestry-stratified meta-analysis and Tractor local-ancestry GWAS",
             size=12, italic=True, color=MUTED)

    # left branch: trans-ancestry
    left_x = 0.55; right_x = 6.9; branch_w = 5.85; branch_y = 1.7; branch_h = 4.5
    for x, title, color, steps, note in [
        (left_x, "Aim 3.a  ·  Trans-ancestry", TEAL,
         ["1.  GRAF-anc groups →  stratified SAIGE GWAS  (4 groups)",
          "2.  Trans-ancestry meta-analysis  (n ≈ 5,357)",
          "3.  Fine-mapping with MGflashfm  (multi-trait, multi-ethnic)"],
         "Uses shared trans-ancestry effects across REDIAL + COG + St. Jude"),
        (right_x, "Aim 3.b  ·  Local-ancestry-informed", ORANGE,
         ["1.  RFMix v1 on 1,660 Latin American individuals",
          "2.  Tractor GWAS  →  per-ancestry β for AMR / EUR / AFR",
          "3.  Functionally-prioritized fine-mapping (SuSiE + PolyFun)"],
         "Detects Amerindigenous-enriched risk variants missed by global-ancestry GWAS"),
    ]:
        card = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(x), Inches(branch_y), Inches(branch_w), Inches(branch_h))
        card.fill.solid(); card.fill.fore_color.rgb = PILL
        card.line.color.rgb = color; card.line.width = Pt(1.5)
        add_text(s, x + 0.25, branch_y + 0.2, branch_w - 0.5, 0.5,
                 title, size=15, bold=True, color=color)
        add_text(s, x + 0.25, branch_y + 0.9, branch_w - 0.5, 2.8,
                 "\n\n".join(steps), size=12, color=INK)
        add_text(s, x + 0.25, branch_y + branch_h - 0.7, branch_w - 0.5, 0.5,
                 note, size=10.5, italic=True, color=MUTED)

    add_text(s, 0.55, 6.4, 12.2, 0.35,
             "Aim 2's LAI-accuracy work directly picks the reference panel used by Aim 3.b's RFMix step.",
             size=11.5, italic=True, color=INK, align=PP_ALIGN.CENTER)
    add_footer(s, 10, 12, "Aim 3 · methods")
    set_notes(s,
        "Aim 3 splits into two GWAS tracks. 3.a is the trans-ancestry track: "
        "GRAF-anc first partitions the combined REDIAL + COG + St. Jude "
        "cohort into 4 ancestry groups (African-American, European-American, "
        "Latin American 1, Latin American 2), stratified SAIGE GWAS per "
        "group, then trans-ancestry meta-analysis at n≈5,357, followed by "
        "MGflashfm fine-mapping which handles multi-trait multi-ethnic LD "
        "differences. 3.b is the local-ancestry track: RFMix v1 (the exact "
        "tool this Aim 2 work is benchmarking) on the 1,660 self-reported "
        "Latino individuals, then Tractor to fit per-ancestry β for AMR, "
        "EUR, and AFR tracts jointly. Fine-mapping downstream uses SuSiE + "
        "PolyFun. The point of Aim 2 is to make sure the RFMix step in 3.b "
        "is using the best reference panel — HGDPMXB by our current data.")
    return s

# ---------------------------------------------------------------- SLIDE 11 (NEW): timeline
def slide_timeline():
    s = prs.slides.add_slide(BLANK)
    add_headline(s, "Timeline to graduation — July 2027 target.")

    # simple Gantt-like row list
    milestones = [
        ("Aim 1", "REDIAL manuscript revision + submission (Leukemia)",  "Apr–Jun 2026",  ACCENT),
        ("Aim 2", "22-chr LAI sweep + full-coverage figures",             "Aug 2026",      TEAL),
        ("Aim 2", "Nextflow pipeline benchmarking on TOPMed / Michigan / AoU", "Sep–Dec 2026", TEAL),
        ("Aim 2", "Pipeline manuscript drafting",                         "Jan–Mar 2027",  TEAL),
        ("Aim 3", "GRAF-anc stratification + PC calculation",             "Oct–Nov 2026",  ORANGE),
        ("Aim 3", "SAIGE trans-ancestry GWAS + MGflashfm fine-mapping",   "Dec 2026 – Feb 2027", ORANGE),
        ("Aim 3", "RFMix + Tractor local-ancestry GWAS on REDIAL/COG",    "Feb – Apr 2027", ORANGE),
        ("Aim 3", "PolyFun-SuSiE fine-mapping + candidate follow-up",     "Apr – May 2027", ORANGE),
        ("All",   "Thesis writing + defense",                             "May – Jul 2027", NAVY),
    ]

    y = 1.55; row_h = 0.52
    for aim, task, when, color in milestones:
        pill = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(0.55), Inches(y),
                                  Inches(0.95), Inches(row_h - 0.07))
        pill.fill.solid(); pill.fill.fore_color.rgb = color; pill.line.fill.background()
        tf = pill.text_frame
        tf.margin_top = Emu(0); tf.margin_bottom = Emu(0)
        p = tf.paragraphs[0]; p.text = aim; p.alignment = PP_ALIGN.CENTER
        for rr in p.runs:
            rr.font.size = Pt(11); rr.font.bold = True
            rr.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        add_text(s, 1.65, y + 0.05, 7.6, row_h,
                 task, size=12, color=INK)
        add_text(s, 9.35, y + 0.05, 3.6, row_h,
                 when, size=12, bold=True, color=color, align=PP_ALIGN.RIGHT)
        y += row_h

    add_footer(s, 11, 12, "Timeline")
    set_notes(s,
        "Timeline broken down by aim. Aim 1 wraps this quarter with the "
        "REDIAL manuscript. Aim 2 finishes the 22-chr LAI sweep this month "
        "(currently 126/198 done on the WGS side, backfill queued for HOMOG "
        "/ PEL_EAS panels and full chip-density coverage), then benchmarks "
        "the Nextflow pipeline on all three imputation servers (TOPMed, "
        "Michigan, AoU AnVIL) through end of 2026. Pipeline manuscript "
        "targets Q1 2027. Aim 3 sequences behind Aim 2 — GRAF-anc + PCs in "
        "the fall, SAIGE + fine-mapping over winter, then Tractor local-"
        "ancestry GWAS in spring 2027 (which needs Aim 2's reference panel "
        "decision locked in). Thesis writing spring 2027, defense July 2027.")
    return s

# ---------------------------------------------------------------- SLIDE 12 (NEW): acknowledgements
def slide_acks():
    s = prs.slides.add_slide(BLANK)
    add_headline(s, "Acknowledgements")

    cols = [
        ("Lupo Lab & TCH EpiCenter",
         "Dr. Philip J. Lupo  (Emory, thesis advisor)\n"
         "Dr. Karen Rabin  (UCSF)\n"
         "Dr. Melissa A. Richard\n"
         "Dr. Austin Brown\n"
         "Dr. Jeremy Schraw\n"
         "Dr. Michael Scheurer"),
        ("Atkinson Lab  (BCM · local advisor)",
         "Dr. Elizabeth Atkinson  (local advisor)\n"
         "Nirav Shah  ·  Jessica Honorato-Mauer\n"
         "Grace Tietz  ·  Hatoon Al Ali\n"
         "Helen Lin  ·  Pragati Kore\n"
         "Erik Stricker  ·  Aishi Ayyanathan\n"
         "Astrid Manuel  ·  Shalini Dhamodharan"),
        ("Collaborators + funding",
         "St. Jude Children's Research Hospital\n"
         "  Dr. Jun Yang, Zenhua Li\n"
         "Children's Oncology Group\n"
         "MX Biobank  (Dr. Andrés Moreno-Estrada)\n"
         "\n"
         "AIM-AHEAD Consortium\n"
         "Robert & Janice McNair Foundation\n"
         "BCM MSTP · G&G Graduate Program"),
    ]

    col_w = 4.15; col_h = 5.0; gap = 0.15
    x0 = (13.333 - (3 * col_w + 2 * gap)) / 2
    y0 = 1.55
    for i, (hd, body) in enumerate(cols):
        x = x0 + i * (col_w + gap)
        add_text(s, x, y0, col_w, 0.5,
                 hd, size=13, bold=True, color=NAVY)
        add_text(s, x, y0 + 0.6, col_w, col_h,
                 body, size=11, color=INK)

    add_text(s, 0.55, 6.85, 12.2, 0.35,
             "Thank you — questions?",
             size=16, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
    add_footer(s, 12, 12, "Acknowledgements")
    set_notes(s,
        "Standard acks — matches the Spring TAC deck acknowledgements slide. "
        "Ready for questions.")
    return s

# ---------------------------------------------------------------- build
if __name__ == "__main__":
    slide_title()                # 1
    slide_aims_overview()        # 2  (new)
    slide_aim1_redial()          # 3  (new)
    slide_motivation()           # 4
    slide_methods()              # 5
    slide_ancestry_overview()    # 6
    slide_wgs_result()           # 7  AMR zoom
    slide_chip_result()          # 8
    slide_pipeline()             # 9
    slide_aim3_methods()         # 10 (new)
    slide_timeline()             # 11 (new)
    slide_acks()                 # 12 (new)
    prs.save(OUT)
    n_slides = len(prs.slides)
    sz = OUT.stat().st_size
    print(f"wrote {OUT} ({sz:,} bytes, {n_slides} slides)")
