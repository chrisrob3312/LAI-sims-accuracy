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

# palette lifted from the Spring TAC deck feel
NAVY   = RGBColor(0x1F, 0x3A, 0x5F)
INK    = RGBColor(0x1A, 0x1A, 0x1A)
MUTED  = RGBColor(0x66, 0x66, 0x66)
ACCENT = RGBColor(0xC9, 0x4A, 0x53)   # MXB red
BLUE   = RGBColor(0x3E, 0x6B, 0xB0)
GREEN  = RGBColor(0x4C, 0x9F, 0x70)
PILL   = RGBColor(0xEE, 0xE9, 0xE0)   # cream pill for tags

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
             "Building an ancestry-informed genotype-preprocessing pipeline\n"
             "and benchmarking local ancestry inference for Latino cohorts",
             size=30, bold=True, color=NAVY)
    add_text(s, 0.7, 3.7, 12, 0.5,
             "Christina Magyar  ·  GS4 Thesis Committee Update",
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

# ---------------------------------------------------------------- SLIDE 2
def slide_motivation():
    s = prs.slides.add_slide(BLANK)
    add_headline(s,
        "Public Amerindigenous references are ~7× smaller than other continents.")
    add_text(s, 0.55, 1.15, 12.2, 0.4,
             "HGDP + 1000 Genomes reference haplotypes, by continental group",
             size=12, italic=True, color=MUTED)
    add_picture(s, FIGDIR / "fig_amr_disparity.png", x=0.4, y=1.65, w=12.5)
    add_text(s, 0.55, 6.65, 12.2, 0.45,
             "Adding 50 MX Biobank WGS samples raises AMR to 138 haps — still leaves "
             "Central America, Southern Cone, Amazon, Caribbean underrepresented.",
             size=11, color=INK)
    add_footer(s, 2, 6, "Motivation")
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

    add_footer(s, 3, 6, "Methods")
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
    add_footer(s, 4, 6, "Result 1 · WGS")
    set_notes(s,
        "This is the headline result. On the Mexican-like cohort (NATMXB track, "
        "blue), adding 25 MXB samples to the HGDP baseline lifts NAT recall from "
        "0.931 to 0.942. Adding 50 gets the same benefit (0.942) — 25 is enough. "
        "On the Brazilian-like cohort (NAT track, red), the picture flips: the "
        "MXB-augmented panel actually LOSES ~2% NAT recall because Brazilian "
        "NAT ancestry is closer to Amazonian sources than to Mexican; PEL (Peru) "
        "is the best panel there. The homogeneous-only NAT_HOMOG panel matches "
        "HGDPMXB on the Mexican cohort while being donor-safe. Punchline for "
        "the committee: there is no single best panel — the answer depends on "
        "the ancestry composition of the cohort being called.")
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
    add_footer(s, 5, 6, "Result 2 · Chip")
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

    add_footer(s, 6, 6, "Pipeline update")
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

# ---------------------------------------------------------------- build
if __name__ == "__main__":
    slide_title()
    slide_motivation()
    slide_methods()
    slide_wgs_result()
    slide_chip_result()
    slide_pipeline()
    prs.save(OUT)
    n_slides = len(prs.slides)
    sz = OUT.stat().st_size
    print(f"wrote {OUT} ({sz:,} bytes, {n_slides} slides)")
