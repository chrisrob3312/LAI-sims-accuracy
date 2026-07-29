#!/usr/bin/env Rscript
# =============================================================================
# accuracy_v2.R -- LAI accuracy evaluation for the (panel x track x chr) grid
#
# Reads:
#   Truth (2b):   ${SIM_DIR}/${track}.Brasa${chr}.hanc              (per-site x per-hap ancestry, 0/1/2)
#                 ${SIM_DIR}/1kg_hgdp.chr${chr}.snp                 (admix-simu .snp: snp_id chr cM phys_pos ref alt)
#                 ${SIM_DIR}/${track}.Brasa.chr${chr}.sample        (30 admixed sample IDs)
#   RFMix (3b):   ${WORKDIR}/${track}.${panel}.gen${GEN}_chr${chr}.recoded    (per-site x per-hap, 0/1/2)
#                 ${WORKDIR}/${track}.${panel}.gen${GEN}_chr${chr}_chr${chr}.map (phys_pos gen_pos snp_id)
#                 ${WORKDIR}/${track}.${panel}.gen${GEN}_chr${chr}.sample      (per-hap sample IDs)
#
# Emits under ${OUTDIR}:
#   accuracy_long.tsv         one row per (track, panel, chr, hap, ancestry)
#   accuracy_summary.tsv      mean/SE per (track, panel, ancestry) genome-wide
#   accuracy_per_chr.tsv      mean per (track, panel, ancestry, chr) for the chr line plots
#   panel_vs_baseline.tsv     paired Wilcoxon comparisons vs HGDP-only baseline
#   plots/                    PDFs for each figure below
#
# Ancestry class encoding (matches 3b's Viterbi recoding sed -e 's/1/0/g' ...):
#   0 = NAT, 1 = EUR, 2 = AFR
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(purrr)
})

# --- Config --------------------------------------------------------------
PROJECT_ROOT   <- Sys.getenv("PROJECT_ROOT",
                             "/storage/atkinson/home/magyar/Projects/01_REDIAL_Projects/01_LAI_Accuracy_MXBiobank")
ADMIX_POP      <- Sys.getenv("ADMIX_POP", "Brasa")
GEN            <- Sys.getenv("GEN", "12")
SIM_DIR        <- file.path(PROJECT_ROOT, "02_simulations", ADMIX_POP, paste0("gen", GEN))
RFMIX_ROOT     <- Sys.getenv("RFMIX_ROOT", file.path(PROJECT_ROOT, "03_rfmix"))
WORKDIR        <- Sys.getenv("WORKDIR", RFMIX_ROOT)  # 3b's per-(pop, gen) dir
OUTDIR         <- Sys.getenv("OUTDIR", file.path(PROJECT_ROOT, "04_accuracy"))
# When OUTDIR points inside a git repo (e.g. results/accuracy_pilot), the
# .tsv summary and PDF plots become committable review artifacts. Big files
# (accuracy_long.tsv) still stay under OUTDIR only.
BASELINE_PANEL <- Sys.getenv("BASELINE_PANEL", "NAT_HGDP")

CHRS   <- 1:22
TRACKS <- c("NAT", "NATMXB")
PANELS <- strsplit(Sys.getenv("PANELS_TO_RUN",
                              "NAT_HGDP NAT_PEL NAT_PEL_EAS NAT_HGDPMXB NAT_HGDPMXB_FULL"),
                   " +")[[1]]

# Invalid combo per PIPELINE.md: NATMXB sim + HGDPMXB_FULL panel double-dips donors.
is_valid_combo <- function(track, panel) !(track == "NATMXB" && panel == "NAT_HGDPMXB_FULL")

ANC_LABELS <- c(`0` = "NAT", `1` = "EUR", `2` = "AFR")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "plots"), showWarnings = FALSE, recursive = TRUE)

# --- Readers -------------------------------------------------------------
# Truth .hanc: 2b writes ONE LINE PER HAPLOTYPE with one character per SNP
# (no separators). So an n_haps x n_sites matrix — we transpose to
# n_sites x n_haps to match RFMix's row=SNP orientation.
read_hanc <- function(path) {
  L <- readLines(path)
  if (length(L) == 0) stop("empty hanc: ", path)
  m <- do.call(rbind, lapply(strsplit(L, "", fixed = TRUE), as.integer))
  storage.mode(m) <- "integer"
  t(m)   # -> n_sites x n_haps
}

# Truth positions: admix-simu .snp is 5-col "snp_id chr cM phys_pos allele".
# (2b's insert-map.pl output layout; column 4 is the physical hg38 position.)
read_truth_pos <- function(chr) {
  f <- file.path(SIM_DIR, sprintf("1kg_hgdp.chr%d.snp", chr))
  d <- read.table(f, header = FALSE, sep = "", stringsAsFactors = FALSE,
                  colClasses = c("character","character","numeric","integer","character","character")[1:6])
  as.integer(d$V4)
}

# RFMix ancestry + positions from a .Lat3 file. 3b deletes the intermediate
# .recoded, so read directly from .Lat3 which is guaranteed to exist for any
# combo whose .done marker was written.
#
# .Lat3 format (from 3b line 401): each row is a SNP site, space-separated:
#   col 1  = physical position (hg38, integer)
#   col 2  = snp_id (usually chr:pos)
#   col 3+ = per-haplotype recoded ancestry (0/1/2)
read_lat3 <- function(path) {
  d <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE,
                  colClasses = "character")
  positions <- as.integer(d[[1]])
  m <- as.matrix(sapply(d[, -c(1, 2), drop = FALSE], as.integer))
  storage.mode(m) <- "integer"
  list(pos = positions, mat = m)
}

# 2b doesn't write per-track admixed sample IDs, so synth "sim_<h>" hap labels.
# We keep the arg signature so callers stay stable.
read_admixed_samples <- function(chr, track, n_admixed) {
  sprintf("sim_%03d", seq_len(n_admixed / 2))
}

# .classes has one integer per haplotype in RFMix's column order:
#   ref haps first (class 1, 2, 3), admixed haps last (class 0).
# Returns the column indices INTO the .Lat3 hap columns (which start at 3 in
# the raw .Lat3 file) where class == 0.
read_admixed_hap_cols <- function(chr, track, panel) {
  f <- file.path(WORKDIR, sprintf("%s.%s.gen%s_chr%d.classes",
                                  track, panel, GEN, chr))
  cls <- scan(f, what = integer(), quiet = TRUE)
  which(cls == 0L)
}

# --- Core: per-hap per-ancestry accuracy for one (track, panel, chr) ----
score_one <- function(track, panel, chr) {
  hanc_f <- file.path(SIM_DIR, sprintf("%s.%s%d.hanc",   track, ADMIX_POP, chr))
  lat3_f <- file.path(WORKDIR, sprintf("%s.%s.gen%s_chr%d.Lat3",
                                       track, panel, GEN, chr))
  done_f <- file.path(WORKDIR, sprintf("%s.%s.gen%s_chr%d.done",
                                       track, panel, GEN, chr))
  # Silently skip anything without a .done marker so partial runs (e.g. 55/198)
  # produce clean output for whatever finished.
  if (!file.exists(hanc_f) || !file.exists(lat3_f) || !file.exists(done_f))
    return(NULL)

  truth <- read_hanc(hanc_f)                 # n_sites_truth x n_haps
  tpos  <- read_truth_pos(chr)               # length n_sites_truth
  if (length(tpos) != nrow(truth))
    stop(sprintf("[chr%d] truth pos vs hanc rows mismatch: %d vs %d",
                 chr, length(tpos), nrow(truth)))

  lat3 <- read_lat3(lat3_f)
  rec_all <- lat3$mat                        # n_sites x n_ALL_haps (ref + admixed)
  rpos    <- lat3$pos

  # RFMix writes ancestry for every hap (ref + admixed) in Viterbi.txt, which
  # is what 3b pastes into .Lat3. Only the class-0 columns (admixed) match
  # truth.
  adm_cols <- read_admixed_hap_cols(chr, track, panel)
  rec <- rec_all[, adm_cols, drop = FALSE]
  if (ncol(truth) != ncol(rec))
    stop(sprintf("[%s/%s chr%d] admixed-hap count mismatch: truth %d vs .classes 0-count %d (n_all=%d)",
                 track, panel, chr, ncol(truth), ncol(rec), ncol(rec_all)))

  # Intersect on physical position (both are 1-based, hg38).
  shared <- intersect(tpos, rpos)
  if (length(shared) < 100) {
    warning(sprintf("[%s/%s chr%d] only %d shared sites; skipping",
                    track, panel, chr, length(shared)))
    return(NULL)
  }
  ti <- match(shared, tpos)                  # rows into truth
  ri <- match(shared, rpos)                  # rows into rec
  T  <- truth[ti, , drop = FALSE]
  R  <- rec  [ri, , drop = FALSE]

  n_haps <- ncol(T)
  samps  <- read_admixed_samples(chr, track, n_haps)
  # Each simulated diploid contributes 2 haps: sample name replicated once per hap.
  hap_sample <- rep(samps, each = 2)[seq_len(n_haps)]

  # Long-form per hap x per ancestry class.
  # For each class a (only sites where truth == a): concordance = mean(R == T)
  out <- vector("list", n_haps * length(ANC_LABELS))
  k <- 0L
  for (h in seq_len(n_haps)) {
    tv <- T[, h]; rv <- R[, h]
    for (a in as.integer(names(ANC_LABELS))) {
      mask <- (tv == a)
      n_a  <- sum(mask)
      if (n_a == 0L) next
      hits <- sum(rv[mask] == a)
      k <- k + 1L
      out[[k]] <- data.frame(
        track     = track,
        panel     = panel,
        chr       = chr,
        hap       = h,
        sample    = hap_sample[h],
        ancestry  = ANC_LABELS[[as.character(a)]],
        n_sites   = as.integer(n_a),
        n_correct = as.integer(hits),
        concordance = hits / n_a,
        stringsAsFactors = FALSE)
    }
  }
  bind_rows(out[seq_len(k)])
}

# --- Main sweep ----------------------------------------------------------
message("[accuracy_v2] sweeping (track x panel x chr) grid")
grid <- expand.grid(track = TRACKS, panel = PANELS, chr = CHRS,
                    stringsAsFactors = FALSE) |>
  filter(mapply(is_valid_combo, track, panel))

results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  message(sprintf("  [%d/%d] %s / %s chr%d", i, nrow(grid), g$track, g$panel, g$chr))
  results[[i]] <- tryCatch(score_one(g$track, g$panel, g$chr),
                           error = function(e) { warning(e$message); NULL })
}
long <- bind_rows(results)
if (nrow(long) == 0) stop("No accuracy rows produced -- check paths and file existence.")

write_tsv(long, file.path(OUTDIR, "accuracy_long.tsv"))
message(sprintf("[accuracy_v2] wrote %d long rows", nrow(long)))

# --- Summaries -----------------------------------------------------------
# Genome-wide per (track, panel, ancestry): weighted mean across chrs & haps,
# weights = n_sites in truth for that ancestry class.
weighted_summary <- long |>
  group_by(track, panel, ancestry) |>
  summarise(
    weighted_concordance = sum(n_correct) / sum(n_sites),
    unweighted_mean      = mean(concordance),
    se                   = sd(concordance) / sqrt(n()),
    n_haps_chr           = n(),
    total_sites          = sum(n_sites),
    .groups = "drop")
write_tsv(weighted_summary, file.path(OUTDIR, "accuracy_summary.tsv"))

per_chr <- long |>
  group_by(track, panel, ancestry, chr) |>
  summarise(concordance = sum(n_correct) / sum(n_sites), .groups = "drop")
write_tsv(per_chr, file.path(OUTDIR, "accuracy_per_chr.tsv"))

# --- Paired comparisons vs baseline (Wilcoxon signed-rank, Bonferroni) ---
message("[accuracy_v2] paired Wilcoxon vs baseline: ", BASELINE_PANEL)
wide <- long |>
  select(track, panel, chr, hap, sample, ancestry, concordance) |>
  pivot_wider(names_from = panel, values_from = concordance)

non_baseline <- setdiff(PANELS, BASELINE_PANEL)
comparisons <- expand.grid(track = TRACKS, panel = non_baseline,
                           ancestry = unname(ANC_LABELS),
                           stringsAsFactors = FALSE) |>
  filter(mapply(is_valid_combo, track, panel))

comp_rows <- vector("list", nrow(comparisons))
for (i in seq_len(nrow(comparisons))) {
  c_ <- comparisons[i, ]
  sub <- wide |> filter(track == c_$track, ancestry == c_$ancestry) |>
    select(chr, hap, all_of(BASELINE_PANEL), all_of(c_$panel)) |>
    drop_na()
  if (nrow(sub) < 10) {
    comp_rows[[i]] <- data.frame(track = c_$track, panel = c_$panel,
                                 ancestry = c_$ancestry, n = nrow(sub),
                                 base_mean = NA, alt_mean = NA, delta = NA,
                                 wilcox_p = NA, stringsAsFactors = FALSE)
    next
  }
  base <- sub[[BASELINE_PANEL]]; alt <- sub[[c_$panel]]
  w <- suppressWarnings(wilcox.test(alt, base, paired = TRUE))
  comp_rows[[i]] <- data.frame(
    track = c_$track, panel = c_$panel, ancestry = c_$ancestry,
    n = nrow(sub), base_mean = mean(base), alt_mean = mean(alt),
    delta = mean(alt - base), wilcox_p = w$p.value, stringsAsFactors = FALSE)
}
comp <- bind_rows(comp_rows) |>
  mutate(bonferroni_p = pmin(1, wilcox_p * sum(!is.na(wilcox_p))),
         signif = case_when(is.na(bonferroni_p) ~ "",
                            bonferroni_p < 0.001 ~ "***",
                            bonferroni_p < 0.01  ~ "**",
                            bonferroni_p < 0.05  ~ "*",
                            TRUE ~ "n.s."))
write_tsv(comp, file.path(OUTDIR, "panel_vs_baseline.tsv"))

# --- Plots ---------------------------------------------------------------
theme_lai <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# 1. Boxplot per panel per ancestry, split by track.
p1 <- ggplot(long, aes(x = panel, y = concordance, fill = panel)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  facet_grid(track ~ ancestry) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  labs(x = NULL, y = "Per-hap concordance",
       title = "LAI concordance by panel, split by sim track and ancestry") +
  theme_lai + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                    legend.position = "none")
ggsave(file.path(OUTDIR, "plots", "01_boxplot_panel_ancestry.pdf"),
       p1, width = 10, height = 6)

# 2. Mean +/- SE bar plot per panel per ancestry.
p2 <- ggplot(weighted_summary,
             aes(x = panel, y = weighted_concordance, fill = panel)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = weighted_concordance - se,
                    ymax = weighted_concordance + se), width = 0.3) +
  facet_grid(track ~ ancestry) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  labs(x = NULL, y = "Weighted concordance",
       title = "Weighted mean LAI concordance (site-weighted)") +
  theme_lai + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                    legend.position = "none")
ggsave(file.path(OUTDIR, "plots", "02_barplot_weighted_mean.pdf"),
       p2, width = 10, height = 6)

# 3. Per-chr line plot per panel per ancestry.
p3 <- ggplot(per_chr, aes(x = chr, y = concordance,
                          colour = panel, group = panel)) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  facet_grid(track ~ ancestry) +
  scale_x_continuous(breaks = seq(2, 22, 4)) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  labs(x = "Chromosome", y = "Concordance",
       title = "Per-chromosome concordance trajectory by panel") +
  theme_lai
ggsave(file.path(OUTDIR, "plots", "03_per_chr_lines.pdf"),
       p3, width = 12, height = 6)

# 4. Heatmap: (panel x chr) with concordance.
p4 <- ggplot(per_chr, aes(x = factor(chr), y = panel, fill = concordance)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0.5, 1), option = "mako", direction = -1,
                       labels = scales::percent_format(1),
                       oob = scales::squish) +
  facet_grid(track ~ ancestry) +
  labs(x = "Chromosome", y = NULL, fill = "Concordance",
       title = "Per-(panel, chr) concordance heatmap") +
  theme_lai
ggsave(file.path(OUTDIR, "plots", "04_heatmap_panel_chr.pdf"),
       p4, width = 12, height = 6)

# 5. Panel-vs-baseline paired scatter (points per (chr, hap)) for HGDPMXB panels.
mxb_panels <- intersect(PANELS, c("NAT_HGDPMXB", "NAT_HGDPMXB_FULL"))
if (length(mxb_panels) && BASELINE_PANEL %in% PANELS) {
  paired_long <- long |>
    select(track, panel, chr, hap, ancestry, concordance) |>
    pivot_wider(names_from = panel, values_from = concordance) |>
    pivot_longer(all_of(mxb_panels), names_to = "mxb_panel",
                 values_to = "mxb_concordance") |>
    rename(baseline_concordance = all_of(BASELINE_PANEL)) |>
    drop_na(baseline_concordance, mxb_concordance) |>
    filter(mapply(is_valid_combo, track, mxb_panel))

  p5 <- ggplot(paired_long, aes(x = baseline_concordance,
                                y = mxb_concordance)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2) +
    geom_point(alpha = 0.25, size = 0.7) +
    facet_grid(track + ancestry ~ mxb_panel) +
    coord_equal(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
    labs(x = paste("Baseline (", BASELINE_PANEL, ") concordance", sep = ""),
         y = "HGDPMXB panel concordance",
         title = "HGDPMXB panels vs. baseline, paired by (chr, hap)") +
    theme_lai
  ggsave(file.path(OUTDIR, "plots", "05_mxb_vs_baseline_scatter.pdf"),
         p5, width = 9, height = 10)
}

# 6. Panel deltas vs baseline as forest plot with Wilcoxon significance.
if (nrow(comp) > 0) {
  p6 <- comp |>
    mutate(label = sprintf("%s (%s)  n=%d  %s", panel, ancestry, n, signif)) |>
    ggplot(aes(x = delta, y = label, colour = track)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    geom_point(size = 2) +
    facet_wrap(~ track, scales = "free_y", ncol = 1) +
    scale_x_continuous(labels = scales::percent_format(0.1)) +
    labs(x = paste0("Concordance delta vs ", BASELINE_PANEL,
                    " (positive = alt panel better)"),
         y = NULL, colour = "Sim track",
         title = "Panel vs baseline deltas (paired chr x hap, Bonferroni-adjusted)") +
    theme_lai + theme(legend.position = "none")
  ggsave(file.path(OUTDIR, "plots", "06_delta_forest.pdf"),
         p6, width = 10, height = 7)
}

message("[accuracy_v2] DONE. Outputs under ", OUTDIR)
message("  Tables : accuracy_long.tsv, accuracy_summary.tsv, accuracy_per_chr.tsv, panel_vs_baseline.tsv")
message("  Plots  : plots/01..06_*.pdf")
