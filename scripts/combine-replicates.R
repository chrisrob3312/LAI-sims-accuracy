#!/usr/bin/env Rscript
# ============================================================================
# combine-replicates.R
#
# Combine accuracy_summary.tsv across independent admix-simu REPLICATES and
# report across-replicate uncertainty.
#
# WHY THIS IS THE RIGOROUS ESTIMATOR, and what it fixes:
#   Within a single simulation run, every uncertainty estimator is limited by
#   the same two dependencies:
#     (1) sites inside one ancestry tract are near-duplicates of each other;
#     (2) every simulated haplotype is a mosaic of the SAME donor pool, so an
#         atypical donor haplotype correlates errors across every block and
#         every tract that drew from it.
#   The block bootstrap in accuracy_v2.R handles (1) and explicitly cannot
#   handle (2) -- it resamples blocks within one realisation, so the donor
#   pool is held fixed in every draw.
#   Replicate simulations re-draw BOTH: new recombination breakpoints and a
#   new sampling of donor haplotypes into mosaics. The spread of the pooled
#   metric ACROSS replicates therefore captures the full variance that a
#   within-run resample structurally cannot see.
#
#   Consequence to expect, and to report honestly: the across-replicate
#   interval should be WIDER than the block-bootstrap interval. If it is not,
#   the replicates are not independent -- check that each used a different
#   seed and that they did not silently reuse a cached simulation.
#
# Usage:
#   Rscript scripts/combine-replicates.R <rep1_dir> <rep2_dir> ... [--out DIR]
#     where each <repN_dir> contains an accuracy_summary.tsv
#
# Output (into --out, default the first replicate's parent):
#   accuracy_replicates_long.tsv     one row per replicate x track x panel x ancestry
#   accuracy_replicates_summary.tsv  mean, SD, SE, t-CI across replicates,
#                                    side by side with the block-bootstrap CI
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
out_i <- which(args == "--out")
OUTDIR <- if (length(out_i) == 1 && length(args) > out_i) args[out_i + 1] else NA_character_
dirs <- args[!seq_along(args) %in% c(out_i, out_i + 1)]
if (length(dirs) < 2) {
  stop("need at least 2 replicate directories; got ", length(dirs),
       "\nusage: Rscript scripts/combine-replicates.R <rep1> <rep2> ... [--out DIR]")
}
if (is.na(OUTDIR)) OUTDIR <- dirname(normalizePath(dirs[1], mustWork = FALSE))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- read every replicate -------------------------------------------------
read_rep <- function(d, i) {
  f <- file.path(d, "accuracy_summary.tsv")
  if (!file.exists(f)) {
    warning("missing ", f, " -- replicate skipped"); return(NULL)
  }
  suppressMessages(read_tsv(f, show_col_types = FALSE)) |>
    mutate(replicate = i, rep_dir = d, .before = 1)
}
reps <- bind_rows(lapply(seq_along(dirs), function(i) read_rep(dirs[i], i)))
if (nrow(reps) == 0) stop("no replicate summaries could be read")

n_rep <- length(unique(reps$replicate))
message(sprintf("[replicates] read %d replicate(s), %d rows", n_rep, nrow(reps)))
if (n_rep < length(dirs))
  message(sprintf("[replicates] WARNING: %d of %d directories had no summary",
                  length(dirs) - n_rep, length(dirs)))

write_tsv(reps, file.path(OUTDIR, "accuracy_replicates_long.tsv"))

# --- across-replicate statistics -----------------------------------------
# The replicate is the independent unit here, so n is the replicate count --
# small, which is why a t-interval is used rather than a normal one. With 5
# replicates t(.975, 4) = 2.776, appreciably wider than 1.96; using the normal
# quantile at this n would understate the interval by ~40%.
metric_cols <- c("weighted_recall", "weighted_precision", "weighted_f1")
present <- intersect(metric_cols, names(reps))
if (length(present) == 0) stop("no metric columns found in the replicate summaries")

summ <- reps |>
  group_by(track, panel, ancestry) |>
  summarise(
    n_reps = n(),
    across(all_of(present),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    # Carry one replicate's block-bootstrap CI through for comparison. These
    # are within-run intervals and are expected to be narrower.
    block_recall_lo = first(na.omit(c(recall_lo, NA_real_))),
    block_recall_hi = first(na.omit(c(recall_hi, NA_real_))),
    n_blocks        = first(na.omit(c(n_blocks, NA_integer_))),
    .groups = "drop")

tcrit <- function(n) if (n >= 2) stats::qt(0.975, df = n - 1) else NA_real_

for (m in present) {
  mu <- summ[[paste0(m, "_mean")]]; s <- summ[[paste0(m, "_sd")]]
  se <- s / sqrt(summ$n_reps)
  tc <- vapply(summ$n_reps, tcrit, numeric(1))
  summ[[paste0(m, "_se")]] <- se
  summ[[paste0(m, "_lo")]] <- mu - tc * se
  summ[[paste0(m, "_hi")]] <- mu + tc * se
}

# Width comparison: across-replicate vs within-run block bootstrap.
summ <- summ |>
  mutate(rep_ci_width   = weighted_recall_hi - weighted_recall_lo,
         block_ci_width = block_recall_hi - block_recall_lo,
         width_ratio    = rep_ci_width / block_ci_width)

write_tsv(summ, file.path(OUTDIR, "accuracy_replicates_summary.tsv"))

# --- report ---------------------------------------------------------------
message(sprintf("[replicates] wrote %s", file.path(OUTDIR, "accuracy_replicates_summary.tsv")))
if (n_rep < 3)
  message("[replicates] NOTE: with <3 replicates the t-interval is very wide. ",
          "5 is the practical minimum for a reportable error bar.")

ratios <- summ$width_ratio[is.finite(summ$width_ratio)]
if (length(ratios) > 0) {
  message(sprintf("[replicates] CI width, across-replicate / block-bootstrap: median %.2fx (range %.2f-%.2f)",
                  median(ratios), min(ratios), max(ratios)))
  if (median(ratios) < 1) {
    message("[replicates] *** WARNING: across-replicate intervals are NARROWER than the")
    message("[replicates]     within-run block bootstrap. That should not happen if the")
    message("[replicates]     replicates are independent. Check that each replicate used a")
    message("[replicates]     different seed and did not reuse a cached simulation.")
  }
}

# Headline table to stdout, so a run is readable without opening the TSV.
hdr <- summ |>
  filter(ancestry %in% c("NAT", "AMR")) |>
  select(track, panel, ancestry, n_reps,
         mean = weighted_recall_mean, sd = weighted_recall_sd,
         lo = weighted_recall_lo, hi = weighted_recall_hi)
if (nrow(hdr) > 0) {
  message("\n[replicates] Amerindigenous recall across replicates:")
  print(as.data.frame(hdr), row.names = FALSE, digits = 4)
}
