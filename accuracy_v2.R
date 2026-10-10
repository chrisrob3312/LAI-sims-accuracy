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

# Comma-separated chr override, e.g. CHRS=20,21,22 to restrict a pilot run.
CHRS   <- if (nzchar(Sys.getenv("CHRS"))) {
    as.integer(strsplit(Sys.getenv("CHRS"), "[, ]+")[[1]])
} else {
    1:22
}
TRACKS <- c("NAT", "NATMXB")
PANELS <- strsplit(Sys.getenv("PANELS_TO_RUN",
                              "NAT_HGDP NAT_PEL NAT_PEL_EAS NAT_HGDPMXB NAT_HGDPMXB_FULL NAT_HOMOG"),
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
  n_haps  <- length(L)
  n_sites <- nchar(L[1L])
  # Character-by-character parse in one pass: paste all lines, tear into
  # single-char vector, coerce to integer, reshape as (n_sites x n_haps).
  # Bulk vector allocation is ~4x less RAM than rbind-of-split.
  chars <- unlist(strsplit(L, "", fixed = TRUE), use.names = FALSE)
  rm(L); gc(FALSE)
  # Column-major fill: chars[1..n_sites] = hap 1 sites 1..N so
  # M[s, h] = char at site s of hap h with byrow=FALSE (R default).
  m <- matrix(as.integer(chars), nrow = n_sites, ncol = n_haps)
  rm(chars); gc(FALSE)
  m
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
  # Use scan() into a raw integer vector — much less RAM than read.table,
  # which materialises a data.frame + several copies during coercion.
  # We know the column count from the first line: 2 metadata + N haps.
  ncol_total <- length(scan(path, what = integer(), nlines = 1, quiet = TRUE))
  flat <- scan(path, what = integer(), quiet = TRUE)
  m <- matrix(flat, ncol = ncol_total, byrow = TRUE)
  rm(flat); gc(FALSE)
  list(pos = m[, 1L], mat = m[, -c(1L, 2L), drop = FALSE])
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
  # Skip anything without a .done marker so partial runs (e.g. 55/198)
  # produce clean output for whatever finished. Emit a message() (not
  # warning) so it's visible in the run log without polluting the summary
  # -- these are expected during partial pilots but hide silent bugs (e.g.
  # missing .hanc truth files for a chr range) if not surfaced.
  missing <- c(
    hanc = !file.exists(hanc_f),
    lat3 = !file.exists(lat3_f),
    done = !file.exists(done_f))
  if (any(missing)) {
    message(sprintf("      [SKIP %s/%s chr%d] missing: %s",
                    track, panel, chr,
                    paste(names(missing)[missing], collapse = ",")))
    return(NULL)
  }

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
  # For each class a we record:
  #   n_sites      = sites where truth == a  (support for recall)
  #   n_correct    = sites where truth == a AND rec == a  (TP)
  #   concordance  = n_correct / n_sites  = per-class recall (Jessica's TPR)
  #   n_pred       = sites where rec == a   (support for precision)
  #   n_tp         = same as n_correct (for clarity in F1 downstream)
  #   n_fp         = sites where rec == a AND truth != a
  #   n_fn         = sites where truth == a AND rec != a
  # Downstream aggregation computes precision, recall, F1 per (track,panel,ancestry).
  out <- vector("list", n_haps * length(ANC_LABELS))
  k <- 0L
  for (h in seq_len(n_haps)) {
    tv <- T[, h]; rv <- R[, h]
    for (a in as.integer(names(ANC_LABELS))) {
      truth_mask <- (tv == a)
      pred_mask  <- (rv == a)
      n_a   <- sum(truth_mask)
      n_p   <- sum(pred_mask)
      # BUGFIX (precision inflation). This used to be `if (n_a == 0L) next`,
      # which dropped the whole (hap, chr, ancestry) stratum whenever truth
      # carried no sites of that class. Any site PREDICTED as that class on
      # such a stratum is a pure false positive, and skipping the row meant
      # those FPs were never counted in n_fp anywhere -- so genome-wide
      # precision = sum(tp)/(sum(tp)+sum(fp)) was computed against an
      # artificially small denominator. The effect scales with how often an
      # ancestry is ABSENT from a stratum, so it inflates the RAREST ancestry
      # most -- i.e. NAT/AMR, the headline metric of this study.
      # Recall is unaffected: tp and fn both require truth == a, so a
      # zero-truth stratum contributes 0 to each and skipping it was a no-op.
      # Keep the row whenever the stratum is informative in EITHER direction.
      if (n_a == 0L && n_p == 0L) next
      n_tp  <- sum(truth_mask & pred_mask)
      n_fp  <- n_p - n_tp
      n_fn  <- n_a - n_tp
      k <- k + 1L
      out[[k]] <- data.frame(
        track       = track,
        panel       = panel,
        chr         = chr,
        hap         = h,
        sample      = hap_sample[h],
        ancestry    = ANC_LABELS[[as.character(a)]],
        n_sites     = as.integer(n_a),
        n_correct   = as.integer(n_tp),
        # Per-class recall / TPR. Undefined (not zero) when truth carries no
        # sites of this class -- NA keeps it out of macro-averages instead of
        # dragging them toward 0.
        concordance = if (n_a > 0L) n_tp / n_a else NA_real_,
        n_pred      = as.integer(n_p),
        n_fp        = as.integer(n_fp),
        n_fn        = as.integer(n_fn),
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
# Account for EVERY combo. Previously a combo that returned NULL -- whether
# from a caught error, a missing file, or the <100-shared-sites guard (which
# raises warning(), not error(), so it slipped past the error handler) -- just
# vanished from the output with nothing downstream to notice. A silently
# missing panel x chr reads as "scored and fine" in the summary tables.
skipped <- vector("list", nrow(grid)); n_skip <- 0L
last_warn <- NULL; err_msg <- NULL
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  message(sprintf("  [%d/%d] %s / %s chr%d", i, nrow(grid), g$track, g$panel, g$chr))
  last_warn <- NULL; err_msg <- NULL
  res <- withCallingHandlers(
    tryCatch(score_one(g$track, g$panel, g$chr),
             error = function(e) { err_msg <<- conditionMessage(e); NULL }),
    # The <100-shared-sites guard raises warning(), not error(), so it would
    # otherwise slip past tryCatch entirely. Capture it, then muffle so the
    # loop behaves exactly as before.
    warning = function(w) { last_warn <<- conditionMessage(w)
                            invokeRestart("muffleWarning") })
  if (is.null(res) || nrow(res) == 0L) {
    n_skip <- n_skip + 1L
    reason <- if (!is.null(err_msg))   paste("error:", err_msg)
         else if (!is.null(last_warn)) last_warn
         else "no rows (missing input, missing .done, or <100 shared sites)"
    skipped[[n_skip]] <- data.frame(track = g$track, panel = g$panel, chr = g$chr,
                                    reason = reason, stringsAsFactors = FALSE)
    message(sprintf("    [SKIPPED] %s / %s chr%d -- %s", g$track, g$panel, g$chr, reason))
  } else {
    results[[i]] <- res
  }
  # Explicit gc() every iteration -- big chrs (1-3) allocate ~2 GB matrices,
  # R's lazy GC leaves them resident and kills the process by chr 4-5.
  gc(FALSE)
}
long <- bind_rows(results)
if (nrow(long) == 0) stop("No accuracy rows produced -- check paths and file existence.")

# Persist the skip ledger so a partial sweep is visible in the outputs, not
# only in scrollback.
skip_tbl <- if (n_skip > 0L) bind_rows(skipped[seq_len(n_skip)]) else
  data.frame(track = character(0), panel = character(0), chr = integer(0),
             reason = character(0), stringsAsFactors = FALSE)
write_tsv(skip_tbl, file.path(OUTDIR, "accuracy_skipped_combos.tsv"))
message(sprintf("[accuracy_v2] %d/%d combos scored, %d skipped (see accuracy_skipped_combos.tsv)",
                nrow(grid) - n_skip, nrow(grid), n_skip))

write_tsv(long, file.path(OUTDIR, "accuracy_long.tsv"))
message(sprintf("[accuracy_v2] wrote %d long rows", nrow(long)))

# --- Summaries -----------------------------------------------------------
# Genome-wide per (track, panel, ancestry): weighted mean across chrs & haps,
# weights = n_sites in truth for that ancestry class.
# Metrics reported:
#   recall (aka TPR, aka concordance): TP / (TP + FN) = n_correct / n_sites
#   precision:                          TP / (TP + FP) = n_correct / n_pred
#   F1:                                 2 * P * R / (P + R)
#   These are the standard metrics Honorato-Mauer 2025 reports; F1 is preferred
#   for imbalanced ancestries (NAT is 15% of Brasa sites, so recall alone
#   overstates when the model predicts NAT rarely).
weighted_summary <- long |>
  group_by(track, panel, ancestry) |>
  summarise(
    tp = sum(n_correct),
    fp = sum(n_fp),
    fn = sum(n_fn),
    total_sites          = sum(n_sites),
    n_haps_chr           = n(),
    .groups = "drop") |>
  mutate(
    weighted_recall      = tp / (tp + fn),
    weighted_precision   = tp / (tp + fp),
    weighted_f1          = 2 * weighted_precision * weighted_recall /
                            (weighted_precision + weighted_recall),
    weighted_concordance = weighted_recall  # legacy alias for compat with old plots
  ) |>
  # keep same column order as before + append precision/F1
  select(track, panel, ancestry,
         weighted_concordance, weighted_recall, weighted_precision, weighted_f1,
         total_sites, n_haps_chr, tp, fp, fn)

# --- Uncertainty: block bootstrap over (chr x hap) -------------------------
# What is being estimated here is a SITE-LEVEL RATE: how often the caller gets
# a genome position right. It is a methods benchmark, not a biological
# experiment, so the question is "which genome sites, and how often correct".
# That rules out two tempting units:
#
#   NOT the site. Ancestry segments span megabases -- at GEN=12 the expected
#   segment is ~8 cM, so one true ancestry state covers thousands of adjacent
#   SNPs and the caller usually labels the whole run identically. Treating
#   ~10^7 sites as independent overstates the effective sample size by three to
#   four orders of magnitude and yields a meaninglessly tiny SE.
#
#   NOT the individual. admix-simu builds every simulated haplotype as its own
#   mosaic drawn from a SHARED donor pool. The two haplotypes of one simulated
#   "individual" are no more related to each other than to any other simulated
#   haplotype, so clustering by individual imposes structure the generative
#   process does not have.
#
# The real unit is the independent ancestry SEGMENT (~3600 cM / ~8 cM ~= 450
# per haplotype). The available proxy is the (chr x hap) BLOCK: an independent
# recombination realization containing a set of segments. Resampling those
# blocks and recomputing the POOLED metric respects within-segment correlation
# without pretending sites are independent.
#
# This also fixes the original estimator mismatch, which was a separate bug:
# `se = sd(concordance)/sqrt(n())` reported the SD of the unweighted per-block
# RECALL next to a point estimate that is a pooled micro-average over summed
# TP/FP/FN -- the SE of a different quantity than the number it annotated.
#
# RESIDUAL LIMITATION, stated honestly: blocks still share the donor pool, so
# an atypical donor haplotype correlates errors across blocks. This bootstrap
# does not capture that. The rigorous version is REPLICATE SIMULATIONS with
# different admix-simu seeds, taking the SE across replicates; set
# BOOT_UNIT=replicate once replicate runs exist.
BOOT_B <- as.integer(Sys.getenv("BOOT_B", "1000"))
BOOT_SEED <- as.integer(Sys.getenv("BOOT_SEED", "20261010"))

boot_pooled <- function(df, B, seed) {
  # Block = one (chr, hap) recombination realization.
  blk <- paste(df$chr, df$hap, sep = ":")
  ids <- unique(blk)
  n <- length(ids)
  if (n < 2L || B < 1L) {
    return(data.frame(n_blocks = n,
                      recall_lo = NA_real_, recall_hi = NA_real_,
                      precision_lo = NA_real_, precision_hi = NA_real_,
                      f1_lo = NA_real_, f1_hi = NA_real_,
                      recall_se = NA_real_, precision_se = NA_real_, f1_se = NA_real_))
  }
  # Precompute row indices per block so each bootstrap draw is a cheap index
  # gather rather than an O(rows) scan.
  idx_by_id <- split(seq_len(nrow(df)), blk)
  tpv <- df$n_correct; fpv <- df$n_fp; fnv <- df$n_fn
  set.seed(seed)
  # NB: do not name these R / F -- F is a base alias for FALSE and R shadows
  # the inferred-call matrix name used elsewhere in this script.
  bootR <- numeric(B); bootP <- numeric(B); bootF <- numeric(B)
  for (b in seq_len(B)) {
    pick <- sample(ids, n, replace = TRUE)
    idx  <- unlist(idx_by_id[pick], use.names = FALSE)
    tp <- sum(tpv[idx]); fp <- sum(fpv[idx]); fn <- sum(fnv[idx])
    r <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    p <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    bootR[b] <- r; bootP[b] <- p
    bootF[b] <- if (!is.na(p) && !is.na(r) && (p + r) > 0) 2 * p * r / (p + r) else NA_real_
  }
  q <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  qr <- q(bootR); qp <- q(bootP); qf <- q(bootF)
  data.frame(n_blocks = n,
             recall_lo = qr[1], recall_hi = qr[2],
             precision_lo = qp[1], precision_hi = qp[2],
             f1_lo = qf[1], f1_hi = qf[2],
             recall_se = sd(bootR, na.rm = TRUE),
             precision_se = sd(bootP, na.rm = TRUE),
             f1_se = sd(bootF, na.rm = TRUE))
}

message(sprintf("[accuracy_v2] cluster bootstrap over individuals (B=%d)", BOOT_B))
boot_stats <- long |>
  group_by(track, panel, ancestry) |>
  group_modify(~ boot_pooled(.x, BOOT_B, BOOT_SEED)) |>
  ungroup()

# Macro (unweighted) mean of per-(hap, chr) recall. Kept for continuity with
# the old plots but RENAMED so it is never mistaken for the pooled figure: it
# is an average of averages, weighting chr22 equally with chr1 and a hap with
# 50 NAT sites equally with one carrying 50,000. na.rm because recall is NA on
# strata where truth carried none of that ancestry (see the n_a == 0 note).
macro_stats <- long |>
  group_by(track, panel, ancestry) |>
  summarise(
    macro_mean_recall_hapchr = mean(concordance, na.rm = TRUE),
    macro_sd_recall_hapchr   = sd(concordance, na.rm = TRUE),
    n_hapchr_strata          = sum(!is.na(concordance)),
    unweighted_mean          = mean(concordance, na.rm = TRUE),  # legacy alias
    .groups = "drop")

weighted_summary <- weighted_summary |>
  left_join(boot_stats,  by = c("track", "panel", "ancestry")) |>
  left_join(macro_stats, by = c("track", "panel", "ancestry")) |>
  # legacy `se` column: now the bootstrap SE of the POOLED recall, which is
  # the quantity the point estimate actually reports.
  mutate(se = recall_se)
write_tsv(weighted_summary, file.path(OUTDIR, "accuracy_summary.tsv"))

# Per-chr: recall + F1
per_chr <- long |>
  group_by(track, panel, ancestry, chr) |>
  summarise(
    tp = sum(n_correct), fp = sum(n_fp), fn = sum(n_fn),
    n_sites = sum(n_sites),
    .groups = "drop") |>
  mutate(
    recall      = tp / (tp + fn),
    precision   = tp / (tp + fp),
    f1          = 2 * precision * recall / (precision + recall),
    concordance = recall  # legacy alias
  )
write_tsv(per_chr, file.path(OUTDIR, "accuracy_per_chr.tsv"))

# --- Paired comparisons vs baseline (Wilcoxon signed-rank, Bonferroni) ---
# PAIRING UNIT. Pair on the (chr, hap) BLOCK, matching the bootstrap unit
# above: an independent recombination realization. Individual is NOT the unit
# -- simulated haplotypes are separate mosaics from a shared donor pool, so the
# two haps of one simulated individual are no more related than any other pair.
# Pairing removes block-level variation, which is what makes the comparison
# informative: the same block is scored under both panels.
#
# CAVEAT, not fixable by the test: blocks share the donor pool, so they are not
# strictly independent and these p-values are mildly anti-conservative. Treat a
# marginal result as suggestive. The block-bootstrap CI on the paired delta
# (below) is the more honest summary, and replicate simulations with different
# seeds would settle it properly.
message("[accuracy_v2] paired Wilcoxon vs baseline (paired by chr x hap block): ",
        BASELINE_PANEL)

# Per-individual pooled recall is still written out as a descriptive view --
# useful for spotting an outlier simulated individual -- but it is NOT the
# inference unit.
per_individual <- long |>
  group_by(track, panel, ancestry, sample) |>
  summarise(tp = sum(n_correct), fn = sum(n_fn), fp = sum(n_fp), .groups = "drop") |>
  mutate(recall = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_))
write_tsv(per_individual, file.path(OUTDIR, "accuracy_per_individual.tsv"))

wide <- long |>
  select(track, panel, chr, hap, sample, ancestry, concordance) |>
  pivot_wider(names_from = panel, values_from = concordance)

# Only compare panels that actually landed in the wide table -- during a
# partial run (e.g. NAT_PEL_EAS still finishing), pivot_wider will omit any
# panel that has zero completed combos, so all_of(panel) would blow up.
available_panels <- intersect(PANELS, names(wide))
non_baseline <- setdiff(available_panels, BASELINE_PANEL)
if (!(BASELINE_PANEL %in% available_panels)) {
  message("[accuracy_v2] baseline panel ", BASELINE_PANEL,
          " has no done combos yet -- skipping panel_vs_baseline.tsv")
  non_baseline <- character(0)
}
comparisons <- expand.grid(track = TRACKS, panel = non_baseline,
                           ancestry = unname(ANC_LABELS),
                           stringsAsFactors = FALSE) |>
  filter(mapply(is_valid_combo, track, panel))

comp_rows <- vector("list", max(nrow(comparisons), 1L))
if (nrow(comparisons) == 0L) comp_rows <- list()
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
comp <- if (length(comp_rows)) {
  bind_rows(comp_rows) |>
    mutate(bonferroni_p = pmin(1, wilcox_p * sum(!is.na(wilcox_p))),
           signif = case_when(is.na(bonferroni_p) ~ "",
                              bonferroni_p < 0.001 ~ "***",
                              bonferroni_p < 0.01  ~ "**",
                              bonferroni_p < 0.05  ~ "*",
                              TRUE ~ "n.s."))
} else {
  data.frame(track = character(0), panel = character(0),
             ancestry = character(0), n = integer(0),
             base_mean = numeric(0), alt_mean = numeric(0),
             delta = numeric(0), wilcox_p = numeric(0),
             bonferroni_p = numeric(0), signif = character(0))
}
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

# 7. MXB-in-reference vs MXB-not-in-reference, split by sim track.
#    Answers: "does adding MXB to the RFMix reference help, and does it depend
#    on whether the simulated admixed haps contain MXB donors?"
mxb_in_ref <- c("NAT_HGDPMXB", "NAT_HGDPMXB_FULL", "NAT_HOMOG")
long_mxb_flag <- long |>
  mutate(mxb_in_ref = ifelse(panel %in% mxb_in_ref, "MXB in reference", "no MXB in reference"),
         track_lab  = ifelse(track == "NATMXB", "MXB in simulation", "no MXB in simulation"))

p7 <- ggplot(long_mxb_flag,
             aes(x = panel, y = concordance, fill = mxb_in_ref)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  facet_grid(track_lab ~ ancestry) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  scale_fill_manual(values = c("MXB in reference" = "#4C9F70",
                               "no MXB in reference" = "#B95E82")) +
  labs(x = NULL, y = "Per-hap concordance", fill = NULL,
       title = "Effect of MXB in reference panel, split by whether MXB is in simulation",
       subtitle = "Rows: sim track (with/without MXB donors)  |  Columns: ancestry class") +
  theme_lai + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "plots", "07_mxb_in_ref_by_sim_track.pdf"),
       p7, width = 12, height = 7)

# 8. Mexican-cohort framing: NAT-ancestry accuracy only, ranked panel bar with
#    both sim tracks side-by-side. This is the "for my Mexican-admixed GWAS
#    cohort, which reference panel gives me the best NAT ancestry calls" plot.
nat_only <- weighted_summary |>
  filter(ancestry == "NAT") |>
  mutate(track_lab = ifelse(track == "NATMXB",
                            "MXB donors in simulated admixed",
                            "HGDP-only donors in simulated admixed"),
         mxb_in_ref = ifelse(panel %in% mxb_in_ref, "MXB in ref", "no MXB in ref"))

p8 <- ggplot(nat_only,
             aes(x = reorder(panel, weighted_concordance),
                 y = weighted_concordance, fill = mxb_in_ref)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_errorbar(aes(ymin = weighted_concordance - se,
                    ymax = weighted_concordance + se), width = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * weighted_concordance)),
            hjust = -0.15, size = 3.2) +
  coord_flip(ylim = c(0, 1.05)) +
  facet_wrap(~ track_lab, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(1)) +
  scale_fill_manual(values = c("MXB in ref" = "#4C9F70",
                               "no MXB in ref" = "#B95E82")) +
  labs(x = NULL, y = "Native-American ancestry concordance",
       fill = NULL,
       title = "NAT-ancestry concordance by reference panel",
       subtitle = "Split by whether MXB donors are in the simulated admixed sample; error bars = SE across haps") +
  theme_lai + theme(legend.position = "top")
ggsave(file.path(OUTDIR, "plots", "08_nat_ancestry_by_panel_ranked.pdf"),
       p8, width = 10, height = 8)

message("[accuracy_v2] DONE. Outputs under ", OUTDIR)
message("  Tables : accuracy_long.tsv, accuracy_summary.tsv, accuracy_per_chr.tsv, panel_vs_baseline.tsv")
message("  Plots  : plots/01..08_*.pdf")
