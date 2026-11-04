#!/usr/bin/env Rscript
# qc_truth_set.R  <pca.eigenvec>  <samples.tsv (sample_id, cohort)>  <outdir>
# Plots PC1-PC4 of the simulated truth set coloured by cohort and writes
# per-cohort PC centroids. Homogeneous cohorts should form tight clusters;
# admixed cohorts should lie between their source clusters in proportion
# to the designed mix (e.g. ADM_MEX ~ half-way EUR->AMR, slightly toward AFR).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: qc_truth_set.R pca.eigenvec samples.tsv outdir")
ev  <- read.table(args[1], header = TRUE, comment.char = "", check.names = FALSE)
names(ev)[1:2] <- c("FID", "IID")
smp <- read.table(args[2], header = FALSE, sep = "\t", col.names = c("IID", "cohort"))
d <- merge(ev, smp, by = "IID")
d$cohort <- factor(d$cohort)
cent <- aggregate(d[, grep("^PC", names(d))[1:4]], list(cohort = d$cohort), mean)
write.table(cent, file.path(args[3], "pca_centroids.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pal <- grDevices::hcl.colors(nlevels(d$cohort), "Dark 3")
pdf(file.path(args[3], "pca_by_cohort.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (pcs in list(c("PC1", "PC2"), c("PC3", "PC4"))) {
  plot(d[[pcs[1]]], d[[pcs[2]]], col = pal[d$cohort], pch = 16, cex = 0.6,
       xlab = pcs[1], ylab = pcs[2], main = "Simulated benchmark truth set")
  text(cent[[pcs[1]]], cent[[pcs[2]]], labels = cent$cohort, cex = 0.6, font = 2)
}
dev.off()
cat("wrote", file.path(args[3], "pca_by_cohort.pdf"), "and pca_centroids.tsv\n")
