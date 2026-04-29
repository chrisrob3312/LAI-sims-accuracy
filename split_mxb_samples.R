#!/usr/bin/env Rscript
# Stratified 25/25 split of the 50 MX Biobank WGS samples by inferred genetic
# cluster. Half become additional NAT haplotype donors for ADMIX-SIMU; half
# join the RFMix reference (panel #4: HGDP-NAT + MXB).
#
# Allocation per cluster (sums to 25 SIMU + 25 RFMIX):
#   Nahua             18 -> 9 / 9
#   Zapotec/Mazatec   10 -> 5 / 5
#   Maya               7 -> 3 / 4
#   Totonac            6 -> 3 / 3
#   Tzotzil            4 -> 2 / 2
#   Tarahumara         3 -> 2 / 1
#   Huichol            2 -> 1 / 1
#
# The split within each cluster is deterministic (lexicographic Sample_ID
# sort, first n_simu -> SIMU, remainder -> RFMIX) so it is reproducible from
# any language without seed-RNG concerns.

popinfo_path <- "reference_ids/MXB50genomes_popinfo.tsv"
hgdp_rfmix   <- "reference_ids/amr_rfmix.txt"
out_simu     <- "reference_ids/mxb_simu.txt"
out_rfmix    <- "reference_ids/mxb_rfmix.txt"
out_combined <- "reference_ids/amr_hgdpmxb_rfmix.txt"

popinfo <- read.table(popinfo_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)

allocation <- list(
  "Nahua"           = c(simu = 9, rfmix = 9),
  "Zapotec/Mazatec" = c(simu = 5, rfmix = 5),
  "Maya"            = c(simu = 3, rfmix = 4),
  "Totonac"         = c(simu = 3, rfmix = 3),
  "Tzotzil"         = c(simu = 2, rfmix = 2),
  "Tarahumara"      = c(simu = 2, rfmix = 1),
  "Huichol"         = c(simu = 1, rfmix = 1)
)

simu_ids  <- character()
rfmix_ids <- character()

for (cluster in names(allocation)) {
  ids <- sort(popinfo$Sample_ID[popinfo$inferred_genetic_cluster == cluster])
  n_simu  <- allocation[[cluster]]["simu"]
  n_rfmix <- allocation[[cluster]]["rfmix"]
  stopifnot(length(ids) == n_simu + n_rfmix)
  simu_ids  <- c(simu_ids,  ids[seq_len(n_simu)])
  rfmix_ids <- c(rfmix_ids, ids[(n_simu + 1):(n_simu + n_rfmix)])
}

stopifnot(length(simu_ids)  == 25)
stopifnot(length(rfmix_ids) == 25)
stopifnot(length(intersect(simu_ids, rfmix_ids)) == 0)
stopifnot(length(union(simu_ids, rfmix_ids))     == 50)

writeLines(sort(simu_ids),  out_simu)
writeLines(sort(rfmix_ids), out_rfmix)

# Panel #4 RFMix reference = HGDP-NAT (rfmix half) + MXB (rfmix half).
hgdp <- readLines(hgdp_rfmix)
writeLines(c(hgdp, sort(rfmix_ids)), out_combined)

cat(sprintf("SIMU      %2d -> %s\n",  length(simu_ids),  out_simu))
cat(sprintf("RFMIX     %2d -> %s\n",  length(rfmix_ids), out_rfmix))
cat(sprintf("COMBINED  %2d -> %s\n",
            length(c(hgdp, rfmix_ids)), out_combined))
