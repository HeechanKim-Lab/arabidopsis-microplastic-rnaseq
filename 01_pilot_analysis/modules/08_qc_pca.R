#!/usr/bin/env Rscript
# ==============================================================================
# Module 08: log-CPM values and PCA of the six libraries (log 2026-04-08_01)
#
# Input : results/intermediate/dge_filtered_tmm.rds (module 07), phenodata_synced.rds (module 06)
# Output: assets/001_pilot_qc_pca.png
#         results/intermediate/logCPM_6samples.rds
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

dge       <- readRDS(file.path(inter_dir, "dge_filtered_tmm.rds"))
phenodata <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))

# Log2-CPM from the TMM-normalised counts; the prior count of 2 avoids log(0)
logCPM <- cpm(dge, log = TRUE, prior.count = 2)

# PCA with samples as rows; centred, not scaled (prcomp defaults)
pca_res <- prcomp(t(logCPM))

pca_df <- data.frame(
    Sample    = rownames(pca_res$x),
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    Condition = phenodata$condition
)
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))
cat("[08] Variance explained (%):", paste(percentVar, collapse = ", "), "\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
    geom_point(size = 4) +
    geom_text(vjust = 1.5) +
    theme_minimal() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot: Arabidopsis Root Samples")
ggsave(file.path(asset_dir, "001_pilot_qc_pca.png"), p_pca, width = 7, height = 5, dpi = 200, bg = "white")

saveRDS(logCPM, file.path(inter_dir, "logCPM_6samples.rds"))
