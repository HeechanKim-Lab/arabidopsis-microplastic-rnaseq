#!/usr/bin/env Rscript
# ==============================================================================
# Module: 08_qc_pca.R
# Description: Log-CPM transformation, SVD, and ordination plotting (PC1 vs PC2)
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
})

base_dir  <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

dge <- readRDS(file.path(inter_dir, "dge_filtered_tmm.rds"))

logCPM  <- cpm(dge, log = TRUE, prior.count = 2)
pca_res <- prcomp(t(logCPM), scale. = TRUE)
pct_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- data.frame(
    Sample    = rownames(pca_res$x),
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    Condition = dge$samples$group
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
    geom_point(size = 4.5, alpha = 0.85) +
    geom_text(vjust = 1.8, size = 4) +
    scale_color_manual(values = c("WT" = "#377EB8", "Cond" = "#E41A1C")) +
    theme_minimal(base_size = 12) +
    labs(
        title = "PCA Plot: Arabidopsis Root Samples (Pilot N=6)",
        x     = paste0("PC1: ", pct_var[1], "% variance"),
        y     = paste0("PC2: ", pct_var[2], "% variance")
    ) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(face = "bold", size = 12)
    )

out_png <- file.path(asset_dir, "001_pilot_qc_pca.png")
ggsave(out_png, plot = p_pca, width = 8, height = 6, dpi = 300)
cat(">>> [08] Saved PCA diagnostic asset to:", out_png, "\n")

saveRDS(logCPM, file.path(inter_dir, "logCPM_6samples.rds"))