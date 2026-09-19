#!/usr/bin/env Rscript
# ==============================================================================
# Module: 10_visualization.R
# Description: Publication Volcano plot and Z-score clustering heatmap
# ==============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
    library(dplyr)
})

base_dir  <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")

res_lrt         <- readRDS(file.path(inter_dir, "res_lrt.rds"))
final_sig_genes <- readRDS(file.path(inter_dir, "final_sig_genes.rds"))
logCPM_6        <- readRDS(file.path(inter_dir, "logCPM_6samples.rds"))
phenodata       <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))

# 1. Publication Volcano Plot
volcano_df <- res_lrt
volcano_df$sig <- ifelse(volcano_df$FDR < 0.05 & abs(volcano_df$logFC) > 1.0, "Significant", "Not Sig")

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(PValue), color = sig)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("Not Sig" = "grey", "Significant" = "red")) +
    theme_minimal(base_size = 12) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
        title = "Volcano Plot: Conditioned vs WT (LRT Method)",
        x     = "Log2 Fold Change",
        y     = "-log10(P-value)"
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(face = "bold", size = 12)
    )

ggsave(file.path(asset_dir, "004_pilot_lrt_volcano.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)

# 2. Z-Score Heatmap of Top 50 DEGs
top50_ids      <- rownames(final_sig_genes)[1:min(50, nrow(final_sig_genes))]
heatmap_matrix <- logCPM_6[top50_ids, , drop = FALSE]
heatmap_scaled <- t(scale(t(heatmap_matrix)))

png(file.path(asset_dir, "005_pilot_top50_heatmap.png"), width = 1200, height = 1500, res = 150)
pheatmap(
    heatmap_scaled,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = FALSE,
    annotation_col = phenodata[, "condition", drop = FALSE],
    color          = colorRampPalette(c("blue", "white", "red"))(100),
    border_color   = NA,
    main           = "Top 50 Differentially Expressed Transcripts (Z-score)"
)
dev.off()

cat(">>> [10] Visualization assets successfully exported to assets/\n")