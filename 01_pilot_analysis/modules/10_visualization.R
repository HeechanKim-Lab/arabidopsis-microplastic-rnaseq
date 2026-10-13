#!/usr/bin/env Rscript
# ==============================================================================
# Module 10: volcano plot and heatmap of the LRT genes (log 2026-04-09)
#
# Input : results/intermediate/res_lrt.rds, final_sig_genes.rds (module 09)
#         results/intermediate/logCPM_6samples.rds (module 08), phenodata_synced.rds (module 06)
# Output: assets/004_pilot_lrt_volcano.png, assets/005_pilot_top50_heatmap.png
# ==============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
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

res_lrt         <- readRDS(file.path(inter_dir, "res_lrt.rds"))
final_sig_genes <- readRDS(file.path(inter_dir, "final_sig_genes.rds"))
logCPM          <- readRDS(file.path(inter_dir, "logCPM_6samples.rds"))
phenodata       <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))

# Volcano plot (figure 004). log2FC is WT relative to Cond (see module 09).
plot_data <- res_lrt
plot_data$sig <- ifelse(plot_data$FDR < 0.05 & abs(plot_data$logFC) > 1, "Significant", "Not Sig")

p_volcano <- ggplot(plot_data, aes(x = logFC, y = -log10(PValue), color = sig)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggtitle("Volcano Plot: Conditioned vs WT (LRT Method)") +
    xlab("Log2 Fold Change") +
    ylab("-log10(P-value)")
ggsave(file.path(asset_dir, "004_pilot_lrt_volcano.png"), p_volcano, width = 7, height = 5, dpi = 200, bg = "white")

# Heatmap of the 50 genes with the smallest p-values (figure 005); rows scaled to Z-scores
top50_ids      <- head(rownames(final_sig_genes), 50)
heatmap_matrix <- logCPM[top50_ids, , drop = FALSE]
heatmap_scaled <- t(scale(t(heatmap_matrix)))

pheatmap(heatmap_scaled,
         cluster_rows   = TRUE,
         cluster_cols   = TRUE,
         show_rownames  = FALSE,
         annotation_col = phenodata[, "condition", drop = FALSE],
         main           = "Top 50 Differentially Expressed Transcripts (Z-score)",
         color          = colorRampPalette(c("blue", "white", "red"))(100),
         border_color   = NA,
         filename       = file.path(asset_dir, "005_pilot_top50_heatmap.png"),
         width = 6, height = 8)

cat("[10] Figures 004 and 005 written to", asset_dir, "\n")
