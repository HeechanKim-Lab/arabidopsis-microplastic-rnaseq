#!/usr/bin/env Rscript
# ==============================================================================
# Module 03: volcano plot and heatmaps of the Root01 interaction genes (log 2026-04-13_01)
#
# Input : results/intermediate/lrt_root01_list.rds (module 02), TMM_root01_log.rds (module 01)
#         data/Arabidopsis_plastics_2023.RData (root01_group_02 for the column annotation)
# Output: assets/028_root00_lrt_volcano.png (despite the file name, this is Root01),
#         assets/029_root01_interaction_heatmap.png, assets/030_root01_four_contrast_heatmap.png
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
    library(EnhancedVolcano)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
data_dir  <- file.path(base_dir, "data")
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

find_input <- function(name) {
    for (p in c(file.path(data_dir, name), file.path(getwd(), name))) if (file.exists(p)) return(p)
    stop("Input not found: ", name, " (expected in ", data_dir, ")")
}
load(find_input("Arabidopsis_plastics_2023.RData"))
TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
lrt_root01     <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

# ---- Volcano plot of the interaction contrast (figure 028) -------------------
# The y-axis shows -log10(FDR); the axis label is the package default ("-Log10 P").
res01_inter_df <- as.data.frame(topTags(lrt_root01$res01_inter, n = Inf))
n_inter <- sum(res01_inter_df$FDR < 0.05)

p_volcano <- EnhancedVolcano(res01_inter_df,
                             lab = rownames(res01_inter_df),
                             x = "logFC",
                             y = "FDR",
                             title = paste0("Plastic x Time Incubation (", n_inter, " Genes)"),
                             subtitle = "LRT Method",
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             pointSize = 3.0,
                             labSize = 4.0,
                             col = c("black", "black", "black", "red3"),
                             legendLabels = c("Not sig.", "LogFC", "FDR", "FDR & LogFC"),
                             caption = paste0("Total Genes: ", nrow(res01_inter_df)))
ggsave(file.path(asset_dir, "028_root00_lrt_volcano.png"), p_volcano,
       width = 9.5, height = 6, dpi = 200, bg = "white")

# ---- Heatmap of the interaction genes (figure 029) ---------------------------
sig_inter_genes <- rownames(res01_inter_df[res01_inter_df$FDR < 0.05, ])
heatmap_data    <- TMM_root01_log[sig_inter_genes, ]
ann_col <- data.frame(Group = root01_group_02)
rownames(ann_col) <- colnames(heatmap_data)

pheatmap(heatmap_data,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         annotation_col = ann_col,
         main = "Incubation Signature: Plastic Response Over Time",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = file.path(asset_dir, "029_root01_interaction_heatmap.png"),
         width = 7, height = 8)

# ---- Log2FC heatmap, union of four contrasts (figure 030) --------------------
# 1 = c0_d0, 2 = c5_d0, 3 = c0_d120, 4 = c5_d120
tab <- function(res) topTags(res, n = Inf)$table
sig_genes_list <- unique(c(
    rownames(tab(lrt_root01$res01_p0)[tab(lrt_root01$res01_p0)$FDR < 0.05, ]),
    rownames(tab(lrt_root01$res01_p120)[tab(lrt_root01$res01_p120)$FDR < 0.05, ]),
    rownames(tab(lrt_root01$res01_aging_c)[tab(lrt_root01$res01_aging_c)$FDR < 0.05, ]),
    rownames(tab(lrt_root01$res01_aging_p)[tab(lrt_root01$res01_aging_p)$FDR < 0.05, ])))

heatmap_matrix <- cbind(
    "4 vs. 2" = tab(lrt_root01$res01_aging_p)[sig_genes_list, "logFC"],
    "4 vs. 3" = tab(lrt_root01$res01_p120)[sig_genes_list, "logFC"],
    "2 vs. 1" = tab(lrt_root01$res01_p0)[sig_genes_list, "logFC"],
    "3 vs. 1" = tab(lrt_root01$res01_aging_c)[sig_genes_list, "logFC"])
rownames(heatmap_matrix) <- sig_genes_list

pheatmap(heatmap_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-4, 4, length.out = 101),
         main = NA,
         border_color = NA,
         fontsize_col = 14,
         angle_col = 0,
         filename = file.path(asset_dir, "030_root01_four_contrast_heatmap.png"),
         width = 6, height = 8)

cat("[03] Figures 028-030 written to", asset_dir, "\n")
