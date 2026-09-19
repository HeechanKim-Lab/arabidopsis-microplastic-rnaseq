#!/usr/bin/env Rscript
# ==============================================================================
# Module: 03_interaction_visualization.R
# Description: Volcano and heatmap rendering for interaction and contrast suites
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")

TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
lrt_root01     <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

rdata_file <- file.path(base_dir, "data", "Arabidopsis_plastics_2023.RData")
if (!file.exists(rdata_file)) rdata_file <- "Arabidopsis_plastics_2023.RData"
load(rdata_file)

# 1. Publication Volcano Plot (Matching 06_volcanoplot.jpg format)
res01_inter_df <- as.data.frame(topTags(lrt_root01$res01_inter, n = Inf))
res01_inter_df$Gene <- rownames(res01_inter_df)
res01_inter_df$sig_status <- case_when(
    res01_inter_df$FDR < 0.05 & abs(res01_inter_df$logFC) >= 1.0 ~ "|log2 FC| > 1 and p.adjust < 0.05",
    abs(res01_inter_df$logFC) >= 1.0                              ~ "|log2 FC| > 1",
    TRUE                                                          ~ "Not significant"
)

n_up   <- sum(res01_inter_df$sig_status == "|log2 FC| > 1 and p.adjust < 0.05" & res01_inter_df$logFC > 0)
n_down <- sum(res01_inter_df$sig_status == "|log2 FC| > 1 and p.adjust < 0.05" & res01_inter_df$logFC < 0)

p_volcano <- ggplot(res01_inter_df, aes(x = logFC, y = -log10(PValue), color = sig_status)) +
    geom_point(alpha = 0.45, size = 1.8) +
    scale_color_manual(values = c(
        "Not significant" = "gray50",
        "|log2 FC| > 1" = "#2CA02C",
        "|log2 FC| > 1 and p.adjust < 0.05" = "#E41A1C"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
    annotate("text", x = 4.0, y = max(-log10(res01_inter_df$PValue)) * 0.95, 
             label = as.character(n_up), fontface = "bold", size = 6) +
    annotate("text", x = -4.0, y = max(-log10(res01_inter_df$PValue)) * 0.95, 
             label = as.character(n_down), fontface = "bold", size = 6) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        axis.title = element_text(face = "bold", size = 12)
    ) +
    labs(x = "Log2 fold change", y = "-Log10 adjusted p-value")

ggsave(file.path(assets_dir, "028_root00_lrt_volcano.png"), plot = p_volcano, width = 7, height = 7, dpi = 300)

# 2. Heatmap of 146 Interaction Genes
sig_inter_genes <- rownames(res01_inter_df[res01_inter_df$FDR < 0.05, ])
heatmap_inter_data <- TMM_root01_log[sig_inter_genes, ]
ann_col <- data.frame(Group = root01_group_02)
rownames(ann_col) <- colnames(heatmap_inter_data)

png(file.path(assets_dir, "029_root01_interaction_heatmap.png"), width = 1200, height = 1500, res = 150)
pheatmap(heatmap_inter_data, 
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         annotation_col = ann_col,
         main = "Incubation Signature: Plastic Response Over Time",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# 3. 4-Contrast LogFC Heatmap (Matching 05_heatmap.png)
sig_union <- unique(c(
    rownames(topTags(lrt_root01$res01_p0, n = Inf)$table[topTags(lrt_root01$res01_p0, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(lrt_root01$res01_p120, n = Inf)$table[topTags(lrt_root01$res01_p120, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(lrt_root01$res01_aging_c, n = Inf)$table[topTags(lrt_root01$res01_aging_c, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(lrt_root01$res01_aging_p, n = Inf)$table[topTags(lrt_root01$res01_aging_p, n = Inf)$table$FDR < 0.05, ])
))

mat_4c <- cbind(
    "4 vs. 2" = topTags(lrt_root01$res01_aging_p, n = Inf)$table[sig_union, "logFC"],
    "4 vs. 3" = topTags(lrt_root01$res01_p120, n = Inf)$table[sig_union, "logFC"],
    "2 vs. 1" = topTags(lrt_root01$res01_p0, n = Inf)$table[sig_union, "logFC"],
    "3 vs. 1" = topTags(lrt_root01$res01_aging_c, n = Inf)$table[sig_union, "logFC"]
)
rownames(mat_4c) <- sig_union

png(file.path(assets_dir, "030_root01_four_contrast_heatmap.png"), width = 1200, height = 1200, res = 150)
pheatmap(mat_4c,
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
         angle_col = 0)
dev.off()

cat(">>> [03] Volcano and Heatmap visualizations exported to assets/.\n")