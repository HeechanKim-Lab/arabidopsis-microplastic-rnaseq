#!/usr/bin/env Rscript
# ==============================================================================
# Script: run_full_experiment.R
# Project: arabidopsis-microplastic-rnaseq (02_full_experiment)
# Description: Monolithic production pipeline executing PCA ordination, LRT 
#              factorial modeling, PERMANOVA/ANOVA, Venn calculations, 
#              functional enrichment, and targeted pathway heatmaps.
# ==============================================================================

# --- 1. Dependencies and Environment Setup ---
suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(pheatmap)
    library(gridExtra)
    library(clusterProfiler)
    library(org.At.tair.db)
    library(vegan)
    library(multcompView)
    library(readxl)
})

# Establish repository path anchors
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
    base_dir <- normalizePath(file.path(dirname(script_path), ".."))
} else {
    base_dir <- getwd()
}

data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
diff_dir    <- file.path(results_dir, "differential")
tables_dir  <- file.path(results_dir, "tables")
assets_dir  <- file.path(base_dir, "assets")

dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

# Ingest RData and Reference Annotation Table
rdata_file <- file.path(data_dir, "Arabidopsis_plastics_2023.RData")
anno_file  <- file.path(data_dir, "Arabidopsis_genes_description.xlsx")

if (!file.exists(rdata_file)) {
    # Fallback to local working directory structure if not in data/
    if (file.exists("Arabidopsis_plastics_2023.RData")) {
        rdata_file <- "Arabidopsis_plastics_2023.RData"
    } else {
        stop("FATAL: Arabidopsis_plastics_2023.RData not found.")
    }
}
load(rdata_file)

if (!file.exists(anno_file)) {
    if (file.exists("Arabidopsis_genes_description.xlsx")) {
        anno_file <- "Arabidopsis_genes_description.xlsx"
    } else {
        stop("FATAL: Arabidopsis_genes_description.xlsx not found.")
    }
}
anno_db <- read_excel(anno_file)


# ==============================================================================
# PHASE 1: Principal Component Analysis (Root00 & Root01)
# ==============================================================================
cat("=== [PHASE 1] Unsupervised Principal Component Analysis ===\n")

# 1. Global 27-Sample Dataset (Root00)
TMM_root00_log <- cpm(root00_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root00_pca <- prcomp(t(TMM_root00_log), center = TRUE, scale. = TRUE)
sum_pca00      <- summary(TMM_root00_pca)
pc1_var00      <- round(sum_pca00$importance[2, 1] * 100, 1)
pc2_var00      <- round(sum_pca00$importance[2, 2] * 100, 1)
pc3_var00      <- round(sum_pca00$importance[2, 3] * 100, 1)

# Export Root00 PC1 vs PC2
png(file.path(assets_dir, "020_root00_full_matrix_pca.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root00_pca$x[, 1], TMM_root00_pca$x[, 2],
     pch = root00_pca_points, col = TMM_root00_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_var00, "%)"),
     ylab = paste0("PC2 (", pc2_var00, "%)"),
     main = "PCA: Full 27-Sample Root Dataset")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

# Export Root00 PC1 vs PC3
png(file.path(assets_dir, "021_root00_pca_pc1_vs_pc3.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root00_pca$x[, 1], TMM_root00_pca$x[, 3],
     pch = root00_pca_points, col = TMM_root00_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_var00, "%)"),
     ylab = paste0("PC3 (", pc3_var00, "%)"),
     main = "PCA: PC1 vs. PC3 (Root00)")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

# 2. Characteristic 12-Sample Subset (Root01)
TMM_root01_log <- cpm(root01_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root01_pca <- prcomp(t(TMM_root01_log), center = TRUE, scale. = TRUE)
sum_pca01      <- summary(TMM_root01_pca)
pc1_var01      <- round(sum_pca01$importance[2, 1] * 100, 2)
pc2_var01      <- round(sum_pca01$importance[2, 2] * 100, 2)

# Export Root01 Dose-Response Diagnostic PCA
png(file.path(assets_dir, "022_root01_dose_response_pca.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root01_pca$x[, 1], TMM_root01_pca$x[, 2],
     pch = root01_pca_points, col = TMM_root01_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", round(pc1_var01, 1), "%)"),
     ylab = paste0("PC2 (", round(pc2_var01, 1), "%)"),
     main = "PCA: Characteristic 12-Sample Subset")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

# Export Publication 2D Monochrome PCA
sample_groups <- gsub("-.*", "", rownames(TMM_root01_pca$x))
pub_pch <- ifelse(sample_groups == "R1", 1,
           ifelse(sample_groups == "R3", 16,
           ifelse(sample_groups == "R7", 0, 15)))

png(file.path(assets_dir, "023_root01_publication_pca_2d.png"), width = 1200, height = 1200, res = 150)
par(mar = c(5, 5, 4, 2), font.main = 2)
plot(TMM_root01_pca$x[, 1], TMM_root01_pca$x[, 2],
     pch = pub_pch, col = "black", cex = 2.2, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_var01, " %)"),
     ylab = paste0("PC2 (", pc2_var01, " %)"),
     cex.lab = 1.2, las = 1, asp = 1)
legend("topright",
       legend = c("Condition 1", "Condition 2", "Condition 3", "Condition 4"),
       pch = c(1, 16, 0, 15), pt.cex = 1.6, cex = 0.9, bty = "n")
dev.off()


# ==============================================================================
# PHASE 2: Differential Expression Modeling via LRT
# ==============================================================================
cat("\n=== [PHASE 2] LRT Differential Expression Modeling ===\n")

# 1. Root00 LRT Modeling (27 Samples)
root00_fit_LRT <- glmFit(root00_y4, root00_design_01)
contrasts_root00 <- makeContrasts(
    Plastic_5p_d0   = c5_d0   - c0_d0,
    Plastic_5p_d30  = c5_d30  - c0_d30,
    Plastic_5p_d120 = c5_d120 - c0_d120,
    Plastic_1p_d120 = c1_d120 - c0_d120,
    Incubation_5p   = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    levels          = root00_design_01
)

res00_0d    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d0"])
res00_30d   <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d30"])
res00_120d  <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d120"])
res00_1p    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_1p_d120"])
res00_inter <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Incubation_5p"])

# 2. Root01 LRT Modeling (12 Samples)
root01_fit_LRT <- glmFit(root01_y4, root01_design_01)
contrasts_root01 <- makeContrasts(
    p5_d0           = c5_d0   - c0_d0,
    p5_d120         = c5_d120 - c0_d120,
    d120_p0         = c0_d120 - c0_d0,
    d120_p5         = c5_d120 - c5_d0,
    plastic_5p_120d = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    incubation      = (c5_d120 - c0_d120) - (c5_d0 - c0_d0),
    levels          = root01_design_01
)

res01_p0      <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d0"])
res01_p120    <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d120"])
res01_aging_c <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p0"])
res01_aging_p <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p5"])
res01_inter   <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "plastic_5p_120d"])

# 3. Export Summary Barplot Across Subsets
summary_counts <- data.frame(
    Dataset = c(rep("Full Dataset (Root00)", 5), rep("Subset (Root01)", 6)),
    Contrast = c(
        "D0 Plastic Effect (5%)", "D30 Plastic Effect (5%)", "D120 Plastic Effect (5%)", 
        "D120 Plastic Effect (1%)", "Overall Incubation",
        "D0 Plastic Effect", "D120 Plastic Effect", 
        "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)", 
        "Plastic x Time Interaction", "Incubation Effect"
    ),
    Count = c(
        sum(topTags(res00_0d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_30d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_120d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_1p, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_inter, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_p0, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_p120, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_aging_c, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_aging_p, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_inter, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_inter, n=Inf)$table$FDR < 0.05)
    )
)

summary_counts$Contrast <- factor(summary_counts$Contrast, levels = rev(c(
    "D0 Plastic Effect", "D0 Plastic Effect (5%)", 
    "D30 Plastic Effect (5%)", 
    "D120 Plastic Effect (1%)", "D120 Plastic Effect (5%)", "D120 Plastic Effect", 
    "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)", 
    "Plastic x Time Interaction", "Overall Incubation", "Incubation Effect"
)))

p_bar <- ggplot(summary_counts, aes(x = Contrast, y = Count, fill = Dataset)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.3, size = 4.2, fontface = "bold") +
    coord_flip() +
    facet_wrap(~ Dataset, scales = "free") +
    scale_fill_manual(values = c("Full Dataset (Root00)" = "#B3E2CD", "Subset (Root01)" = "#FDCDAC")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Significant DEGs Yield (FDR < 0.05)", x = "", y = "Number of Significant Genes")

ggsave(file.path(assets_dir, "027_deg_contrasts_summary.png"), plot = p_bar, width = 10, height = 6, dpi = 300)


# ==============================================================================
# PHASE 3: Visualizing the Interaction Term (Volcano & Heatmap)
# ==============================================================================
cat("\n=== [PHASE 3] Factorial Interaction Visualizations ===\n")

# 1. Publication Volcano Plot (146 Genes)
res01_inter_df <- as.data.frame(topTags(res01_inter, n = Inf))
res01_inter_df$Gene <- rownames(res01_inter_df)
res01_inter_df$sig_status <- case_when(
    res01_inter_df$FDR < 0.05 & abs(res01_inter_df$logFC) >= 1.0 ~ "|log2 FC| > 1 and p.adjust < 0.05",
    abs(res01_inter_df$logFC) >= 1.0                              ~ "|log2 FC| > 1",
    TRUE                                                          ~ "Not significant"
)

n_up   <- sum(res01_inter_df$sig_status == "|log2 FC| > 1 and p.adjust < 0.05" & res01_inter_df$logFC > 0)
n_down <- sum(res01_inter_df$sig_status == "|log2 FC| > 1 and p.adjust < 0.05" & res01_inter_df$logFC < 0)

p_volc <- ggplot(res01_inter_df, aes(x = logFC, y = -log10(PValue), color = sig_status)) +
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
        plot.title = element_blank(),
        axis.title = element_text(face = "bold", size = 12)
    ) +
    labs(x = "Log2 fold change", y = "-Log10 adjusted p-value")

ggsave(file.path(assets_dir, "028_root00_lrt_volcano.png"), plot = p_volc, width = 7, height = 7, dpi = 300)

# 2. Heatmap of the 146 Interaction Genes
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

# 3. 4-Contrast LogFC Heatmap
sig_genes_union <- unique(c(
    rownames(topTags(res01_p0, n = Inf)$table[topTags(res01_p0, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(res01_p120, n = Inf)$table[topTags(res01_p120, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(res01_aging_c, n = Inf)$table[topTags(res01_aging_c, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(res01_aging_p, n = Inf)$table[topTags(res01_aging_p, n = Inf)$table$FDR < 0.05, ])
))

mat_4contrast <- cbind(
    "4 vs. 2" = topTags(res01_aging_p, n = Inf)$table[sig_genes_union, "logFC"],
    "4 vs. 3" = topTags(res01_p120, n = Inf)$table[sig_genes_union, "logFC"],
    "2 vs. 1" = topTags(res01_p0, n = Inf)$table[sig_genes_union, "logFC"],
    "3 vs. 1" = topTags(res01_aging_c, n = Inf)$table[sig_genes_union, "logFC"]
)
rownames(mat_4contrast) <- sig_genes_union

png(file.path(assets_dir, "030_root01_four_contrast_heatmap.png"), width = 1200, height = 1200, res = 150)
pheatmap(mat_4contrast,
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


# ==============================================================================
# PHASE 4: Functional Enrichment Analysis (GO & KEGG)
# ==============================================================================
cat("\n=== [PHASE 4] Functional Pathway Deconvolution ===\n")

# 1. 146 Interaction Genes GO & KEGG
go_inter_BP <- enrichGO(gene = sig_inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
go_inter_MF <- enrichGO(gene = sig_inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "MF", pvalueCutoff = 0.05)
go_inter_CC <- enrichGO(gene = sig_inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "CC", pvalueCutoff = 0.05)

ggsave(file.path(assets_dir, "035_interaction_go_bp.png"), dotplot(go_inter_BP, showCategory = 15), width = 8, height = 7, dpi = 300)
ggsave(file.path(assets_dir, "036_interaction_go_mf.png"), dotplot(go_inter_MF, showCategory = 15), width = 8, height = 7, dpi = 300)
ggsave(file.path(assets_dir, "037_interaction_go_cc.png"), dotplot(go_inter_CC, showCategory = 15), width = 8, height = 7, dpi = 300)

kegg_inter <- enrichKEGG(gene = sig_inter_genes, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
if (!is.null(kegg_inter) && nrow(as.data.frame(kegg_inter)) > 0) {
    kegg_summary <- as.data.frame(table(kegg_inter$ID))
    colnames(kegg_summary) <- c("Pathway", "GeneCount")
    kegg_descriptions <- c(
        "ath01100" = "Metabolic pathways",
        "ath01110" = "Biosynthesis of secondary metabolites",
        "ath04120" = "Ubiquitin mediated proteolysis",
        "ath04141" = "Protein processing in endoplasmic reticulum",
        "ath04075" = "Plant hormone signal transduction",
        "ath01200" = "Carbon metabolism",
        "ath00940" = "Phenylpropanoid biosynthesis",
        "ath00270" = "Cysteine and methionine metabolism",
        "ath00030" = "Pentose phosphate pathway",
        "ath00010" = "Glycolysis / Gluconeogenesis"
    )
    kegg_summary$Description <- kegg_descriptions[as.character(kegg_summary$Pathway)]
    kegg_summary <- kegg_summary[!is.na(kegg_summary$Description), ]
    kegg_summary <- kegg_summary[order(-kegg_summary$GeneCount), ]
    
    p_kegg <- ggplot(head(kegg_summary, 10), aes(x = reorder(Description, GeneCount), y = GeneCount)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        theme_minimal(base_size = 12) +
        labs(title = "Top 10 KEGG Pathways (Incubation Genes)",
             subtitle = paste0(nrow(as.data.frame(kegg_inter)), "/146 genes mapped"),
             x = "Pathway Name", y = "Number of Genes")
    ggsave(file.path(assets_dir, "038_interaction_kegg_top10.png"), plot = p_kegg, width = 8, height = 6, dpi = 300)
}

# 2. Plastic Aging (5%) GO Enrichment
res01_aging_p_df  <- as.data.frame(topTags(res01_aging_p, n = Inf))
sig_aging_p_genes <- rownames(res01_aging_p_df[res01_aging_p_df$FDR < 0.05, ])
go_aging_p_BP     <- enrichGO(gene = sig_aging_p_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
ggsave(file.path(assets_dir, "034_aging_5p_go_bp.png"), dotplot(go_aging_p_BP, showCategory = 20), width = 8, height = 7, dpi = 300)

# 3. Formal GO Overrepresentation Tables (Cond 4 vs 2)
up_4vs2   <- rownames(res01_aging_p_df[res01_aging_p_df$logFC > 1  & res01_aging_p_df$PValue < 0.01 & res01_aging_p_df$FDR < 0.05, ])
down_4vs2 <- rownames(res01_aging_p_df[res01_aging_p_df$logFC < -1 & res01_aging_p_df$PValue < 0.01 & res01_aging_p_df$FDR < 0.05, ])

run_go_table <- function(genes) {
    res <- enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    as.data.frame(res) %>%
        mutate(
            obs = as.numeric(sub("/.*", "", GeneRatio)),
            test = as.numeric(sub(".*/", "", GeneRatio)),
            exp = as.numeric(sub("/.*", "", BgRatio)),
            ref = as.numeric(sub(".*/", "", BgRatio)),
            Enr = round((obs / test) / (exp / ref), 2),
            "No. gene" = paste0(obs, "/", exp)
        ) %>%
        dplyr::select(Description, `No. gene`, Enr, p.adjust) %>%
        head(5)
}

table_4vs2_up   <- run_go_table(up_4vs2)
table_4vs2_down <- run_go_table(down_4vs2)
write.csv(table_4vs2_up,   file.path(tables_dir, "08_go_bp_4vs2_up.csv"), row.names = FALSE)
write.csv(table_4vs2_down, file.path(tables_dir, "08_go_bp_4vs2_down.csv"), row.names = FALSE)


# ==============================================================================
# PHASE 5: Multivariate PERMANOVA & Latent ANOVA
# ==============================================================================
cat("\n=== [PHASE 5] Multivariate Statistical Hypothesis Testing ===\n")

meta_root01 <- data.frame(
    Day     = factor(rep(c("0d", "120d"), each = 6)),
    Plastic = factor(rep(c("0p", "5p", "0p", "5p"), each = 3))
)
rownames(meta_root01) <- colnames(TMM_root01_log)

# 1. PERMANOVA via adonis2
set.seed(42)
perma_results <- adonis2(t(TMM_root01_log) ~ Day * Plastic, 
                         data = meta_root01, 
                         method = "euclidean",
                         by = "terms",
                         permutations = 999)

sink(file.path(tables_dir, "02_anova_result.txt"))
cat("################ PERMANOVA ################\n")
print(perma_results)

# 2. Latent PC 2-Way ANOVAs
pca_scores <- as.data.frame(TMM_root01_pca$x[, 1:3])
anova_data <- cbind(pca_scores, meta_root01)

cat("\n################ 2-WAY ANOVA ################\n")
cat("--- 2-way ANOVA Results for PC1 ---\n")
print(summary(aov(PC1 ~ Day * Plastic, data = anova_data)))
cat("\n--- 2-way ANOVA Results for PC2 ---\n")
print(summary(aov(PC2 ~ Day * Plastic, data = anova_data)))
cat("\n--- 2-way ANOVA Results for PC3 ---\n")
print(summary(aov(PC3 ~ Day * Plastic, data = anova_data)))
sink()

# 3. Latent PC Boxplots with Tukey HSD Letters
plot_df <- as.data.frame(TMM_root01_pca$x[, 1:3])
plot_df$Group <- gsub("-.*", "", rownames(plot_df))
plot_df$Condition <- factor(ifelse(plot_df$Group == "R1", "1",
                            ifelse(plot_df$Group == "R3", "2",
                            ifelse(plot_df$Group == "R7", "3", "4"))))

make_publication_boxplot <- function(df, pc_col, ylabel, out_png) {
    fit    <- aov(df[[pc_col]] ~ df$Condition)
    tukey  <- TukeyHSD(fit)
    let    <- multcompLetters4(fit, tukey)
    let_df <- data.frame(
        Condition = names(let$`df$Condition`$Letters),
        Letter    = let$`df$Condition`$Letters
    )
    pos_df <- df %>% 
        group_by(Condition) %>% 
        summarise(y_max = max(get(pc_col)) + (max(df[[pc_col]]) - min(df[[pc_col]])) * 0.1) %>%
        left_join(let_df, by = "Condition")

    p <- ggplot(df, aes_string(x = "Condition", y = pc_col)) +
        geom_boxplot(outlier.shape = NA, width = 0.5, color = "black", linewidth = 0.8) +
        geom_jitter(color = "darkgrey", width = 0.15, size = 2, alpha = 0.8) +
        geom_text(data = pos_df, aes(x = Condition, y = y_max, label = Letter), 
                  size = 5, fontface = "bold", vjust = 0) +
        theme_bw(base_size = 14) +
        theme(
            panel.grid = element_blank(),
            axis.title = element_text(face = "bold", size = 14),
            axis.text  = element_text(color = "black", size = 12)
        ) +
        labs(x = "Conditions", y = ylabel)
    ggsave(out_png, plot = p, width = 7, height = 7, dpi = 300)
}

make_publication_boxplot(plot_df, "PC1", "PC1 score", file.path(assets_dir, "024_root01_anova_pc1.png"))
make_publication_boxplot(plot_df, "PC2", "PC2 score", file.path(assets_dir, "025_root01_anova_pc2.png"))
make_publication_boxplot(plot_df, "PC3", "PC3 score", file.path(assets_dir, "026_root01_anova_pc3.png"))


# ==============================================================================
# PHASE 6: Targeted Pathway Modules (ROS & Photosynthesis/Plastid)
# ==============================================================================
cat("\n=== [PHASE 6] Targeted Physiological Pathway Mining ===\n")

# Combine all significant LRT loci across Root01
all_lrt_sig <- unique(c(
    rownames(topTags(res01_p120, n = Inf)$table[topTags(res01_p120, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(res01_inter, n = Inf)$table[topTags(res01_inter, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(res01_aging_p, n = Inf)$table[topTags(res01_aging_p, n = Inf)$table$FDR < 0.05, ])
))

anno_sig <- anno_db %>% filter(Gene_locus %in% all_lrt_sig)

# 1. Oxidative Stress / ROS Heatmap (12 Genes)
ox_keywords <- "hydrogen peroxide|oxidative|reactive oxygen|peroxidase"
ox_genes    <- anno_sig %>% filter(grepl(ox_keywords, Description, ignore.case = TRUE))

if (nrow(ox_genes) > 1) {
    ox_matrix <- TMM_root01_log[ox_genes$Gene_locus, ]
    rownames(ox_matrix) <- ifelse(is.na(ox_genes$Symbol) | ox_genes$Symbol == "", 
                                  ox_genes$Gene_locus, ox_genes$Symbol)
    
    png(file.path(assets_dir, "039_targeted_ros_heatmap.png"), width = 1200, height = 900, res = 150)
    pheatmap(ox_matrix, 
             scale = "row", 
             annotation_col = ann_col,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Oxidative Stress Responsiveness (", nrow(ox_genes), " Genes)"),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize_row = 9)
    dev.off()
}

# 2. Photosynthesis & Plastid Shifts Heatmap (7 Genes)
photo_keywords <- "photosynthe|chloroplast|photosystem|light-harvesting|thylakoid"
photo_genes    <- anno_sig %>% filter(grepl(photo_keywords, Description, ignore.case = TRUE))

if (nrow(photo_genes) > 1) {
    photo_matrix <- TMM_root01_log[photo_genes$Gene_locus, ]
    rownames(photo_matrix) <- ifelse(is.na(photo_genes$Symbol) | photo_genes$Symbol == "", 
                                     photo_genes$Gene_locus, photo_genes$Symbol)
    
    png(file.path(assets_dir, "040_targeted_photosynthesis_heatmap.png"), width = 1200, height = 900, res = 150)
    pheatmap(photo_matrix, 
             scale = "row", 
             annotation_col = ann_col,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Photosynthesis & Plastid Shifts (", nrow(photo_genes), " Genes)"),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize_row = 9)
    dev.off()
}


# ==============================================================================
# PHASE 7: Master Annotated DEG Registry Export
# ==============================================================================
cat("\n=== [PHASE 7] Master Annotation Spreadsheet Generation ===\n")

list_00_p120    <- rownames(topTags(res00_120d, n = Inf)$table[topTags(res00_120d, n = Inf)$table$FDR < 0.05, ])
list_00_inter   <- rownames(topTags(res00_inter, n = Inf)$table[topTags(res00_inter, n = Inf)$table$FDR < 0.05, ])
list_01_p120    <- rownames(topTags(res01_p120, n = Inf)$table[topTags(res01_p120, n = Inf)$table$FDR < 0.05, ])
list_01_inter   <- rownames(topTags(res01_inter, n = Inf)$table[topTags(res01_inter, n = Inf)$table$FDR < 0.05, ])
list_01_aging_p <- rownames(topTags(res01_aging_p, n = Inf)$table[topTags(res01_aging_p, n = Inf)$table$FDR < 0.05, ])

all_sig_genes <- unique(c(list_00_p120, list_00_inter, list_01_p120, list_01_inter, list_01_aging_p))
master_report <- data.frame(Gene_locus = all_sig_genes)

master_report$root00_Plastic_D120 <- ifelse(master_report$Gene_locus %in% list_00_p120, "Yes", "")
master_report$root00_Incubation   <- ifelse(master_report$Gene_locus %in% list_00_inter, "Yes", "")
master_report$root00_Intersection <- ifelse(master_report$root00_Plastic_D120 == "Yes" & 
                                            master_report$root00_Incubation == "Yes", "Yes", "")

master_report$root01_Plastic_D120 <- ifelse(master_report$Gene_locus %in% list_01_p120, "Yes", "")
master_report$root01_Incubation   <- ifelse(master_report$Gene_locus %in% list_01_inter, "Yes", "")
master_report$root01_Aging_5p     <- ifelse(master_report$Gene_locus %in% list_01_aging_p, "Yes", "")
master_report$root01_Triple_Intersection <- ifelse(
    master_report$root01_Plastic_D120 == "Yes" & 
    master_report$root01_Incubation == "Yes" & 
    master_report$root01_Aging_5p == "Yes", "Yes", ""
)

master_final <- left_join(master_report, anno_db, by = "Gene_locus") %>%
    dplyr::select(Gene_locus, Symbol, Description, everything())

write.csv(master_final, file.path(tables_dir, "Final_Core_Plastic_Response_Annotated.csv"), row.names = FALSE)
cat(">>> [PHASE 7] Exported master annotated report with", nrow(master_final), "loci.\n")
cat("=== Full Experiment Pipeline Completed Successfully ===\n")