#!/usr/bin/env Rscript
# ==============================================================================
# Module: 05_multivariate_permanova.R
# Description: PERMANOVA dissimilarity testing, 2-way ANOVA, and boxplot generation
# ==============================================================================

suppressPackageStartupMessages({
    library(vegan)
    library(ggplot2)
    library(dplyr)
    library(multcompView)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")
tables_dir <- file.path(base_dir, "results", "tables")

TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
TMM_root01_pca <- readRDS(file.path(inter_dir, "TMM_root01_pca.rds"))

meta_root01 <- data.frame(
    Day     = factor(rep(c("0d", "120d"), each = 6)),
    Plastic = factor(rep(c("0p", "5p", "0p", "5p"), each = 3))
)
rownames(meta_root01) <- colnames(TMM_root01_log)

cat(">>> [05] Executing PERMANOVA (adonis2) with 999 permutations...\n")
set.seed(42)
perma_results <- adonis2(t(TMM_root01_log) ~ Day * Plastic, 
                         data = meta_root01, 
                         method = "euclidean",
                         by = "terms",
                         permutations = 999)

sink(file.path(tables_dir, "02_anova_result.txt"))
cat("################ PERMANOVA ################\n")
print(perma_results)

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

cat(">>> [05] Generating 2-way ANOVA publication boxplots...\n")
plot_df <- as.data.frame(TMM_root01_pca$x[, 1:3])
plot_df$Group <- gsub("-.*", "", rownames(plot_df))
plot_df$Condition <- factor(ifelse(plot_df$Group == "R1", "1",
                            ifelse(plot_df$Group == "R3", "2",
                            ifelse(plot_df$Group == "R7", "3", "4"))))

make_publication_boxplot <- function(df, pc_col, ylabel, out_png) {
    fit    <- aov(df[[pc_col]] ~ df$Condition)     tukey  <- TukeyHSD(fit)     let    <- multcompLetters4(fit, tukey)     let_df <- data.frame(         Condition = names(let$`df$Condition`$Letters),         Letter    = let$`df$Condition`$Letters
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
        theme_bw() +
        theme(
            panel.grid   = element_blank(),
            panel.border = element_rect(linewidth = 1.2),
            axis.title   = element_text(size = 14, face = "bold"),
            axis.text    = element_text(size = 12, color = "black"),
            axis.ticks   = element_line(linewidth = 0.8)
        ) +
        labs(x = "Conditions", y = ylabel)
    ggsave(out_png, plot = p, width = 7, height = 7, dpi = 300)
}

make_publication_boxplot(plot_df, "PC1", "PC1 score", file.path(assets_dir, "024_root01_anova_pc1.png"))
make_publication_boxplot(plot_df, "PC2", "PC2 score", file.path(assets_dir, "025_root01_anova_pc2.png"))
make_publication_boxplot(plot_df, "PC3", "PC3 score", file.path(assets_dir, "026_root01_anova_pc3.png"))

cat(">>> [05] PERMANOVA & ANOVA modeling complete.\n")