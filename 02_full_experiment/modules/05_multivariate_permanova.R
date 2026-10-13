#!/usr/bin/env Rscript
# ==============================================================================
# Module 05: PERMANOVA, two-way ANOVA on PC scores and PCA figure for Root01 (log 2026-04-15)
#
# Input : results/intermediate/TMM_root01_log.rds, TMM_root01_pca.rds (module 01)
# Output: assets/023_root01_publication_pca_2d.png
#         assets/024_root01_anova_pc1.png, assets/025_root01_anova_pc2.png, assets/026_root01_anova_pc3.png
#         results/tables/permanova_anova_root01.txt
#
# Conditions: 1 = c0_d0 (R1), 2 = c5_d0 (R3), 3 = c0_d120 (R7), 4 = c5_d120 (R9).
# ==============================================================================

suppressPackageStartupMessages({
    library(vegan)
    library(ggplot2)
    library(dplyr)
    library(multcompView)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir   <- get_base_dir()
inter_dir  <- file.path(base_dir, "results", "intermediate")
tables_dir <- file.path(base_dir, "results", "tables")
asset_dir  <- file.path(base_dir, "assets")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
TMM_root01_pca <- readRDS(file.path(inter_dir, "TMM_root01_pca.rds"))

# The metadata below are assigned by column position: R1, R3, R7, R9, three libraries each
stopifnot(identical(gsub("-.*", "", colnames(TMM_root01_log)), rep(c("R1", "R3", "R7", "R9"), each = 3)))

# ---- PCA, black and white (figure 023) ----------------------------------------
pca_sum <- summary(TMM_root01_pca)
pc1_var <- round(pca_sum$importance[2, 1] * 100, 2)
pc2_var <- round(pca_sum$importance[2, 2] * 100, 2)
sample_groups <- gsub("-.*", "", rownames(TMM_root01_pca$x))
pub_pch <- ifelse(sample_groups == "R1", 1,
           ifelse(sample_groups == "R3", 16,
           ifelse(sample_groups == "R7", 0, 15)))

png(file.path(asset_dir, "023_root01_publication_pca_2d.png"), width = 1200, height = 1200, res = 150)
par(mar = c(5, 5, 4, 2), font.main = 2)
plot(TMM_root01_pca$x[, 1], TMM_root01_pca$x[, 2],
     pch = pub_pch, col = "black", cex = 2.2, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_var, " %)"),
     ylab = paste0("PC2 (", pc2_var, " %)"),
     cex.lab = 1.2, las = 1, asp = 1)
legend("topright",
       legend = c("Condition 1", "Condition 2", "Condition 3", "Condition 4"),
       pch = c(1, 16, 0, 15), pt.cex = 1.6, cex = 0.85, bty = "n")
invisible(dev.off())

# ---- PERMANOVA on Euclidean distances of TMM log-CPM ---------------------------
meta_root01 <- data.frame(
    Day     = factor(rep(c("0d", "120d"), each = 6)),
    Plastic = factor(rep(c("0p", "5p", "0p", "5p"), each = 3)))
rownames(meta_root01) <- colnames(TMM_root01_log)

# 999 permutations (default). No seed was set in the original run, so p-values can
# differ slightly between runs.
perma_results <- adonis2(t(TMM_root01_log) ~ Day * Plastic,
                         data = meta_root01, method = "euclidean", by = "terms")

# ---- Two-way ANOVA on the PC1-PC3 scores ---------------------------------------
anova_data <- cbind(as.data.frame(TMM_root01_pca$x[, 1:3]), meta_root01)
anova_pc1 <- aov(PC1 ~ Day * Plastic, data = anova_data)
anova_pc2 <- aov(PC2 ~ Day * Plastic, data = anova_data)
anova_pc3 <- aov(PC3 ~ Day * Plastic, data = anova_data)

report <- c("PERMANOVA (adonis2, Euclidean, ~ Day * Plastic, terms added sequentially)",
            capture.output(print(perma_results)),
            "", "Two-way ANOVA, PC1", capture.output(summary(anova_pc1)),
            "", "Two-way ANOVA, PC2", capture.output(summary(anova_pc2)),
            "", "Two-way ANOVA, PC3", capture.output(summary(anova_pc3)))
writeLines(report, file.path(tables_dir, "permanova_anova_root01.txt"))
cat(report, sep = "\n")

# ---- Boxplots with Tukey letters (figures 024-026) ---------------------------
plot_df <- as.data.frame(TMM_root01_pca$x[, 1:3])
plot_df$Group <- gsub("-.*", "", rownames(plot_df))
plot_df$Condition <- factor(ifelse(plot_df$Group == "R1", "1",
                            ifelse(plot_df$Group == "R3", "2",
                            ifelse(plot_df$Group == "R7", "3", "4"))))

make_publication_boxplot <- function(df, pc_col, ylabel) {
    # One-way ANOVA across the four conditions, Tukey HSD, compact letter display
    fit   <- aov(df[[pc_col]] ~ df$Condition)
    tukey <- TukeyHSD(fit)
    letters_list <- multcompLetters4(fit, tukey)
    letters_df <- data.frame(
        Condition = names(letters_list$`df$Condition`$Letters),
        Letter    = letters_list$`df$Condition`$Letters)

    # Letters are placed 10% of the data range above the highest point of each condition
    pos_df <- df %>%
        group_by(Condition) %>%
        summarise(y_max = max(.data[[pc_col]]) + (max(df[[pc_col]]) - min(df[[pc_col]])) * 0.1) %>%
        left_join(letters_df, by = "Condition")

    ggplot(df, aes(x = Condition, y = .data[[pc_col]])) +
        geom_boxplot(outlier.shape = NA, width = 0.5, color = "black", linewidth = 0.8) +
        geom_jitter(color = "darkgrey", width = 0.15, size = 2, alpha = 0.8) +
        geom_text(data = pos_df, aes(x = Condition, y = y_max, label = Letter),
                  size = 5, fontface = "bold", vjust = 0) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(linewidth = 1.2),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, color = "black"),
            axis.ticks = element_line(linewidth = 0.8)) +
        labs(x = "Conditions", y = ylabel)
}

for (k in 1:3) {
    p <- make_publication_boxplot(plot_df, paste0("PC", k), paste0("PC", k, " score"))
    ggsave(file.path(asset_dir, sprintf("%03d_root01_anova_pc%d.png", 23 + k, k)), p,
           width = 6, height = 6, dpi = 200, bg = "white")
}

cat("[05] Figures 023-026 written to", asset_dir, "\n")
