#!/usr/bin/env Rscript
# ==============================================================================
# Module: 01_full_matrix_and_subset_pca.R
# Description: Unsupervised PCA across Root00 (N=27) and Root01 (N=12) cohorts
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
data_dir   <- file.path(base_dir, "data")
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")

dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

rdata_file <- file.path(data_dir, "Arabidopsis_plastics_2023.RData")
if (!file.exists(rdata_file)) rdata_file <- "Arabidopsis_plastics_2023.RData"
load(rdata_file)

cat(">>> [01] Running PCA for Root00 (Full 27-sample matrix)...\n")
TMM_root00_log <- cpm(root00_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root00_pca <- prcomp(t(TMM_root00_log), center = TRUE, scale. = TRUE)
sum_pca00      <- summary(TMM_root00_pca)
pc1_v00        <- round(sum_pca00$importance[2, 1] * 100, 1)
pc2_v00        <- round(sum_pca00$importance[2, 2] * 100, 1)
pc3_v00        <- round(sum_pca00$importance[2, 3] * 100, 1)

# Root00 PC1 vs PC2
png(file.path(assets_dir, "020_root00_full_matrix_pca.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root00_pca$x[, 1], TMM_root00_pca$x[, 2],
     pch = root00_pca_points, col = TMM_root00_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_v00, "%)"),
     ylab = paste0("PC2 (", pc2_v00, "%)"),
     main = "PCA: Full 27-Sample Root Dataset")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

# Root00 PC1 vs PC3
png(file.path(assets_dir, "021_root00_pca_pc1_vs_pc3.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root00_pca$x[, 1], TMM_root00_pca$x[, 3],
     pch = root00_pca_points, col = TMM_root00_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_v00, "%)"),
     ylab = paste0("PC3 (", pc3_v00, "%)"),
     main = "PCA: PC1 vs. PC3 (root00)")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

cat(">>> [01] Running PCA for Root01 (Characteristic 12-sample subset)...\n")
TMM_root01_log <- cpm(root01_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root01_pca <- prcomp(t(TMM_root01_log), center = TRUE, scale. = TRUE)
sum_pca01      <- summary(TMM_root01_pca)
pc1_v01        <- round(sum_pca01$importance[2, 1] * 100, 2)
pc2_v01        <- round(sum_pca01$importance[2, 2] * 100, 2)

# Root01 Dose-Response Diagnostic Plot
png(file.path(assets_dir, "022_root01_dose_response_pca.png"), width = 1200, height = 1000, res = 150)
plot(TMM_root01_pca$x[, 1], TMM_root01_pca$x[, 2],
     pch = root01_pca_points, col = TMM_root01_subtype.col.vec01,
     cex = 1.6, lwd = 1.5,
     xlab = paste0("PC1 (", round(pc1_v01, 1), "%)"),
     ylab = paste0("PC2 (", round(pc2_v01, 1), "%)"),
     main = "PCA: Characteristic 12-Sample Subset")
legend("topright", legend = c("0% Plastic", "1% Plastic", "5% Plastic"),
       col = c("black", "red", "blue"), pch = 16, bty = "n")
dev.off()

# Publication Monochrome 2D Plot
sample_groups <- gsub("-.*", "", rownames(TMM_root01_pca$x))
pub_pch <- ifelse(sample_groups == "R1", 1,
           ifelse(sample_groups == "R3", 16,
           ifelse(sample_groups == "R7", 0, 15)))

png(file.path(assets_dir, "023_root01_publication_pca_2d.png"), width = 1200, height = 1200, res = 150)
par(mar = c(5, 5, 4, 2), font.main = 2)
plot(TMM_root01_pca$x[, 1], TMM_root01_pca$x[, 2],
     pch = pub_pch, col = "black", cex = 2.2, lwd = 1.5,
     xlab = paste0("PC1 (", pc1_v01, " %)"),
     ylab = paste0("PC2 (", pc2_v01, " %)"),
     cex.lab = 1.2, las = 1, asp = 1)
legend("topright",
       legend = c("Condition 1", "Condition 2", "Condition 3", "Condition 4"),
       pch = c(1, 16, 0, 15), pt.cex = 1.6, cex = 0.9, bty = "n")
dev.off()

# Serialize intermediate structures
saveRDS(TMM_root00_log, file.path(inter_dir, "TMM_root00_log.rds"))
saveRDS(TMM_root00_pca, file.path(inter_dir, "TMM_root00_pca.rds"))
saveRDS(TMM_root01_log, file.path(inter_dir, "TMM_root01_log.rds"))
saveRDS(TMM_root01_pca, file.path(inter_dir, "TMM_root01_pca.rds"))

cat(">>> [01] PCA processing completed and serialized.\n")