#!/usr/bin/env Rscript
# ==============================================================================
# Module: 09_differential_analysis_lrt.R
# Description: Outlier excision, QLF vs LRT evaluation, and final DEG extraction
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(reshape2)
    library(ggplot2)
})

base_dir  <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
diff_dir  <- file.path(base_dir, "results", "differential")
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)

counts    <- readRDS(file.path(inter_dir, "counts_synced.rds"))
phenodata <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))
dge_6     <- readRDS(file.path(inter_dir, "dge_filtered_tmm.rds"))
logCPM_6  <- readRDS(file.path(inter_dir, "logCPM_6samples.rds"))

# 1. 6-Sample QLF Diagnostic
design_6 <- model.matrix(~factor(phenodata$condition))
dge_6    <- estimateDisp(dge_6, design_6)
fit_6    <- glmQLFit(dge_6, design_6)
qlf_6    <- glmQLFTest(fit_6, coef = 2)
res_6    <- topTags(qlf_6, n = Inf)$table
cat(">>> [09] 6-sample QLF Significant DEGs (FDR < 0.05):", sum(res_6$FDR < 0.05), "\n")

# 2. Expression Density Plot
melted_counts <- melt(logCPM_6)
colnames(melted_counts) <- c("Gene", "Sample", "Log2CPM")
melted_counts$Condition <- ifelse(grepl("WT", melted_counts$Sample), "WT", "Cond")

p_density <- ggplot(melted_counts, aes(x = Log2CPM, color = Sample, linetype = Condition)) +
    geom_density(linewidth = 0.9) +
    theme_minimal(base_size = 12) +
    labs(
        title = "Per-Sample Log2-CPM Density Distribution (N=6)",
        x     = "Log2 Counts Per Million",
        y     = "Density"
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(asset_dir, "002_pilot_expression_density.png"), plot = p_density, width = 8, height = 6, dpi = 300)

# 3. Cleaned 5-Sample Cohort Initialization (Excluding R1-1_WT)
samples_clean   <- setdiff(colnames(counts), "R1-1_WT")
counts_clean    <- counts[, samples_clean]
phenodata_clean <- phenodata[samples_clean, , drop = FALSE]

group_clean     <- factor(phenodata_clean$condition)
dge_clean       <- DGEList(counts = counts_clean, group = group_clean)
dge_clean       <- dge_clean[filterByExpr(dge_clean), , keep.lib.sizes = FALSE]
dge_clean       <- calcNormFactors(dge_clean)

# 4. Cleaned 5-Sample QLF Fit & P-value Histograms
design_clean <- model.matrix(~group_clean)
dge_clean    <- estimateDisp(dge_clean, design_clean)
fit_clean    <- glmQLFit(dge_clean, design_clean)
qlf_clean    <- glmQLFTest(fit_clean, coef = 2)
res_clean    <- topTags(qlf_clean, n = Inf)$table
cat(">>> [09] 5-sample QLF Significant DEGs (FDR < 0.05):", sum(res_clean$FDR < 0.05), "\n")

png(file.path(asset_dir, "003_pilot_pvalue_histogram_6vs5.png"), width = 1200, height = 600, res = 150)
par(mfrow = c(1, 2))
hist(res_6$PValue, breaks = 50, col = "salmon", main = "Original (6 Samples)", xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
hist(res_clean$PValue, breaks = 50, col = "skyblue", main = "Cleaned (5 Samples)", xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
dev.off()

# 5. Final Production Likelihood Ratio Test (LRT)
fit_lrt  <- glmFit(dge_clean, design_clean)
lrt_test <- glmLRT(fit_lrt, coef = 2)
res_lrt  <- topTags(lrt_test, n = Inf)$table

final_sig_genes <- res_lrt %>% filter(FDR < 0.05 & abs(logFC) > 1.0)
cat(">>> [09] Final Significant DEGs via LRT (FDR < 0.05, |logFC| > 1):", nrow(final_sig_genes), "\n")

write.csv(res_lrt, file.path(diff_dir, "pilot_lrt_all_genes.csv"))
write.csv(final_sig_genes, file.path(diff_dir, "pilot_lrt_significant_degs.csv"))
saveRDS(res_lrt, file.path(inter_dir, "res_lrt.rds"))
saveRDS(final_sig_genes, file.path(inter_dir, "final_sig_genes.rds"))