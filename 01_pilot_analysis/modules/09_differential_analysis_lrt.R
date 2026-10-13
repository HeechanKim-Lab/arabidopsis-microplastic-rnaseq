#!/usr/bin/env Rscript
# ==============================================================================
# Module 09: QLF and LRT models (log 2026-04-08_02)
#
# Input : results/intermediate/counts_synced.rds, phenodata_synced.rds (module 06)
#         results/intermediate/dge_filtered_tmm.rds (module 07)
# Output: assets/002_pilot_expression_density.png, assets/003_pilot_pvalue_histogram_6vs5.png
#         results/differential/pilot_lrt_all_genes.csv, pilot_lrt_significant_genes.csv
#         results/intermediate/res_lrt.rds, final_sig_genes.rds
#
# The gene list used downstream comes from the LRT on all six libraries (section E).
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(reshape2)
    library(ggplot2)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
inter_dir <- file.path(base_dir, "results", "intermediate")
diff_dir  <- file.path(base_dir, "results", "differential")
asset_dir <- file.path(base_dir, "assets")
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

counts    <- readRDS(file.path(inter_dir, "counts_synced.rds"))
phenodata <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))
dge       <- readRDS(file.path(inter_dir, "dge_filtered_tmm.rds"))

# Fit one QLF model on a subset of libraries (filtering and TMM are redone for the subset)
fit_qlf_subset <- function(samples) {
    group_sub <- factor(phenodata[samples, "condition"])
    dge_sub   <- DGEList(counts = counts[, samples], group = group_sub)
    dge_sub   <- dge_sub[filterByExpr(dge_sub), , keep.lib.sizes = FALSE]
    dge_sub   <- calcNormFactors(dge_sub)
    design_sub <- model.matrix(~group_sub)
    dge_sub   <- estimateDisp(dge_sub, design_sub)
    fit_sub   <- glmQLFit(dge_sub, design_sub)
    glmQLFTest(fit_sub, coef = 2)
}

# ---- A. Six libraries, QLF ---------------------------------------------------
# The reference level of group is "Cond" (alphabetical order), so coefficient 2 estimates
# WT - Cond: a positive log2FC means higher expression in WT. The column name below is the
# one used in the original run.
group  <- factor(phenodata$condition)
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "Conditioned_vs_WT")

dge   <- estimateDisp(dge, design)
fit   <- glmQLFit(dge, design)
qlf   <- glmQLFTest(fit, coef = 2)
res_6 <- topTags(qlf, n = Inf)$table
cat("\n[09] A. Six libraries, QLF (FDR < 0.05)\n")
print(summary(decideTests(qlf)))

# ---- B. Six libraries, robust QLF ---------------------------------------------
dge_robust <- estimateDisp(dge, design, robust = TRUE)
fit_robust <- glmQLFit(dge_robust, design, robust = TRUE)
qlf_robust <- glmQLFTest(fit_robust, coef = 2)
cat("\n[09] B. Six libraries, robust QLF\n")
print(summary(decideTests(qlf_robust)))

# ---- C. Five libraries (R1-1_WT removed), QLF ---------------------------------
qlf_clean <- fit_qlf_subset(setdiff(colnames(counts), "R1-1_WT"))
res_clean <- topTags(qlf_clean, n = Inf)$table
cat("\n[09] C. Five libraries (R1-1_WT removed), QLF; genes tested:", nrow(res_clean), "\n")
print(summary(decideTests(qlf_clean)))
cat("Genes with P < 0.01: five libraries", sum(res_clean$PValue < 0.01),
    "; six libraries", sum(res_6$PValue < 0.01), "\n")

# ---- D. P-value histograms (figure 003) and expression distributions (figure 002)
png(file.path(asset_dir, "003_pilot_pvalue_histogram_6vs5.png"), width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2))
hist(res_6$PValue, breaks = 50, col = "salmon", main = "Original (6 Samples)",
     xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
hist(res_clean$PValue, breaks = 50, col = "skyblue", main = "Cleaned (5 Samples)",
     xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
invisible(dev.off())

# Log2-CPM of all genes from the raw counts (library sizes only, no filtering or TMM)
logCPM_all <- cpm(as.matrix(counts), log = TRUE, prior.count = 2)
melted_counts <- melt(logCPM_all)
colnames(melted_counts) <- c("Gene", "Sample", "Log2CPM")
melted_counts$Condition <- ifelse(grepl("WT", melted_counts$Sample), "WT", "Cond")
p_density <- ggplot(melted_counts, aes(x = Log2CPM, color = Sample, linetype = Condition)) +
    geom_density(linewidth = 1) +
    theme_minimal() +
    ggtitle("Per-Sample Expression Distribution (All 6 Samples)") +
    xlab("Log2 Counts Per Million") +
    ylab("Density")
ggsave(file.path(asset_dir, "002_pilot_expression_density.png"), p_density,
       width = 8, height = 5, dpi = 200, bg = "white")

# ---- E. Six libraries, LRT (gene list used for figures 004-013) -------------
fit_lrt  <- glmFit(dge, design)
lrt_test <- glmLRT(fit_lrt, coef = 2)
cat("\n[09] E. Six libraries, LRT (FDR < 0.05)\n")
print(summary(decideTests(lrt_test)))

res_lrt <- topTags(lrt_test, n = Inf)$table
final_sig_genes <- res_lrt %>% filter(FDR < 0.05 & abs(logFC) > 1)
cat("LRT genes with FDR < 0.05 and |log2FC| > 1:", nrow(final_sig_genes), "\n")

# ---- F. Four libraries (R1-1_WT and R3-9_Cond removed), QLF --------------------
qlf_final <- fit_qlf_subset(c("R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond"))
cat("\n[09] F. Four libraries, QLF; genes tested:", nrow(qlf_final$table), "\n")
print(summary(decideTests(qlf_final)))

# ---- Outputs -----------------------------------------------------------------
write.csv(res_lrt, file.path(diff_dir, "pilot_lrt_all_genes.csv"))
write.csv(final_sig_genes, file.path(diff_dir, "pilot_lrt_significant_genes.csv"))
saveRDS(res_lrt, file.path(inter_dir, "res_lrt.rds"))
saveRDS(final_sig_genes, file.path(inter_dir, "final_sig_genes.rds"))
