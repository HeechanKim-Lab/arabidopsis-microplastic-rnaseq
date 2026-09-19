#!/usr/bin/env Rscript
# ==============================================================================
# Script: run_downstream_pilot.R
# Project: arabidopsis-microplastic-rnaseq (01_pilot_analysis)
# Description: End-to-end downstream statistical modeling and enrichment pipeline
# ==============================================================================

# --- 1. Dependencies and Environment Initialization ---
suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(pheatmap)
    library(clusterProfiler)
    library(org.At.tair.db)
})

# Setup directory anchors relative to script execution
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
    base_dir <- normalizePath(file.path(dirname(script_path), ".."))
} else {
    base_dir <- getwd()
}

results_dir  <- file.path(base_dir, "results")
matrices_dir <- file.path(results_dir, "matrices")
ballgown_dir <- file.path(results_dir, "ballgown")
assets_dir   <- file.path(base_dir, "assets")
diff_dir     <- file.path(results_dir, "differential")

dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== [PHASE 1] FPKM to TPM Linear Transformation ===\n")
samples <- c("R1-1_WT", "R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond", "R3-9_Cond")

get_tpm_safe <- function(s) {
    file_path <- file.path(ballgown_dir, s, "t_data.ctab")
    if (!file.exists(file_path)) {
        stop(paste("Error: File not found:", file_path))
    }
    df <- read.table(file_path, header = TRUE, sep = "\t")
    tpm_vals <- (df$FPKM / sum(df$FPKM)) * 1e6
    res <- data.frame(t_name = df$t_name, gene_id = df$gene_id, TPM = tpm_vals)
    colnames(res)[3] <- s
    return(res)
}

tpm_list <- lapply(samples, get_tpm_safe)
final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by = c("t_name", "gene_id")), tpm_list)
tpm_out_file <- file.path(matrices_dir, "arabidopsis_transcript_tpm_matrix.csv")
write.csv(final_tpm_matrix, tpm_out_file, row.names = FALSE)
cat("TPM matrix compiled and written to:", tpm_out_file, "\n")


cat("\n=== [PHASE 2] Data Ingestion & Metadata Synchronization ===\n")
counts_file <- file.path(matrices_dir, "gene_count_matrix.csv")
pheno_file  <- file.path(matrices_dir, "phenodata.csv")

counts    <- read.csv(counts_file, row.names = 1, check.names = FALSE)
phenodata <- read.csv(pheno_file, row.names = 1, check.names = FALSE)

# Harmonize column order to metadata row sequence
counts <- counts[, rownames(phenodata)]
stopifnot(all(colnames(counts) == rownames(phenodata)))
cat("Counts and metadata synchronized across", ncol(counts), "libraries.\n")


cat("\n=== [PHASE 3] Pre-filtering & TMM Normalization ===\n")
group <- factor(phenodata$condition)
dge   <- DGEList(counts = counts, group = group)

keep  <- filterByExpr(dge)
dge   <- dge[keep, , keep.lib.sizes = FALSE]
dge   <- calcNormFactors(dge)
cat("Filtered feature space retained:", nrow(dge), "genes out of", nrow(counts), "\n")


cat("\n=== [PHASE 4] Quality Control: Unsupervised PCA ===\n")
logCPM  <- cpm(dge, log = TRUE, prior.count = 2)
pca_res <- prcomp(t(logCPM), scale. = TRUE)
pct_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- data.frame(
    Sample    = rownames(pca_res$x),
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    Condition = phenodata$condition
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text(vjust = 1.8, size = 3.8) +
    scale_color_manual(values = c("WT" = "#377EB8", "Cond" = "#E41A1C")) +
    theme_minimal(base_size = 12) +
    labs(
        title = "PCA Plot: Arabidopsis Root Pilot (N=6)",
        x     = paste0("PC1: ", pct_var[1], "% variance"),
        y     = paste0("PC2: ", pct_var[2], "% variance")
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(assets_dir, "001_pilot_qc_pca.png"), plot = p_pca, width = 8, height = 6, dpi = 300)


cat("\n=== [PHASE 5] Outlier Diagnostic Diagnostics (QLF vs LRT) ===\n")
# 1. Baseline 6-sample Quasi-Likelihood Model Fit
design_6 <- model.matrix(~group)
dge_6    <- estimateDisp(dge, design_6)
fit_6    <- glmQLFit(dge_6, design_6)
qlf_6    <- glmQLFTest(fit_6, coef = 2)
res_6    <- topTags(qlf_6, n = Inf)$table
cat("Original 6-sample QLF Significant DEGs (FDR < 0.05):", sum(res_6$FDR < 0.05), "\n")

# 2. Density Plot Export
melted_counts <- melt(logCPM)
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

ggsave(file.path(assets_dir, "002_pilot_expression_density.png"), plot = p_density, width = 8, height = 6, dpi = 300)

# 3. Cleaned 5-sample QLF Model Fit (Excluding R1-1_WT)
samples_clean   <- setdiff(colnames(counts), "R1-1_WT")
counts_clean    <- counts[, samples_clean]
phenodata_clean <- phenodata[samples_clean, , drop = FALSE]

group_clean     <- factor(phenodata_clean$condition)
dge_clean       <- DGEList(counts = counts_clean, group = group_clean)
dge_clean       <- dge_clean[filterByExpr(dge_clean), , keep.lib.sizes = FALSE]
dge_clean       <- calcNormFactors(dge_clean)

design_clean    <- model.matrix(~group_clean)
dge_clean       <- estimateDisp(dge_clean, design_clean)
fit_clean       <- glmQLFit(dge_clean, design_clean)
qlf_clean       <- glmQLFTest(fit_clean, coef = 2)
res_clean       <- topTags(qlf_clean, n = Inf)$table
cat("Cleaned 5-sample QLF Significant DEGs (FDR < 0.05):", sum(res_clean$FDR < 0.05), "\n")

# 4. P-Value Comparison Histogram Export
png(file.path(assets_dir, "003_pilot_pvalue_histogram_6vs5.png"), width = 1200, height = 600, res = 150)
par(mfrow = c(1, 2))
hist(res_6$PValue, breaks = 50, col = "salmon", main = "Original (6 Samples)", xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
hist(res_clean$PValue, breaks = 50, col = "skyblue", main = "Cleaned (5 Samples)", xlab = "P-value", ylim = c(0, 1600))
abline(v = 0.05, col = "red", lty = 2)
dev.off()


cat("\n=== [PHASE 6] Final LRT Production Modeling ===\n")
# Fit GLM under Likelihood Ratio Test formulation
fit_lrt  <- glmFit(dge_clean, design_clean)
lrt_test <- glmLRT(fit_lrt, coef = 2)
res_lrt  <- topTags(lrt_test, n = Inf)$table

final_sig_genes <- res_lrt %>% filter(FDR < 0.05 & abs(logFC) > 1.0)
cat("Final Significant DEGs via LRT (FDR < 0.05, |logFC| > 1):", nrow(final_sig_genes), "\n")

write.csv(res_lrt, file.path(diff_dir, "pilot_lrt_all_genes.csv"))
write.csv(final_sig_genes, file.path(diff_dir, "pilot_lrt_significant_degs.csv"))


cat("\n=== [PHASE 7] Publication Visualizations (Volcano & Heatmap) ===\n")
# 1. Volcano Plot
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
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(assets_dir, "004_pilot_lrt_volcano.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)

# 2. Z-score Heatmap (Top 50 DEGs across all 6 original libraries)
top50_ids <- rownames(final_sig_genes)[1:min(50, nrow(final_sig_genes))]
heatmap_matrix <- logCPM[top50_ids, , drop = FALSE]
heatmap_scaled <- t(scale(t(heatmap_matrix)))

png(file.path(assets_dir, "005_pilot_top50_heatmap.png"), width = 1200, height = 1500, res = 150)
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


cat("\n=== [PHASE 8] Identifier Mapping & Functional Enrichment (GO & KEGG) ===\n")
# 1. Rosetta Stone ID Dictionary from Ballgown ctab
ctab <- read.table(file.path(ballgown_dir, "R1-1_WT", "t_data.ctab"), header = TRUE, sep = "\t")
id_map <- ctab %>%
    dplyr::select(gene_id, t_name) %>%
    mutate(tair_id = gsub("\\..*$", "", t_name)) %>%
    filter(grepl("^AT", tair_id)) %>%
    distinct(gene_id, tair_id)

final_sig_genes$gene_id <- rownames(final_sig_genes)
final_mapped <- merge(final_sig_genes, id_map, by = "gene_id")
cat("Mapped", nrow(final_mapped), "transcripts to TAIR locus IDs.\n")

up_ids   <- unique(final_mapped$tair_id[final_mapped$logFC > 1.0])
down_ids <- unique(final_mapped$tair_id[final_mapped$logFC < -1.0])

# 2. Gene Ontology Enrichment (BP, MF, CC)
run_and_save_go <- function(genes, ont_type, direction, filename) {
    res <- enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont_type, pvalueCutoff = 0.05)
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        p <- dotplot(res, showCategory = 15) + 
            ggtitle(paste0(direction, ": Biological Ontology (", ont_type, ")"))
        ggsave(file.path(assets_dir, filename), plot = p, width = 8, height = 6, dpi = 300)
    }
}

run_and_save_go(up_ids,   "BP", "UP",   "006_pilot_go_bp_up.png")
run_and_save_go(down_ids, "BP", "DOWN", "007_pilot_go_bp_down.png")
run_and_save_go(up_ids,   "MF", "UP",   "008_pilot_go_mf_up.png")
run_and_save_go(down_ids, "MF", "DOWN", "009_pilot_go_mf_down.png")
run_and_save_go(up_ids,   "CC", "UP",   "010_pilot_go_cc_up.png")
run_and_save_go(down_ids, "CC", "DOWN", "011_pilot_go_cc_down.png")

# 3. KEGG Pathway Deconvolution
ids_mapped <- bitr(final_mapped$tair_id, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)

up_entrez   <- unique(ids_mapped$ENTREZID[ids_mapped$TAIR %in% up_ids])
down_entrez <- unique(ids_mapped$ENTREZID[ids_mapped$TAIR %in% down_ids])

kegg_up   <- enrichKEGG(gene = up_entrez,   organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_entrez, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)

if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
    p_kegg_up <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG: Up-regulated (Activated Circuits)")
    ggsave(file.path(assets_dir, "012_pilot_kegg_pathways_up.png"), plot = p_kegg_up, width = 8, height = 6, dpi = 300)
}

if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
    p_kegg_down <- dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG: Down-regulated (Suppressed Circuits)")
    ggsave(file.path(assets_dir, "013_pilot_kegg_pathways_down.png"), plot = p_kegg_down, width = 8, height = 6, dpi = 300)
}

cat("=== Downstream R Pipeline Completed Successfully ===\n")