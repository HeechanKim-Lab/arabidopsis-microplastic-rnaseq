#!/usr/bin/env Rscript
# ==============================================================================
# Module: 07_filtering_and_tmm.R
# Description: DGEList initialization, empirical pre-filtering, and TMM scaling
# ==============================================================================

suppressPackageStartupMessages(library(edgeR))

base_dir  <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir <- file.path(base_dir, "results", "intermediate")

counts    <- readRDS(file.path(inter_dir, "counts_synced.rds"))
phenodata <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))

group <- factor(phenodata$condition)
dge   <- DGEList(counts = counts, group = group)

# Empirical pre-filtering
keep  <- filterByExpr(dge)
dge   <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization
dge   <- calcNormFactors(dge)

cat(">>> [07] Filtering retained:", nrow(dge), "genes out of", nrow(counts), "\n")
cat(">>> [07] TMM factors calculated:\n")
print(dge$samples)

saveRDS(dge, file.path(inter_dir, "dge_filtered_tmm.rds"))