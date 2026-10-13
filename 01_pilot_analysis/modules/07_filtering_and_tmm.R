#!/usr/bin/env Rscript
# ==============================================================================
# Module 07: edgeR filtering and TMM normalisation (log 2026-04-07)
#
# Input : results/intermediate/counts_synced.rds, phenodata_synced.rds (module 06)
# Output: results/intermediate/dge_filtered_tmm.rds
# ==============================================================================

suppressPackageStartupMessages(library(edgeR))

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
inter_dir <- file.path(base_dir, "results", "intermediate")

counts    <- readRDS(file.path(inter_dir, "counts_synced.rds"))
phenodata <- readRDS(file.path(inter_dir, "phenodata_synced.rds"))

group <- factor(phenodata$condition)
dge   <- DGEList(counts = counts, group = group)

# Keep genes with enough reads in enough libraries (uses the group factor of the DGEList)
keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalisation factors
dge <- calcNormFactors(dge)

cat("[07] Genes kept:", nrow(dge), "of", nrow(counts), "\n")
print(dge$samples)
saveRDS(dge, file.path(inter_dir, "dge_filtered_tmm.rds"))
