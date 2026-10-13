#!/usr/bin/env Rscript
# ==============================================================================
# Module 06: gene count matrix and sample table (log 2026-04-06)
#
# Input : results/matrices/gene_count_matrix.csv, results/matrices/phenodata.csv (module 04)
# Output: results/intermediate/counts_synced.rds, results/intermediate/phenodata_synced.rds
# ==============================================================================

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir     <- get_base_dir()
matrices_dir <- file.path(base_dir, "results", "matrices")
inter_dir    <- file.path(base_dir, "results", "intermediate")
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

# check.names = FALSE keeps the hyphens in the sample names
counts    <- read.csv(file.path(matrices_dir, "gene_count_matrix.csv"), row.names = 1, check.names = FALSE)
phenodata <- read.csv(file.path(matrices_dir, "phenodata.csv"), row.names = 1, check.names = FALSE)

# prepDE.py3 writes the libraries in alphabetical order; edgeR matches samples by position,
# so the columns are put in the order of phenodata.csv
counts <- counts[, rownames(phenodata)]
stopifnot(all(colnames(counts) == rownames(phenodata)))

cat("[06]", nrow(counts), "genes x", ncol(counts), "libraries; columns aligned with phenodata.csv\n")
saveRDS(counts, file.path(inter_dir, "counts_synced.rds"))
saveRDS(phenodata, file.path(inter_dir, "phenodata_synced.rds"))
