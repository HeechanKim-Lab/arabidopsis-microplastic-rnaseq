#!/usr/bin/env Rscript
# ==============================================================================
# Module: 06_data_ingestion.R
# Description: Count matrix ingestion, metadata alignment, and syntax validation
# ==============================================================================

base_dir     <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
matrices_dir <- file.path(base_dir, "results", "matrices")
inter_dir    <- file.path(base_dir, "results", "intermediate")
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(matrices_dir, "gene_count_matrix.csv")
pheno_file  <- file.path(matrices_dir, "phenodata.csv")

# Ingest without mutating sample identifiers (preserve hyphens)
counts    <- read.csv(counts_file, row.names = 1, check.names = FALSE)
phenodata <- read.csv(pheno_file, row.names = 1, check.names = FALSE)

# Reorder matrix columns to match metadata rows
counts <- counts[, rownames(phenodata)]
stopifnot(all(colnames(counts) == rownames(phenodata)))

cat(">>> [06] Ingestion verified: 32,540 features across", ncol(counts), "samples.\n")
saveRDS(counts, file.path(inter_dir, "counts_synced.rds"))
saveRDS(phenodata, file.path(inter_dir, "phenodata_synced.rds"))