#!/usr/bin/env Rscript
# ==============================================================================
# Module: 05_fpkm_to_tpm.R
# Description: Linear transformation of Ballgown FPKM estimates to TPM matrix
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

base_dir     <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
ballgown_dir <- file.path(base_dir, "results", "ballgown")
matrices_dir <- file.path(base_dir, "results", "matrices")
samples      <- c("R1-1_WT", "R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond", "R3-9_Cond")

get_tpm_safe <- function(s) {
    file_path <- file.path(ballgown_dir, s, "t_data.ctab")
    if (!file.exists(file_path)) stop(paste("File missing:", file_path))
    
    df <- read.table(file_path, header = TRUE, sep = "\t")
    tpm_vals <- (df$FPKM / sum(df$FPKM)) * 1e6
    res <- data.frame(t_name = df$t_name, gene_id = df$gene_id, TPM = tpm_vals)
    colnames(res)[3] <- s
    return(res)
}

cat(">>> [05] Extracting and transforming FPKM to TPM across 6 samples...\n")
tpm_list <- lapply(samples, get_tpm_safe)
final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by = c("t_name", "gene_id")), tpm_list)

out_file <- file.path(matrices_dir, "arabidopsis_transcript_tpm_matrix.csv")
write.csv(final_tpm_matrix, out_file, row.names = FALSE)
cat(">>> [05] Saved TPM matrix to:", out_file, "\n")