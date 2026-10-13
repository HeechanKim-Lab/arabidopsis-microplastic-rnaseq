#!/usr/bin/env Rscript
# ==============================================================================
# Module 05: transcript TPM table from StringTie FPKM (log 2026-03-31)
#
# Input : results/ballgown/<sample>/t_data.ctab (module 04)
# Output: results/matrices/arabidopsis_transcript_tpm_matrix.csv
# ==============================================================================

# Stage folder (01_pilot_analysis/): from the script path, or the working directory when sourced
get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir     <- get_base_dir()
ballgown_dir <- file.path(base_dir, "results", "ballgown")
matrices_dir <- file.path(base_dir, "results", "matrices")
dir.create(matrices_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("R1-1_WT", "R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond", "R3-9_Cond")

get_tpm_safe <- function(s) {
    file_path <- file.path(ballgown_dir, s, "t_data.ctab")
    df <- read.table(file_path, header = TRUE, sep = "\t")
    # TPM: FPKM rescaled to sum to one million within the library
    tpm_vals <- (df$FPKM / sum(df$FPKM)) * 1e6
    # t_name is unique per transcript; gene_id is kept to show the parent gene
    res <- data.frame(t_name = df$t_name, gene_id = df$gene_id, TPM = tpm_vals)
    colnames(res)[3] <- s
    res
}

tpm_list <- lapply(samples, get_tpm_safe)
# Join on both keys; joining on gene_id alone would duplicate rows for genes with several transcripts
final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by = c("t_name", "gene_id")), tpm_list)

out_file <- file.path(matrices_dir, "arabidopsis_transcript_tpm_matrix.csv")
write.csv(final_tpm_matrix, out_file, row.names = FALSE)
cat("[05] TPM table:", nrow(final_tpm_matrix), "transcripts ->", out_file, "\n")
