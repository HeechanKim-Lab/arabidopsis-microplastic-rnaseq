#!/usr/bin/env Rscript
# ==============================================================================
# Pilot downstream analysis: runs modules 05-11 in order, each in its own R process.
#
# Input : results/ballgown/ and results/matrices/ from pipeline/run_upstream_tuxedo.sh
# Output: figures 001-013 in assets/, tables in results/differential/ and results/enrichment/,
#         intermediate objects in results/intermediate/
# ==============================================================================

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
base_dir <- if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
rscript  <- file.path(R.home("bin"), "Rscript")

modules <- c(
    "05_fpkm_to_tpm.R",
    "06_data_ingestion.R",
    "07_filtering_and_tmm.R",
    "08_qc_pca.R",
    "09_differential_analysis_lrt.R",
    "10_visualization.R",
    "11_functional_enrichment.R"
)

for (m in modules) {
    cat("\n=== ", m, " ===\n", sep = "")
    status <- system2(rscript, shQuote(file.path(base_dir, "modules", m)))
    if (status != 0) stop("Module failed: ", m)
}
cat("\n=== Pilot downstream analysis finished ===\n")
