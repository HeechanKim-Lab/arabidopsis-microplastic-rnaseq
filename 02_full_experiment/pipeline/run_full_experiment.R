#!/usr/bin/env Rscript
# ==============================================================================
# Full experiment: runs modules 01-07 in order, each in its own R process.
#
# Input (not in the repository; put them in data/):
#   data/Arabidopsis_plastics_2023.RData, data/Arabidopsis_genes_description.xlsx
# Output: figures 020-040 in assets/ (except 032, drawn outside R), tables in results/tables/
#         and results/enrichment/, intermediate objects in results/intermediate/
# ==============================================================================

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
base_dir <- if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
rscript  <- file.path(R.home("bin"), "Rscript")

modules <- c(
    "01_full_matrix_and_subset_pca.R",
    "02_lrt_factorial_contrasts.R",
    "03_interaction_visualization.R",
    "04_functional_enrichment.R",
    "05_multivariate_permanova.R",
    "06_targeted_pathways.R",
    "07_master_annotation_registry.R"
)

for (m in modules) {
    cat("\n=== ", m, " ===\n", sep = "")
    status <- system2(rscript, shQuote(file.path(base_dir, "modules", m)))
    if (status != 0) stop("Module failed: ", m)
}
cat("\n=== Full experiment finished ===\n")
