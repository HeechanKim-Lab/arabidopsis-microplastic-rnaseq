#!/usr/bin/env Rscript
# ==============================================================================
# Module 07: annotated table of all significant genes and their gene sets
#
# Input : results/intermediate/lrt_root00_list.rds, lrt_root01_list.rds (module 02)
#         data/Arabidopsis_genes_description.xlsx (lab table: Gene_locus, Symbol, Description)
# Output: results/tables/Final_Core_Plastic_Response_Annotated.csv
#
# Gene sets (FDR < 0.05): Root00 5% vs 0% at day 120 and interaction; Root01 5% vs 0% at
# day 120, interaction, and day 120 vs day 0 at 5%.
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(readxl)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir   <- get_base_dir()
data_dir   <- file.path(base_dir, "data")
inter_dir  <- file.path(base_dir, "results", "intermediate")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

find_input <- function(name) {
    for (p in c(file.path(data_dir, name), file.path(getwd(), name))) if (file.exists(p)) return(p)
    stop("Input not found: ", name, " (expected in ", data_dir, ")")
}
anno_file  <- read_excel(find_input("Arabidopsis_genes_description.xlsx"))
lrt_root00 <- readRDS(file.path(inter_dir, "lrt_root00_list.rds"))
lrt_root01 <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

sig_genes <- function(res) {
    tab <- topTags(res, n = Inf)$table
    rownames(tab[tab$FDR < 0.05, ])
}
list_00_p120  <- sig_genes(lrt_root00$res00_120d)
list_00_inter <- sig_genes(lrt_root00$res00_inter)
list_01_p120  <- sig_genes(lrt_root01$res01_p120)
list_01_inter <- sig_genes(lrt_root01$res01_inter)
list_aging_p  <- sig_genes(lrt_root01$res01_aging_p)

all_sig_genes  <- unique(c(list_00_p120, list_00_inter, list_01_p120, list_01_inter, list_aging_p))
master_results <- data.frame(Gene_locus = all_sig_genes)

# Membership flags ("Yes" or empty)
master_results$root00_Plastic_D120 <- ifelse(master_results$Gene_locus %in% list_00_p120, "Yes", "")
master_results$root00_Incubation   <- ifelse(master_results$Gene_locus %in% list_00_inter, "Yes", "")
master_results$root00_Intersection <- ifelse(master_results$root00_Plastic_D120 == "Yes" &
                                             master_results$root00_Incubation == "Yes", "Yes", "")
master_results$root01_Plastic_D120 <- ifelse(master_results$Gene_locus %in% list_01_p120, "Yes", "")
master_results$root01_Incubation   <- ifelse(master_results$Gene_locus %in% list_01_inter, "Yes", "")
master_results$root01_Aging_5p     <- ifelse(master_results$Gene_locus %in% list_aging_p, "Yes", "")
master_results$root01_Triple_Intersection <- ifelse(
    master_results$root01_Plastic_D120 == "Yes" &
    master_results$root01_Incubation == "Yes" &
    master_results$root01_Aging_5p == "Yes", "Yes", "")

# Add symbols and descriptions from the lab annotation table
final_report <- left_join(master_results, anno_file, by = "Gene_locus") %>%
    dplyr::select(Gene_locus, Symbol, Description, everything())

out_csv <- file.path(tables_dir, "Final_Core_Plastic_Response_Annotated.csv")
write.csv(final_report, out_csv, row.names = FALSE)
cat("[07]", nrow(final_report), "genes ->", out_csv, "\n")
