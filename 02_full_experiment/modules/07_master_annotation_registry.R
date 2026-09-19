#!/usr/bin/env Rscript
# ==============================================================================
# Module: 07_master_annotation_registry.R
# Description: Master DEG membership registry and TAIR annotation catalog generation
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(readxl)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
data_dir   <- file.path(base_dir, "data")
inter_dir  <- file.path(base_dir, "results", "intermediate")
tables_dir <- file.path(base_dir, "results", "tables")

lrt_root00 <- readRDS(file.path(inter_dir, "lrt_root00_list.rds"))
lrt_root01 <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

anno_file <- file.path(data_dir, "Arabidopsis_genes_description.xlsx")
if (!file.exists(anno_file)) anno_file <- "Arabidopsis_genes_description.xlsx"
anno_db <- read_excel(anno_file)

list_00_p120    <- rownames(topTags(lrt_root00$res00_120d, n = Inf)$table[topTags(lrt_root00$res00_120d, n = Inf)$table$FDR < 0.05, ])
list_00_inter   <- rownames(topTags(lrt_root00$res00_inter, n = Inf)$table[topTags(lrt_root00$res00_inter, n = Inf)$table$FDR < 0.05, ])
list_01_p120    <- rownames(topTags(lrt_root01$res01_p120, n = Inf)$table[topTags(lrt_root01$res01_p120, n = Inf)$table$FDR < 0.05, ])
list_01_inter   <- rownames(topTags(lrt_root01$res01_inter, n = Inf)$table[topTags(lrt_root01$res01_inter, n = Inf)$table$FDR < 0.05, ])
list_01_aging_p <- rownames(topTags(lrt_root01$res01_aging_p, n = Inf)$table[topTags(lrt_root01$res01_aging_p, n = Inf)$table$FDR < 0.05, ])

all_sig_genes <- unique(c(list_00_p120, list_00_inter, list_01_p120, list_01_inter, list_01_aging_p))
cat(">>> [07] Compiling membership flags across", length(all_sig_genes), "unique significant loci...\n")

master_report <- data.frame(Gene_locus = all_sig_genes)

# Root00 Memberships
master_report$root00_Plastic_D120 <- ifelse(master_report$Gene_locus %in% list_00_p120, "Yes", "")
master_report$root00_Incubation   <- ifelse(master_report$Gene_locus %in% list_00_inter, "Yes", "")
master_report$root00_Intersection <- ifelse(master_report$root00_Plastic_D120 == "Yes" & 
                                            master_report$root00_Incubation == "Yes", "Yes", "")

# Root01 Memberships
master_report$root01_Plastic_D120 <- ifelse(master_report$Gene_locus %in% list_01_p120, "Yes", "")
master_report$root01_Incubation   <- ifelse(master_report$Gene_locus %in% list_01_inter, "Yes", "")
master_report$root01_Aging_5p     <- ifelse(master_report$Gene_locus %in% list_01_aging_p, "Yes", "")
master_report$root01_Triple_Intersection <- ifelse(
    master_report$root01_Plastic_D120 == "Yes" & 
    master_report$root01_Incubation == "Yes" & 
    master_report$root01_Aging_5p == "Yes", "Yes", ""
)

# Join with functional annotation
master_final <- left_join(master_report, anno_db, by = "Gene_locus") %>%
    dplyr::select(Gene_locus, Symbol, Description, everything())

out_csv <- file.path(tables_dir, "Final_Core_Plastic_Response_Annotated.csv")
write.csv(master_final, out_csv, row.names = FALSE)
cat(">>> [07] Master annotation table written to:", out_csv, "\n")