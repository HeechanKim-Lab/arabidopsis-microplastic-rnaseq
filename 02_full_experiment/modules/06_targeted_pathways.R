#!/usr/bin/env Rscript
# ==============================================================================
# Module: 06_targeted_pathways.R
# Description: Targeted regex mining and heatmap generation for ROS and Plastid loci
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(pheatmap)
    library(readxl)
    library(dplyr)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
data_dir   <- file.path(base_dir, "data")
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")

TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
lrt_root01     <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

rdata_file <- file.path(data_dir, "Arabidopsis_plastics_2023.RData")
if (!file.exists(rdata_file)) rdata_file <- "Arabidopsis_plastics_2023.RData"
load(rdata_file)

anno_file <- file.path(data_dir, "Arabidopsis_genes_description.xlsx")
if (!file.exists(anno_file)) anno_file <- "Arabidopsis_genes_description.xlsx"
anno_db <- read_excel(anno_file)

# Extract non-redundant union of significant LRT loci across Root01
all_lrt_sig <- unique(c(
    rownames(topTags(lrt_root01$res01_p120, n = Inf)$table[topTags(lrt_root01$res01_p120, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(lrt_root01$res01_inter, n = Inf)$table[topTags(lrt_root01$res01_inter, n = Inf)$table$FDR < 0.05, ]),
    rownames(topTags(lrt_root01$res01_aging_p, n = Inf)$table[topTags(lrt_root01$res01_aging_p, n = Inf)$table$FDR < 0.05, ])
))

anno_sig <- anno_db %>% filter(Gene_locus %in% all_lrt_sig)
ann_col  <- data.frame(Group = root01_group_02)
rownames(ann_col) <- colnames(TMM_root01_log)

cat(">>> [06] Mining Oxidative Stress / ROS genes...\n")
ox_keywords <- "hydrogen peroxide|oxidative|reactive oxygen|peroxidase"
ox_genes    <- anno_sig %>% filter(grepl(ox_keywords, Description, ignore.case = TRUE))

if (nrow(ox_genes) > 1) {
    ox_matrix <- TMM_root01_log[ox_genes$Gene_locus, ]
    rownames(ox_matrix) <- ifelse(is.na(ox_genes$Symbol) | ox_genes$Symbol == "", 
                                  ox_genes$Gene_locus, ox_genes$Symbol)
    
    png(file.path(assets_dir, "039_targeted_ros_heatmap.png"), width = 1200, height = 900, res = 150)
    pheatmap(ox_matrix, 
             scale = "row", 
             annotation_col = ann_col,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Oxidative Stress Responsiveness (", nrow(ox_genes), " Genes)"),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize_row = 9)
    dev.off()
}

cat(">>> [06] Mining Photosynthesis / Plastid shift genes...\n")
photo_keywords <- "photosynthe|chloroplast|photosystem|light-harvesting|thylakoid"
photo_genes    <- anno_sig %>% filter(grepl(photo_keywords, Description, ignore.case = TRUE))

if (nrow(photo_genes) > 1) {
    photo_matrix <- TMM_root01_log[photo_genes$Gene_locus, ]
    rownames(photo_matrix) <- ifelse(is.na(photo_genes$Symbol) | photo_genes$Symbol == "", 
                                     photo_genes$Gene_locus, photo_genes$Symbol)
    
    png(file.path(assets_dir, "040_targeted_photosynthesis_heatmap.png"), width = 1200, height = 900, res = 150)
    pheatmap(photo_matrix, 
             scale = "row", 
             annotation_col = ann_col,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Photosynthesis & Plastid Shifts (", nrow(photo_genes), " Genes)"),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize_row = 9)
    dev.off()
}

cat(">>> [06] Targeted pathway heatmaps generated successfully.\n")