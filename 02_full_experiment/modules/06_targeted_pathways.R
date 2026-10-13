#!/usr/bin/env Rscript
# ==============================================================================
# Module 06: keyword-selected heatmaps, oxidative stress and photosynthesis/plastid
#            (log 2026-04-17)
#
# Input : results/intermediate/lrt_root01_list.rds (module 02), TMM_root01_log.rds (module 01)
#         data/Arabidopsis_genes_description.xlsx (lab table: Gene_locus, Symbol, Description)
#         data/Arabidopsis_plastics_2023.RData (root01_group_02 for the column annotation)
# Output: assets/039_pathway_ros_heatmap.png, assets/040_pathway_photosynthesis_heatmap.png
#
# The genes are chosen by keyword from the significant genes; this is not an enrichment test.
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(pheatmap)
    library(readxl)
    library(dplyr)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
data_dir  <- file.path(base_dir, "data")
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

find_input <- function(name) {
    for (p in c(file.path(data_dir, name), file.path(getwd(), name))) if (file.exists(p)) return(p)
    stop("Input not found: ", name, " (expected in ", data_dir, ")")
}
load(find_input("Arabidopsis_plastics_2023.RData"))
anno_file      <- read_excel(find_input("Arabidopsis_genes_description.xlsx"))
TMM_root01_log <- readRDS(file.path(inter_dir, "TMM_root01_log.rds"))
lrt_root01     <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

# Gene pool: FDR < 0.05 in 5% vs 0% at day 120, the interaction, or day 120 vs day 0 at 5%
df_p120    <- as.data.frame(topTags(lrt_root01$res01_p120, n = Inf))
df_inter   <- as.data.frame(topTags(lrt_root01$res01_inter, n = Inf))
df_aging_p <- as.data.frame(topTags(lrt_root01$res01_aging_p, n = Inf))
all_lrt_sig <- unique(c(
    rownames(df_p120[df_p120$FDR < 0.05, ]),
    rownames(df_inter[df_inter$FDR < 0.05, ]),
    rownames(df_aging_p[df_aging_p$FDR < 0.05, ])))
anno_sig <- anno_file %>% filter(Gene_locus %in% all_lrt_sig)

ann_col <- data.frame(Group = root01_group_02)
rownames(ann_col) <- colnames(TMM_root01_log)

keyword_heatmap <- function(keywords, title, out_png) {
    genes <- anno_sig %>% filter(grepl(keywords, Description, ignore.case = TRUE))
    cat("[06]", title, ":", nrow(genes), "genes\n")
    if (nrow(genes) < 2) {
        message("Fewer than two genes; no heatmap for ", title)
        return(invisible(NULL))
    }
    mat <- TMM_root01_log[genes$Gene_locus, ]
    # Row labels: gene symbol when available, otherwise the TAIR locus
    rownames(mat) <- ifelse(is.na(genes$Symbol) | genes$Symbol == "", genes$Gene_locus, genes$Symbol)
    pheatmap(mat, scale = "row", annotation_col = ann_col,
             main = paste(title, "(", nrow(genes), " Genes)"),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize_row = 8,
             filename = out_png, width = 8, height = 6)
}

keyword_heatmap("hydrogen peroxide|oxidative|reactive oxygen|peroxidase",
                "Oxidative Stress Responsiveness",
                file.path(asset_dir, "039_pathway_ros_heatmap.png"))
keyword_heatmap("photosynthe|chloroplast|photosystem|light-harvesting|thylakoid",
                "Photosynthesis & Plastid Shifts",
                file.path(asset_dir, "040_pathway_photosynthesis_heatmap.png"))
