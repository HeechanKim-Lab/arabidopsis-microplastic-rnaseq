#!/usr/bin/env Rscript
# ==============================================================================
# Module 11: StringTie-to-TAIR ID mapping, GO and KEGG enrichment (log 2026-04-11)
#
# Input : results/intermediate/final_sig_genes.rds (module 09)
#         results/ballgown/R1-1_WT/t_data.ctab (module 04)
# Output: assets/006_pilot_go_bp_up.png ... assets/013_pilot_kegg_pathways_down.png
#         results/enrichment/*.csv
#
# Gene sets are split by the sign of log2FC, which is WT relative to Cond:
# "UP" = higher in WT, "DOWN" = higher in Cond.
# KEGG uses TAIR locus IDs, which are the KEGG gene IDs for Arabidopsis (keyType = "kegg").
# enrichKEGG downloads the KEGG data, so it needs an internet connection.
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.At.tair.db)
    library(dplyr)
    library(ggplot2)
})

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir     <- get_base_dir()
inter_dir    <- file.path(base_dir, "results", "intermediate")
ballgown_dir <- file.path(base_dir, "results", "ballgown")
enrich_dir   <- file.path(base_dir, "results", "enrichment")
asset_dir    <- file.path(base_dir, "assets")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

final_sig_genes <- readRDS(file.path(inter_dir, "final_sig_genes.rds"))

# 1. Map StringTie gene IDs (MSTRG.*) and TAIR loci to TAIR IDs using t_data.ctab
ctab <- read.table(file.path(ballgown_dir, "R1-1_WT", "t_data.ctab"), header = TRUE, sep = "\t")
id_map <- ctab %>%
    dplyr::select(gene_id, t_name) %>%
    mutate(tair_id = gsub("\\..*$", "", t_name)) %>%
    filter(grepl("^AT", tair_id)) %>%
    distinct(gene_id, tair_id)

final_sig_genes$gene_id <- rownames(final_sig_genes)
final_mapped <- merge(final_sig_genes, id_map, by = "gene_id")
cat("[11] Rows with a TAIR ID:", nrow(final_mapped), "for", nrow(final_sig_genes), "genes\n")

up_ids_final   <- unique(final_mapped$tair_id[final_mapped$logFC > 1])    # higher in WT
down_ids_final <- unique(final_mapped$tair_id[final_mapped$logFC < -1])   # higher in Cond

save_result <- function(res, title, png_name, csv_name) {
    if (is.null(res)) {
        message("No result for ", png_name, " (no gene could be mapped); figure not written")
        return(invisible(NULL))
    }
    write.csv(as.data.frame(res), file.path(enrich_dir, csv_name), row.names = FALSE)
    p <- dotplot(res, showCategory = 15) + ggtitle(title)
    ggsave(file.path(asset_dir, png_name), p, width = 7, height = 6, dpi = 200, bg = "white")
}

# 2. GO enrichment (BP, MF, CC); background: all annotated genes in org.At.tair.db
run_go <- function(genes, ont) {
    enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont, pvalueCutoff = 0.05)
}
save_result(run_go(up_ids_final, "BP"),   "UP: Biological Process (BP)",   "006_pilot_go_bp_up.png",   "go_bp_up.csv")
save_result(run_go(down_ids_final, "BP"), "DOWN: Biological Process (BP)", "007_pilot_go_bp_down.png", "go_bp_down.csv")
save_result(run_go(up_ids_final, "MF"),   "UP: Molecular Function (MF)",   "008_pilot_go_mf_up.png",   "go_mf_up.csv")
save_result(run_go(down_ids_final, "MF"), "DOWN: Molecular Function (MF)", "009_pilot_go_mf_down.png", "go_mf_down.csv")
save_result(run_go(up_ids_final, "CC"),   "UP: Cellular Component (CC)",   "010_pilot_go_cc_up.png",   "go_cc_up.csv")
save_result(run_go(down_ids_final, "CC"), "DOWN: Cellular Component (CC)", "011_pilot_go_cc_down.png", "go_cc_down.csv")

# 3. KEGG enrichment with TAIR IDs
kegg_up   <- enrichKEGG(gene = up_ids_final,   organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_ids_final, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
save_result(kegg_up,   "KEGG: Up-regulated (Activated Circuits)",    "012_pilot_kegg_pathways_up.png",   "kegg_up.csv")
save_result(kegg_down, "KEGG: Down-regulated (Suppressed Circuits)", "013_pilot_kegg_pathways_down.png", "kegg_down.csv")

cat("[11] Figures 006-013 written to", asset_dir, "\n")
