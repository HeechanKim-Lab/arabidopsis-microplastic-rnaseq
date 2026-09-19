#!/usr/bin/env Rscript
# ==============================================================================
# Module: 11_functional_enrichment.R
# Description: Locus ID mapping, clusterProfiler GO (BP/MF/CC) and KEGG pathways
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.At.tair.db)
    library(dplyr)
    library(ggplot2)
})

base_dir     <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir    <- file.path(base_dir, "results", "intermediate")
ballgown_dir <- file.path(base_dir, "results", "ballgown")
asset_dir    <- file.path(base_dir, "assets")

final_sig_genes <- readRDS(file.path(inter_dir, "final_sig_genes.rds"))

# 1. Structural ID Dictionary Mapping
ctab <- read.table(file.path(ballgown_dir, "R1-1_WT", "t_data.ctab"), header = TRUE, sep = "\t")
id_map <- ctab %>%
    dplyr::select(gene_id, t_name) %>%
    mutate(tair_id = gsub("\\..*$", "", t_name)) %>%
    filter(grepl("^AT", tair_id)) %>%
    distinct(gene_id, tair_id)

final_sig_genes$gene_id <- rownames(final_sig_genes)
final_mapped <- merge(final_sig_genes, id_map, by = "gene_id")
cat(">>> [11] Successfully mapped", nrow(final_mapped), "transcripts to TAIR locus IDs.\n")

up_ids   <- unique(final_mapped$tair_id[final_mapped$logFC > 1.0])
down_ids <- unique(final_mapped$tair_id[final_mapped$logFC < -1.0])

# 2. Gene Ontology Enrichment (BP, MF, CC)
run_and_save_go <- function(genes, ont_type, direction, filename) {
    res <- enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont_type, pvalueCutoff = 0.05)
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        p <- dotplot(res, showCategory = 15) + 
            ggtitle(paste0(direction, ": Biological Ontology (", ont_type, ")"))
        ggsave(file.path(asset_dir, filename), plot = p, width = 8, height = 6, dpi = 300)
    }
}

run_and_save_go(up_ids,   "BP", "UP",   "006_pilot_go_bp_up.png")
run_and_save_go(down_ids, "BP", "DOWN", "007_pilot_go_bp_down.png")
run_and_save_go(up_ids,   "MF", "UP",   "008_pilot_go_mf_up.png")
run_and_save_go(down_ids, "MF", "DOWN", "009_pilot_go_mf_down.png")
run_and_save_go(up_ids,   "CC", "UP",   "010_pilot_go_cc_up.png")
run_and_save_go(down_ids, "CC", "DOWN", "011_pilot_go_cc_down.png")

# 3. KEGG Pathway Deconvolution
ids_mapped <- bitr(final_mapped$tair_id, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)

up_entrez   <- unique(ids_mapped$ENTREZID[ids_mapped$TAIR %in% up_ids])
down_entrez <- unique(ids_mapped$ENTREZID[ids_mapped$TAIR %in% down_ids])

kegg_up   <- enrichKEGG(gene = up_entrez,   organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_entrez, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)

if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
    p_kegg_up <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG: Up-regulated (Activated Circuits)")
    ggsave(file.path(asset_dir, "012_pilot_kegg_pathways_up.png"), plot = p_kegg_up, width = 8, height = 6, dpi = 300)
}

if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
    p_kegg_down <- dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG: Down-regulated (Suppressed Circuits)")
    ggsave(file.path(asset_dir, "013_pilot_kegg_pathways_down.png"), plot = p_kegg_down, width = 8, height = 6, dpi = 300)
}

cat(">>> [11] Functional enrichment analysis completed successfully.\n")