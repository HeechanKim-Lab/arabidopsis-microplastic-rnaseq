#!/usr/bin/env Rscript
# ==============================================================================
# Module 04: GO and KEGG analysis of the interaction genes and of day 120 vs day 0
#            in 5% soil (log 2026-04-13_02)
#
# Input : results/intermediate/lrt_root01_list.rds (module 02)
# Output: assets/034_aging_5p_go_bp.png, assets/035_interaction_go_bp.png,
#         assets/036_interaction_go_mf.png, assets/037_interaction_go_cc.png,
#         assets/038_interaction_kegg_top10.png
#         results/enrichment/*.csv
#
# Gene sets: all genes with FDR < 0.05, not split by direction. KEGG uses TAIR locus IDs
# (keyType = "kegg"); enrichKEGG and bitr_kegg download KEGG data and need an internet connection.
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
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
base_dir   <- get_base_dir()
inter_dir  <- file.path(base_dir, "results", "intermediate")
enrich_dir <- file.path(base_dir, "results", "enrichment")
asset_dir  <- file.path(base_dir, "assets")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

lrt_root01 <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

save_dotplot <- function(res, title, png_name, csv_name, n = 15) {
    if (is.null(res)) {
        message("No result for ", png_name, " (no gene could be mapped); figure not written")
        return(invisible(NULL))
    }
    write.csv(as.data.frame(res), file.path(enrich_dir, csv_name), row.names = FALSE)
    p <- dotplot(res, showCategory = n) + ggtitle(title)
    ggsave(file.path(asset_dir, png_name), p, width = 7, height = 5.5, dpi = 200, bg = "white")
}
run_go <- function(genes, ont) {
    enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont, pvalueCutoff = 0.05)
}

# ---- Interaction genes (figures 035-037) --------------------------------------
res01_inter_df  <- as.data.frame(topTags(lrt_root01$res01_inter, n = Inf))
sig_inter_genes <- rownames(res01_inter_df[res01_inter_df$FDR < 0.05, ])
cat("[04] Interaction genes:", length(sig_inter_genes), "\n")

save_dotplot(run_go(sig_inter_genes, "BP"), "Incubation: Biological Process (BP)",
             "035_interaction_go_bp.png", "interaction_go_bp.csv")
save_dotplot(run_go(sig_inter_genes, "MF"), "Incubation: Molecular Function (MF)",
             "036_interaction_go_mf.png", "interaction_go_mf.csv")
save_dotplot(run_go(sig_inter_genes, "CC"), "Incubation: Cellular Component (CC)",
             "037_interaction_go_cc.png", "interaction_go_cc.csv")
# GO:0180055 = response to metal ion starvation (the term name is missing in figure 035)

# ---- Day 120 vs day 0 in 5% soil (figure 034) ---------------------------------
res01_aging_p_df  <- as.data.frame(topTags(lrt_root01$res01_aging_p, n = Inf))
sig_aging_p_genes <- rownames(res01_aging_p_df[res01_aging_p_df$FDR < 0.05, ])
cat("[04] Day 120 vs day 0 at 5%:", length(sig_aging_p_genes), "genes\n")
save_dotplot(run_go(sig_aging_p_genes, "BP"), "Aging Effect (5% Plastic): Total Transcriptomic Shift",
             "034_aging_5p_go_bp.png", "aging_5p_go_bp.csv", n = 20)

# ---- KEGG for the interaction genes -----------------------------------------
kegg_inter <- enrichKEGG(gene = sig_inter_genes, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
if (!is.null(kegg_inter)) {
    write.csv(as.data.frame(kegg_inter), file.path(enrich_dir, "interaction_kegg_enrichment.csv"), row.names = FALSE)
}

# ---- KEGG pathway membership of the interaction genes (figure 038) -------------
# Counts of genes per pathway, not an enrichment test. The original session used a table
# kegg_check with the columns kegg (gene) and Path (pathway); bitr_kegg() returns this table.
kegg_check <- bitr_kegg(sig_inter_genes, fromType = "kegg", toType = "Path", organism = "ath")
write.csv(kegg_check, file.path(enrich_dir, "interaction_kegg_membership.csv"), row.names = FALSE)

kegg_summary <- as.data.frame(table(kegg_check$Path))
colnames(kegg_summary) <- c("Pathway", "GeneCount")
kegg_summary <- kegg_summary[order(-kegg_summary$GeneCount), ]

kegg_descriptions <- c(
    "ath01100" = "Metabolic pathways",
    "ath01110" = "Biosynthesis of secondary metabolites",
    "ath04120" = "Ubiquitin mediated proteolysis",
    "ath04141" = "Protein processing in endoplasmic reticulum",
    "ath04075" = "Plant hormone signal transduction",
    "ath01200" = "Carbon metabolism",
    "ath00940" = "Phenylpropanoid biosynthesis",
    "ath00270" = "Cysteine and methionine metabolism",
    "ath00030" = "Pentose phosphate pathway",
    "ath00010" = "Glycolysis / Gluconeogenesis")
kegg_summary$Description <- kegg_descriptions[as.character(kegg_summary$Pathway)]
kegg_summary$Description <- ifelse(is.na(kegg_summary$Description),
                                   as.character(kegg_summary$Pathway), kegg_summary$Description)

n_kegg <- length(unique(kegg_check$kegg))
p_kegg <- ggplot(kegg_summary[1:min(10, nrow(kegg_summary)), ],
                 aes(x = reorder(Description, GeneCount), y = GeneCount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 10 KEGG Pathways (Incubation Genes)",
         subtitle = paste0(n_kegg, "/", length(sig_inter_genes), " genes mapped"),
         x = "Pathway Name",
         y = "Number of Genes")
ggsave(file.path(asset_dir, "038_interaction_kegg_top10.png"), p_kegg,
       width = 7.5, height = 5.5, dpi = 200, bg = "white")

# ---- Mapping rates --------------------------------------------------------------
go_check <- bitr(sig_inter_genes, fromType = "TAIR", toType = "GO", OrgDb = org.At.tair.db)
cat("KEGG Mapping Rate: ", n_kegg, "out of", length(sig_inter_genes), "\n")
cat("GO Mapping Rate:   ", length(unique(go_check$TAIR)), "out of", length(sig_inter_genes), "\n")

# ---- GO BP tables for day 120 vs day 0 at 5% (|log2FC| > 1, P < 0.01, FDR < 0.05) --
up_4vs2   <- rownames(res01_aging_p_df[res01_aging_p_df$logFC > 1  & res01_aging_p_df$PValue < 0.01 &
                                        res01_aging_p_df$FDR < 0.05, ])
down_4vs2 <- rownames(res01_aging_p_df[res01_aging_p_df$logFC < -1 & res01_aging_p_df$PValue < 0.01 &
                                        res01_aging_p_df$FDR < 0.05, ])

run_go_table <- function(genes) {
    go_res <- run_go(genes, "BP")
    if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0) return(NULL)
    # Enrichment = (observed / tested) / (expected / reference), from GeneRatio and BgRatio
    as.data.frame(go_res) %>%
        mutate(obs  = as.numeric(sub("/.*", "", GeneRatio)),
               test = as.numeric(sub(".*/", "", GeneRatio)),
               exp  = as.numeric(sub("/.*", "", BgRatio)),
               ref  = as.numeric(sub(".*/", "", BgRatio)),
               Enr  = round((obs / test) / (exp / ref), 2),
               "No. gene" = paste0(obs, "/", exp)) %>%
        dplyr::select(Description, `No. gene`, Enr, p.adjust) %>%
        head(5)
}
table_up   <- run_go_table(up_4vs2)
table_down <- run_go_table(down_4vs2)
cat("\n--- GO BP, day 120 vs day 0 at 5%, up ---\n");   print(table_up)
cat("\n--- GO BP, day 120 vs day 0 at 5%, down ---\n"); print(table_down)
if (!is.null(table_up))   write.csv(table_up,   file.path(enrich_dir, "go_bp_4vs2_up_top5.csv"),   row.names = FALSE)
if (!is.null(table_down)) write.csv(table_down, file.path(enrich_dir, "go_bp_4vs2_down_top5.csv"), row.names = FALSE)

cat("[04] Figures 034-038 written to", asset_dir, "\n")
