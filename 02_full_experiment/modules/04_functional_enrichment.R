#!/usr/bin/env Rscript
# ==============================================================================
# Module: 04_functional_enrichment.R
# Description: Functional overrepresentation analysis across interaction and aging suites
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.At.tair.db)
    library(dplyr)
    library(ggplot2)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")
tables_dir <- file.path(base_dir, "results", "tables")

lrt_root01 <- readRDS(file.path(inter_dir, "lrt_root01_list.rds"))

# 1. Enrichment for 146 Interaction Genes
res_inter_df <- as.data.frame(topTags(lrt_root01$res01_inter, n = Inf))
inter_genes  <- rownames(res_inter_df[res_inter_df$FDR < 0.05, ])

go_inter_BP <- enrichGO(gene = inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
go_inter_MF <- enrichGO(gene = inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "MF", pvalueCutoff = 0.05)
go_inter_CC <- enrichGO(gene = inter_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "CC", pvalueCutoff = 0.05)

ggsave(file.path(assets_dir, "035_interaction_go_bp.png"), dotplot(go_inter_BP, showCategory = 15), width = 8, height = 7, dpi = 300)
ggsave(file.path(assets_dir, "036_interaction_go_mf.png"), dotplot(go_inter_MF, showCategory = 15), width = 8, height = 7, dpi = 300)
ggsave(file.path(assets_dir, "037_interaction_go_cc.png"), dotplot(go_inter_CC, showCategory = 15), width = 8, height = 7, dpi = 300)

kegg_inter <- enrichKEGG(gene = inter_genes, organism = "ath", keyType = "kegg", pvalueCutoff = 0.05)
if (!is.null(kegg_inter) && nrow(as.data.frame(kegg_inter)) > 0) {
    kegg_summary <- as.data.frame(table(kegg_inter$ID))
    colnames(kegg_summary) <- c("Pathway", "GeneCount")
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
        "ath00010" = "Glycolysis / Gluconeogenesis"
    )
    kegg_summary$Description <- kegg_descriptions[as.character(kegg_summary$Pathway)]
    kegg_summary <- kegg_summary[!is.na(kegg_summary$Description), ]
    kegg_summary <- kegg_summary[order(-kegg_summary$GeneCount), ]
    
    p_kegg <- ggplot(head(kegg_summary, 10), aes(x = reorder(Description, GeneCount), y = GeneCount)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        theme_minimal(base_size = 12) +
        labs(title = "Top 10 KEGG Pathways (Incubation Genes)",
             subtitle = paste0(nrow(as.data.frame(kegg_inter)), "/146 genes mapped"),
             x = "Pathway Name", y = "Number of Genes")
    ggsave(file.path(assets_dir, "038_interaction_kegg_top10.png"), plot = p_kegg, width = 8, height = 6, dpi = 300)
}

# 2. Plastic Aging (5%) GO Enrichment
res_aging_p_df  <- as.data.frame(topTags(lrt_root01$res01_aging_p, n = Inf))
sig_aging_genes <- rownames(res_aging_p_df[res_aging_p_df$FDR < 0.05, ])
go_aging_BP     <- enrichGO(gene = sig_aging_genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
ggsave(file.path(assets_dir, "034_aging_5p_go_bp.png"), dotplot(go_aging_BP, showCategory = 20), width = 8, height = 7, dpi = 300)

# 3. Formal GO BP Tables for Condition 4 vs 2 (Aging Plastic)
up_4vs2   <- rownames(res_aging_p_df[res_aging_p_df$logFC > 1  & res_aging_p_df$PValue < 0.01 & res_aging_p_df$FDR < 0.05, ]) down_4vs2 <- rownames(res_aging_p_df[res_aging_p_df$logFC < -1 & res_aging_p_df$PValue < 0.01 & res_aging_p_df$FDR < 0.05, ])

run_go_table <- function(genes) {
    res <- enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    as.data.frame(res) %>%
        mutate(
            obs = as.numeric(sub("/.*", "", GeneRatio)),
            test = as.numeric(sub(".*/", "", GeneRatio)),
            exp = as.numeric(sub("/.*", "", BgRatio)),
            ref = as.numeric(sub(".*/", "", BgRatio)),
            Enr = round((obs / test) / (exp / ref), 2),
            "No. gene" = paste0(obs, "/", exp)
        ) %>%
        dplyr::select(Description, `No. gene`, Enr, p.adjust) %>%
        head(5)
}

t_up   <- run_go_table(up_4vs2)
t_down <- run_go_table(down_4vs2)
write.csv(t_up,   file.path(tables_dir, "08_go_bp_4vs2_up.csv"), row.names = FALSE)
write.csv(t_down, file.path(tables_dir, "08_go_bp_4vs2_down.csv"), row.names = FALSE)

cat(">>> [04] Functional enrichment analysis completed and tables exported.\n")