#!/usr/bin/env Rscript
# ==============================================================================
# Module 02: likelihood ratio tests for the factorial contrasts (log 2026-04-11_02)
#
# Input : data/Arabidopsis_plastics_2023.RData (root00_y4, root01_y4 with dispersions,
#         root00_design_01, root01_design_01 with design ~0 + group, root01_y3)
# Output: assets/027_deg_contrasts_summary.png
#         assets/031_venn_root00_static_vs_incubation.png (four gene sets: Root00 and Root01,
#         day-120 plastic effect and interaction)
#         assets/033_venn_root00_vs_root01_total.png
#         results/tables/root01_venn_counts.txt (numbers for figure 032, which was drawn outside R)
#         results/tables/lrt_gene_counts_fdr005.csv, results/tables/lrt_<contrast>.csv
#         results/intermediate/lrt_root00_list.rds, lrt_root01_list.rds
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
    library(ggvenn)
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
asset_dir  <- file.path(base_dir, "assets")
for (d in c(inter_dir, tables_dir, asset_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

find_input <- function(name) {
    for (p in c(file.path(data_dir, name), file.path(getwd(), name))) if (file.exists(p)) return(p)
    stop("Input not found: ", name, " (expected in ", data_dir, ")")
}
load(find_input("Arabidopsis_plastics_2023.RData"))

# ---- Root00 (27 libraries) ----------------------------------------------------
root00_fit_LRT <- glmFit(root00_y4, root00_design_01)
contrasts_root00 <- makeContrasts(
    Plastic_5p_d0   = c5_d0   - c0_d0,
    Plastic_5p_d30  = c5_d30  - c0_d30,
    Plastic_5p_d120 = c5_d120 - c0_d120,
    Plastic_1p_d120 = c1_d120 - c0_d120,
    Incubation_5p   = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    levels = root00_design_01)
res00_0d    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d0"])
res00_30d   <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d30"])
res00_120d  <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d120"])
res00_1p    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_1p_d120"])
res00_inter <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Incubation_5p"])

# ---- Root01 (12 libraries) ----------------------------------------------------
root01_fit_LRT <- glmFit(root01_y4, root01_design_01)
contrasts_root01 <- makeContrasts(
    p5_d0   = c5_d0   - c0_d0,
    p5_d120 = c5_d120 - c0_d120,
    d120_p0 = c0_d120 - c0_d0,
    d120_p5 = c5_d120 - c5_d0,
    plastic_5p_120d = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    incubation      = (c5_d120 - c0_d120) - (c5_d0 - c0_d0),   # same contrast as plastic_5p_120d
    levels = root01_design_01)
res01_p0      <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d0"])
res01_p120    <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d120"])
res01_aging_c <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p0"])
res01_aging_p <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p5"])
res01_inter   <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "plastic_5p_120d"])
res01_inc     <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "incubation"])

lrt_root00 <- list(res00_0d = res00_0d, res00_30d = res00_30d, res00_120d = res00_120d,
                   res00_1p = res00_1p, res00_inter = res00_inter)
lrt_root01 <- list(res01_p0 = res01_p0, res01_p120 = res01_p120, res01_aging_c = res01_aging_c,
                   res01_aging_p = res01_aging_p, res01_inter = res01_inter, res01_inc = res01_inc)
for (nm in names(c(lrt_root00, lrt_root01))) {
    write.csv(topTags(c(lrt_root00, lrt_root01)[[nm]], n = Inf)$table,
              file.path(tables_dir, paste0("lrt_", nm, ".csv")))
}

# ---- Genes at FDR < 0.05 per contrast (figure 027) ----------------------------
n_sig <- function(res) sum(topTags(res, n = Inf)$table$FDR < 0.05)
summary_data <- data.frame(
    Dataset = c(rep("Full Dataset (Root00)", 5), rep("Subset (Root01)", 6)),
    Contrast = c(
        "D0 Plastic Effect (5%)", "D30 Plastic Effect (5%)", "D120 Plastic Effect (5%)",
        "D120 Plastic Effect (1%)", "Overall Incubation",
        "D0 Plastic Effect", "D120 Plastic Effect",
        "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)",
        "Plastic x Time Interaction", "Incubation Effect"),
    Count = c(n_sig(res00_0d), n_sig(res00_30d), n_sig(res00_120d), n_sig(res00_1p), n_sig(res00_inter),
              n_sig(res01_p0), n_sig(res01_p120), n_sig(res01_aging_c), n_sig(res01_aging_p),
              n_sig(res01_inter), n_sig(res01_inc)))
cat("[02] Genes at FDR < 0.05\n")
print(summary_data, row.names = FALSE)
write.csv(summary_data, file.path(tables_dir, "lrt_gene_counts_fdr005.csv"), row.names = FALSE)

summary_data$Contrast <- factor(summary_data$Contrast, levels = rev(c(
    "D0 Plastic Effect", "D0 Plastic Effect (5%)",
    "D30 Plastic Effect (5%)",
    "D120 Plastic Effect (1%)", "D120 Plastic Effect (5%)", "D120 Plastic Effect",
    "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)",
    "Plastic x Time Interaction", "Overall Incubation", "Incubation Effect")))

p_summary <- ggplot(summary_data, aes(x = Contrast, y = Count, fill = Dataset)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.3, size = 4.5, fontface = "bold") +
    coord_flip() +
    facet_wrap(~ Dataset, scales = "free") +
    scale_fill_manual(values = c("Full Dataset (Root00)" = "#B3E2CD", "Subset (Root01)" = "#FDCDAC")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0"),
        plot.title = element_text(size = 15, face = "bold"),
        panel.grid.major.y = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Summary of Significant Differentially Expressed Genes (FDR < 0.05)",
         subtitle = "Comparing Contrast Yields Across Data Subsets",
         x = "", y = "Number of Significant Genes")
ggsave(file.path(asset_dir, "027_deg_contrasts_summary.png"), p_summary,
       width = 11, height = 4.7, dpi = 200, bg = "white")

# ---- Venn diagrams of the FDR < 0.05 gene sets (figures 031 and 033) ----------
sig_genes <- function(res) {
    tab <- topTags(res, n = Inf)$table
    rownames(tab[tab$FDR < 0.05, ])
}
list_00_p120    <- sig_genes(res00_120d)
list_00_inter   <- sig_genes(res00_inter)
list_01_p120    <- sig_genes(res01_p120)
list_01_inter   <- sig_genes(res01_inter)
list_01_aging_p <- sig_genes(res01_aging_p)

venn_list_4way <- list(
    "root00: Plastic"    = list_00_p120,
    "root00: Incubation" = list_00_inter,
    "root01: Plastic"    = list_01_p120,
    "root01: Incubation" = list_01_inter)
p_031 <- ggvenn(venn_list_4way,
                fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                stroke_size = 0.4, set_name_size = 3) +
    labs(title = "LRT: Full(root00) vs. Subset(root01)", subtitle = "Plastic and Incubation")
ggsave(file.path(asset_dir, "031_venn_root00_static_vs_incubation.png"), p_031,
       width = 7, height = 5, dpi = 200, bg = "white")

venn_final <- list(
    "root00 Total Impact" = unique(c(list_00_p120, list_00_inter)),
    "root01 Total Impact" = unique(c(list_01_p120, list_01_inter, list_01_aging_p)))
p_033 <- ggvenn(venn_final, fill_color = c("#B3E2CD", "#FDCDAC"),
                stroke_size = 0.5, set_name_size = 5, show_percentage = TRUE) +
    labs(title = "Total System Impact: Full Dataset vs. Characteristic Subset",
         subtitle = "Comparing the Union of All Significant LRT Discoveries")
ggsave(file.path(asset_dir, "033_venn_root00_vs_root01_total.png"), p_033,
       width = 7, height = 5, dpi = 200, bg = "white")

# ---- Numbers for figure 032 (Root01; |log2FC| > 1, P < 0.01, FDR < 0.05) ------
# Conditions: 1 = c0_d0, 2 = c5_d0, 3 = c0_d120, 4 = c5_d120
get_sig <- function(res_obj) {
    tab  <- as.data.frame(topTags(res_obj, n = Inf))
    up   <- rownames(tab[tab$logFC > 1  & tab$PValue < 0.01 & tab$FDR < 0.05, ])
    down <- rownames(tab[tab$logFC < -1 & tab$PValue < 0.01 & tab$FDR < 0.05, ])
    list(up = up, down = down, all = c(up, down))
}
set_A <- get_sig(res01_aging_c)   # 3 vs 1
set_B <- get_sig(res01_p120)      # 4 vs 3
set_C <- get_sig(res01_aging_p)   # 4 vs 2
total_expressed <- nrow(root01_y3$counts)

count_dir <- function(gene_list) {
    u <- sum(gene_list %in% unique(c(set_A$up, set_B$up, set_C$up)))
    d <- sum(gene_list %in% unique(c(set_A$down, set_B$down, set_C$down)))
    paste0(u, " up, ", d, " down")
}
A_only     <- setdiff(set_A$all, union(set_B$all, set_C$all))
B_only     <- setdiff(set_B$all, union(set_A$all, set_C$all))
C_only     <- setdiff(set_C$all, union(set_A$all, set_B$all))
AB_inter   <- setdiff(intersect(set_A$all, set_B$all), set_C$all)
BC_inter   <- setdiff(intersect(set_B$all, set_C$all), set_A$all)
AC_inter   <- setdiff(intersect(set_A$all, set_C$all), set_B$all)
ABC_center <- intersect(intersect(set_A$all, set_B$all), set_C$all)

venn_lines <- "Circle totals: genes [n, % of expressed genes] (up, down)"
for (i in list(list("3 vs. 1", set_A), list("4 vs. 3", set_B), list("4 vs. 2", set_C))) {
    total <- length(i[[2]]$all)
    perc  <- round((total / total_expressed) * 100, 2)
    venn_lines <- c(venn_lines, paste0(i[[1]], ": [", total, ", ", perc, "%] (",
                                       length(i[[2]]$up), " up, ", length(i[[2]]$down), " down)"))
}
venn_lines <- c(venn_lines, "", "Venn sections",
                paste("4 vs 3 only:          ", count_dir(B_only)),
                paste("3 vs 1 only:          ", count_dir(A_only)),
                paste("4 vs 2 only:          ", count_dir(C_only)),
                paste("3 vs 1 and 4 vs 3:    ", count_dir(AB_inter)),
                paste("4 vs 3 and 4 vs 2:    ", count_dir(BC_inter)),
                paste("3 vs 1 and 4 vs 2:    ", count_dir(AC_inter)),
                paste("All three:            ", count_dir(ABC_center)))
writeLines(venn_lines, file.path(tables_dir, "root01_venn_counts.txt"))
cat("\n", venn_lines, sep = "\n")

saveRDS(lrt_root00, file.path(inter_dir, "lrt_root00_list.rds"))
saveRDS(lrt_root01, file.path(inter_dir, "lrt_root01_list.rds"))
cat("[02] Figures 027, 031 and 033 written to", asset_dir, "\n")
