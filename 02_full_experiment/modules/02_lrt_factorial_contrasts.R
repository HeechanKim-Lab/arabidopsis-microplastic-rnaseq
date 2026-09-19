#!/usr/bin/env Rscript
# ==============================================================================
# Module: 02_lrt_factorial_contrasts.R
# Description: GLM LRT contrast modeling, DEG extraction, and Venn statistics
# ==============================================================================

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(ggplot2)
})

base_dir   <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]), ".."))
data_dir   <- file.path(base_dir, "data")
inter_dir  <- file.path(base_dir, "results", "intermediate")
assets_dir <- file.path(base_dir, "assets")
tables_dir <- file.path(base_dir, "results", "tables")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

rdata_file <- file.path(data_dir, "Arabidopsis_plastics_2023.RData")
if (!file.exists(rdata_file)) rdata_file <- "Arabidopsis_plastics_2023.RData"
load(rdata_file)

cat(">>> [02] Fitting GLM LRT for Root00 (N=27)...\n")
root00_fit_LRT <- glmFit(root00_y4, root00_design_01)
contrasts_root00 <- makeContrasts(
    Plastic_5p_d0   = c5_d0   - c0_d0,
    Plastic_5p_d30  = c5_d30  - c0_d30,
    Plastic_5p_d120 = c5_d120 - c0_d120,
    Plastic_1p_d120 = c1_d120 - c0_d120,
    Incubation_5p   = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    levels          = root00_design_01
)

res00_0d    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d0"])
res00_30d   <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d30"])
res00_120d  <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_5p_d120"])
res00_1p    <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Plastic_1p_d120"])
res00_inter <- glmLRT(root00_fit_LRT, contrast = contrasts_root00[, "Incubation_5p"])

cat(">>> [02] Fitting GLM LRT for Root01 (N=12)...\n")
root01_fit_LRT <- glmFit(root01_y4, root01_design_01)
contrasts_root01 <- makeContrasts(
    p5_d0           = c5_d0   - c0_d0,
    p5_d120         = c5_d120 - c0_d120,
    d120_p0         = c0_d120 - c0_d0,
    d120_p5         = c5_d120 - c5_d0,
    plastic_5p_120d = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
    incubation      = (c5_d120 - c0_d120) - (c5_d0 - c0_d0),
    levels          = root01_design_01
)

res01_p0      <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d0"])
res01_p120    <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d120"])
res01_aging_c <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p0"])
res01_aging_p <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p5"])
res01_inter   <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "plastic_5p_120d"])

# Export Summary Counts Plot
summary_data <- data.frame(
    Dataset = c(rep("Full Dataset (Root00)", 5), rep("Subset (Root01)", 6)),
    Contrast = c(
        "D0 Plastic Effect (5%)", "D30 Plastic Effect (5%)", "D120 Plastic Effect (5%)", 
        "D120 Plastic Effect (1%)", "Overall Incubation",
        "D0 Plastic Effect", "D120 Plastic Effect", 
        "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)", 
        "Plastic x Time Interaction", "Incubation Effect"
    ),
    Count = c(
        sum(topTags(res00_0d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_30d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_120d, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_1p, n=Inf)$table$FDR < 0.05),
        sum(topTags(res00_inter, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_p0, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_p120, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_aging_c, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_aging_p, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_inter, n=Inf)$table$FDR < 0.05),
        sum(topTags(res01_inter, n=Inf)$table$FDR < 0.05)
    )
)

summary_data$Contrast <- factor(summary_data$Contrast, levels = rev(c(
    "D0 Plastic Effect", "D0 Plastic Effect (5%)", 
    "D30 Plastic Effect (5%)", 
    "D120 Plastic Effect (1%)", "D120 Plastic Effect (5%)", "D120 Plastic Effect", 
    "Aging Effect (0% Plastic)", "Aging Effect (5% Plastic)", 
    "Plastic x Time Interaction", "Overall Incubation", "Incubation Effect"
)))

p_bar <- ggplot(summary_data, aes(x = Contrast, y = Count, fill = Dataset)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.3, size = 4.2, fontface = "bold") +
    coord_flip() +
    facet_wrap(~ Dataset, scales = "free") +
    scale_fill_manual(values = c("Full Dataset (Root00)" = "#B3E2CD", "Subset (Root01)" = "#FDCDAC")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Summary of Significant Differentially Expressed Genes (FDR < 0.05)", x = "", y = "Number of Significant Genes")

ggsave(file.path(assets_dir, "027_deg_contrasts_summary.png"), plot = p_bar, width = 10, height = 6, dpi = 300)

# Calculate Venn Partition Values from strict manuscript thresholds
get_sig <- function(res_obj) {
    tab  <- as.data.frame(topTags(res_obj, n = Inf))
    up   <- rownames(tab[tab$logFC > 1  & tab$PValue < 0.01 & tab$FDR < 0.05, ])
    down <- rownames(tab[tab$logFC < -1 & tab$PValue < 0.01 & tab$FDR < 0.05, ])
    return(list(up = up, down = down, all = c(up, down)))
}

set_A <- get_sig(res01_aging_c)
set_B <- get_sig(res01_p120)
set_C <- get_sig(res01_aging_p)
total_expressed <- nrow(root01_y3$counts)

count_dir <- function(gene_list) {
    u <- sum(gene_list %in% unique(c(set_A$up, set_B$up, set_C$up)))
    d <- sum(gene_list %in% unique(c(set_A$down, set_B$down, set_C$down)))
    return(paste0(u, " ↑, ", d, " ↓"))
}

A_only <- setdiff(set_A$all, union(set_B$all, set_C$all))
B_only <- setdiff(set_B$all, union(set_A$all, set_C$all))
C_only <- setdiff(set_C$all, union(set_A$all, set_B$all))
AB_inter   <- setdiff(intersect(set_A$all, set_B$all), set_C$all)
BC_inter   <- setdiff(intersect(set_B$all, set_C$all), set_A$all)
AC_inter   <- setdiff(intersect(set_A$all, set_C$all), set_B$all)
ABC_center <- intersect(intersect(set_A$all, set_B$all), set_C$all)

sink(file.path(tables_dir, "07_venndiagram_data.txt"))
cat("--- Circle Totals (Brackets & Parentheses) ---\n")
for (i in list(list("3 vs. 1", set_A), list("4 vs. 3", set_B), list("4 vs. 2", set_C))) {
    tot <- length(i[[2]]$all)
    prc <- round((tot / total_expressed) * 100, 2)
    cat(i[[1]], ": [", tot, ", ", prc, "%] (", length(i[[2]]$up), " ↑, ", length(i[[2]]$down), " ↓)\n", sep = "")
}
cat("\n--- Venn Section Values (Type these in manually) ---\n")
cat("Top (4 vs 3) Only:        ", count_dir(B_only), "\n")
cat("Left (3 vs 1) Only:       ", count_dir(A_only), "\n")
cat("Right (4 vs 2) Only:      ", count_dir(C_only), "\n")
cat("Top-Left Inter (A & B):   ", count_dir(AB_inter), "\n")
cat("Top-Right Inter (B & C):  ", count_dir(BC_inter), "\n")
cat("Bottom Inter (A & C):     ", count_dir(AC_inter), "\n")
cat("Center Triple Inter:      ", count_dir(ABC_center), "\n")
sink()

# Serialize intermediate LRT tables
saveRDS(list(res00_0d=res00_0d, res00_30d=res00_30d, res00_120d=res00_120d, res00_1p=res00_1p, res00_inter=res00_inter),
        file.path(inter_dir, "lrt_root00_list.rds"))
saveRDS(list(res01_p0=res01_p0, res01_p120=res01_p120, res01_aging_c=res01_aging_c, res01_aging_p=res01_aging_p, res01_inter=res01_inter),
        file.path(inter_dir, "lrt_root01_list.rds"))

cat(">>> [02] LRT contrast modeling complete. Outputs saved.\n")