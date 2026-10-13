#!/usr/bin/env Rscript
# ==============================================================================
# Module 01: PCA of Root00 (27 libraries) and Root01 (12 libraries) (log 2026-04-11_01)
#
# Input : data/Arabidopsis_plastics_2023.RData (lab file; root00_y3, root01_y3 and the
#         plotting symbols and colours of the lab script)
# Output: assets/020_root00_full_matrix_pca.png, assets/021_root00_pca_pc1_vs_pc3.png,
#         assets/022_root01_dose_response_pca.png
#         results/intermediate/TMM_root00_log.rds, TMM_root00_pca.rds,
#         TMM_root01_log.rds, TMM_root01_pca.rds
#
# Colour and symbol shape show the pre-incubation time; symbol fill shows the plastic level.
# The legend lists the conditions with their symbols. The committed figures 020-022 have a
# legend that labels the colours as plastic levels and axis percentages typed by hand.
# ==============================================================================

suppressPackageStartupMessages(library(edgeR))

get_base_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 1) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
}
base_dir  <- get_base_dir()
data_dir  <- file.path(base_dir, "data")
inter_dir <- file.path(base_dir, "results", "intermediate")
asset_dir <- file.path(base_dir, "assets")
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

find_input <- function(name) {
    for (p in c(file.path(data_dir, name), file.path(getwd(), name))) if (file.exists(p)) return(p)
    stop("Input not found: ", name, " (expected in ", data_dir, ")")
}
load(find_input("Arabidopsis_plastics_2023.RData"))

# One legend entry per condition, with the symbol and colour used for its libraries
plot_pca <- function(pca, x_pc, y_pc, pch, col, groups, main, out_png) {
    var_pct <- round(summary(pca)$importance[2, ] * 100, 1)
    first   <- !duplicated(groups)
    png(out_png, width = 1500, height = 1000, res = 150)
    par(mar = c(5, 4, 4, 8))
    plot(pca$x[, x_pc], pca$x[, y_pc],
         pch = pch, col = col, cex = 1.5, lwd = 1.5,
         xlab = paste0("PC", x_pc, " (", var_pct[x_pc], "%)"),
         ylab = paste0("PC", y_pc, " (", var_pct[y_pc], "%)"),
         main = main)
    # Legend in the right margin, outside the plotting area
    usr <- par("usr")
    legend(x = usr[2] + 0.02 * diff(usr[1:2]), y = usr[4], xjust = 0, xpd = TRUE,
           legend = groups[first], pch = pch[first], col = col[first],
           pt.lwd = 1.5, bty = "n", cex = 0.9)
    invisible(dev.off())
}

# ---- Root00 (27 libraries) ----------------------------------------------------
TMM_root00_log <- cpm(root00_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root00_pca <- prcomp(t(TMM_root00_log), center = TRUE, scale = TRUE)
cat("[01] Root00 variance explained:\n")
print(summary(TMM_root00_pca)$importance[, 1:5])

groups00 <- as.character(root00_y3$samples$group)
plot_pca(TMM_root00_pca, 1, 2, root00_pca_points, TMM_root00_subtype.col.vec01, groups00,
         "PCA: Full 27-Sample Root Dataset", file.path(asset_dir, "020_root00_full_matrix_pca.png"))
plot_pca(TMM_root00_pca, 1, 3, root00_pca_points, TMM_root00_subtype.col.vec01, groups00,
         "PCA: PC1 vs. PC3 (root00)", file.path(asset_dir, "021_root00_pca_pc1_vs_pc3.png"))

# ---- Root01 (12 libraries) ----------------------------------------------------
TMM_root01_log <- cpm(root01_y3, normalized.lib.sizes = TRUE, log = TRUE)
TMM_root01_pca <- prcomp(t(TMM_root01_log), center = TRUE, scale = TRUE)
cat("[01] Root01 variance explained:\n")
print(summary(TMM_root01_pca)$importance[, 1:5])

groups01 <- as.character(root01_y3$samples$group)
plot_pca(TMM_root01_pca, 1, 2, root01_pca_points, TMM_root01_subtype.col.vec01, groups01,
         "PCA: Characteristic 12-Sample Subset", file.path(asset_dir, "022_root01_dose_response_pca.png"))

saveRDS(TMM_root00_log, file.path(inter_dir, "TMM_root00_log.rds"))
saveRDS(TMM_root00_pca, file.path(inter_dir, "TMM_root00_pca.rds"))
saveRDS(TMM_root01_log, file.path(inter_dir, "TMM_root01_log.rds"))
saveRDS(TMM_root01_pca, file.path(inter_dir, "TMM_root01_pca.rds"))
cat("[01] Figures 020-022 written to", asset_dir, "\n")
