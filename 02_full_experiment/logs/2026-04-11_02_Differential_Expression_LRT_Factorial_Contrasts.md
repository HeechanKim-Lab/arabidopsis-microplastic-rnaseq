# Likelihood Ratio Tests for the Factorial Contrasts

**Date:** 2026-04-11  
**Stage:** 02_full_experiment  
**Script:** `modules/02_lrt_factorial_contrasts.R`  
**Input:** `root00_y4`, `root01_y4` (dispersions already estimated in the lab script), design matrices `root00_design_01` and `root01_design_01` (`~0 + group`)  
**Output:** `assets/027_deg_contrasts_summary.png`, `assets/031_venn_root00_static_vs_incubation.png`, `assets/032_venn_root01_factorial.png`, `assets/033_venn_root00_vs_root01_total.png`

---

## 1. Methods

* **Design:** one coefficient per condition (`~0 + group`), so each comparison is written as a contrast between condition means.
* **Test:** edgeR GLM likelihood ratio test (`glmFit`, `glmLRT`). With three replicates per condition the LRT is less conservative than the quasi-likelihood F-test, so the gene counts below should be read with that in mind.
* **Significance:** the counts in Section 3A use FDR < 0.05 only (no fold-change filter). The up/down counts in Section 3B use |log2FC| > 1, P < 0.01 and FDR < 0.05.

| Contrast | Definition | Datasets |
| :--- | :--- | :--- |
| 5% vs 0% at day 0 | `c5_d0 - c0_d0` | Root00, Root01 |
| 5% vs 0% at day 30 | `c5_d30 - c0_d30` | Root00 |
| 5% vs 0% at day 120 | `c5_d120 - c0_d120` | Root00, Root01 |
| 1% vs 0% at day 120 | `c1_d120 - c0_d120` | Root00 |
| Day 120 vs day 0 at 0% | `c0_d120 - c0_d0` | Root01 |
| Day 120 vs day 0 at 5% | `c5_d120 - c5_d0` | Root01 |
| Plastic × incubation interaction | `(c5_d120 - c5_d0) - (c0_d120 - c0_d0)` | Root00 ("Overall Incubation"), Root01 |

In Root01 the script also tests `incubation = (c5_d120 - c0_d120) - (c5_d0 - c0_d0)`. This is algebraically the same contrast as the interaction, so it gives the same 146 genes. Figure 027 shows the number twice.

---

## 2. Code (R)

```R
# Root00 (27 libraries)
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
sum(topTags(res00_120d, n = Inf)$table$FDR < 0.05)   # and likewise for each contrast

# Root01 (12 libraries)
root01_fit_LRT <- glmFit(root01_y4, root01_design_01)
contrasts_root01 <- makeContrasts(
  p5_d0   = c5_d0   - c0_d0,
  p5_d120 = c5_d120 - c0_d120,
  d120_p0 = c0_d120 - c0_d0,
  d120_p5 = c5_d120 - c5_d0,
  plastic_5p_120d = (c5_d120 - c5_d0) - (c0_d120 - c0_d0),
  incubation      = (c5_d120 - c0_d120) - (c5_d0 - c0_d0),
  levels = root01_design_01)
res01_p0      <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d0"])
res01_p120    <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "p5_d120"])
res01_aging_c <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p0"])
res01_aging_p <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "d120_p5"])
res01_inter   <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "plastic_5p_120d"])
res01_inc     <- glmLRT(root01_fit_LRT, contrast = contrasts_root01[, "incubation"])

# Up/down counts used for the Venn diagram
get_sig <- function(res_obj) {
  tab  <- as.data.frame(topTags(res_obj, n = Inf))
  up   <- rownames(tab[tab$logFC > 1  & tab$PValue < 0.01 & tab$FDR < 0.05, ])
  down <- rownames(tab[tab$logFC < -1 & tab$PValue < 0.01 & tab$FDR < 0.05, ])
  return(list(up = up, down = down, all = c(up, down)))
}
set_A <- get_sig(res01_aging_c)   # 3 vs 1
set_B <- get_sig(res01_p120)      # 4 vs 3
set_C <- get_sig(res01_aging_p)   # 4 vs 2
total_expressed <- nrow(root01_y3$counts)   # 17,560 genes
```

---

## 3. Results

### A. Genes at FDR < 0.05 (console output)

| Contrast | Root00 (27 libraries) | Root01 (12 libraries) |
| :--- | ---: | ---: |
| 5% vs 0% at day 0 | 0 | 0 |
| 5% vs 0% at day 30 | 0 | — |
| 5% vs 0% at day 120 | 48 | 138 |
| 1% vs 0% at day 120 | 0 | — |
| Day 120 vs day 0 at 0% | — | 49 |
| Day 120 vs day 0 at 5% | — | 1,240 |
| Plastic × incubation interaction | 99 | 146 |

![Number of genes per contrast (FDR < 0.05)](../assets/027_deg_contrasts_summary.png)

### B. Up/down counts with |log2FC| > 1, P < 0.01 and FDR < 0.05 (Root01)

Conditions are numbered as in the figures: 1 = `c0_d0`, 2 = `c5_d0`, 3 = `c0_d120`, 4 = `c5_d120`.

| Comparison | Genes | % of 17,560 | Up | Down |
| :--- | ---: | ---: | ---: | ---: |
| 3 vs 1 (day 120 vs day 0 at 0%) | 49 | 0.28% | 41 | 8 |
| 4 vs 3 (5% vs 0% at day 120) | 136 | 0.77% | 43 | 93 |
| 4 vs 2 (day 120 vs day 0 at 5%) | 1,123 | 6.40% | 540 | 583 |

Overlaps (up / down):

| Region | Genes |
| :--- | :--- |
| 4 vs 3 only | 12 / 17 |
| 3 vs 1 only | 16 / 4 |
| 4 vs 2 only | 488 / 507 |
| 3 vs 1 and 4 vs 3 | 4 / 4 |
| 4 vs 3 and 4 vs 2 | 30 / 73 |
| 3 vs 1 and 4 vs 2 | 22 / 3 |
| All three | 0 / 0 |

![Venn diagram of the three Root01 comparisons](../assets/032_venn_root01_factorial.png)

Figure 032 was drawn outside R from the `get_sig()` output above. The numbers in it match that output.

### C. Overlap between datasets

![Four-set Venn: Root00 and Root01, plastic (day 120) and interaction gene sets](../assets/031_venn_root00_static_vs_incubation.png)

Figure 031 compares the four FDR < 0.05 sets: Root00 5% vs 0% at day 120 (48), Root00 interaction (99), Root01 5% vs 0% at day 120 (138) and Root01 interaction (146). Together they contain 231 genes, and 27 genes are in all four sets.

![Union of significant genes, Root00 vs Root01](../assets/033_venn_root00_vs_root01_total.png)

Figure 033 compares the union of the Root00 sets (day-120 plastic effect and interaction; 114 genes) with the union of the Root01 sets (day-120 plastic effect, interaction, and day 120 vs day 0 at 5%; 1,286 genes). 106 genes are shared, 8 are only in Root00 and 1,180 are only in Root01.

---

## 4. Interpretation

* No gene differed between 5% and 0% plastic at day 0 or day 30, or between 1% and 0% plastic at day 120.
* After 120 days of soil pre-incubation, 5% plastic changed the expression of 48 genes in the 27-library model and 138 in the 12-library model. The difference reflects the different dispersion estimates of the two models.
* Comparing day 120 with day 0, 1,240 genes changed in 5% soil and 49 in 0% soil. The interaction contrast, which tests this difference formally, gave 146 genes (Root01) and 99 genes (Root00).
* These lists come from LRTs with three replicates per condition and have not been validated by another method.

Next step: volcano plot and heatmap of the interaction genes (`modules/03_interaction_visualization.R`).
