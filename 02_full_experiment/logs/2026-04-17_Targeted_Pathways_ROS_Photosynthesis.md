# Keyword-Selected Heatmaps: Oxidative Stress and Photosynthesis/Plastid Genes (Root01)

**Date:** 2026-04-17  
**Stage:** 02_full_experiment  
**Script:** `modules/06_targeted_pathways.R` (see the note in Section 5)  
**Input:** Root01 LRT results (`res01_p120`, `res01_inter`, `res01_aging_p`), `TMM_root01_log`, `Arabidopsis_genes_description.xlsx` (lab annotation table with `Gene_locus`, `Symbol` and `Description`)  
**Output:** `assets/039_pathway_ros_heatmap.png`, `assets/040_pathway_photosynthesis_heatmap.png`

---

## 1. Methods

1. **Gene pool:** the union of genes with FDR < 0.05 in three Root01 contrasts: 5% vs 0% at day 120 (138 genes), the interaction (146) and day 120 vs day 0 at 5% (1,240).
2. **Keyword selection:** genes in the pool whose `Description` in the annotation table matched (case-insensitive):
   * oxidative stress: `hydrogen peroxide|oxidative|reactive oxygen|peroxidase` (12 genes);
   * photosynthesis/plastid: `photosynthe|chloroplast|photosystem|light-harvesting|thylakoid` (7 genes).
3. **Heatmaps:** TMM log-CPM values of the selected genes, scaled by row, clustered with `pheatmap` defaults (Euclidean distance, complete linkage). Rows are labelled with the gene symbol when one is available, otherwise with the TAIR locus.

This is a descriptive selection by keyword, not an enrichment test.

---

## 2. Code (R)

```R
df_p120    <- as.data.frame(topTags(res01_p120, n = Inf))
df_inter   <- as.data.frame(topTags(res01_inter, n = Inf))
df_aging_p <- as.data.frame(topTags(res01_aging_p, n = Inf))
all_lrt_sig <- unique(c(
  rownames(df_p120[df_p120$FDR < 0.05, ]),
  rownames(df_inter[df_inter$FDR < 0.05, ]),
  rownames(df_aging_p[df_aging_p$FDR < 0.05, ])))

anno_file <- read_excel("Arabidopsis_genes_description.xlsx")
anno_sig  <- anno_file %>% filter(Gene_locus %in% all_lrt_sig)

# Oxidative stress genes
ox_keywords <- "hydrogen peroxide|oxidative|reactive oxygen|peroxidase"
ox_genes  <- anno_sig %>% filter(grepl(ox_keywords, Description, ignore.case = TRUE))
ox_matrix <- TMM_root01_log[ox_genes$Gene_locus, ]
rownames(ox_matrix) <- ifelse(is.na(ox_genes$Symbol) | ox_genes$Symbol == "",
                              ox_genes$Gene_locus, ox_genes$Symbol)
pheatmap(ox_matrix, scale = "row", annotation_col = ann_col,
         main = paste("Oxidative Stress Responsiveness (", nrow(ox_genes), " Genes)"),
         color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_row = 8)

# Photosynthesis / plastid genes
photo_keywords <- "photosynthe|chloroplast|photosystem|light-harvesting|thylakoid"
photo_genes  <- anno_sig %>% filter(grepl(photo_keywords, Description, ignore.case = TRUE))
photo_matrix <- TMM_root01_log[photo_genes$Gene_locus, ]
rownames(photo_matrix) <- ifelse(is.na(photo_genes$Symbol) | photo_genes$Symbol == "",
                                 photo_genes$Gene_locus, photo_genes$Symbol)
pheatmap(photo_matrix, scale = "row", annotation_col = ann_col,
         main = paste("Photosynthesis & Plastid Shifts (", nrow(photo_genes), " Genes)"),
         color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_row = 8)
```

---

## 3. Results

### A. Oxidative stress (12 genes)

![Oxidative stress genes, row Z-scores](../assets/039_pathway_ros_heatmap.png)

Genes in the figure: APX1, TPX2, DOX1, OXS3, AT1G71695, AT5G14110, AT2G18980, AT1G30870, AT3G49960, AT5G19880, AT3G28200 and AT5G14130.

* The `c5_d120` libraries (R9-1, R9-2, R9-3) cluster together with one `c0_d120` library (R7-3).
* TPX2, AT1G71695, DOX1, AT5G14110 and OXS3 are lower in the `c5_d120` libraries than in most other libraries.
* APX1 and AT5G14130 are higher in two of the three `c5_d120` libraries (R9-2 and R9-3). AT5G19880 and AT3G28200 are higher in all three.
* AT1G30870 and AT3G49960 are higher in R9-3, R7-3 and R7-1 and lower in R1-1 and R3-1.

### B. Photosynthesis/plastid (7 genes)

![Photosynthesis and plastid genes, row Z-scores](../assets/040_pathway_photosynthesis_heatmap.png)

Genes in the figure: CAO, PSB27, AT2G28660, LHCB4.3, CT-BMY, AT1G45545 and AT5G16730.

* The three `c5_d120` libraries form one cluster.
* CAO, PSB27 and AT2G28660 are lower in all three `c5_d120` libraries.
* AT5G16730 is higher in all three `c5_d120` libraries. LHCB4.3, CT-BMY and AT1G45545 are higher in two of the three.

---

## 4. Interpretation and limits

* Within these keyword-selected sets, the `c5_d120` libraries differ most from the others. This fits the DE results, where most genes changed between day 120 and day 0 in 5% soil.
* The genes were chosen by keyword from a list of significant genes and are shown for description only. The heatmaps do not test whether oxidative-stress or plastid pathways as a whole are affected.

---

## 5. Note on the script

`modules/06_targeted_pathways.R` saves the figures as `039_targeted_ros_heatmap.png` and `040_targeted_photosynthesis_heatmap.png`. The committed files are named `039_pathway_ros_heatmap.png` and `040_pathway_photosynthesis_heatmap.png`.

---

## 6. Summary of the full experiment

* **PCA:** day-120 libraries tend to separate from day-0 libraries; no grouping by plastic level (2026-04-11 log).
* **DE (LRT, FDR < 0.05):** no plastic effect at day 0 or day 30. At day 120, 48 genes (Root00) and 138 genes (Root01) differed between 5% and 0%. Day 120 vs day 0 changed 1,240 genes in 5% soil and 49 in 0% soil. The interaction gave 99 genes (Root00) and 146 genes (Root01).
* **PERMANOVA (Root01):** Day R² = 0.195, p = 0.010; Plastic p = 0.722; Day × Plastic p = 0.184 (2026-04-15 log).
* **Enrichment:** both gene sets are enriched for responses to nutrient levels and nitrogen compounds. The interaction genes are also enriched for metal-ion responses and anion transporters (2026-04-13 log).
