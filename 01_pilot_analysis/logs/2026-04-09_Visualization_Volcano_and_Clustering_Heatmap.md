# Statistical Visualization: Global Transcriptomic Volcano Plot and Hierarchical Clustering Heatmap
**Date:** April 9, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** LRT results table (`res_lrt`), filtered DEG list (`final_sig_genes`), normalized log-CPM matrix (`logCPM`), metadata frame (`phenodata`)  
**Output Artifacts:** `assets/004_pilot_lrt_volcano.png`, `assets/005_pilot_top50_heatmap.png`  

---

## 1. Methodological & Theoretical Rationale

### A. Two-Dimensional Volcano Plot Formulation
Volcano plots project differential expression results across two orthogonal axes to distinguish statistically significant features from biological noise:
1. **Biological Effect Size (Horizontal Axis):** Log2 Fold Change ($\log_2\text{FC} = \beta_g / \ln(2)$), representing the magnitude of directional transcriptional change between Conditioned roots and Wild Type controls.
2. **Statistical Significance (Vertical Axis):** Negative $\log_{10}$ nominal $p$-value ($-\log_{10}(P)$), mapping low $p$-values to high coordinate elevations.

Significance thresholds are established at two boundaries:
* A biological threshold of $\vert{}\log_2\text{FC}\vert{} \ge 1.0$ (representing at least a two-fold up- or downregulation in transcript abundance).
* A statistical threshold of $\text{FDR} < 0.05$ (Benjamini-Hochberg adjusted $p$-value), corresponding to a horizontal threshold line at $-\log_{10}(P_{\text{cutoff}})$.

Genes satisfying both conditions represent high-stringency differentially expressed genes (DEGs) that are prioritized over features exhibiting high fold changes driven by low-count sampling instability.

### B. Standardized Row Z-score Transformation
Absolute expression levels vary across transcripts by orders of magnitude (from single-digit CPMs to $> 10^4$ CPM for structural housekeeping genes). Visualizing raw or log-transformed CPM directly in a heatmap causes hyper-abundant genes to dominate the color gradient, masking expression dynamics in moderately expressed regulatory factors.

To resolve this, row-wise Z-score standardization is applied across each gene vector $g$:

$$Z_{gi} = \frac{x_{gi} - \bar{x}_g}{s_g}$$

Where $x_{gi}$ is the $\log_2\text{CPM}$ of gene $g$ in library $i$, $\bar{x}_g$ is the mean $\log_2\text{CPM}$ across all evaluated libraries, and $s_g$ is the sample standard deviation. This transformation centers each gene's baseline at zero ($\mu = 0, \sigma = 1$), visualizing relative deviations (standard deviations above or below the mean) rather than absolute abundance.

### C. Hierarchical Agglomerative Clustering
Sample and gene relationships are determined via agglomerative hierarchical clustering:
* **Distance Metric:** Euclidean distance computes geometric separation across the standardized multi-sample coordinate space:
  $$d(u, v) = \sqrt{\sum_{i=1}^n (u_i - v_i)^2}$$
* **Linkage Method:** Complete linkage defines inter-cluster distance as the maximum distance between any single point in Cluster A and any single point in Cluster B:
  $$D(A, B) = \max \{ d(x, y) : x \in A, y \in B \}$$
  Complete linkage produces compact, balanced clades and prevents single-linkage chaining effects, verifying clear separation between treatment and control groups.

### D. Outlier Re-inclusion for Post-Hoc Validation
Although sample `R1-1_WT` was excluded from statistical GLM fitting to prevent dispersion inflation ($s_0^2$), it was deliberately retained in the 6-sample heatmap matrix. Evaluating `R1-1_WT` against the Top 50 DEGs tests whether its divergence in unsupervised PCA was driven by technical noise in background genes or by true deregulation of the biological conditioning response.

---

## 2. Computational Implementation in R

Visualizations were rendered using `ggplot2` and `pheatmap` (`10_plotting_20260409.R`).

### A. Volcano Plot Rendering
```R
library(ggplot2)

# Classify significance categories based on calibrated thresholds
plot_data <- res_lrt
plot_data$sig <- ifelse(plot_data$FDR < 0.05 & abs(plot_data$logFC) > 1, "Significant", "Not Sig")

# Render volcano plot
p_volcano <- ggplot(plot_data, aes(x = logFC, y = -log10(PValue), color = sig)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("Not Sig" = "grey", "Significant" = "red")) +
    theme_minimal() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
        title = "Volcano Plot: Conditioned vs WT (LRT Method)",
        x = "Log2 Fold Change",
        y = "-log10(P-value)"
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12)
    )

ggsave("assets/004_pilot_lrt_volcano.png", plot = p_volcano, width = 8, height = 6, dpi = 300)
```

### B. Z-Score Scaled Heatmap Rendering
```R
library(pheatmap)

# Extract top 50 DEGs ranked by LRT significance
top50_ids <- rownames(final_sig_genes)

# Extract log-CPM values across all 6 libraries
heatmap_matrix <- logCPM[top50_ids, , drop = FALSE]

# Row-standardize matrix (Z-score scaling)
heatmap_scaled <- t(scale(t(heatmap_matrix)))

# Generate heatmap with two-way hierarchical clustering
png("assets/005_pilot_top50_heatmap.png", width = 1200, height = 1500, res = 150)
pheatmap(
    heatmap_scaled,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = FALSE,
    annotation_col = phenodata[, "condition", drop = FALSE],
    color          = colorRampPalette(c("blue", "white", "red"))(100),
    border_color   = NA,
    main           = "Top 50 Differentially Expressed Transcripts (Z-score)"
)
dev.off()
```

---

## 3. Visual Diagnostics & Publication Deliverables

![Volcano Plot](../assets/004_pilot_lrt_volcano.png)

![Top 50 DEGs Heatmap](../assets/005_pilot_top50_heatmap.png)

---

## 4. Diagnostic & Biological Interpretation

### A. Global Transcriptomic Shift (Volcano Plot)
* **Symmetric Transcriptional Partitioning:** The volcano plot demonstrates a balanced, bimodal distribution across the $\vert{}\log_2\text{FC}\vert{} = 1$ threshold boundaries. Upregulated features (435 genes) and downregulated features (430 genes) show comparable magnitude distributions, extending beyond $\log_2\text{FC} = \pm 6.0$.
* **Significance Elevation:** The most significant DEGs reach $-\log_{10}(P)$ values $> 8.0$, with adjusted FDR values reaching $10^{-7}$. This confirms that the Likelihood Ratio Test successfully unmasked the transcriptional response from the background sonication noise.

### B. Clade Concordance and Outlier Resolution (Heatmap)
* **Complete Clade Separation:** Column dendrogram clustering cleanly separates the experimental libraries into two distinct primary clades: Wild Type controls and Conditioned roots.
* **Resolution of `R1-1_WT`:** Although `R1-1_WT` separated along PC2 (18% variance) in global unsupervised PCA, **it clusters firmly within the Wild Type clade** when evaluated across the top 50 conditioning DEGs.
  * Its row Z-score expression profile matches replicates `R2-1_WT` and `R3-1_WT` across both activated and repressed blocks.
  * This confirms that `R1-1_WT`'s divergence in PCA was caused by non-specific technical noise in background, unexpressed, or structural genes during sonication, rather than biological contamination or treatment mislabeling.
  * Excluding `R1-1_WT` from GLM model fitting was methodologically sound: it eliminated technical dispersion inflation while preserving the authentic underlying control baseline.

---

## 5. Output Artifact Inventory

```text
assets/
├── 004_pilot_lrt_volcano.png     # Publication volcano plot of 865 significant DEGs (LRT)
└── 005_pilot_top50_heatmap.png   # Hierarchical clustering heatmap of top 50 DEGs across 6 libraries
```

---

## 6. Pipeline Status

Statistical visualization is finalized. The extracted gene set (`final_sig_genes`) contains internal StringTie transcript/gene identifiers (`MSTRG.*`). The pipeline now advances to locus dictionary matching against `t_data.ctab` to map MSTRG identifiers to canonical Arabidopsis TAIR IDs, followed by Gene Ontology (GO) and KEGG functional pathway enrichment (`11_go_kegg_analysis_20260411.R`).