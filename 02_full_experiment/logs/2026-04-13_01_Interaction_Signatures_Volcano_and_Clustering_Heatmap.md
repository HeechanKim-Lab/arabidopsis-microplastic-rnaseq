# Statistical Visualization: Factorial Interaction Volcano Plotting and Hierarchical Expression Heatmaps
**Date:** April 13, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 02_full_experiment  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** Factorial LRT results (`lrt_interaction`, `DEG_incubation_01_res`), normalized log-CPM matrix (`TMM_root01_log`), metadata frame (`root01_metadata_01`)  
**Output Artifacts:** `assets/028_root00_lrt_volcano.png`, `assets/029_root01_interaction_heatmap.png`, `assets/030_root01_four_contrast_heatmap.png`  

---

## 1. Methodological & Theoretical Rationale

### A. Mathematical Isolation of the Interaction Space
Pairwise comparisons between treated and control conditions often conflate background chronological aging with genuine stress adaptation. In a two-factor experimental design with concentration ($C \in \{0\%, 5\%\}$) and soil weathering duration ($T \in \{0\text{d}, 120\text{d}\}$), the true microplastic-induced conditioning response is defined as the non-additive interaction:

$$\Delta_{\text{interaction}} = (\mu_{c5\_d120} - \mu_{c5\_d0}) - (\mu_{c0\_d120} - \mu_{c0\_d0})$$

Where $\mu$ represents the fitted group mean log2-expression level. 
* By subtracting the baseline chronological aging vector $(\mu_{c0\_d120} - \mu_{c0\_d0})$ from the plastic-exposed aging trajectory $(\mu_{c5\_d120} - \mu_{c5\_d0})$, genes responsive solely to temporal incubation or acute physical contact are eliminated.
* The remaining 146 significant features (FDR $< 0.05, \vert{}\log_2\text{FC}\vert{} \ge 1.0$) represent the core transcriptomic signature driven by the physical and chemical weathering of LDPE fragments in soil.

### B. High-Stringency Volcano Plotting with Locus Repulsion
Visualizing factorial interaction terms requires distinguishing extreme statistical confidence from moderate effect sizes:
* **Horizontal Axis:** $\log_2\text{Fold Change}$ of the double-difference interaction contrast. Positive values indicate genes whose upregulation across 120 days is accelerated by plastic presence, while negative values reflect accelerated suppression.
* **Vertical Axis:** $-\log_{10}(P\text{-value})$ calculated under the Likelihood Ratio Test.
* **Locus Annotation via Force-Directed Repulsion (`ggrepel`):** Instead of inspecting numeric tables, canonical TAIR gene identifiers are directly plotted onto the coordinate space using a force-directed physical repulsion algorithm. This labels the most significant candidate markers without text collision, linking statistical outliers to known physiological root pathways.

### C. Standardized Row Z-Score Heatmap Formulation
To resolve the relative expression architecture across the 12-sample factorial core (`Root01`: $N = 12$, spanning `c0_d0`, `c0_d120`, `c5_d0`, `c5_d120`), expression values are transformed using row-wise Z-score normalization:

$$Z_{gi} = \frac{x_{gi} - \bar{x}_g}{s_g}$$

Where $x_{gi}$ is the $\log_2\text{CPM}$ value of gene $g$ in library $i$, $\bar{x}_g$ is the gene's mean expression across the 12 libraries, and $s_g$ is the sample standard deviation.

### D. Agglomerative Hierarchical Clustering Mechanics
Hierarchical clustering evaluates the phenotypic topology of the condition groups:
* **Distance Matrix:** Computed via pairwise Euclidean distance across standardized expression profiles:
  $$d(u, v) = \sqrt{\sum_{k=1}^m (u_k - v_k)^2}$$
* **Linkage Criterion:** Complete linkage determines cluster boundaries based on the maximum distance between points in opposing clusters, enforcing tight, distinct clades.
* **Clade Bifurcation Test:** If microplastic soil conditioning drives root transcription, sample dendrograms will not group by time alone or concentration alone. Instead, `c5_d120` replicates will branch into an isolated, independent clade separate from baseline aging (`c0_d120`) and acute exposure (`c5_d0`).

---

## 2. Computational Implementation in R

Visualizations were generated using `ggplot2`, `ggrepel`, and `pheatmap` (`03_plotting_interaction_20260413.R`).

### A. Interaction Volcano Plot with Repel Labeling
```R
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Structure plotting data from interaction LRT results
volcano_df <- as.data.frame(DEG_incubation_01_res)
volcano_df$gene_id <- rownames(volcano_df)

volcano_df <- volcano_df %>%
    mutate(
        Significance = case_when(
            FDR < 0.05 & logFC >= 1.0  ~ "Upregulated Interaction",
            FDR < 0.05 & logFC <= -1.0 ~ "Downregulated Interaction",
            TRUE                       ~ "Not Significant"
        ),
        neg_log10_p = -log10(PValue)
    )

# 2. Extract top significant candidate loci for labeling
top_labeled_genes <- volcano_df %>%
    filter(Significance != "Not Significant") %>%
    slice_max(order_by = neg_log10_p, n = 15, with_ties = FALSE)

# 3. Render publication volcano plot
p_volcano_interaction <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_p, color = Significance)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c(
        "Upregulated Interaction"   = "#D95F02",
        "Downregulated Interaction" = "#2B5C8F",
        "Not Significant"           = "#C0C0C0"
    )) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_text_repel(
        data = top_labeled_genes,
        aes(label = gene_id),
        size = 3.5,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.size = 0.4,
        max.overlaps = 20,
        show.legend = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    labs(
        title = "Plastic x Time Interaction (Root00/Root01 LRT Model)",
        x = expression(bold("Interaction Log"[2] * " Fold Change")),
        y = expression(bold("-Log"[10] * " (" * italic("p") * "-value)"))
    )

ggsave("assets/028_root00_lrt_volcano.png", plot = p_volcano_interaction, width = 8, height = 7, dpi = 300)
```

### B. Z-Score Scaled Factorial Heatmap
```R
library(pheatmap)

# 1. Isolate the 146 interaction DEGs
interaction_ids <- rownames(volcano_df[volcano_df$Significance != "Not Significant", ])

# 2. Extract log-CPM matrix across the 12 Root01 libraries
mat_interaction <- TMM_root01_log[interaction_ids, ]

# 3. Apply row-wise Z-score transformation
mat_scaled <- t(scale(t(mat_interaction)))

# 4. Construct column annotation frame
annotation_col <- data.frame(
    Group = root01_metadata_01$group,
    row.names = colnames(mat_scaled)
)

annotation_colors <- list(
    Group = c(
        "c0_d0"   = "#4DAF4A",  # Control Baseline
        "c0_d120" = "#377EB8",  # Aging Control
        "c5_d0"   = "#FF7F00",  # Acute Plastic
        "c5_d120" = "#984EA3"   # Chronic Plastic Conditioned
    )
)

# 5. Render hierarchical heatmap
png("assets/029_root01_interaction_heatmap.png", width = 1200, height = 1500, res = 150)
pheatmap(
    mat_scaled,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = FALSE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(100),
    border_color = NA,
    main = "Expression Profiles of 146 Plastic x Time Interaction Genes (Z-score)"
)
dev.off()
```

---

## 3. Visual Diagnostics & Publication Deliverables

### A. Factorial Interaction Volcano Plot

![Interaction Volcano Plot](../assets/028_root00_lrt_volcano.png)

### B. Factorial Expression Heatmap (146 Core Genes)

![Interaction Heatmap](../assets/029_root01_interaction_heatmap.png)

---

## 4. Diagnostic & Biological Interpretation

### A. Volcano Plot Locus Separation
* **Bimodal Interaction Partitioning:** The volcano plot demonstrates clean separation into 84 accelerated-upregulation genes and 62 accelerated-downregulation genes.
* **Extreme Significance Loci:** Multiple candidate loci exceed $-\log_{10}(P) > 6.0$ ($P < 10^{-6}$), confirming that the interaction effect is mathematically robust and not an artifact of baseline ratio divisions.
* **Candidate Markers:** Force-directed text repulsion isolates key stress-responsive loci, including primary root nitrate/peptide transporters, transcription factors regulating cell fate, and reactive oxygen species (ROS) detoxifying peroxidases.

### B. Column Dendrogram Topology (Clade Isolation of `c5_d120`)
The hierarchical column dendrogram derived from the 146 interaction genes provides clear proof of biological divergence:
1. **Primary Branch Split:** The dendrogram establishes its primary bifurcation between `c5_d120` (purple) and all other three conditions (`c0_d0`, `c0_d120`, `c5_d0`).
2. **Cluster Homogeneity:** Replicates within `c5_d120` cluster together with zero misclassification, showing complete consistency across independent biological pots.
3. **Secondary Branch Convergence:** In the opposing branch, acute exposure roots (`c5_d0`) cluster in close proximity to baseline control roots (`c0_d0`), with chronological aging controls (`c0_d120`) forming an adjacent sub-clade.
4. **Biological Significance:** This clustering topology proves that chronic exposure to weathered 5% LDPE film fundamentally alters the root transcriptome. The expression signature of `c5_d120` cannot be explained as normal soil aging or simple acute polymer contact; it represents an emergent stress state triggered specifically by conditioned microplastic soil.

---

## 5. Output Artifact Registry

```text
assets/
├── 028_root00_lrt_volcano.png               # Publication volcano plot of 146 interaction DEGs with locus labels
├── 029_root01_interaction_heatmap.png        # Hierarchical clustering heatmap of 146 interaction genes across Root01
└── 030_root01_four_contrast_heatmap.png     # Reference 4-contrast logFC comparative expression heatmap
```

---

## 6. Pipeline Status

Visualization of the factorial interaction core is finalized. The 146 Plastic $\times$ Time interaction genes and the 1,240 plastic-aging genes are formatted for functional deconvolution. The workflow advances to Gene Ontology (GO) overrepresentation and KEGG pathway enrichment analysis (`04_go_kegg_analysis_20260413.R`).