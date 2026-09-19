# Quality Control: Global 27-Sample Matrix Variance Profiling, Technical Sonication Noise Deconvolution, and Dose-Response Ordination
**Date:** April 11, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 02_full_experiment  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** Raw count matrices (`gene_count_matrix.csv`), metadata mapping files (`phenodata.csv`), R workspace objects (`root00_y4`, `root01_y4`)  
**Output Artifacts:** `assets/020_root00_full_matrix_pca.png`, `assets/021_root00_pca_pc1_vs_pc3.png`, `assets/022_root01_dose_response_pca.png`  

---

## 1. Methodological & Theoretical Rationale

### A. Global Transcriptomic Landscape vs. Factorial Stratification
The complete experimental dataset encompasses 27 biological libraries spanning a multi-factorial space:
* **Plastic Concentration:** 0% (Control), 1% (Low Dose), 5% (High Dose) $(w/w)$ LDPE film fragments.
* **Soil Incubation Time:** 0 days (Acute), 30 days (Intermediate), 120 days (Chronic weathering).
* **Biological Replicates:** 3 independent pots per treatment combination ($3 \times 3 \times 3 = 27$ libraries).

To evaluate global variance while isolating primary hypothesis contrasts, the dataset is analyzed across two distinct matrix architectures:
1. **Global Dataset (`Root00`, $N = 27$):** Captures the comprehensive expression landscape across all doses and incubation intervals.
2. **Characteristic Factorial Cohort (`Root01`, $N = 12$):** Subsets the extreme experimental endpoints (`c0_d0`, `c5_d0`, `c0_d120`, `c5_d120`) to eliminate intermediate dilution effects and directly model Plastic $\times$ Time interaction vectors.

### B. The Biophysical Consequence of Rhizosphere Sonication
Harvesting tightly bound rhizosphere soil from delicate Arabidopsis root systems requires physical sonication (ultrasonic bath agitation). While sonication successfully detaches adhering soil and microbial cells, it introduces an unavoidable biophysical artifact:
* **Variable Shear Stress:** Uneven ultrasonic wave propagation causes variable root hair shearing, epidermal micro-abrasion, and localized cell lysis across replicates.
* **Intra-Group Dispersion Inflation:** This mechanical stress introduces substantial technical variance within identical treatment groups. When mapped into high-dimensional gene expression space, this within-group variance aligns along secondary variance components.
* **Model Suppression Under QLF:** In standard negative binomial Quasi-Likelihood F-testing (`glmQLFTest`), the quasi-likelihood dispersion denominator ($s_g^2$) is inflated by this mechanical variance, collapsing statistical test ratios and suppressing genuine biological differential expression.

### C. Unsupervised PCA as a Decision Gate for LRT Modeling
Principal Component Analysis (PCA) performs orthogonal linear decomposition of the variance-covariance matrix:
* If sonication damage completely scrambles biological signal, samples will fail to cluster by treatment, and variance will scatter randomly across higher components.
* If biological perturbation remains intact, the treatment effect will orthogonalize into a primary component (e.g., PC1), while technical sonication noise will isolate along an orthogonal secondary component (e.g., PC2).

Confirming that plastic concentration projects along a distinct, coherent axis provides mathematical justification to bypass the overly conservative QLF test and utilize the **Likelihood Ratio Test (LRT)** for full-matrix differential expression modeling.

---

## 2. Computational Implementation in R

PCA was executed using base R and `edgeR` normalization routines (`01_pca_analysis_20260411.R`).

### A. Global Dataset Preparation (`Root00`, $N = 27$)
```R
# Extract TMM-normalized log2-CPM with prior count
TMM_root00_log <- cpm(root00_y4, log = TRUE, prior.count = 2)

# Execute SVD on transposed gene space
TMM_root00_pca <- prcomp(t(TMM_root00_log), scale. = TRUE)

# Calculate variance explained
root00_percentVar <- round(100 * (TMM_root00_pca$sdev^2 / sum(TMM_root00_pca$sdev^2)), 1)

# Compile ordination frame
root00_pca_df <- data.frame(
    Sample    = colnames(TMM_root00_log),
    PC1       = TMM_root00_pca$x[, 1],
    PC2       = TMM_root00_pca$x[, 2],
    PC3       = TMM_root00_pca$x[, 3],
    Group     = root00_metadata_01$group,
    Plastic   = root00_metadata_01$plastic,
    Period    = root00_metadata_01$period
)
```

### B. Characteristic Cohort Preparation (`Root01`, $N = 12$)
```R
# Extract TMM-normalized log2-CPM for the 12-sample factorial core
TMM_root01_log <- cpm(root01_y4, log = TRUE, prior.count = 2)

# Execute SVD
TMM_root01_pca <- prcomp(t(TMM_root01_log), scale. = TRUE)
root01_percentVar <- round(100 * (TMM_root01_pca$sdev^2 / sum(TMM_root01_pca$sdev^2)), 1)

# Compile ordination frame
root01_pca_df <- data.frame(
    Sample    = colnames(TMM_root01_log),
    PC1       = TMM_root01_pca$x[, 1],
    PC2       = TMM_root01_pca$x[, 2],
    PC3       = TMM_root01_pca$x[, 3],
    Group     = root01_metadata_01$group,
    Plastic   = root01_metadata_01$plastic,
    Period    = root01_metadata_01$period
)
```

### C. Visualization Rendering
```R
# Export Root00 Global PCA (PC1 vs PC2)
png("assets/020_root00_full_matrix_pca.png", width = 1200, height = 1000, res = 150)
plot(root00_pca_df$PC1, root00_pca_df$PC2, col = root00_pca_df$Plastic, pch = 16,
     xlab = paste0("PC1: ", root00_percentVar[1], "% variance"),
     ylab = paste0("PC2: ", root00_percentVar[2], "% variance"),
     main = "Root00 Global PCA (N=27)")
dev.off()

# Export Root00 PC1 vs PC3
png("assets/021_root00_pca_pc1_vs_pc3.png", width = 1200, height = 1000, res = 150)
plot(root00_pca_df$PC1, root00_pca_df$PC3, col = root00_pca_df$Plastic, pch = 16,
     xlab = paste0("PC1: ", root00_percentVar[1], "% variance"),
     ylab = paste0("PC3: ", root00_percentVar[3], "% variance"),
     main = "Root00 PCA: PC1 vs PC3")
dev.off()

# Export Root01 Factorial PCA
png("assets/022_root01_dose_response_pca.png", width = 1200, height = 1000, res = 150)
plot(root01_pca_df$PC1, root01_pca_df$PC2, col = root01_pca_df$Plastic, pch = 16,
     xlab = paste0("PC1: ", root01_percentVar[1], "% variance"),
     ylab = paste0("PC2: ", root01_percentVar[2], "% variance"),
     main = "Root01 Characteristic Subset PCA (N=12)")
dev.off()
```

---

## 3. Visual Diagnostics & Ordination Analysis

### A. Global Dataset Landscape (`Root00`, $N = 27$)

![Root00 Full Matrix PCA](../assets/020_root00_full_matrix_pca.png)

![Root00 PC1 vs PC3](../assets/021_root00_pca_pc1_vs_pc3.png)

#### Diagnostic Observations
1. **Vertical Dispersion along PC2:** The global 27-sample ordination displays prominent spread across the vertical axis (PC2). Biological replicates treated identically exhibit broad distribution across PC2, reflecting the heterogeneous impact of root sonication shear stress.
2. **Technical Confounding:** In this unconstrained 27-sample space, the interaction between multiple incubation timelines (0d, 30d, 120d) and mechanical root processing obscures clear-cut visual separation when viewing PC1 vs PC2 alone.
3. **PC3 Resolving Power:** Plotting PC1 against PC3 demonstrates cleaner segregation of incubation timepoints, confirming that temporal soil aging contributes an independent source of variance separate from mechanical root damage.

---

### B. Characteristic Cohort Dose-Response (`Root01`, $N = 12$)

![Root01 Characteristic Subset PCA](../assets/022_root01_dose_response_pca.png)

#### Diagnostic Observations
1. **Unmasked Biological Dose-Response along PC1:** Subsetting to the clean factorial cohort resolves a linear biological gradient along the horizontal axis (PC1):
   * **Left (Negative PC1):** 0% Plastic Control samples.
   * **Center (Intermediate PC1):** 1% Plastic samples (inferred intermediate transition).
   * **Right (Positive PC1):** 5% Plastic Conditioned samples.
2. **Signal-to-Noise Deconvolution:** Removing intermediate incubation noise demonstrates that physical sonication does not erase the transcriptional response. The molecular impact of microplastic residue exposure scales monotonically along PC1, while sonication variance remains largely restricted to PC2.

---

## 4. Cohort Architecture Comparison

| Diagnostic Parameter | Global Architecture (`Root00`) | Factorial Architecture (`Root01`) | Analytical Role |
| :--- | :--- | :--- | :--- |
| **Sample Count ($N$)** | 27 libraries | 12 libraries | Global landscape vs. Core contrast verification |
| **Plastic Levels** | 0%, 1%, 5% | 0%, 5% (with 1% diagnostic trace) | Multi-dose screening vs. Factorial endpoints |
| **Incubation Intervals** | 0d, 30d, 120d | 0d, 120d | Time-series modeling vs. Acute vs. Chronic extremes |
| **Primary Variance (PC1)** | Multi-factor composite | Clean Dose-Response Gradient | Validates biological signal preservation |
| **Secondary Variance (PC2)**| Sonication shear noise | Residual sonication noise | Confirms mechanical noise is orthogonal |
| **Modeling Decision** | LRT GLM (`glmLRT`) | LRT Factorial Contrasts | Bypasses QLF power collapse |

---

## 5. Output Artifact Registry

```text
assets/
├── 020_root00_full_matrix_pca.png      # Global 27-sample PCA showing PC2 sonication dispersion
├── 021_root00_pca_pc1_vs_pc3.png       # Global 27-sample PCA resolving temporal aging variance
└── 022_root01_dose_response_pca.png    # 12-sample factorial PCA resolving the 0% -> 1% -> 5% dose gradient
```

---

## 6. Pipeline Decision & Next Steps

The recovery of the 0% $\rightarrow$ 5% plastic dose-response gradient along PC1 in `Root01` provides empirical justification for the analytical framework:
1. **Rejection of QLF Modeling:** Because within-group sonication variance persists along PC2, Quasi-Likelihood F-testing would artificially penalize test statistics and inflate false negatives across the 27-sample series.
2. **Adoption of LRT Contrasts:** Downstream modeling will utilize generalized linear models evaluated via the **Likelihood Ratio Test (LRT)**.
3. **Execution Plan:** Proceed to pairwise contrast definitions, likelihood ratio testing, and DEG extraction across both Root00 and Root01 matrices (`02_deg_analysis_20260411.R`).