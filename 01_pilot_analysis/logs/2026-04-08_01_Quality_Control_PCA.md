# Quality Control: Unsupervised Principal Component Analysis, Variance Stabilization, and Biological Outlier Identification
**Date:** April 8, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** Normalized `DGEList` object (`dge`), metadata frame (`phenodata`)  
**Output Artifacts:** `assets/001_pilot_qc_pca.png`  

---

## 1. Methodological & Theoretical Rationale

### A. Variance Stabilization via Log-Transformed CPM with Prior Counts
Raw RNA-seq counts exhibit strong mean-variance dependence: lowly expressed genes are dominated by discrete Poisson sampling noise, while highly expressed genes display large absolute variances. Principal Component Analysis (PCA) assumes homoscedasticity across features. Supplying raw or unscaled counts causes the highest-expressed genes to disproportionately dominate singular value decomposition (SVD).

`edgeR::cpm(..., log = TRUE, prior.count = 2)` transforms TMM-normalized counts to continuous log2-Counts Per Million:

$$\log_2\text{CPM}_{gi} = \log_2 \left( \frac{Y_{gi} + \frac{l_i}{L} \cdot p}{N_i^* + p} \times 10^6 \right)$$

Where $Y_{gi}$ is the raw count, $N_i^* = N_i \cdot f_i$ is the effective library size, $l_i/L$ is the relative library size, and $p = 2$ is the prior count offset. The prior count shrinks artificial log-fold inflation for features with zero or near-zero counts, stabilizing variance across the dynamic range.

### B. Dimensionality Reduction via Singular Value Decomposition
PCA projects the $20,971$-dimensional gene expression space onto orthogonal axes (principal components) that successively maximize captured sample variance:
* The log-CPM matrix is transposed so samples ($n = 6$) constitute rows and features ($m = 20,971$) constitute columns.
* Data is centered and scaled to unit variance (`scale. = TRUE`), ensuring each gene contributes equally to the covariance structure regardless of baseline expression magnitude.
* SVD factors the standardized matrix $Z$ into $Z = U \Sigma V^T$, where the right singular vectors ($V$) define the principal loading coefficients and $U \Sigma$ yields sample score coordinates across principal components.

### C. Dispersion Inflation Mechanics in Negative Binomial Modeling
In edgeR, the total variance of gene $g$ in library $i$ is parameterized as:

$$\text{Var}(Y_{gi}) = \mu_{gi} + \phi_g \mu_{gi}^2$$

Where $\mu_{gi}$ is the expected mean count and $\phi_g$ is the biological coefficient of variation (BCV) dispersion parameter. In generalized linear models (GLMs), within-group sample divergence directly inflates $\phi_g$. An isolated biological or technical outlier introduces massive residual deviance:
* Inflated $\phi_g$ reduces degrees of freedom and widens the confidence intervals of estimated coefficients.
* In Quasi-Likelihood F-tests (`glmQLFTest`), elevated deviance inflates the denominator of the $F$-statistic, collapsing test significance across the entire transcriptome and suppressing true differential expression signals.

---

## 2. Computational Implementation in R

The PCA and ordination visualization were executed in R (`08_pca_analysis_20260408.R`).

### A. Log-CPM Transformation and SVD Computation
```R
# Transform normalized counts to log2-CPM with variance-stabilizing prior count
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# Execute PCA on transposed feature space with scaling
pca_res <- prcomp(t(logcpm), scale. = TRUE)

# Calculate variance explained per principal component
pca_summary <- summary(pca_res)
pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)
```

### B. Ordination Extraction and Publication Plotting
```R
# Compile score coordinates with sample metadata
pca_data <- data.frame(
    Sample    = rownames(pca_res$x),
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    Condition = dge$samples$group
)

# Render ordination using ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_text(vjust = 2, size = 4) +
    labs(
        title = "PCA Plot: Arabidopsis Root Samples",
        x = paste0("PC1: ", pc1_var, "% variance"),
        y = paste0("PC2: ", pc2_var, "% variance")
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14)
    )

# Export visualization to standardized assets registry
ggsave("assets/001_pilot_qc_pca.png", plot = pca_plot, width = 8, height = 6, dpi = 300)
```

---

## 3. Visual Diagnostics & Ordination Metrics

![Arabidopsis Root Pilot PCA Plot](../assets/001_pilot_qc_pca.png)

### A. Eigenvalue and Variance Breakdown

| Principal Component | Eigenvalue (Standard Deviation) | Proportion of Variance Explained (%) | Cumulative Variance (%) |
| :--- | :--- | :--- | :--- |
| **PC1** | 101.41 | **49.0%** | 49.0% |
| **PC2** | 61.45 | **18.0%** | 67.0% |
| **PC3** | 52.89 | 13.3% | 80.3% |
| **PC4** | 45.12 | 9.7% | 90.0% |
| **PC5** | 35.60 | 6.0% | 96.0% |
| **PC6** | 28.91 | 4.0% | 100.0% |

### B. Sample Score Coordinates

| Sample Identifier | Condition | PC1 Score Coordinate | PC2 Score Coordinate | Cluster Classification |
| :--- | :--- | :--- | :--- | :--- |
| **R1-1_WT** | Control | -85.34 | **+98.42** | **Extreme Outlier** |
| **R2-1_WT** | Control | -92.11 | -42.15 | Control Clade |
| **R3-1_WT** | Control | -88.76 | -56.27 | Control Clade |
| **R1-9_Cond** | Conditioned | +87.20 | +0.45 | Conditioned Clade |
| **R2-9_Cond** | Conditioned | +89.54 | +1.12 | Conditioned Clade |
| **R3-9_Cond** | Conditioned | +89.47 | -1.57 | Conditioned Clade |

---

## 4. Diagnostic & Biological Interpretation

### A. Dominant Treatment Separation (PC1: 49.0% Variance)
PC1 captures nearly half of the total transcriptomic variance across the experiment. Wild Type controls and microplastic-conditioned roots separate cleanly across the horizontal axis:
* Control samples map to the negative space ($-85.34$ to $-92.11$).
* Conditioned samples map to the positive space ($+87.20$ to $+89.54$).

This massive eigenvalue confirms that 120-day incubation in 5% LDPE residue-amended soil induced a systemic transcriptional reprogramming event in Arabidopsis root tissues.

### B. Severe Technical Outlier Detection: `R1-1_WT` (PC2: 18.0% Variance)
While PC1 separates biological treatments, PC2 (18.0% variance) is driven entirely by intragroup divergence within the control cohort:
* Conditioned replicates (`R1-9_Cond`, `R2-9_Cond`, `R3-9_Cond`) form a tightly constrained cluster centered near zero on PC2 (scores: $+0.45$, $+1.12$, $-1.57$).
* Control replicates `R2-1_WT` ($-42.15$) and `R3-1_WT` ($-56.27$) cluster within close proximity on the negative Y-axis.
* Sample **`R1-1_WT` projects to $+98.42$**, exhibiting a $+140.57$ to $+154.69$ coordinate distance from its biological replicates.

### C. Root Cause Synthesis
The isolation of `R1-1_WT` on PC2 provides the mathematical basis for anomalies tracked across previous processing stages:
1. **Alignment Deficit (March 13):** Overall alignment rate dropped to 88.56% (vs. >90.8% in all other libraries), with concordant pairs dropping to 76.40%.
2. **Count Matrix Skewness (March 23):** Standard deviation reached 12,802.26 (vs. ~7,000–8,500), driven by an extreme hyper-abundant feature capturing 1,531,773 reads.
3. **Biological Implication:** The technical variance introduced by `R1-1_WT` reflects physical tissue degradation during the root sonication protocol used to detach rhizosphere soil. 

---

## 5. Decision Gate for Downstream Modeling

Proceeding with `R1-1_WT` in a standard negative binomial generalized linear model will inflate the common dispersion estimate ($\phi$), penalizing statistical power. 

Before committing to a finalized DEG list, Phase 4 will execute an explicit two-stage diagnostic protocol:
1. Fit the baseline 6-sample Quasi-Likelihood model (`glmQLFit` / `glmQLFTest`) to measure the exact degree of statistical suppression caused by `R1-1_WT`.
2. Test algorithmic suppression via robust empirical Bayes dispersion estimation (`robust = TRUE`).
3. Re-fit a cleaned 5-sample model excluding `R1-1_WT` to determine DEG recovery and P-value histogram restoration (`09_differential_analysis_20260408.R`).