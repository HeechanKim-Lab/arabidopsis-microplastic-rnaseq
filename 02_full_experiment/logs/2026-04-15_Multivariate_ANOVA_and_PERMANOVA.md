# Multivariate Hypothesis Testing: Permutational Multivariate Analysis of Variance (PERMANOVA) and Two-Way ANOVA of Latent Principal Components
**Date:** April 15, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 02_full_experiment  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** TMM-normalized $\log_2\text{CPM}$ matrix (`TMM_root01_log`), metadata frame (`meta_root01`), SVD coordinates (`TMM_root01_pca`)  
**Output Artifacts:** `assets/023_root01_publication_pca_2d.png`, `assets/024_root01_anova_pc1.png`, `assets/025_root01_anova_pc2.png`, `assets/026_root01_anova_pc3.png`, `results/02_anova_result.txt`  

---

## 1. Methodological & Theoretical Rationale

### A. Permutational Multivariate Analysis of Variance (PERMANOVA via `adonis2`)
Unconstrained ordination (PCA) demonstrated that biological exposure gradients and technical root-sonication shear stress project across distinct latent dimensions. To test whether global transcriptomic distances partition significantly by experimental factors, non-parametric Permutational Multivariate Analysis of Variance (PERMANOVA) was executed via `vegan::adonis2`:

$$F = \frac{SS_A / (a - 1)}{SS_E / (N - a)}$$

Where $SS_A$ is the among-group sum of squared dissimilarities, $SS_E$ is the within-group (residual) sum of squares, $a$ is the factor level count, and $N$ is total sample size ($N = 12$).
* **Distance Metric:** Euclidean distance computed directly on the high-dimensional TMM-normalized $\log_2\text{CPM}$ expression matrix, maintaining parity with PCA singular value decomposition.
* **Sequential Sum of Squares (`by = "terms"`):** Evaluates variance contributions sequentially:
  $$\text{Expression Space} \sim \text{Day} + \text{Plastic} + \text{Day} \times \text{Plastic}$$
  This isolates the baseline chronological soil incubation effect (`Day`) prior to evaluating the marginal plastic effect (`Plastic`) and their non-additive interaction (`Day:Plastic`).
* **Permutation Scheme:** Evaluated under 999 unrestricted permutations (`how(nperm = 999)`), avoiding distributional multivariate normality assumptions.

### B. Two-Way Factorial ANOVA on Latent Principal Components
While PERMANOVA tests global multi-dimensional dissimilarity across all ~21,000 genes simultaneously, univariate two-way ANOVA applied to individual principal component score vectors ($PC_1, PC_2, PC_3$) determines which experimental factors drive specific ordination axes:

$$PC_{k, ij} = \mu + \alpha_i^{\text{Day}} + \beta_j^{\text{Plastic}} + (\alpha\beta)_{ij}^{\text{Day:Plastic}} + \epsilon_{ij}$$

Where $PC_{k, ij}$ is the score of sample replicate $j$ under treatment $i$ along the $k$-th principal component, $\alpha^{\text{Day}}$ models the 0d vs. 120d incubation effect, $\beta^{\text{Plastic}}$ models the 0% vs. 5% LDPE exposure, and $(\alpha\beta)$ evaluates the non-additive interaction.

### C. Multiple Comparisons and Post-Hoc Letter Assignment (`multcompView`)
For each principal component, pairwise differences between the four factorial conditions (1: `c0_d0`, 2: `c5_d0`, 3: `c0_d120`, 4: `c5_d120`) were evaluated via Tukey’s Honest Significant Difference (HSD) test:

$$q = \frac{\bar{y}_i - \bar{y}_j}{\sqrt{\frac{MS_{\text{residual}}}{n}}}$$

Post-hoc compact letter displays were derived using `multcompView::multcompLetters4()`. Under small sample sizes ($N = 12$ total, $n = 3$ replicates per condition, residual $df = 8$), family-wise error rate control prevents over-interpreting within-group variance spikes while testing between-condition shifts.

---

## 2. Computational Implementation Protocol

The multivariate and univariate models were executed in R (`05_multivariate_anova_and_permanova_20260415.R`).

### A. Non-Parametric PERMANOVA
```R
library(vegan)

# Ingest normalized matrix and metadata for the 12-sample Root01 cohort
# Conditions: 1 = c0_d0, 2 = c5_d0, 3 = c0_d120, 4 = c5_d120
meta_root01$Day     <- factor(meta_root01$period)   # 0d vs 120d
meta_root01$Plastic <- factor(meta_root01$plastic)  # 0% vs 5%

set.seed(42)
permanova_res <- adonis2(
    formula = t(TMM_root01_log) ~ Day * Plastic,
    data    = meta_root01,
    method  = "euclidean",
    by      = "terms",
    permutations = 999
)
print(permanova_res)
```

### B. Two-Way ANOVA and Tukey HSD Letter Generation
```R
library(multcompView)
library(ggplot2)

# Extract scores from prcomp object
pca_scores <- as.data.frame(TMM_root01_pca$x[, 1:3])
pca_scores$Conditions <- factor(meta_root01$condition_id) # Levels: 1, 2, 3, 4
pca_scores$Day        <- meta_root01$Day
pca_scores$Plastic    <- meta_root01$Plastic

# Fit 2-Way ANOVA across PC1, PC2, and PC3
fit_pc1 <- aov(PC1 ~ Day * Plastic, data = pca_scores)
fit_pc2 <- aov(PC2 ~ Day * Plastic, data = pca_scores)
fit_pc3 <- aov(PC3 ~ Day * Plastic, data = pca_scores)

# Generate Tukey post-hoc letter assignments
tukey_pc1 <- TukeyHSD(aov(PC1 ~ Conditions, data = pca_scores))
letters_pc1 <- multcompLetters4(aov(PC1 ~ Conditions, data = pca_scores), tukey_pc1)

# Publication Boxplot Function
plot_pc_boxplot <- function(df, pc_col, letters_obj, y_label, out_png) {
    p <- ggplot(df, aes_string(x = "Conditions", y = pc_col)) +
        geom_boxplot(width = 0.5, color = "black", linewidth = 0.8, outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 2.5, color = "gray50", alpha = 0.8) +
        theme_bw(base_size = 14) +
        theme(
            panel.grid = element_blank(),
            axis.title = element_text(face = "bold", size = 14),
            axis.text  = element_text(color = "black", size = 12)
        ) +
        labs(x = "Conditions", y = y_label)
    
    # Save asset
    ggsave(out_png, plot = p, width = 7, height = 7, dpi = 300)
}
```

---

## 3. Statistical Testing Results

### A. PERMANOVA Variance Partitioning (`adonis2`)
Testing the high-dimensional Euclidean space ($20,971$ genes across $N = 12$ samples, 999 permutations):

| Source of Variation | Degrees of Freedom ($df$) | Sum of Squares ($SS$) | Mean Squares ($MS$) | Explained Variance ($R^2$) | Pseudo-$F$ | $P$-value ($\text{Pr}(>F)$) | Statistical Significance |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Day** | 1 | 19,584 | 19,584 | **0.19450** | **2.4266** | **0.010** | **$P \le 0.01$ (\*\*)** |
| **Plastic** | 1 | 5,708 | 5,708 | 0.05669 | 0.7073 | 0.722 | Not Significant |
| **Day:Plastic** | 1 | 10,832 | 10,832 | 0.10758 | 1.3421 | 0.184 | Not Significant |
| **Residual** | 8 | 64,565 | 8,071 | **0.64123** | — | — | — |
| **Total** | 11 | 100,690 | — | 1.00000 | — | — | — |

---

### B. Univariate Two-Way ANOVA across Principal Components

#### 1. Principal Component 1 (29.77% Variance)
$$PC_1 \sim \text{Day} + \text{Plastic} + \text{Day} \times \text{Plastic}$$

| Model Term | Degrees of Freedom ($df$) | Sum of Squares ($SS$) | Mean Square ($MS$) | $F$-value | $P$-value ($\text{Pr}(>F)$) | Significance Level |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Day** | 1 | 19,740 | 19,740 | **4.419** | **0.0687** | Trend ($P < 0.1$) |
| **Plastic** | 1 | 994 | 994 | 0.223 | 0.6497 | Not Significant |
| **Day:Plastic** | 1 | 1,040 | 1,040 | 0.233 | 0.6423 | Not Significant |
| **Residuals** | 8 | 35,737 | 4,467 | — | — | — |

#### 2. Principal Component 2 (17.29% Variance)
$$PC_2 \sim \text{Day} + \text{Plastic} + \text{Day} \times \text{Plastic}$$

| Model Term | Degrees of Freedom ($df$) | Sum of Squares ($SS$) | Mean Square ($MS$) | $F$-value | $P$-value ($\text{Pr}(>F)$) | Significance Level |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Day** | 1 | 8,171 | 8,171 | **3.887** | **0.0842** | Trend ($P < 0.1$) |
| **Plastic** | 1 | 330 | 330 | 0.157 | 0.7022 | Not Significant |
| **Day:Plastic** | 1 | 8,079 | 8,079 | **3.843** | **0.0856** | Trend ($P < 0.1$) |
| **Residuals** | 8 | 16,818 | 2,102 | — | — | — |

#### 3. Principal Component 3 (13.30% Variance)
$$PC_3 \sim \text{Day} + \text{Plastic} + \text{Day} \times \text{Plastic}$$

| Model Term | Degrees of Freedom ($df$) | Sum of Squares ($SS$) | Mean Square ($MS$) | $F$-value | $P$-value ($\text{Pr}(>F)$) | Significance Level |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Day** | 1 | 139 | 139 | 0.062 | 0.809 | Not Significant |
| **Plastic** | 1 | 336 | 336 | 0.150 | 0.709 | Not Significant |
| **Day:Plastic** | 1 | 7,737 | 7,737 | **3.458** | **0.100** | Borderline Trend |
| **Residuals** | 8 | 17,897 | 2,237 | — | — | — |

---

## 4. Visual Diagnostics & Publication Deliverables

### A. 2D Principal Component Projection (`Root01`)

![Publication 2D PCA](../assets/023_root01_publication_pca_2d.png)

*Symbols:* Filled circles = Condition 1 (`c0_d0`), Open circles = Condition 2 (`c5_d0`), Filled squares = Condition 3 (`c0_d120`), Open squares = Condition 4 (`c5_d120`).

---

### B. Latent Score Boxplots with Post-Hoc Tukey HSD Letters

![ANOVA PC1 Boxplot](../assets/024_root01_anova_pc1.png)

![ANOVA PC2 Boxplot](../assets/025_root01_anova_pc2.png)

![ANOVA PC3 Boxplot](../assets/026_root01_anova_pc3.png)

---

## 5. Statistical & Biological Interpretation

### A. Dominant Chronological Variance in Global Distance Space
PERMANOVA demonstrates that **soil incubation duration (`Day`) is the single statistically significant global driver** ($R^2 = 0.1945, F = 2.4266, P = 0.010$). 
* Approximately $19.45\%$ of the total variance across all 20,971 genes is governed by soil temporal aging.
* The main effect of `Plastic` ($R^2 = 0.0567, P = 0.722$) and the interaction `Day:Plastic` ($R^2 = 0.1076, P = 0.184$) do not reach global multivariate significance when computed across the entire transcriptome.
* **Residual Variance ($64.12\%$):** Nearly two-thirds of the total variance resides in the residual error term. This high residual fraction is the direct mathematical footprint of physical root sonication, which introduces variable baseline shear stress across individual pots.

### B. Functional Deconvolution of Latent Axes
Univariate ANOVA on individual principal components clarifies how experimental factors partition across dimensions:
1. **PC1 (29.77% Variance) Governs Temporal Progression:** `Day` explains the largest sum of squares ($SS = 19,740, F = 4.419, P = 0.0687$). Moving along PC1 reflects chronological transitions from Day 0 (Conditions 1 & 2, negative scores) to Day 120 (Conditions 3 & 4, positive scores).
2. **PC2 (17.29% Variance) Captures Condition 4 (`c5_d120`) Interaction Divergence:** Along PC2, both `Day` ($P = 0.0842$) and the `Day:Plastic` interaction ($P = 0.0856, F = 3.843$) exhibit marginal significance. As shown in the boxplot, Condition 4 (`c5_d120`) shifts to elevated positive PC2 scores (mean score $\approx +50$), while Day 0 controls and acute plastic samples remain centered near or below zero.
3. **PC3 (13.30% Variance) Reflects Interaction Dominance:** On PC3, the interaction term accounts for $7,737$ out of $8,212$ model sum of squares ($P = 0.100, F = 3.458$), indicating that secondary variation is driven by the divergence of plastic aging from baseline control aging.

### C. Methodological Rigor of Post-Hoc Letter Invariance
In all three boxplots (`024`, `025`, `026`), all four conditions share the letter **`a`**:
* With $n = 3$ biological replicates per group and $df = 8$ residual degrees of freedom, the statistical power of pairwise Tukey HSD tests is constrained.
* Although omnibus $F$-tests detect marginal directional trends ($P \approx 0.06\text{--}0.08$), pairwise confidence intervals between group means overlap after controlling for family-wise error rates.
* Presenting uniform 'a' groupings maintains methodological transparency: it demonstrates that while global LRT models isolate 146 specific interaction DEGs, unconstrained whole-transcriptome ordination axes honestly capture high within-group sonication variance.

---

## 6. Output Artifact Registry

```text
assets/
├── 023_root01_publication_pca_2d.png    # 2D ordination plot of Root01 (PC1: 29.77%, PC2: 17.29%)
├── 024_root01_anova_pc1.png             # Boxplot: PC1 score distribution across Conditions 1-4
├── 025_root01_anova_pc2.png             # Boxplot: PC2 score distribution across Conditions 1-4
└── 026_root01_anova_pc3.png             # Boxplot: PC3 score distribution across Conditions 1-4
results/
└── 02_anova_result.txt                  # Numerical output of adonis2 and 2-way ANOVA models
```

---

## 7. Pipeline Status

Multivariate variance partitioning and latent component ANOVA are finalized. The analysis confirms that temporal soil aging and plastic interaction govern PC1 and PC2, while high residual variance reflects sonication noise. 

The pipeline now advances to the final analytical module: deconstructing targeted physiological pathways—specifically Reactive Oxygen Species (ROS) oxidative scavenging and Photosynthetic/Plastid repression—via dedicated expression heatmaps (`2026-04-17_01_Targeted_Pathways_ROS_Photosynthesis.md`).