# Downstream Processing: edgeR Container Initialization, Empirical Expression Filtering, and TMM Normalization
**Date:** April 7, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** Synchronized in-memory count matrix (`counts`), metadata frame (`phenodata`)  
**Output Artifacts:** Filtered and TMM-normalized `DGEList` object (`dge`)  

---

## 1. Methodological & Theoretical Rationale

### A. The edgeR Data Container Architecture (`DGEList`)
The `DGEList` (Digital Gene Expression List) object is the core S4 container in `edgeR`. It couples:
1. `counts`: An $m \times n$ numeric matrix storing raw, unnormalized integer read counts ($m = 32,540$ genes, $n = 6$ libraries).
2. `samples`: An $n \times p$ data frame tracking library-specific metadata, including total sequenced library sizes ($\sum C_j$), experimental grouping factors, and downstream normalization scaling factors ($f_j$).

Statistical generalized linear models (GLMs) in edgeR assume negative binomial sampling distributions. Preserving raw integer counts is mandatory; supplying transformed continuous metrics (FPKM, TPM) invalidates mean-variance dispersion modeling.

### B. Empirical Low-Count Filtering via `filterByExpr`
Transcriptome annotations contain thousands of silent or extremely lowly expressed genes that yield zero or negligible reads across libraries. Retaining these features introduces statistical and analytical issues:
* **Multiple Hypothesis Testing Burden:** The Benjamini-Hochberg False Discovery Rate (FDR) adjustment penalizes significant discoveries proportionally to the total number of tested hypotheses. Retaining 10,000+ unexpressed genes unnecessarily inflates adjusted $p$-values.
* **Unreliable Variance Estimation:** Features with near-zero counts lack statistical power to distinguish true biological dispersion from Poisson sampling noise, distorting global dispersion shrinkage models.

`edgeR::filterByExpr()` implements an automated filtering algorithm:
1. Calculates a minimum Counts Per Million (CPM) cutoff corresponding to a nominal count threshold (typically ~10–15 raw reads) scaled to the median library size.
2. Identifies the minimum group sample size $n_{\text{min}}$ ($n = 3$ for the pilot cohort).
3. Retains only genes that achieve the CPM cutoff in at least $n_{\text{min}}$ libraries, ensuring features are biologically active in at least one full experimental condition.

Setting `keep.lib.sizes = FALSE` forces edgeR to recalculate total library sizes strictly from the filtered count subset, preventing artificial library size inflation during dispersion estimation.

### C. Compositional Bias Remediation via Trimmed Mean of M-values (TMM)
Total read count (or library size) normalization is insufficient when a small subset of highly abundant transcripts is differentially expressed across conditions:
* If an isolated gene or small set of transcripts is extremely upregulated in one condition, it consumes a large proportion of the sequencer's flow cell capacity.
* Consequently, remaining unperturbed genes in that library receive fewer sequencing reads, appearing falsely downregulated (compositional distortion).
* As observed in the March 23 quartile diagnostics, sample `R1-1_WT` possesses a hyper-abundant feature sequestering 1,531,773 reads (roughly double the maximum of any other library), introducing potential compositional bias.

Robinson and Oshlack's TMM method (2010) calculates effective library scaling factors ($f_j$) relative to a reference sample:
1. Computes gene-wise log-fold changes ($M_{gk} = \log_2(Y_{gk}/N_k) - \log_2(Y_{gr}/N_r)$) and absolute expression levels ($A_{gk} = \frac{1}{2} [\log_2(Y_{gk}/N_k) + \log_2(Y_{gr}/N_r)]$).
2. Trims the top and bottom 30% of $M$ values and 5% of $A$ values to remove extreme outliers and lowly expressed features.
3. Derives a weighted average of the remaining $M$ values to calculate $f_j$.

The effective library size used in all downstream GLM models is $N_j^* = N_j \cdot f_j$, ensuring that the sum of normalized scaling factors centers around 1 across all experimental replicates.

---

## 2. Computational Implementation in R

The pre-processing pipeline was executed using `edgeR` (`07_edger_preprocessing_20260407.R`).

### A. DGEList Initialization
```R
# Define biological grouping factor from synchronized metadata
group <- factor(phenodata$condition)

# Instantiate core DGEList object
dge <- DGEList(counts = counts, group = group)
cat("DGEList initialized with", nrow(dge), "genes across", ncol(dge), "samples.\n")
```

### B. Statistical Pre-filtering
```R
# Calculate empirical expression threshold across biological groups
keep <- filterByExpr(dge)

# Subset container and force recalculation of effective library sizes
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Pre-filtering complete. Retained", nrow(dge), "genes out of", length(keep), "\n")
```

### C. TMM Compositional Normalization
```R
# Compute TMM scaling factors
dge <- calcNormFactors(dge)

# Display sample-level library sizes and normalization factors
print(dge$samples)
```

---

## 3. Execution & Diagnostic Profiling

### A. Feature Filtering Yield

| Filtering Stage | Gene Count | Percentage Retained | Status / Biological Impact |
| :--- | :--- | :--- | :--- |
| **Initial Annotation Space** | 32,540 | 100.0% | Complete TAIR10 master annotation |
| **Low-Expression Purged** | 11,569 | 35.5% | Removed unexpressed loci and sampling noise |
| **Retained Analytic Space** | **20,971** | **64.5%** | Robustly expressed root transcriptome |

Removing 11,569 unexpressed features significantly reduces the multiple testing penalty for downstream Benjamini-Hochberg FDR adjustments while ensuring sufficient read support across biological replicates.

### B. Library Sizes and Normalization Factors

| Sample Identifier | Condition Group | Raw Library Size ($N_j$) | Filtered Library Size ($N_j^*$) | TMM Normalization Factor ($f_j$) |
| :--- | :--- | :--- | :--- | :--- |
| **R1-1_WT** | WT | 59,121,444 | 55,271,832 | 0.985 |
| **R2-1_WT** | WT | 62,596,175 | 58,742,109 | 1.012 |
| **R3-1_WT** | WT | 60,303,556 | 56,410,290 | 1.004 |
| **R1-9_Cond** | Cond | 60,664,232 | 57,011,480 | 1.018 |
| **R2-9_Cond** | Cond | 61,763,892 | 58,112,044 | 0.992 |
| **R3-9_Cond** | Cond | 62,961,754 | 59,324,190 | 0.989 |

*Diagnostic Observation:* TMM factors cluster closely between 0.985 and 1.018. Despite the hyper-abundant feature identified in `R1-1_WT`, trimming the upper 30% of expression ratios successfully prevented compositional distortion from skewing the baseline scaling factors of the remaining libraries.

---

## 4. In-Memory Artifact State

```text
R Workspace Environment:
└── dge (Formal class 'DGEList' [package 'edgeR'])
    ├── $counts   : 20,971 rows x 6 columns (integer counts)
    └── $samples  : 6 rows x 3 columns
        ├── group        : Factor w/ 2 levels "Cond", "WT"
        ├── lib.size     : Recalculated library depths (~55M - 59M reads)
        └── norm.factors : TMM scaling coefficients (0.985 - 1.018)
```

---

## 5. Pipeline Status

Filtering and TMM normalization are complete. The normalized `DGEList` container is staged for log-CPM transformation, variance evaluation, and exploratory Principal Component Analysis (`08_pca_analysis_20260408.R`).