# edgeR: DGEList, Filtering and TMM Normalisation

**Date:** 2026-04-07  
**Stage:** 01_pilot_analysis  
**Script:** `modules/07_filtering_and_tmm.R`  
**Input:** `counts`, `phenodata`  
**Output:** filtered, TMM-normalised `DGEList` (`dge`)

---

## 1. Method

* A `DGEList` was built from the raw gene counts with `group = factor(phenodata$condition)`. The factor levels are `Cond` and `WT` (alphabetical order), so `Cond` is the reference level in the models that follow.
* `filterByExpr` keeps genes whose CPM is above a cutoff (equivalent to about 10 counts at the median library size) in at least as many libraries as the smallest group (3 here). Library sizes were recalculated after filtering (`keep.lib.sizes = FALSE`).
* `calcNormFactors` computes TMM scaling factors from a weighted, trimmed mean of log-ratios between libraries (by default 30% of M values and 5% of A values are trimmed).

---

## 2. Code (R)

```R
group <- factor(phenodata$condition)
dge   <- DGEList(counts = counts, group = group)

keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Pre-filtering complete. Genes remaining:", nrow(dge), "out of", nrow(counts), "\n")
# Pre-filtering complete. Genes remaining: 20971 out of 32540

dge <- calcNormFactors(dge)
```

---

## 3. Results

* **Genes:** 20,971 of 32,540 kept (64.4%); 11,569 removed.
* **Total counts per library before filtering** (column sums of `gene_count_matrix.csv`, from the means in the 2026-03-23 summary):

| Library | Group | Total counts (millions) |
| :--- | :--- | ---: |
| R1-1_WT | WT | 59.1 |
| R2-1_WT | WT | 62.6 |
| R3-1_WT | WT | 60.3 |
| R1-9_Cond | Cond | 60.7 |
| R2-9_Cond | Cond | 61.8 |
| R3-9_Cond | Cond | 63.0 |

The TMM factors were stored in `dge$samples$norm.factors` but were not written to these notes.

Next step: log-CPM values and PCA (`modules/08_qc_pca.R`).
