# Downstream Processing: Count Matrix Ingestion, Syntactic Identifier Preservation, and Metadata Harmonization
**Date:** April 6, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `results/matrices/gene_count_matrix.csv`, `results/matrices/phenodata.csv`  
**Output Artifacts:** Synchronized in-memory R data structures (`counts`, `phenodata`)  

---

## 1. Methodological & Theoretical Rationale

### A. Syntactic Identifier Preservation (`check.names = FALSE`)
Base R's `read.csv()` enforces syntactically valid variable names by default (`check.names = TRUE`). During ingestion, characters outside alphanumeric standards (such as hyphens `-`) are silently coerced into periods (`.`):
* An identifier such as `R1-1_WT` is mutated into `R1.1_WT`.
* This transformation breaks relational integrity across the pipeline: sample IDs defined in the metadata (`phenodata.csv`), upstream FASTQ/BAM filenames, and Ballgown directories fail exact string matching.
* Enforcing `check.names = FALSE` preserves literal hyphenated string keys across data frames.

### B. Strict Index Concordance in Matrix-Metadata Architectures
Linear modeling frameworks (`edgeR`, `DESeq2`) formulate differential expression by coupling a discrete count matrix $Y_{m \times n}$ (genes $m$, libraries $n$) to a model design matrix $X_{n \times p}$ derived from sample metadata:
* **The Index Invariant:** The $j$-th column vector of the count matrix must correspond to the $j$-th row record of the metadata data frame ($1 \le j \le n$).
* **Silent Failure Risk:** Linear algebra operations in R do not validate sample names automatically. If column names in $Y$ do not match the row order of $X$, the generalized linear model (GLM) will map biological treatments to incorrect sample counts without throwing a runtime error. This mislabels condition vectors and invalidates downstream hypothesis testing.

### C. Lexicographical Sorting Discrepancy
`prepDE.py3` scans the `ballgown/` directory and writes count matrix columns in filesystem/lexicographical sorting order:
* Count columns sort by directory string: `R1-1_WT`, `R1-9_Cond`, `R2-1_WT`, `R2-9_Cond`, `R3-1_WT`, `R3-9_Cond`.
* Phenotype metadata (`phenodata.csv`) is grouped categorically by biological factor: Wild Type controls first (`R1-1_WT`, `R2-1_WT`, `R3-1_WT`), followed by Conditioned replicates (`R1-9_Cond`, `R2-9_Cond`, `R3-9_Cond`).

Re-indexing the count matrix against the exact row indices of `phenodata` harmonizes column-row correspondence prior to object instantiation.

---

## 2. Computational Implementation in R

The ingestion and alignment validation were executed in R (`06_data_loading_20260406.R`).

### A. Ingestion with Name Checking Disabled
```R
# Import raw gene-level count matrix
counts <- read.csv("gene_count_matrix.csv", row.names = 1, check.names = FALSE)

# Import sample metadata
phenodata <- read.csv("phenodata.csv", row.names = 1, check.names = FALSE)
```

### B. Column Re-indexing and Structural Harmonization
```R
# Reorder count matrix columns to match metadata row order exactly
counts <- counts[, rownames(phenodata)]
```

### C. Programmatic Alignment Verification
```R
# Assert strict identity across matrix columns and metadata rows
if (all(colnames(counts) == rownames(phenodata))) {
    cat("PASSED: Matrix columns and metadata rows are strictly synchronized.\n")
} else {
    stop("ERROR: Dimension or identifier mismatch detected between counts and metadata.\n")
}
```

---

## 3. Structural Diagnostics & Alignment Audit

### Column and Row Sequence Comparison

| Array Index ($j$) | Raw Count Columns (`gene_count_matrix.csv`) | Metadata Rows (`phenodata.csv`) | Harmonized Count Columns (`counts[, rownames(phenodata)]`) | Experimental Group |
| :--- | :--- | :--- | :--- | :--- |
| **[, 1]** | `R1-1_WT` | `R1-1_WT` | `R1-1_WT` | WT |
| **[, 2]** | `R1-9_Cond` | `R2-1_WT` | `R2-1_WT` | WT |
| **[, 3]** | `R2-1_WT` | `R3-1_WT` | `R3-1_WT` | WT |
| **[, 4]** | `R2-9_Cond` | `R1-9_Cond` | `R1-9_Cond` | Cond |
| **[, 5]** | `R3-1_WT` | `R2-9_Cond` | `R2-9_Cond` | Cond |
| **[, 6]** | `R3-9_Cond` | `R3-9_Cond` | `R3-9_Cond` | Cond |

### Verification Metrics
* **Total Transcribed Loci ($m$):** 32,540 rows.
* **Biological Replicates ($n$):** 6 columns (3 WT controls, 3 Conditioned samples).
* **Missing Value Audit:** Zero `NA` or `NaN` cells across the $32,540 \times 6$ matrix space (`sum(is.na(counts)) == 0`).
* **Identity Audit:** `all(colnames(counts) == rownames(phenodata))` evaluated to `TRUE`.

---

## 4. Pipeline Status

Input matrices are verified, unmutated, and synchronized. The environment is prepared for edgeR `DGEList` data structure initialization, empirical pre-filtering (`filterByExpr`), and Trimmed Mean of M-values (TMM) compositional normalization (`07_edger_preprocessing_20260407.R`).