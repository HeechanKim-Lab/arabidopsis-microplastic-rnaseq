# Loading the Count Matrix and Sample Table into R

**Date:** 2026-04-06  
**Stage:** 01_pilot_analysis  
**Script:** `modules/06_data_ingestion.R`  
**Input:** `gene_count_matrix.csv`, `phenodata.csv`  
**Output:** R objects `counts` and `phenodata`

---

## 1. Points checked

* `read.csv()` changes hyphens in column names to dots by default (`R1-1_WT` becomes `R1.1_WT`). `check.names = FALSE` keeps the original sample names.
* edgeR matches samples to the design by position, so the columns of the count matrix must be in the same order as the rows of the sample table. `prepDE.py3` writes the columns in alphabetical order, which differs from the order in `phenodata.csv`.

---

## 2. Code (R)

```R
counts    <- read.csv("gene_count_matrix.csv", row.names = 1, check.names = FALSE)
phenodata <- read.csv("phenodata.csv", row.names = 1, check.names = FALSE)

print(colnames(counts))
# "R1-1_WT" "R1-9_Cond" "R2-1_WT" "R2-9_Cond" "R3-1_WT" "R3-9_Cond"
print(rownames(phenodata))
# "R1-1_WT" "R2-1_WT" "R3-1_WT" "R1-9_Cond" "R2-9_Cond" "R3-9_Cond"

# Reorder the count columns to match the sample table
counts <- counts[, rownames(phenodata)]

all(colnames(counts) == rownames(phenodata))   # TRUE
```

---

## 3. Result

| Position | Columns as read | Sample table rows | Columns after reordering | Group |
| :--- | :--- | :--- | :--- | :--- |
| 1 | R1-1_WT | R1-1_WT | R1-1_WT | WT |
| 2 | R1-9_Cond | R2-1_WT | R2-1_WT | WT |
| 3 | R2-1_WT | R3-1_WT | R3-1_WT | WT |
| 4 | R2-9_Cond | R1-9_Cond | R1-9_Cond | Cond |
| 5 | R3-1_WT | R2-9_Cond | R2-9_Cond | Cond |
| 6 | R3-9_Cond | R3-9_Cond | R3-9_Cond | Cond |

The matrix has 32,540 genes and 6 libraries (3 WT, 3 Cond).

Next step: filtering and TMM normalisation in edgeR (`modules/07_filtering_and_tmm.R`).
