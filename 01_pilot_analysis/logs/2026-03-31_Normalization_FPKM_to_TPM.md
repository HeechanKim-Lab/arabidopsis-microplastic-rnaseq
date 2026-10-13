# Transcript TPM Table from StringTie FPKM

**Date:** 2026-03-31  
**Stage:** 01_pilot_analysis  
**Script:** `modules/05_fpkm_to_tpm.R`  
**Input:** `ballgown/<sample>/t_data.ctab`  
**Output:** `arabidopsis_transcript_tpm_matrix.csv`

---

## 1. Method

StringTie reports FPKM per transcript in `t_data.ctab`. Within one library, TPM is FPKM rescaled to sum to one million:

$$TPM_i = \frac{FPKM_i}{\sum_j FPKM_j} \times 10^6$$

The TPM table is used only for descriptive purposes. Differential expression used the count matrix from `prepDE.py3` in edgeR.

---

## 2. Code (R)

```R
setwd("~/Desktop/at_pr/at_data")
file.exists("ballgown/R1-1_WT/t_data.ctab")   # check the path

samples <- c("R1-1_WT", "R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond", "R3-9_Cond")

get_tpm_safe <- function(s) {
    file_path <- paste0("ballgown/", s, "/t_data.ctab")
    df <- read.table(file_path, header = TRUE, sep = "\t")
    tpm_vals <- (df$FPKM / sum(df$FPKM)) * 1e6
    # t_name is unique per transcript; gene_id is kept to show the parent gene
    res <- data.frame(t_name = df$t_name, gene_id = df$gene_id, TPM = tpm_vals)
    colnames(res)[3] <- s
    return(res)
}

tpm_list <- lapply(samples, get_tpm_safe)

# Join on both keys; joining on gene_id alone would duplicate rows for genes with several transcripts
final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by = c("t_name", "gene_id")), tpm_list)

write.csv(final_tpm_matrix, "arabidopsis_transcript_tpm_matrix.csv", row.names = FALSE)
```

---

## 3. Output

`arabidopsis_transcript_tpm_matrix.csv`: one row per transcript of `at_merged.gtf` (52,665), with `t_name`, `gene_id` and one TPM column per library.

Next step: loading the gene count matrix into edgeR (`modules/06_data_ingestion.R`).
