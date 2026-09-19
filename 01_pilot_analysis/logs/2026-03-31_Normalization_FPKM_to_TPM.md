# Downstream Processing: Transcript-Level Normalization and FPKM-to-TPM Linear Transformation
**Date:** March 31, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `results/ballgown/*/t_data.ctab`  
**Output Artifacts:** `arabidopsis_transcript_tpm_matrix.csv`  

---

## 1. Methodological & Theoretical Rationale

### A. Technical Biases in RNA-seq Quantification
Direct comparisons of raw sequence fragment counts across genes or biological replicates are invalid due to two systemic technical confounders:
1. **Feature Length Bias:** Longer transcripts produce more sequencing fragments per unit molar concentration than shorter transcripts at identical expression levels.
2. **Sequencing Depth (Library Size) Bias:** Libraries sequenced to higher depths yield proportionally higher raw counts across all features, independent of biological upregulation.

### B. Mathematical Formulations: RPK, FPKM, and TPM

#### Reads Per Kilobase (RPK)
RPK normalizes raw counts strictly for transcript length:

$$RPK_i = \frac{C_i}{L_i / 10^3}$$

Where $C_i$ is the raw read (or fragment) count for transcript $i$, and $L_i$ is the transcript length in base pairs. RPK serves as the unscaled intermediate rate for both FPKM and TPM.

#### Fragments Per Kilobase Million (FPKM)
FPKM normalizes for library depth *after* normalizing for length:

$$FPKM_i = \frac{RPK_i}{\sum_{j} C_j / 10^6} = \frac{C_i}{(L_i / 10^3) \cdot (\sum_{j} C_j / 10^6)}$$

*The Compositional Inconsistency:* The denominator scales against total read counts ($\sum C$). Because the distribution of transcript lengths varies across samples, the sum of FPKM values ($\sum FPKM$) is non-constant across libraries. Consequently, an identical FPKM value in two samples does not represent an identical proportion of the cellular transcript pool, distorting cross-sample comparisons.

#### Transcripts Per Million (TPM)
TPM enforces length normalization *prior* to scaling against total library abundance:

$$TPM_i = \frac{RPK_i}{\sum_{j} RPK_j} \times 10^6$$

*The Proportional Invariance:* Because the scaling denominator is the sum of length-normalized rates ($\sum RPK$), the sum of TPM values across all features in every sample is constrained to an invariant constant:

$$\sum_{i=1}^{M} TPM_i = 1,000,000$$

TPM constitutes a true relative molar concentration metric. A TPM value of $10$ indicates that for every one million transcripts present in the sequenced pool, exactly 10 originate from transcript $i$, ensuring cross-replicate proportionality.

### C. Direct Linear Transformation of StringTie FPKM to TPM
`StringTie` outputs abundance estimates in FPKM within `t_data.ctab`. Re-parsing raw alignment BAMs to compute RPK is computationally redundant because FPKM is related to RPK by a sample-specific scalar constant:

$$k_{\text{sample}} = \frac{10^6}{\sum_j C_j} \implies FPKM_i = RPK_i \cdot k_{\text{sample}}$$

Substituting this identity into the TPM definition:

$$TPM_i = \frac{\frac{FPKM_i}{k_{\text{sample}}}}{\sum_j \left(\frac{FPKM_j}{k_{\text{sample}}}\right)} \times 10^6 = \left( \frac{FPKM_i}{\sum_j FPKM_{\text{sample}}} \right) \times 10^6$$

Dividing each transcript's FPKM by the sum of all FPKM values within that specific library and multiplying by $10^6$ produces exact TPM values in linear time ($O(N)$).

---

## 2. Computational Implementation in R

The pipeline was implemented in R (`05_fpkm_to_tpm_20260331.R`) to parse Ballgown transcript tables, apply the linear transformation, and assemble a consolidated cross-replicate expression matrix.

### A. Modular Transformation Function
The function `get_tpm_safe` verifies file paths, reads the tab-delimited `t_data.ctab` table, and applies the TPM transformation:

```R
get_tpm_safe <- function(s) {
    file_path <- paste0("ballgown/", s, "/t_data.ctab")
    
    if (!file.exists(file_path)) {
        stop(paste("Error: File not found:", file_path))
    }
    
    df <- read.table(file_path, header = TRUE, sep = "\t")
    
    # Apply linear FPKM-to-TPM transformation
    tpm_vals <- (df$FPKM / sum(df$FPKM)) * 1e6
    
    # Use unique transcript ID (t_name) as primary key; preserve gene_id context
    res <- data.frame(t_name = df$t_name, gene_id = df$gene_id, TPM = tpm_vals)
    colnames(res)[3] <- s
    return(res)
}
```

### B. Vectorized Batch Execution and Matrix Consolidation
Processing was vectorized across all 6 experimental libraries using `lapply`, followed by recursive 1-to-1 table merges:

```R
samples <- c("R1-1_WT", "R2-1_WT", "R3-1_WT", "R1-9_Cond", "R2-9_Cond", "R3-9_Cond")

# Vectorized extraction
tpm_list <- lapply(samples, get_tpm_safe)

# Linear-time recursive merge via unique transcript identifiers
final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by = c("t_name", "gene_id")), tpm_list)

# Export consolidated matrix
write.csv(final_tpm_matrix, "arabidopsis_transcript_tpm_matrix.csv", row.names = FALSE)
```

---

## 3. Structural Validation & Verification

### A. Key Invariance and Relational Integrity
* **Primary Key Collision Prevention:** `gene_id` values in `StringTie` contain multiple transcript isoforms (`t_name`). Merging tables on `gene_id` alone creates many-to-many row expansion errors. Binding composite keys `c("t_name", "gene_id")` ensures strict 1-to-1 relational joins without row duplication.
* **Molar Sum Verification:** Column sums across each sample in `final_tpm_matrix` evaluate to $1.0 \times 10^6 \pm 10^{-5}$ (within floating-point rounding margins), confirming scaling parity across all six biological replicates.

### B. Metric Application Boundary
* **Normalized TPM Matrix (`arabidopsis_transcript_tpm_matrix.csv`):** Reserved for exploratory data analysis (EDA), transcript-level proportion tracking, and comparative expression visualizations.
* **Raw Integer Count Matrix (`gene_count_matrix.csv`):** Preserved without transformation as input for `edgeR` differential expression modeling, where library size variance is parameterized through empirical negative binomial dispersion estimation.

---

## 4. Output Artifact Inventory

```text
results/
└── matrices/
    └── arabidopsis_transcript_tpm_matrix.csv   # Consolidated transcript TPM matrix (52,665 features x 6 samples)
```

---

## 5. Pipeline Status

Transcript-level TPM normalization is complete. Upstream Tuxedo quantification is finalized. The workflow advances to R environment ingestion, metadata alignment verification, and edgeR DGEList initialization (`06_data_loading_20260406.R`).