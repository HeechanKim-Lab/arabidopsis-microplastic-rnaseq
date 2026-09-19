# Upstream Processing: StringTie Abundance Estimation, Ballgown Table Generation, and Raw Count Matrix Extraction
**Date:** March 23, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `results/bam/*.bam`, `results/gtf/at_merged.gtf`  
**Output Artifacts:** `results/ballgown/*/*.ctab`, `results/matrices/gene_count_matrix.csv`, `results/matrices/transcript_count_matrix.csv`, `results/matrices/phenodata.csv`  

---

## 1. Methodological & Theoretical Rationale

### A. Restrictive Abundance Estimation (`-e`)
In downstream differential expression workflows, transcript structures must remain fixed across all biological replicates to ensure that count matrices are directly comparable. 
* Executing `StringTie` with `-e` (expression estimation only) restricts read assignment strictly to the genomic intervals defined in the master reference transcriptome (`at_merged.gtf`).
* De novo transfrag assembly is disabled. Spliced alignments falling outside annotated loci or alternative junctions not captured during the initial merge phase are excluded, preventing library-specific transcript fragmentation from distorting cross-sample counts.

### B. Ballgown Table Emission (`-B`)
The `-B` flag directs `StringTie` to write tabular abundance files compatible with the Ballgown R framework. For each library, five relational tables are generated:
* `e_data.ctab`: Exon-level coordinates, genomic strands, and average sequence coverage.
* `i_data.ctab`: Intron-level coordinates and junction read counts.
* `t_data.ctab`: Transcript-level abundance estimates (raw coverage, FPKM) linked to explicit `t_name` and `gene_id` identifiers.
* `e2t.ctab` & `i2t.ctab`: Relational mapping dictionaries connecting exons and introns to parent transcript identifiers.

These tables preserve both normalized metrics (FPKM) for proportion transformations (TPM) and the structural mapping dictionary needed to convert StringTie internal transfrag IDs (`MSTRG`) back to canonical TAIR locus identifiers.

### C. Expected Raw Count Derivation via `prepDE.py3`
Statistical modeling packages (`edgeR`, `DESeq2`) are formulated around discrete probability distributions (negative binomial) that explicitly model sampling variance, library size differences, and overdispersion. Normalized metrics such as FPKM, RPKM, or TPM violate these distribution assumptions and cannot be supplied directly to edgeR.

The Johns Hopkins Center for Computational Biology (CCB) utility `prepDE.py3` computes expected raw integer read counts ($N$) from StringTie transcript coverage ($C$) and length ($L$):

$$N = \frac{C \cdot L}{R}$$

Where $R$ represents the nominal sequencing read length (100 bp for this NovaSeq library). Gene-level matrices are produced by summing expected counts across all transcript models assigned to the same `gene_id`.

---

## 2. Execution Protocol

The quantification and matrix compilation pipeline was automated using a Bash driver (`04_quantification_20260323.sh`).

### A. StringTie Abundance Estimation
Quantification was executed in parallel across 10 threads per library:

```bash
SAMPLES=(
    "R1-1_WT"
    "R2-1_WT"
    "R3-1_WT"
    "R1-9_Cond"
    "R2-9_Cond"
    "R3-9_Cond"
)

for NAME in "${SAMPLES[@]}"; do
    mkdir -p results/ballgown/${NAME}
    stringtie -e -B -p 10 \
        -G results/gtf/at_merged.gtf \
        -o results/ballgown/${NAME}/${NAME}.gtf \
        results/bam/${NAME}.bam
done
```

### B. Count Matrix Generation
Expected integer read counts were extracted from the Ballgown table hierarchy:

```bash
cd results/matrices
curl -O [http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3)
python3 prepDE.py3 -i ../ballgown
```

### C. Sample Phenotype Metadata
A metadata mapping table (`phenodata.csv`) was generated to link library identifiers to their biological condition:

```bash
printf "ids,condition\nR1-1_WT,WT\nR2-1_WT,WT\nR3-1_WT,WT\nR1-9_Cond,Cond\nR2-9_Cond,Cond\nR3-9_Cond,Cond\n" > phenodata.csv
```

---

## 3. Count Distribution Diagnostics & Outlier Anomaly

Distributional properties of the raw gene count matrix (`gene_count_matrix.csv`, 32,540 total genomic loci) were computed using Python:

```bash
python3 -c "import pandas as pd; df = pd.read_csv('gene_count_matrix.csv', index_col=0); print(df.describe())"
```

### Summary Statistics

| Statistic | R1-1_WT | R1-9_Cond | R2-1_WT | R2-9_Cond | R3-1_WT | R3-9_Cond |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Count (Genes)** | 32,540 | 32,540 | 32,540 | 32,540 | 32,540 | 32,540 |
| **Mean Count** | 1,816.89 | 1,864.30 | 1,923.67 | 1,898.09 | 1,853.21 | 1,934.90 |
| **Standard Deviation** | **12,802.26** | 7,093.06 | 8,556.71 | 8,355.40 | 7,894.84 | 7,836.69 |
| **Minimum** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| **25th Percentile** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| **Median (50th)** | **104.00** | 170.00 | 152.50 | 140.00 | 152.00 | 135.00 |
| **75th Percentile** | 1,093.25 | 1,391.25 | 1,379.00 | 1,273.25 | 1,360.00 | 1,311.00 |
| **Maximum Count** | **1,531,773** | 405,114 | 783,371 | 598,885 | 735,200 | 346,520 |

### Technical Observations & Outlier Profile

1. **High Positive Skewness:** Across all libraries, the median expression level (~104–170 counts) is more than an order of magnitude lower than the mean (~1,816–1,935 counts). At least 25% of annotated features possess zero mapped reads. This sparsity is typical in eukaryotic transcriptomes and confirms the necessity of pre-filtering lowly expressed genes (`filterByExpr`) to minimize multiple testing burdens.
2. **Library Distortion in `R1-1_WT`:**
   * **Excessive Dispersion:** The standard deviation of read counts in `R1-1_WT` reaches **12,802.26**, compared to 7,093–8,556 across the other five replicates.
   * **Depressed Median:** The median expression in `R1-1_WT` is **104.00**, noticeably lower than all other samples (135.00–170.00).
   * **Hyper-Abundant Outlier Feature:** The maximum gene count in `R1-1_WT` reaches **1,531,773 reads**, roughly double the maximum observed in any other library (e.g., `R2-1_WT` max = 783,371; `R3-9_Cond` max = 346,520). A massive fraction of the sequencing library is monopolized by an isolated outlier gene, artificially deflating relative count representations of the remaining genes.
3. **Implication for Downstream edgeR Modeling:** This distributional distortion confirms the upstream alignment anomaly documented on March 13 (88.56% alignment rate). Compositional bias of this magnitude will require TMM normalization scaling factors and rigorous PCA diagnostics to determine if `R1-1_WT` compromises biological dispersion estimation.

---

## 4. Output Artifact Inventory

```text
results/
├── ballgown/
│   ├── R1-1_WT/      # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
│   ├── R2-1_WT/      # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
│   ├── R3-1_WT/      # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
│   ├── R1-9_Cond/    # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
│   ├── R2-9_Cond/    # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
│   └── R3-9_Cond/    # e_data.ctab, e2t.ctab, i_data.ctab, i2t.ctab, t_data.ctab
└── matrices/
    ├── prepDE.py3                   # CCB count matrix extraction utility
    ├── transcript_count_matrix.csv  # Raw integer counts per transcript model
    ├── gene_count_matrix.csv        # Raw integer counts per gene locus (32,540 loci)
    └── phenodata.csv                # Sample condition metadata mapping
```

---

## 5. Pipeline Status

Upstream command-line processing (HISAT2 -> StringTie -> prepDE) is complete. The count matrices, metadata, and Ballgown abundance files are ready for transfer to the R statistical environment for TPM conversion (`05_fpkm_to_tpm_20260331.R`) and edgeR modeling (`06_data_loading_20260406.R`).