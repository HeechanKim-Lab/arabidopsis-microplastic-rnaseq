# Quantification and Count Matrices

**Date:** 2026-03-23  
**Stage:** 01_pilot_analysis  
**Script:** `modules/04_quantification.sh`  
**Input:** sorted BAM files, `at_merged.gtf`  
**Output:** `ballgown/<sample>/*.ctab`, `gene_count_matrix.csv`, `transcript_count_matrix.csv`, `phenodata.csv`

---

## 1. Method

* `stringtie -e` estimates abundances only for the transcripts in `at_merged.gtf` (no new assembly), so all six libraries are quantified against the same transcript models.
* `-B` writes Ballgown tables for each library (`e_data.ctab`, `i_data.ctab`, `t_data.ctab`, `e2t.ctab`, `i2t.ctab`). `t_data.ctab` holds coverage and FPKM per transcript and links StringTie gene IDs (`MSTRG.*`) to TAIR transcript names (`t_name`).
* `prepDE.py3` converts coverage to read counts (count = coverage × transcript length / read length) and sums transcripts per gene, giving a gene count matrix for edgeR.

---

## 2. Commands

### A. Abundance estimation (10 threads)

```bash
stringtie -e -B -p 10 -G at_merged.gtf -o ballgown/R1-1_WT/R1-1_WT.gtf R1-1_WT.bam
# same command for R2-1_WT, R3-1_WT, R1-9_Cond, R2-9_Cond and R3-9_Cond

head -n 10 ballgown/R1-1_WT/t_data.ctab   # check FPKM values
```

### B. Count matrices

```bash
curl -O http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
python3 prepDE.py3 -i ballgown      # writes gene_count_matrix.csv and transcript_count_matrix.csv
head transcript_count_matrix.csv
```

`prepDE.py3` was run without `-l`, so it used its default average read length of 75 bp. The reads are 2 × 100 bp (see README), so the derived counts are about 1.33 times the read counts. The factor is the same for every gene.

### C. Sample table

```bash
printf "ids,condition\nR1-1_WT,WT\nR2-1_WT,WT\nR3-1_WT,WT\nR1-9_Cond,Cond\nR2-9_Cond,Cond\nR3-9_Cond,Cond" > phenodata.csv
```

---

## 3. Gene count summary

```bash
python3 -c "import pandas as pd; df = pd.read_csv('gene_count_matrix.csv', index_col=0); print(df.describe())"
```

| Statistic | R1-1_WT | R1-9_Cond | R2-1_WT | R2-9_Cond | R3-1_WT | R3-9_Cond |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Genes | 32,540 | 32,540 | 32,540 | 32,540 | 32,540 | 32,540 |
| Mean | 1,816.89 | 1,864.30 | 1,923.67 | 1,898.09 | 1,853.21 | 1,934.90 |
| SD | 12,802.26 | 7,093.06 | 8,556.71 | 8,355.40 | 7,894.84 | 7,836.69 |
| Minimum | 0 | 0 | 0 | 0 | 0 | 0 |
| 25th percentile | 0 | 0 | 0 | 0 | 0 | 0 |
| Median | 104.00 | 170.00 | 152.50 | 140.00 | 152.00 | 135.00 |
| 75th percentile | 1,093.25 | 1,391.25 | 1,379.00 | 1,273.25 | 1,360.00 | 1,311.00 |
| Maximum | 1,531,773 | 405,114 | 783,371 | 598,885 | 735,200 | 346,520 |

### Observations

* In every library at least 25% of genes have zero counts and the median is far below the mean, so low-count genes were filtered before modelling (`filterByExpr`, 2026-04-07 log).
* `R1-1_WT` differs from the other five libraries: higher SD (12,802 vs 7,093–8,556), lower median (104 vs 135–170) and a higher maximum (1,531,773 vs 346,520–783,371). Together with its lower alignment rate (2026-03-13 log), this was checked again in the PCA.

Next steps: TPM table (`modules/05_fpkm_to_tpm.R`) and loading into edgeR (`modules/06_data_ingestion.R`).
