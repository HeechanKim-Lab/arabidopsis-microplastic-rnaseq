# Arabidopsis Root RNA-seq: LDPE Plastic Residues in Soil

RNA-seq analysis of *Arabidopsis thaliana* (Col-0) roots grown in soil containing low-density polyethylene (LDPE) plastic residues. The repository has two parts: a six-library pilot processed from FASTQ files (HISAT2 → StringTie → edgeR), and a 27-library factorial experiment analysed in edgeR from a gene count matrix.

---

## 1. Overview

LDPE mulching film breaks down into plastic residues (PRs) in agricultural soil. In this system, LDPE residues had little visible effect on seedling fresh weight, survival or flowering time. This project asks whether root gene expression changes with the amount of plastic in the soil (0, 1 or 5% w/w) and with how long the plastic-amended soil was incubated before planting (0, 30 or 120 days).

**Prior work:** [Lee et al. (2025), *Journal of Experimental Botany* 76: 7003–7025](https://academic.oup.com/jxb/article/76/22/7003/8213809) describes this soil and plastic system and its effects on plant growth, the rhizosphere bacterial community and the leaf transcriptome.

---

## 2. Experimental design

* **Plant:** *Arabidopsis thaliana* Col-0.
* **Soil:** silty loam agricultural soil with 1% perlite, kept at 70% water-holding capacity.
* **Plastic:** LDPE mulching film cut into pieces of ≥ 5 mm and surface-sterilised with UV-C; 0, 1 or 5% (w/w).
* **Soil pre-incubation:** 0, 30 or 120 days in the dark at 25 °C and 80% relative humidity.
* **Growth:** 7-day-old seedlings were transplanted and grown for 21–28 days in growth chambers. Roots were sonicated to separate tightly adhering rhizosphere soil.
* **Replication:** 3 biological replicates per condition.

| Condition | Plastic | Pre-incubation | Libraries |
| :--- | :--- | :--- | :--- |
| R1 (`c0_d0`) | 0% | 0 d | R1-1, R1-2, R1-3 |
| R2 (`c1_d0`) | 1% | 0 d | R2-1, R2-2, R2-3 |
| R3 (`c5_d0`) | 5% | 0 d | R3-1, R3-2, R3-3 |
| R4 (`c0_d30`) | 0% | 30 d | R4-1, R4-2, R4-3 |
| R5 (`c1_d30`) | 1% | 30 d | R5-1, R5-2, R5-3 |
| R6 (`c5_d30`) | 5% | 30 d | R6-1, R6-2, R6-3 |
| R7 (`c0_d120`) | 0% | 120 d | R7-1, R7-2, R7-3 |
| R8 (`c1_d120`) | 1% | 120 d | R8-1, R8-2, R8-3 |
| R9 (`c5_d120`) | 5% | 120 d | R9-1, R9-2, R9-3 |

Datasets analysed:

* **Pilot (`01_pilot_analysis`):** 6 libraries. WT: `R1-1_WT`, `R2-1_WT`, `R3-1_WT` (0% plastic control). Cond: `R1-9_Cond`, `R2-9_Cond`, `R3-9_Cond` (5% plastic, 120-day pre-incubated soil).
* **Root00 (`02_full_experiment`):** all 27 libraries.
* **Root01 (`02_full_experiment`):** the 12 libraries of conditions R1, R3, R7 and R9 (0 or 5% plastic × 0 or 120 days). In figures these are labelled Conditions 1–4.

---

## 3. Data

* **Sequencing:** Illumina NovaSeq 6000, paired-end 2 × 100 bp, TruSeq Stranded mRNA (poly-A selection).
* **Pilot:** paired-end FASTQ files, 29.1–30.3 million read pairs per library. Not included in this repository.
* **Full experiment:** a gene-level count matrix for the 27 libraries (TAIR gene IDs) and the R objects derived from it (`Arabidopsis_plastics_2023.RData`) were provided by the lab. The upstream processing of these 27 libraries is not part of this repository.
* **Annotation table:** `Arabidopsis_genes_description.xlsx` (gene symbols and descriptions; lab file).

---

## 4. Methods

### Pilot (`01_pilot_analysis`)

1. **Reference:** TAIR10 genome (Ensembl Plants FASTA; chromosome names changed to `Chr1`–`Chr5`, `ChrM`, `ChrC`) and `TAIR10_GFF3_genes.gff`, converted to GTF with `gffread`.
2. **Alignment:** `hisat2 --dta`.
3. **Assembly:** `samtools sort`, StringTie per library (`-G`), `stringtie --merge`, and `gffcompare` against TAIR10.
4. **Quantification:** `stringtie -e -B` against the merged GTF; gene counts with `prepDE.py3` (32,540 genes); TPM computed from FPKM.
5. **Filtering and normalisation:** `filterByExpr` (20,971 genes kept) and TMM in edgeR. PCA on log-CPM values (prior count 2, `prcomp` without scaling).
6. **Differential expression:** edgeR quasi-likelihood F-tests (QLF) and likelihood ratio tests (LRT); see Key results.
7. **Enrichment:** StringTie gene IDs were mapped to TAIR IDs using `t_data.ctab`. GO enrichment used `enrichGO` (org.At.tair.db; BP, MF and CC). KEGG enrichment used `enrichKEGG` (organism `ath`, TAIR IDs).

### Full experiment (`02_full_experiment`)

1. **Input objects (from the lab):** `filterByExpr` and TMM normalisation per dataset, design `~0 + group`, dispersions from `estimateDisp`. Root01 retains 17,560 genes.
2. **PCA:** `prcomp` on TMM log-CPM values (centred and scaled).
3. **Differential expression:** edgeR GLM likelihood ratio tests (`glmFit`, `glmLRT`) with the contrasts listed below.
4. **Gene sets:** up and down counts used |log2FC| > 1, P < 0.01 and FDR < 0.05.
5. **Enrichment:** GO for the 146 interaction genes and the 1,240 genes of the 5% day-120 vs day-0 contrast; KEGG pathway membership for the 146 interaction genes.
6. **Multivariate tests (Root01):** PERMANOVA (`vegan::adonis2`, Euclidean distance on TMM log-CPM, Day × Plastic, 999 permutations, terms added sequentially) and two-way ANOVA (Day × Plastic) on PC1–PC3 scores. Tukey letters in the boxplots come from a one-way ANOVA across the four conditions.
7. **Targeted heatmaps:** genes significant (FDR < 0.05) in any of three Root01 contrasts whose annotation matched oxidative-stress or photosynthesis/plastid keywords.

| Contrast | Definition |
| :--- | :--- |
| 5% vs 0% at day 0 | `c5_d0 - c0_d0` |
| 5% vs 0% at day 30 (Root00 only) | `c5_d30 - c0_d30` |
| 5% vs 0% at day 120 | `c5_d120 - c0_d120` |
| 1% vs 0% at day 120 (Root00 only) | `c1_d120 - c0_d120` |
| Day 120 vs day 0 at 0% (Root01) | `c0_d120 - c0_d0` |
| Day 120 vs day 0 at 5% (Root01) | `c5_d120 - c5_d0` |
| Plastic × incubation interaction | `(c5_d120 - c5_d0) - (c0_d120 - c0_d0)` |

---

## 5. Key results

### Pilot

| Samples | Test | Genes tested | FDR < 0.05 |
| :--- | :--- | ---: | :--- |
| 6 | QLF | 20,971 | 1 |
| 6 | QLF, `robust = TRUE` | 20,971 | 1 |
| 5 (`R1-1_WT` removed) | QLF | 21,280 | 12 |
| 4 (`R1-1_WT` and `R3-9_Cond` removed) | QLF | 21,132 | 0 |
| 6 | LRT | 20,971 | 932; 865 with \|log2FC\| > 1 |

* The 865-gene LRT list (all six libraries) was used for the volcano plot, the heatmap and the GO/KEGG analysis. With three libraries per group, the LRT is less conservative than the QLF test because it does not account for uncertainty in the dispersion estimates, so this list is exploratory.
* `R1-1_WT` had the lowest alignment rate (88.56%; others 90.88–93.78%) and was separated from the other five libraries along PC1 (49% of variance). It was kept in the LRT model.
* edgeR used `Cond` as the reference level, so a positive log2FC means higher expression in WT. In the GO/KEGG figures, "UP" refers to genes higher in WT and "DOWN" to genes higher in Cond.
* Genes higher in Cond were enriched for cell wall organisation, root hair and root epidermal cell differentiation, peroxidase activity and phenylpropanoid biosynthesis. Genes higher in WT were enriched for responses to nutrient levels, starvation and hypoxia, and for amino-acid catabolism.

### Full experiment (LRT, FDR < 0.05)

| Contrast | Root00 (27 libraries) | Root01 (12 libraries) |
| :--- | ---: | ---: |
| 5% vs 0% at day 0 | 0 | 0 |
| 5% vs 0% at day 30 | 0 | — |
| 5% vs 0% at day 120 | 48 | 138 |
| 1% vs 0% at day 120 | 0 | — |
| Day 120 vs day 0 at 0% | — | 49 |
| Day 120 vs day 0 at 5% | — | 1,240 |
| Plastic × incubation interaction | 99 | 146 |

* With |log2FC| > 1, P < 0.01 and FDR < 0.05 (Root01):
  * 5% vs 0% at day 120: 136 genes (43 up, 93 down).
  * Day 120 vs day 0 at 0%: 49 genes (41 up, 8 down).
  * Day 120 vs day 0 at 5%: 1,123 genes (540 up, 583 down).
* The same 5% vs 0% day-120 contrast gives 48 genes in the 27-library model and 138 in the 12-library model, because the dispersions are estimated from different sets of libraries.
* PERMANOVA (Root01): Day R² = 0.195, p = 0.010; Plastic R² = 0.057, p = 0.722; Day × Plastic R² = 0.108, p = 0.184.
* Two-way ANOVA on PC1–PC3 scores: no term reached p < 0.05 (smallest: Day on PC1, p = 0.069).
* In summary, soil pre-incubation time had the larger effect on root expression. A plastic effect was detected only after 120 days of pre-incubation, and the day-120 vs day-0 difference involved more genes in 5% soil (1,240) than in 0% soil (49).

---

## 6. Limitations

* Three biological replicates per condition. All reported gene lists come from LRTs, which are less conservative than QLF tests at this sample size. No result has been validated by an independent method.
* The legends of PCA figures 020–022 label the colours as plastic levels, but the colours show the pre-incubation time (see the 2026-04-11 full-experiment log).
* The figures were produced in interactive R sessions, and the logs record that code. The scripts in `modules/` and `pipeline/` were reorganised afterwards and have not been re-run end to end. Known differences:
  * the pilot PCA scripts use `scale. = TRUE`, but the figure used unscaled data;
  * the pilot LRT scripts drop `R1-1_WT`, but the reported list used all six libraries;
  * the pilot KEGG scripts use Entrez IDs, but the analysis used TAIR IDs;
  * full-experiment modules 04 and 05 contain lines that do not parse;
  * module 06 writes different file names from the committed figures.

---

## 7. How to run

```bash
# Pilot
cd 01_pilot_analysis
./pipeline/run_upstream_tuxedo.sh      # bash >= 4; uses BSD sed (macOS)
Rscript pipeline/run_downstream_pilot.R

# Full experiment (needs Arabidopsis_plastics_2023.RData and Arabidopsis_genes_description.xlsx)
cd ../02_full_experiment
Rscript pipeline/run_full_experiment.R
```

Each step can also be run on its own from `modules/`. Intermediate objects are saved to and read from `results/intermediate/*.rds`. See Limitations before using the scripts to reproduce the figures.

---

## 8. Repository layout

```text
.
├── README.md
├── 01_pilot_analysis/
│   ├── assets/     # figures 001–013
│   ├── logs/       # analysis logs, 2026-03-10 to 2026-04-11
│   ├── modules/    # 01_reference_and_index.sh … 11_functional_enrichment.R
│   └── pipeline/   # run_upstream_tuxedo.sh, run_downstream_pilot.R
└── 02_full_experiment/
    ├── assets/     # figures 020–040
    ├── logs/       # analysis logs, 2026-04-11 to 2026-04-17
    ├── modules/    # 01_full_matrix_and_subset_pca.R … 07_master_annotation_registry.R
    └── pipeline/   # run_full_experiment.R
```

---

## 9. Software

* **Upstream (conda environment `lab_intel`, macOS):** HISAT2 v2.2.2, StringTie v3.0.3, samtools v1.22.1, gffcompare v0.12.10, gffread, Python 3 with pandas.
* **R:** R v4.4.0; edgeR v4.4.0, clusterProfiler v4.13.3, org.At.tair.db v3.19.1, vegan v2.6-6.1, pheatmap v1.0.12, ggplot2 v3.5.1, multcompView v0.1-10, EnhancedVolcano, ggvenn, plotly, readxl, dplyr, reshape2.

---

## 10. References and annotation

* **Genome:** TAIR10 (`Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`, Ensembl Plants).
* **Annotation:** `TAIR10_GFF3_genes.gff`, converted to `TAIR10_GFF3_genes.gtf` with `gffread`.
* **Gene descriptions:** `Arabidopsis_genes_description.xlsx`.
* **Databases:** `org.At.tair.db` (Bioconductor); KEGG organism code `ath`.
