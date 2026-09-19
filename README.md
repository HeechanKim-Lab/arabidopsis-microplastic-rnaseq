# Arabidopsis Root Transcriptome Pipeline (Microplastic Stress)
Dry-lab computational log and statistical modeling pipeline for profiling transcriptional alterations in *Arabidopsis thaliana* exposed to low-density polyethylene (LDPE) plastic residues across acute and chronic soil incubation regimes.

---

## 1. Project Background & Biological Context

Agricultural mulching films composed of low-density polyethylene (LDPE) fragment into macroscopic plastic residues (PRs, $\ge 5$ mm) and microplastics in arable topsoil. While macroscopic seedling phenotypes (fresh weight, survival rate, flowering time) exhibit high homeostatic buffering against chemically inert LDPE fragments, significant transcriptional perturbation occurs at the root-soil interface.

Physical sonication of roots to isolate tightly adhering rhizosphere soil introduces substantial within-group technical variance. Standard negative binomial Generalized Linear Model (GLM) pipelines utilizing Quasi-Likelihood F-tests (QLF) become overly conservative under sonication-induced dispersion, collapsing statistical power and yielding zero significant differentially expressed genes (DEGs). This pipeline implements a calibrated **Likelihood Ratio Test (LRT)** framework, recovering the authentic transcriptomic stress response across plastic dosage and incubation timelines.

```text
Experimental Design Schema:
Silty Loam Agricultural Soil + 0%, 1%, 5% (w/w) LDPE Residues (≥ 5 mm)
      │
      ├── Controlled Dark Incubation (25°C, 80% RH) for 0 d, 30 d, 120 d
      │
      ▼
Transplantation of 7-day-old Arabidopsis thaliana (Col-0) Seedlings
      │
      ├── Cultivation: 21–28 days in environmental growth chambers
      │
      ▼
Sampling & Computational Modularity:
├── 01_pilot_analysis   : Diagnostic benchmark on isolated subset (LRT proof-of-concept)
└── 02_full_experiment  : Full multi-factorial series (Root00: N=27 global; Root01: N=12 core)
```

---

## 2. Two-Stage Computational Architecture

To resolve the high dispersion introduced during root tissue harvesting, the project was organized into two sequential analytical stages:

### Stage 1: Diagnostic Calibration (`01_pilot_analysis`)
Initial evaluations of root tissue under standard QLF pipelines failed to detect meaningful biological expression changes (yielding 0–1 significant DEGs), suggesting either biological inertness or severe variance inflation. To diagnose this without confounding multi-factorial variables, an isolated 6-sample subset (3 WT vs. 3 Conditioned replicates at 120 days) was allocated as a pilot benchmark. 
* Identified sample `R1-1_WT` as an extreme variance outlier driving PC2 dispersion (18% variance).
* Demonstrated that algorithmic dispersion downweighting (`robust = TRUE`) and sample over-trimming ($N = 4, df = 2$) fail under sonication noise.
* Transitioned from QLF to a calibrated Likelihood Ratio Test (LRT) with outlier excision, unmasking **865 high-confidence DEGs** (FDR $< 0.05$, $\vert{}\log_2\text{FC}\vert{} > 1$). This proved that the conditioning response was biologically intact.

### Stage 2: Production Scaling (`02_full_experiment`)
With statistical viability confirmed, the LRT modeling framework was deployed across the complete experimental series:
* **`Root00` (Global Landscape, $N = 27$):** Screened across 3 plastic concentrations ($0\%, 1\%, 5\%$) and 3 incubation intervals ($0\text{d}, 30\text{d}, 120\text{d}$) to characterize global variance and dose sensitivity.
* **`Root01` (Factorial Core, $N = 12$):** Modeled the 4 extreme condition endpoints (`c0_d0`, `c5_d0`, `c0_d120`, `c5_d120`) to parameterize chronological aging, plastic-accelerated aging (1,240 DEGs), and the non-additive Plastic $\times$ Time interaction (146 DEGs).

---

## 3. Experimental Setup & Cohort Specifications

* **Host Plant:** *Arabidopsis thaliana* L. ecotype Columbia-0 (Col-0)
* **Base Substrate:** Silty loam agricultural field soil mixed with 1% perlite (maintained at 70% water-holding capacity)
* **Contaminant:** Low-density polyethylene (LDPE) agricultural mulching film pulverized to $\ge 5$ mm, UV-C surface-sterilized
* **Treatment Variables:**
  * **Plastic Residue (PR) Concentration:** 0% (Control), 1%, and 5% $(w/w)$
  * **Soil Pre-incubation Duration:** 0 days, 30 days, and 120 days at 25°C, 80% RH
* **Analyzed Cohorts:**
  * **Pilot Dataset (`01_pilot_analysis`):** $N = 6$ root libraries (3 WT vs. 3 Conditioned replicates; `R1-1_WT` diagnosed and excluded as an outlier).
  * **Full Dataset (`02_full_experiment`):**
    * `Root00` (Global Matrix): $N = 27$ libraries covering 9 conditions $\times$ 3 biological replicates across 0d, 30d, and 120d.
    * `Root01` (Characteristic Subset): $N = 12$ libraries isolating the primary factorial endpoints (`c0_d0`, `c5_d0`, `c0_d120`, `c5_d120`).

---

## 4. Sequencing Metadata & Library Architecture

* **Platform:** Illumina NovaSeq 6000
* **Configuration:** Paired-End (PE), $2 \times 100$ bp
* **Library Chemistry:** TruSeq Stranded mRNA Prep Kit (poly-A selection)
* **Reference Genome:** *Arabidopsis thaliana* TAIR10 (Ensembl Plants release)

### Replicate Architecture

```text
Pilot Replicates:
├── R1-1_WT, R2-1_WT, R3-1_WT       : 0% Plastic Control roots
└── R1-9_Cond, R2-9_Cond, R3-9_Cond : 5% Plastic Conditioned roots (120d)

Core Factorial Design (Root01 Subset):
├── c0_d0   (R1-1, R1-2, R1-3) : 0% Plastic, 0-Day incubation (Baseline Control)
├── c5_d0   (R3-1, R3-2, R3-3) : 5% Plastic, 0-Day incubation (Acute Exposure)
├── c0_d120 (R7-1, R7-2, R7-3) : 0% Plastic, 120-Day incubation (Aging Control)
└── c5_d120 (R9-1, R9-2, R9-3) : 5% Plastic, 120-Day incubation (Chronic Conditioned)
```

---

## 5. Repository Structure

The repository maintains an identical architecture across both modules, bifurcated into monolithic drivers (`pipeline/`) and step-by-step unit scripts (`modules/`):

```text
.
├── README.md
├── session_info.txt
├── 01_pilot_analysis/
│   ├── assets/       # Standardized diagnostic assets (001_pilot_qc_pca.png to 013_pilot_kegg_pathways_down.png)
│   ├── logs/         # Markdown logs for upstream Tuxedo and edgeR calibration (2026-03-10 to 2026-04-11)
│   ├── pipeline/     # Monolithic production drivers (run_upstream_tuxedo.sh, run_downstream_pilot.R)
│   └── modules/      # Step-wise reproducible scripts (01_reference_and_index.sh to 11_functional_enrichment.R)
└── 02_full_experiment/
    ├── assets/       # Publication figures (020_root00_full_matrix_pca.png to 040_targeted_photosynthesis_heatmap.png)
    ├── logs/         # Markdown logs for 27-sample matrix and 12-sample subset (2026-04-11 to 2026-04-17)
    ├── pipeline/     # Monolithic production driver (run_full_experiment.R)
    └── modules/      # Step-wise analytical scripts (01_full_matrix_and_subset_pca.R to 07_master_annotation_registry.R)
```

---

## 6. Computational Workflow

```text
                          Raw FASTQ Reads (NovaSeq PE)
                                       │
                    01. HISAT2 Alignment (--dta, TAIR10)
                                       │
                    02. Samtools Coordinate Sorting (BAM)
                                       │
                    03. StringTie Transcript Assembly (-G)
                                       │
                    04. Transcriptome Merge (at_merged.gtf)
                                       │
                    05. StringTie Abundance Estimation (-e -B)
                                       │
          ┌────────────────────────────┴────────────────────────────┐
          ▼                                                         ▼
prepDE.py3 (Gene Count Matrix)                             ballgown/t_data.ctab
          │                                                         │
          ▼                                                         ▼
edgeR DGEList Initialization                               FPKM to TPM Conversion
          │                                                         │
          ├────────────────────────────┬────────────────────────────┘
          ▼                            ▼
[01 Pilot Diagnostic Engine]     [02 Full Experiment Engine]
├── filterByExpr & TMM           ├── 27-Sample Variance Screening
├── PCA QC (R1-1_WT Outlier)     ├── Root01 Factorial Contrasts
├── QLF vs. LRT Benchmark        ├── GLM LRT Model (Interaction)
└── MSTRG-to-TAIR GO/KEGG        ├── PERMANOVA & 2-Way ANOVA
                                 └── ROS & Photosynthesis Heatmaps
```

---

### Step-by-Step Processing Strategy

#### 1. Upstream Read Mapping & Assembly
* Paired-end reads are mapped to the TAIR10 genome assembly using `hisat2` (v2.2.2) with `--dta` enabled to supply splice-site information downstream.
* Alignments are coordinate-sorted via `samtools` (v1.22.1) and individual transcript models are assembled using `stringtie` (v3.0.3).
* Sample GTF models are merged into a non-redundant master transcriptome (`at_merged.gtf`) via `stringtie --merge` and evaluated for novel transcript structures against TAIR10 using `gffcompare` (v0.12.10).

#### 2. Transcript Quantification & Normalization
* Sample-level transcript abundances are extracted using `stringtie -e -B` against the master GTF.
* Raw gene-level read count matrices are compiled using `prepDE.py3` for generalized linear modeling in `edgeR`.
* Ballgown FPKM abundance estimates are converted to Transcripts Per Million (TPM) using the linear transformation:
  $$TPM_i = \left( \frac{FPKM_i}{\sum FPKM_{\text{sample}}} \right) \times 10^6$$
  retaining transcript identity (`t_name`) to eliminate library size bias during cross-replicate expression scaling.

#### 3. Pilot Diagnostics & LRT Calibration (`01_pilot_analysis`)
* Unsupervised PCA on TMM log-CPM values identified sample `R1-1_WT` as an extreme variance outlier along PC2 (18% variance).
* Fitting a standard Quasi-Likelihood F-test (`glmQLFTest`) on all 6 replicates collapsed statistical sensitivity to only 1 significant DEG (FDR $< 0.05$).
* Removing `R1-1_WT` and switching from QLF to the **Likelihood Ratio Test (`glmLRT`)** resolved technical dispersion, extracting **865 significant DEGs** (FDR $< 0.05$, $\vert{}\log_2\text{FC}\vert{} > 1$).

#### 4. Factorial GLM Contrasts & Multivariate Inference (`02_full_experiment`)
* **Dose-Response Verification:** PCA of the characteristic 12-sample cohort (`Root01`) confirmed an uncorrupted biological gradient across PC1 (29.7% variance), separating samples cleanly along 0% $\rightarrow$ 1% $\rightarrow$ 5% plastic.
* **Factorial GLM Contrasts:** Structured contrast vectors under `edgeR`:
  * *Acute Plastic Effect (D0):* `c5_d0 - c0_d0`
  * *Chronic Plastic Effect (D120):* `c5_d120 - c0_d120`
  * *Control Aging Effect:* `c0_d120 - c0_d0`
  * *Plastic Aging Effect:* `c5_d120 - c5_d0`
  * *Plastic $\times$ Time Interaction:* `(c5_d120 - c5_d0) - (c0_d120 - c0_d0)`
* **Multivariate Hypothesis Testing:**
  * Non-parametric multivariate ANOVA (`vegan::adonis2`) evaluating variance partitioning across Treatment, Incubation, and their interaction.
  * Two-way ANOVA with post-hoc Tukey HSD tests (`multcompView`) on PC1, PC2, and PC3 scores to determine axis-specific significance.
* **Targeted Physiological Modules:** Targeted hierarchical heatmaps on critical response pathways:
  * Repression of photosynthesis and light-reaction machinery (`CAO`, `PSB27`).
  * Activation of reactive oxygen species (ROS) scavengers and oxidative stress mediators (`APX1`, class III peroxidases).
  * Perturbation of anion transmembrane transport and starch catabolism (`CT-BMY`).

---

## 7. Execution Protocol

### A. End-to-End Pipeline Drivers
Execute the monolithic drivers to reproduce the complete analysis:

```bash
# 01 Pilot Analysis
cd 01_pilot_analysis
./pipeline/run_upstream_tuxedo.sh
Rscript pipeline/run_downstream_pilot.R

# 02 Full Experiment Analysis
cd ../02_full_experiment
Rscript pipeline/run_full_experiment.R
```

### B. Modular Execution
Individual steps can be executed independently using the scripts located in `01_pilot_analysis/modules/` and `02_full_experiment/modules/`. Intermediate states are automatically saved and loaded via `results/intermediate/*.rds`.

---

## 8. Software Dependencies & Computing Environments

Upstream sequence processing was performed in an isolated Conda environment (`lab_intel`), while downstream statistical modeling and graphical evaluations were conducted within an R session.

### Upstream Processing Environment (`lab_intel`)
* **Operating System:** macOS (Darwin `aarch64` / `x86_64`)
* **HISAT2:** v2.2.2
* **StringTie:** v3.0.3
* **Samtools:** v1.22.1
* **GFFCompare:** v0.12.10
* **Python Runtime:** Python 3 (with `pandas` for matrix summaries)

### Statistical Modeling & Graphics Environment (R Session)
* **Base Platform:** R (v4.4.0)
* **Core Packages:** `edgeR` (v4.4.0), `clusterProfiler` (v4.13.3), `org.At.tair.db` (v3.19.1), `vegan` (v2.6-6.1)
* **Visualization:** `pheatmap` (v1.0.12), `ggplot2` (v3.5.1+), `multcompView` (v0.1-10), `ggrepel`, `plotly`

---

## 9. Genomic References & Annotations

* **Arabidopsis Reference Genome:** TAIR10 assembly (`Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`)
* **Annotation Sets:**
  * Primary Annotation: `TAIR10_GFF3_genes.gff` converted to `TAIR10_GFF3_genes.gtf` via `gffread`
  * Functional Dictionary: Curated locus symbol and descriptive table (`Arabidopsis_genes_description.xlsx`)
* **Organism Database:** `org.At.tair.db` (Bioconductor)
* **Organism KEGG Identifier:** `ath` (*Arabidopsis thaliana*)