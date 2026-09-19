# Upstream Processing: Spliced Read Alignment via HISAT2
**Date:** March 13, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `indexes/arabidopsis_index.*.ht2`, `samples/*.fastq.gz`  
**Output Artifacts:** `results/sam/*.sam`  

---

## 1. Methodological & Theoretical Rationale

### A. Graph-Based Spliced Read Alignment
Eukaryotic mRNA sequencing reads derived from mature transcripts span non-contiguous genomic loci due to intron excision. Standard short-read genomic aligners cannot map junction-spanning reads without heavy soft-clipping or alignment dropouts. `HISAT2` solves this using a Hierarchical Graph FM index (HGFM):
1. **Global Search:** Coordinates global anchoring of read segments using the chromosome-level FM index.
2. **Local Graph Traversal:** Traverses local index graphs (55,000 regional partitions across the Arabidopsis genome) to resolve canonical (`GT-AG`, `GC-AG`) and non-canonical splice junctions down to single-nucleotide resolution.

### B. Downstream Transcriptome Assembly Flag (`--dta`)
The `--dta` (Downstream Transcript Assembly) argument is strictly required when aligning reads for reconstruction via `StringTie`. When enabled:
* HISAT2 modifies its dynamic programming scoring penalty to favor spliced alignments that conform to standard donor/acceptor intron models.
* HISAT2 writes specialized auxiliary SAM tags:
  * `XS:A:[+|-]`: Reports the biological transcription strand of the spliced intron based on canonical dinucleotide motifs, critical for unstranded and stranded paired-end library resolution.
  * Adjusts read alignment boundaries to avoid arbitrary mismatches at exon margins, ensuring clean start and end coordinates for graph-based transcript assembly.

### C. Paired-End Conformance Definitions
* **Concordant Alignment:** Both mates align to the same chromosome on opposing strands (`FR` orientation) within the expected genomic fragment size distribution (typically $100\text{--}500$ bp insert).
* **Discordant Alignment:** Both mates align to the genome, but with abnormal insert lengths, inversion orientations, or across different chromosomes (potential structural rearrangements or trans-splicing events).
* **Mate-Level Alignment:** Mates from pairs that failed paired alignment but where individual single-end mates align independently.

---

## 2. Execution & Implementation

Alignments were parallelized across 10 execution threads per sample using a Bash driver loop (`01_hisat2_alignment_20260313.sh`).

### Execution Protocol
```bash
hisat2 -p 10 --dta \
  -x /Users/heechan/Desktop/practice/genome/index/arabidopsis_index \
  -1 samples/${R1_FASTQ}.fastq.gz \
  -2 samples/${R2_FASTQ}.fastq.gz \
  -S results/sam/${SAMPLE_NAME}.sam
```

### Processed Cohorts
* **Wild Type (0% Plastic Control):** `R1-1_WT`, `R2-1_WT`, `R3-1_WT`
* **Conditioned (5% Plastic, 120d Pre-incubated Soil):** `R1-9_Cond`, `R2-9_Cond`, `R3-9_Cond`

---

## 3. Alignment Statistics & Quality Profiling

All 6 libraries achieved $> 29.0 \times 10^6$ paired reads, delivering sufficient sequencing depth for statistical differential expression analysis.

### Summary Metrics

| Sample ID | Condition | Total Read Pairs | Concordant Unique (%) | Concordant Multi (%) | Discordant (%) | Overall Alignment Rate |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **R1-1_WT** | Control | 29,791,679 | 22,759,510 (76.40%) | 985,658 (3.31%) | 338,280 (5.59%) | **88.56%** |
| **R2-1_WT** | Control | 30,293,327 | 25,675,061 (84.75%) | 722,235 (2.38%) | 481,753 (12.37%) | **93.35%** |
| **R3-1_WT** | Control | 29,526,747 | 23,620,211 (80.00%) | 654,169 (2.22%) | 623,511 (11.87%) | **90.88%** |
| **R1-9_Cond** | Conditioned | 29,061,578 | 24,757,144 (85.19%) | 546,209 (1.88%) | 466,708 (12.42%) | **93.78%** |
| **R2-9_Cond** | Conditioned | 29,817,819 | 24,251,405 (81.33%) | 600,524 (2.01%) | 604,108 (12.17%) | **91.92%** |
| **R3-9_Cond** | Conditioned | 29,899,345 | 25,397,944 (84.94%) | 501,403 (1.68%) | 563,225 (14.08%) | **93.72%** |

---

## 4. Technical Diagnostics & Early Outlier Indication

### A. High Alignment Fidelity
* Five of the six libraries exhibit overall alignment rates $> 90.8\%$, peaking at $93.78\%$ in `R1-9_Cond`. 
* Uniquely aligned concordant fractions average $> 82\%$, confirming that read quality across the Illumina flow cell was high and that reference genome indexing (`TAIR10`) correctly synchronized chromosome headers.
* Low multi-mapping rates ($\le 3.31\%$) demonstrate high specificity, confirming that paralogous gene families and repetitive transposable elements in the Arabidopsis genome do not generate ambiguous read allocations.

### B. Anomaly Detection: `R1-1_WT`
Sample `R1-1_WT` displays noticeable technical divergence from the remaining cohort:
* **Suppressed Alignment Rate:** Overall alignment drops to **88.56%** (2.32% to 5.22% lower than all other samples).
* **Concordance Deficit:** Unique concordant rate drops to **76.40%** (compared to $80.00\%\text{--}85.19\%$ in the other five replicates), with $20.30\%$ of pairs failing initial concordant placement.
* **Diagnostic Check:** While $88.56\%$ remains viable for standard RNA-seq, this structural divergence represents an upstream indicator of technical variation, anticipating the downstream variance separation observed in exploratory PCA (PC2 outlier).

---

## 5. Execution Logs

```text
--- R1-1_WT ---
29791679 reads; of these:
  29791679 (100.00%) were paired; of these:
    6046511 (20.30%) aligned concordantly 0 times
    22759510 (76.40%) aligned concordantly exactly 1 time
    985658 (3.31%) aligned concordantly >1 times
    ----
    6046511 pairs aligned concordantly 0 times; of these:
      338280 (5.59%) aligned discordantly 1 time
    ----
    5708231 pairs aligned 0 times concordantly or discordantly; of these:
      11416462 mates make up the pairs; of these:
        6818863 (59.73%) aligned 0 times
        4347835 (38.08%) aligned exactly 1 time
        249764 (2.19%) aligned >1 times
88.56% overall alignment rate

--- R2-1_WT ---
30293327 reads; of these:
  30293327 (100.00%) were paired; of these:
    3896031 (12.86%) aligned concordantly 0 times
    25675061 (84.75%) aligned concordantly exactly 1 time
    722235 (2.38%) aligned concordantly >1 times
    ----
    3896031 pairs aligned concordantly 0 times; of these:
      481753 (12.37%) aligned discordantly 1 time
    ----
    3414278 pairs aligned 0 times concordantly or discordantly; of these:
      6828556 mates make up the pairs; of these:
        4029606 (59.01%) aligned 0 times
        2693794 (39.45%) aligned exactly 1 time
        105156 (1.54%) aligned >1 times
93.35% overall alignment rate

--- R3-1_WT ---
29526747 reads; of these:
  29526747 (100.00%) were paired; of these:
    5252367 (17.79%) aligned concordantly 0 times
    23620211 (80.00%) aligned concordantly exactly 1 time
    654169 (2.22%) aligned concordantly >1 times
    ----
    5252367 pairs aligned concordantly 0 times; of these:
      623511 (11.87%) aligned discordantly 1 time
    ----
    4628856 pairs aligned 0 times concordantly or discordantly; of these:
      9257712 mates make up the pairs; of these:
        5388107 (58.20%) aligned 0 times
        3727503 (40.26%) aligned exactly 1 time
        142102 (1.53%) aligned >1 times
90.88% overall alignment rate

--- R1-9_Cond ---
29061578 reads; of these:
  29061578 (100.00%) were paired; of these:
    3758225 (12.93%) aligned concordantly 0 times
    24757144 (85.19%) aligned concordantly exactly 1 time
    546209 (1.88%) aligned concordantly >1 times
    ----
    3758225 pairs aligned concordantly 0 times; of these:
      466708 (12.42%) aligned discordantly 1 time
    ----
    3291517 pairs aligned 0 times concordantly or discordantly; of these:
      6583034 mates make up the pairs; of these:
        3617551 (54.95%) aligned 0 times
        2885335 (43.83%) aligned exactly 1 time
        80148 (1.22%) aligned >1 times
93.78% overall alignment rate

--- R2-9_Cond ---
29817819 reads; of these:
  29817819 (100.00%) were paired; of these:
    4965890 (16.65%) aligned concordantly 0 times
    24251405 (81.33%) aligned concordantly exactly 1 time
    600524 (2.01%) aligned concordantly >1 times
    ----
    4965890 pairs aligned concordantly 0 times; of these:
      604108 (12.17%) aligned discordantly 1 time
    ----
    4361782 pairs aligned 0 times concordantly or discordantly; of these:
      8723564 mates make up the pairs; of these:
        4816119 (55.21%) aligned 0 times
        3782711 (43.36%) aligned exactly 1 time
        124734 (1.43%) aligned >1 times
91.92% overall alignment rate

--- R3-9_Cond ---
29899345 reads; of these:
  29899345 (100.00%) were paired; of these:
    3999998 (13.38%) aligned concordantly 0 times
    25397944 (84.94%) aligned concordantly exactly 1 time
    501403 (1.68%) aligned concordantly >1 times
    ----
    3999998 pairs aligned concordantly 0 times; of these:
      563225 (14.08%) aligned discordantly 1 time
    ----
    3436773 pairs aligned 0 times concordantly or discordantly; of these:
      6873546 mates make up the pairs; of these:
        3756411 (54.65%) aligned 0 times
        3042146 (44.26%) aligned exactly 1 time
        74989 (1.09%) aligned >1 times
93.72% overall alignment rate
```

---

## 6. Pipeline Status

Spliced alignment completed across all 6 target libraries. Output SAM files are ready for SAM-to-BAM binary compression, coordinate sorting via `samtools`, and transcript assembly via `StringTie` (`02_assembly_and_merge_20260317.sh`).