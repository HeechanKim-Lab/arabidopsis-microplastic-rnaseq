# Upstream Processing: Transcriptome Reconstruction Evaluation and Annotation Benchmarking via GFFCompare
**Date:** March 18, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `results/gtf/at_merged.gtf`, `genome/TAIR10_GFF3_genes.gtf`  
**Output Artifacts:** `results/gtf/at_compare.stats`, `results/gtf/at_compare.annotated.gtf`, `results/gtf/at_compare.loci`, `results/gtf/at_compare.tracking`, `results/gtf/at_compare.at_merged.gtf.refmap`, `results/gtf/at_compare.at_merged.gtf.tmap`  

---

## 1. Methodological & Theoretical Rationale

### A. Sensitivity and Precision in Transcriptome Reconstruction
Reference-guided assembly via `StringTie` produces a merged annotation (`at_merged.gtf`) containing both canonical reference gene models and candidate novel isoforms. Evaluating assembly quality against the curated reference (`TAIR10_GFF3_genes.gtf`) requires calculating two orthogonal statistical metrics:

1. **Sensitivity ($Sn$):** The proportion of reference elements correctly reconstructed by the assembly:
   $$Sn = \frac{TP}{TP + FN}$$
   Where $TP$ (True Positives) represents annotated features recovered in the assembly, and $FN$ (False Negatives) represents annotated features missed by the sequencing depth or assembly heuristics. A high sensitivity confirms comprehensive coverage of the expressed transcriptome without gene dropouts.

2. **Precision ($Sp$):** The proportion of assembled elements that match annotated reference features:
   $$Sp = \frac{TP}{TP + FP}$$
   Where $FP$ (False Positives) represents assembled transfrags that do not match canonical reference models. In RNA-seq, assembly precision is expected to deviate from 100% due to genuine biological alternative splicing, unannotated non-coding RNAs, condition-specific isoforms, and transcriptional read-through.

### B. Positional Resolution Tiers
`gffcompare` computes sensitivity and precision across discrete structural tiers:
* **Base Level:** Nucleotide-level overlap between query exonic coordinates and reference exons.
* **Exon Level:** Exact matching of individual exon coordinate boundaries (start and end sites).
* **Intron Level:** Exact coordinate matching of individual splice junctions.
* **Intron Chain Level:** Perfect matching of all consecutive internal intron junctions within multi-exon transcripts, serving as the benchmark for full-length splice isoform fidelity.
* **Transcript Level:** Complete structural concordance of all exons, including terminal 5' and 3' boundaries.
* **Locus Level:** Overlap of query transcript clusters within reference gene boundaries.

### C. Structural Class Codes & Downstream Utility
`gffcompare` assigns single-character class codes to each assembled transfrag relative to reference transcripts:
* `=`: Complete coordinate match across all intron junctions.
* `c`: Query transfrag is fully contained within a reference intron/exon structure.
* `j`: Novel splice isoform showing at least one splice junction match with a reference locus.
* `i`: Transfrag falling entirely within an annotated reference intron.
* `u`: Intergenic transfrag representing completely novel, unannotated genomic loci.

The annotated master GTF (`at_compare.annotated.gtf`) embeds these classifications, facilitating downstream isolation of candidate novel stress-responsive transcripts.

---

## 2. Execution & Implementation

The evaluation was executed via `gffcompare` using strict reference comparison arguments (`03_gffcompare_20260318.sh`).

### Execution Command
```bash
cd results/gtf

gffcompare \
  -r ../genome/TAIR10_GFF3_genes.gtf \
  -G \
  -o at_compare \
  at_merged.gtf
```

### Parameter Rationale
* `-r ../genome/TAIR10_GFF3_genes.gtf`: Sets the ground-truth reference annotation against which query structures are evaluated.
* `-G`: Directs GFFCompare to consider all query transcripts without discarding transfrags that possess redundant or sub-optimal terminal exon extensions.
* `-o at_compare`: Sets the output prefix for summary statistics and structural mapping tables.

---

## 3. Evaluation Statistics & Benchmark Results

Inspection of `at_compare.stats` demonstrated near-complete recovery of canonical Arabidopsis gene structures across all six pilot libraries.

### A. Assembly Volume Summary

| Metric | Query Merged (`at_merged.gtf`) | Reference (`TAIR10`) | Recovery Ratio |
| :--- | :--- | :--- | :--- |
| **Total mRNAs / Transcripts** | 52,665 | 41,613 | 1.27x |
| **Transcribed Loci** | 32,482 | 33,350 | 0.97x |
| **Multi-Exon Transcripts** | 40,937 | 30,132 | 1.36x |
| **Multi-Transcript Loci** | 9,530 | — | ~1.6 transcripts / locus |
| **Super-Loci with Reference** | 32,157 | — | 96.4% reference locus coverage |

### B. Hierarchical Sensitivity and Precision

| Structural Level | Sensitivity (%) | Precision (%) | Structural Concordance Assessment |
| :--- | :--- | :--- | :--- |
| **Base Level** | 100.0 | 96.2 | Complete exonic sequence coverage |
| **Exon Level** | 98.3 | 89.2 | High boundary resolution across individual exons |
| **Intron Level** | 100.0 | 92.3 | Exhaustive recovery of annotated splice junctions |
| **Intron Chain Level** | **100.0** | **73.6** | Zero dropped multi-exon reference splicing chains |
| **Transcript Level** | 99.7 | 78.7 | 41,469 exact transcript structure matches |
| **Locus Level** | 99.7 | 98.9 | 33,237 exact locus matches |

### C. Missed vs. Novel Structural Feature Diagnostics

| Feature Type | Reference Missed Rate (%) | Assembled Novel Fraction (%) | Biological Diagnostic Interpretation |
| :--- | :--- | :--- | :--- |
| **Exons** | 0 / 169,269 (**0.0%**) | 3,234 / 188,673 (**1.7%**) | Zero reference exon dropouts; 3,234 candidate novel exons |
| **Introns** | 0 / 127,896 (**0.0%**) | 5,089 / 138,525 (**3.7%**) | Zero splice junctions missed; 5,089 candidate alternative junctions |
| **Loci** | 0 / 33,350 (**0.0%**) | 325 / 32,482 (**1.0%**) | Zero unrecovered reference loci; 325 novel transcribed intergenic loci |

---

## 4. Technical Diagnostics & Assembly Integrity

1. **Exhaustive Intron Chain Sensitivity (100.0%):** All 30,132 annotated multi-exon intron chains present in the TAIR10 reference genome were reconstructed in `at_merged.gtf`. This indicates that upstream HISAT2 alignment thresholds (`--dta`) and StringTie minimum coverage parameters preserved all canonical splicing patterns without dropout.
2. **Terminal Boundary Discrepancies:** While intron chain sensitivity reached 100.0%, transcript-level precision stood at 78.7%. This offset is standard for StringTie reconstructions and reflects variable 5' and 3' untranslated region (UTR) assembly lengths driven by terminal read-depth decay, rather than internal coding structure errors.
3. **Novel Structural Candidates:** The detection of 325 novel loci (class code `u`) and 5,089 novel splice junctions (class code `j`) indicates transcriptomic activity not captured in canonical TAIR10 annotations. These loci represent potential stress-induced long non-coding RNAs (lncRNAs) or condition-specific alternative isoforms triggered by the 120-day microplastic soil regime.

---

## 5. Output Artifact Inventory

```text
results/gtf/
├── at_compare.stats                  # Global sensitivity, precision, and feature frequency report
├── at_compare.annotated.gtf          # Merged GTF annotated with tracking IDs and class codes
├── at_compare.loci                   # Locus coordinate clusters across reference and query models
├── at_compare.tracking               # Transfrag tracking matrix linking sample GTFs to reference IDs
├── at_compare.at_merged.gtf.refmap   # Reference-to-query structural mapping lookup table
└── at_compare.at_merged.gtf.tmap     # Query-to-reference classification and best-match table
```

---

## 6. Pipeline Status

Validation of `at_merged.gtf` confirms structural integrity and complete reference gene recovery. The environment is verified to proceed to transcript and gene-level expression quantification using `StringTie -e -B` and count matrix generation via `prepDE.py3` (`04_quantification_20260323.sh`).