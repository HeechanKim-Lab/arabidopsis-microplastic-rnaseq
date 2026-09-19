# Upstream Processing: Coordinate Sorting, Transcript Assembly, and Master Transcriptome Merging
**Date:** March 17, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `results/sam/*.sam`, `genome/TAIR10_GFF3_genes.gtf`  
**Output Artifacts:** `results/bam/*.bam`, `results/gtf/*.gtf`, `results/gtf/mergelist_at.txt`, `results/gtf/at_merged.gtf`  

---

## 1. Methodological & Theoretical Rationale

### A. SAM-to-BAM Conversion and Coordinate Sorting
Spliced read alignments produced by HISAT2 are emitted as text-based Sequence Alignment/Map (SAM) records ordered arbitrarily by read ingestion order. Downstream transcript assemblers require binary, coordinate-sorted inputs:
* **Storage and I/O Optimization:** The Binary Alignment/Map (BAM) format utilizes Block GZIP Format (BGZF) compression, reducing raw SAM disk footprint by ~75% while permitting random-access decompression.
* **Coordinate Ordering:** `samtools sort` indexes alignments by reference chromosome identifier and leftmost 5' genomic coordinate. Sorting organizes spliced reads linearly along chromosomes, allowing transcript assemblers to process discrete loci sequentially without caching the entire alignment in memory.

### B. Network Flow-Based Transcript Assembly (StringTie)
`StringTie` reconstructs transcript models from spliced read alignments using a maximum flow network optimization strategy:
1. **Splice Graph Construction:** Mapped spliced reads define directed acyclic graphs (splice graphs) where nodes correspond to assembled exons and directed edges correspond to junction reads spanning introns (supported by HISAT2 `XS:A` junction tags).
2. **Maximum Flow Optimization:** A network flow algorithm calculates the minimum set of non-overlapping transcript paths and estimates their relative abundance (expression levels) that parsimoniously accounts for read depth across the splice graph.
3. **Reference Guidance (`-G`):** Supplying `TAIR10_GFF3_genes.gtf` guides graph resolution toward annotated gene models while allowing `StringTie` to resolve novel splice junctions, condition-specific isoforms, and previously unannotated transcripts induced by environmental microplastic stress.

### C. Unified Reference Consolidation (`--merge`)
Transcript assemblies derived from individual replicates (`R1-1_WT` through `R3-9_Cond`) exhibit structural fragmentation: lowly expressed transcripts may assemble partially in one library and fail detection in another.
* **Non-Redundant Master Reference:** `stringtie --merge` ingests the manifest of all replicate GTFs (`mergelist_at.txt`) and consolidates overlapping, redundant, or fragmented transfrags into a single non-redundant reference transcriptome (`at_merged.gtf`).
* **Cross-Replicate Uniformity:** Merging across both Wild Type and Conditioned cohorts establishes a synchronized annotation coordinate system. This uniform framework is strictly required for downstream abundance quantification (Ballgown) and differential expression modeling (edgeR), preventing reference bias against condition-specific transcripts.

---

## 2. Execution & Implementation

The assembly workflow was automated using a sequential shell pipeline (`02_assembly_and_merge_20260317.sh`).

### A. Coordinate Sorting via Samtools
Executed across 6 CPU threads per sample to preserve headroom on the Apple Silicon host:

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
    samtools sort -@ 6 -o results/bam/${NAME}.bam results/sam/${NAME}.sam
done
```

### B. Individual Transcript Assembly
Assembled transcripts for each library against the standardized TAIR10 GTF across 10 threads:

```bash
for NAME in "${SAMPLES[@]}"; do
    stringtie -p 10 \
        -G genome/TAIR10_GFF3_genes.gtf \
        -o results/gtf/${NAME}.gtf \
        -l ${NAME} \
        results/bam/${NAME}.bam
done
```

### C. Master Transcriptome Merging
Compiled the manifest and merged all sample assemblies into a master GTF:

```bash
cd results/gtf
ls R*-*_*.gtf > mergelist_at.txt

stringtie --merge -p 10 \
    -G ../genome/TAIR10_GFF3_genes.gtf \
    -o at_merged.gtf \
    mergelist_at.txt
```

---

## 3. Structural Validation & Assembly Inventory

Integrity checks were performed on intermediate manifest files and the consolidated transcriptome `at_merged.gtf`.

### A. Manifest Verification
Inspection of `mergelist_at.txt` confirmed that all 6 biological replicates across experimental groups were included:
* `R1-1_WT.gtf` (Control Replicate 1)
* `R2-1_WT.gtf` (Control Replicate 2)
* `R3-1_WT.gtf` (Control Replicate 3)
* `R1-9_Cond.gtf` (Conditioned Replicate 1)
* `R2-9_Cond.gtf` (Conditioned Replicate 2)
* `R3-9_Cond.gtf` (Conditioned Replicate 3)

### B. Canonical Locus Validation
Verified the retention and structural preservation of representative canonical genes (e.g., *AT1G01010*):

```bash
grep "AT1G01010" genome/TAIR10_GFF3_genes.gtf
```

### C. Assembly Output Metrics
Total transcript counts in the merged reference were audited using string matching:

```bash
grep -c "transcript" results/gtf/at_merged.gtf
```

| Metric Category | Count / Value | Description |
| :--- | :--- | :--- |
| **Input Assembly Models** | 6 libraries | 3 Wild Type + 3 Conditioned |
| **Consolidated Transcripts** | 52,665 mRNAs | Total assembled isoforms (annotated + candidate novel) |
| **Consolidated Loci** | 32,482 loci | Unique transcribed genomic loci |
| **Multi-Exon Transcripts** | 40,937 transcripts | Spliced transcript structures |
| **Multi-Transcript Loci** | 9,530 loci | Complex loci exhibiting alternative isoform usage (~1.6 transcripts/locus) |

---

## 4. Pipeline Status

Transcript assembly and master reference consolidation are complete. The output master transcriptome (`at_merged.gtf`) is staged for structural sensitivity and precision benchmarking against the TAIR10 reference using `gffcompare` (`03_gffcompare_20260318.sh`).