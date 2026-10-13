# Sorting, Transcript Assembly and Merge

**Date:** 2026-03-17  
**Stage:** 01_pilot_analysis  
**Script:** `modules/03_assembly_and_merge.sh`  
**Input:** one SAM file per library, `genes/TAIR10_GFF3_genes.gtf`  
**Output:** sorted BAM files, one GTF per library, `mergelist_at.txt`, `at_merged.gtf`

---

## 1. Method

* SAM files were converted to coordinate-sorted BAM files with `samtools sort`.
* Each library was assembled with StringTie, guided by the TAIR10 annotation (`-G`).
* The six assemblies were merged with `stringtie --merge` into one transcript set (`at_merged.gtf`), which was then used to quantify every library against the same transcript models.

---

## 2. Commands

### A. Sorting (6 threads)

```bash
samtools sort -@ 6 -o R1-1_WT.bam R1-1_WT.sam
samtools sort -@ 6 -o R2-1_WT.bam R2-1_WT.sam
samtools sort -@ 6 -o R3-1_WT.bam R3-1_WT.sam
samtools sort -@ 6 -o R1-9_Cond.bam R1-9_Cond.sam
samtools sort -@ 6 -o R2-9_Cond.bam R2-9_Cond.sam
samtools sort -@ 6 -o R3-9_Cond.bam R3-9_Cond.sam
```

### B. Assembly per library (10 threads)

```bash
stringtie -p 10 -G genes/TAIR10_GFF3_genes.gtf -o R1-1_WT.gtf -l R1-1_WT R1-1_WT.bam
# same command for R2-1_WT, R3-1_WT, R1-9_Cond, R2-9_Cond and R3-9_Cond
```

### C. Merge

```bash
ls R*-*_*.gtf > mergelist_at.txt
stringtie --merge -p 10 -G genes/TAIR10_GFF3_genes.gtf -o at_merged.gtf mergelist_at.txt
```

`mergelist_at.txt` lists the six assemblies: `R1-1_WT.gtf`, `R2-1_WT.gtf`, `R3-1_WT.gtf`, `R1-9_Cond.gtf`, `R2-9_Cond.gtf` and `R3-9_Cond.gtf`.

A known locus was checked in the annotation as a spot test:

```bash
grep "AT1G01010" TAIR10_GFF3_genes.gtf
```

---

## 3. Merged transcript set

Counts from the `gffcompare` summary (see the 2026-03-18 log):

| Item | Count |
| :--- | ---: |
| Libraries merged | 6 (3 WT, 3 Cond) |
| Transcripts | 52,665 |
| Loci | 32,482 |
| Multi-exon transcripts | 40,937 |
| Loci with more than one transcript | 9,530 (~1.6 transcripts per locus) |

Next step: comparison of `at_merged.gtf` with the TAIR10 annotation using `gffcompare` (also in `modules/03_assembly_and_merge.sh`).
