# Reference Preparation and HISAT2 Index

**Date:** 2026-03-10  
**Stage:** 01_pilot_analysis  
**Script:** `modules/01_reference_and_index.sh`  
**Input:** `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`, `TAIR10_GFF3_genes.gff`, `samples/*.fastq.gz`  
**Output:** `Arabidopsis_thaliana_Chr.fa`, `TAIR10_GFF3_genes.gtf`, `arabidopsis_index.*.ht2`

---

## 1. Purpose

* StringTie and `prepDE.py3` read the `gene_id` and `transcript_id` attributes of a GTF file, so the TAIR10 GFF3 annotation was converted to GTF with `gffread -T`.
* The genome FASTA (Ensembl Plants release) names the chromosomes `1`–`5`, `Mt` and `Pt`, while the TAIR10 annotation uses `Chr1`–`Chr5`, `ChrM` and `ChrC`. The FASTA headers were renamed so that the alignment index and the annotation use the same sequence names.

---

## 2. Commands

### A. Annotation conversion

```bash
gffread TAIR10_GFF3_genes.gff -T -o TAIR10_GFF3_genes.gtf
```

### B. Chromosome names

```bash
# Check the names used in each file
grep ">" Arabidopsis_thaliana.TAIR10.dna.toplevel.fa | head -n 5   # >1 ... >5
cut -f1 TAIR10_GFF3_genes.gff | sort | uniq | head -n 5             # Chr1 ... Chr5

# Rename the FASTA headers (BSD sed syntax on macOS)
sed 's/^>\([1-5]\)/>Chr\1/' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa > Arabidopsis_thaliana_Chr.fa
sed -i '' 's/^>Mt/>ChrM/' Arabidopsis_thaliana_Chr.fa
sed -i '' 's/^>Pt/>ChrC/' Arabidopsis_thaliana_Chr.fa
```

After renaming, the headers read `>Chr1 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF` through `>Chr5 dna:chromosome chromosome:TAIR10:5:1:26975502:1 REF`.

### C. HISAT2 index

```bash
hisat2-build -p 6 Arabidopsis_thaliana_Chr.fa arabidopsis_index
```

The index files (`arabidopsis_index.1.ht2` to `arabidopsis_index.8.ht2`) were read from `indexes/` in the alignment step.

### D. FASTQ integrity

The compressed FASTQ files were tested with `gzip -tv`, for example:

```bash
gzip -tv R2-1_S252_L004_R1_001.fastq.gz
```

---

## 3. Pilot libraries

| Library | FASTQ files (R1 and R2) |
| :--- | :--- |
| R1-1_WT | `R1-1_S0_L009_R*_001.fastq.gz` |
| R2-1_WT | `R2-1_S252_L004_R*_001.fastq.gz` |
| R3-1_WT | `R3-1_S261_L004_R*_001.fastq.gz` |
| R1-9_Cond | `R1-9_S251_L004_R*_001.fastq.gz` |
| R2-9_Cond | `R2-9_S0_L009_R*_001.fastq.gz` |
| R3-9_Cond | `R3-9_S0_L009_R*_001.fastq.gz` |

Next step: read alignment with HISAT2 (`modules/02_hisat2_alignment.sh`).
