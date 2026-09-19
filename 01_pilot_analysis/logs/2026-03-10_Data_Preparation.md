# Upstream Processing: Arabidopsis Reference Genome Standardization and HISAT2 Indexing
**Date:** March 10, 2026  
**Project:** arabidopsis-microplastic-rnaseq  
**Module:** 01_pilot_analysis  
**Target Assays:** *Arabidopsis thaliana* Root mRNA-seq (TAIR10)  
**Input Artifacts:** `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`, `TAIR10_GFF3_genes.gff`, `samples/*.fastq.gz`  
**Output Artifacts:** `Arabidopsis_thaliana_Chr.fa`, `TAIR10_GFF3_genes.gtf`, `indexes/arabidopsis_index.*.ht2`  

---

## 1. Methodological & Theoretical Rationale

### A. Coordinate Header Harmonization
Reference genome assemblies and annotation sets from different repositories often employ conflicting sequence dictionaries:
* The Ensembl Plants FASTA release utilizes bare numeric identifiers (`>1`, `>2`, `>3`, `>4`, `>5`) for nuclear chromosomes and abbreviations (`>Mt`, `>Pt`) for organellar genomes.
* The canonical TAIR10 annotation scheme utilizes prefixed chromosome names (`Chr1`, `Chr2`, `Chr3`, `Chr4`, `Chr5`, `ChrM`, `ChrC`).

Spliced aligners (`HISAT2`) and transcript assemblers (`StringTie`) require identical reference sequence identifiers between the alignment index and the reference annotation GTF. Discrepancies cause silent alignment dropouts or fatal errors during StringTie assembly, where mapped BAM records cannot be assigned to genomic loci in the annotation file.

### B. Annotation Schema Conversion
The Tuxedo protocol relies on the Gene Transfer Format (GTF 2.2) specification. The native annotation file (`TAIR10_GFF3_genes.gff`) uses hierarchical tag-value attributes (`ID=...;Parent=...`) that do not supply explicit `gene_id` and `transcript_id` groupings required by `StringTie` and `prepDE.py3`. Converting GFF3 to GTF using `gffread -T` unpacks parent-child relationships into compliant `gene_id` and `transcript_id` attributes while preserving exon coordinates and coding sequences.

### C. Hierarchical Graph FM Indexing
`hisat2-build` constructs a Hierarchical Graph FM index (HGFM) based on the Burrows-Wheeler Transform (BWT) and Ferragina-Manzini (FM) index. HISAT2 utilizes two index tiers:
1. A global FM index representing the entire Arabidopsis genome.
2. ~55,000 overlapping local FM indices spanning 48–64 kb genomic windows.

This dual-tier index supports rapid spliced alignment across eukaryotic intron-exon junctions with memory consumption under 2 GB for the Arabidopsis genome (~135 Mb).

---

## 2. Execution & Implementation

### A. Annotation File Conversion
The primary GFF3 annotation was converted to GTF format using `gffread`:

```bash
gffread TAIR10_GFF3_genes.gff -T -o TAIR10_GFF3_genes.gtf
```

### B. FASTA Header Normalization
POSIX `sed` stream editor operations were applied to standardize chromosome names to match the GTF annotations. Numeric headers `1` through `5` were prefixed with `Chr`, and organellar genomes `Mt` and `Pt` were mapped to `ChrM` and `ChrC`:

```bash
# Prefix nuclear chromosomes 1-5 with Chr
sed 's/^>\([1-5]\)/>Chr\1/' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa > Arabidopsis_thaliana_Chr.fa

# Normalize mitochondrial and plastid identifiers (BSD sed syntax on macOS)
sed -i '' 's/^>Mt/>ChrM/' Arabidopsis_thaliana_Chr.fa
sed -i '' 's/^>Pt/>ChrC/' Arabidopsis_thaliana_Chr.fa
```

### C. Genome Index Construction
The index was built across 6 computational threads on an Apple Silicon host:

```bash
mkdir -p indexes
hisat2-build -p 6 Arabidopsis_thaliana_Chr.fa indexes/arabidopsis_index
```

The indexing process generated 8 binary index structures (`arabidopsis_index.1.ht2` through `arabidopsis_index.8.ht2`), comprising 4 small-bucket forward indices and 4 mirrored/reversed indices.

---

## 3. Data Integrity & Verification

### A. Archive Integrity Verification
All raw paired-end FASTQ archives from the sequencing facility were validated using `gzip -tv`. This tests the CRC-32 checksum and DEFLATE stream integrity to identify corrupted read pairs prior to compute execution:

```bash
gzip -tv samples/*.fastq.gz
```

### B. Verification Inventory

| File Identifier | Format | Validation Operation | Metric / Checksum Result | Status |
| :--- | :--- | :--- | :--- | :--- |
| `Arabidopsis_thaliana_Chr.fa` | FASTA | `grep "^>" Arabidopsis_thaliana_Chr.fa` | `Chr1, Chr2, Chr3, Chr4, Chr5, ChrM, ChrC` | Validated |
| `TAIR10_GFF3_genes.gtf` | GTF 2.2 | `head -n 5 TAIR10_GFF3_genes.gtf` | Valid `gene_id` & `transcript_id` fields | Validated |
| `arabidopsis_index.*.ht2` | HGFM | `ls -lh indexes/` | 8 `.ht2` files generated (~120 MB total) | Complete |
| `R1-1_S0_L009_R*.fastq.gz` | FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |
| `R2-1_S252_L004_R*.fastq.gz`| FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |
| `R3-1_S261_L004_R*.fastq.gz`| FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |
| `R1-9_S251_L004_R*.fastq.gz`| FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |
| `R2-9_S0_L009_R*.fastq.gz` | FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |
| `R3-9_S0_L009_R*.fastq.gz` | FASTQ.GZ | `gzip -tv` | CRC-32 integrity verified: OK | Uncorrupted |

---

## 4. Pipeline Status

Preprocessing is complete. The synchronized FASTA and GTF reference models are staged, and all 6 paired-end pilot libraries are uncorrupted. The environment is prepared for spliced read alignment using HISAT2 (`01_hisat2_alignment_20260313.sh`).