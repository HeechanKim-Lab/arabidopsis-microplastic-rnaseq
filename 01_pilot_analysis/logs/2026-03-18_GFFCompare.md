# Comparison of the Merged Transcripts with TAIR10 (gffcompare)

**Date:** 2026-03-18  
**Stage:** 01_pilot_analysis  
**Script:** `modules/03_assembly_and_merge.sh` (last step)  
**Input:** `at_merged.gtf`, `genes/TAIR10_GFF3_genes.gtf`  
**Output:** `at_compare.stats`, `at_compare.annotated.gtf`, `at_compare.loci`, `at_compare.tracking`, `at_compare.at_merged.gtf.refmap`, `at_compare.at_merged.gtf.tmap`

---

## 1. Method

`gffcompare` (v0.12.10) compares the merged StringTie transcripts with the TAIR10 annotation and reports sensitivity and precision at the base, exon, intron, intron-chain, transcript and locus levels. It also assigns a class code to each transcript (in `at_compare.annotated.gtf` and the `.tmap` file).

```bash
gffcompare -r genes/TAIR10_GFF3_genes.gtf -G -o at_compare at_merged.gtf
```

* `-r`: reference annotation.
* `-G`: keep all query transcripts (intron-redundant transcripts are not discarded).
* `-o`: output prefix.

---

## 2. Results (`at_compare.stats`)

### A. Transcript and locus counts

| Item | Merged (`at_merged.gtf`) | Reference (TAIR10) |
| :--- | ---: | ---: |
| Transcripts (mRNAs) | 52,665 | 41,613 |
| Loci | 32,482 | 33,350 |
| Multi-exon transcripts | 40,937 | 30,132 |
| Loci with more than one transcript | 9,530 (~1.6 transcripts per locus) | — |
| Super-loci containing reference transcripts | 32,157 | — |

### B. Sensitivity and precision

| Level | Sensitivity (%) | Precision (%) |
| :--- | ---: | ---: |
| Base | 100.0 | 96.2 |
| Exon | 98.3 | 89.2 |
| Intron | 100.0 | 92.3 |
| Intron chain | 100.0 | 73.6 |
| Transcript | 99.7 | 78.7 |
| Locus | 99.7 | 98.9 |

Matching intron chains: 30,132. Matching transcripts: 41,469. Matching loci: 33,237.

### C. Missed and novel features

| Feature | Missed (reference) | Novel (merged set) |
| :--- | ---: | ---: |
| Exons | 0 / 169,269 (0.0%) | 3,234 / 188,673 (1.7%) |
| Introns | 0 / 127,896 (0.0%) | 5,089 / 138,525 (3.7%) |
| Loci | 0 / 33,350 (0.0%) | 325 / 32,482 (1.0%) |

All 52,665 transcripts were written to `at_compare.annotated.gtf` (0 discarded as redundant).

---

## 3. Notes

* `stringtie --merge` was run with `-G`, so the reference transcripts are included in `at_merged.gtf`. The near-100% sensitivity values are therefore expected and do not measure how well the reads reconstructed the annotation.
* Precision below 100% reflects transcripts in the merged set that do not match a reference model. The 325 novel loci and 5,089 novel introns are candidates only; they were not examined further.

Next step: quantification with `stringtie -e -B` and count matrices with `prepDE.py3` (`modules/04_quantification.sh`).
