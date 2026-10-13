# Read Alignment with HISAT2

**Date:** 2026-03-13  
**Stage:** 01_pilot_analysis  
**Script:** `modules/02_hisat2_alignment.sh`  
**Input:** `indexes/arabidopsis_index.*.ht2`, `samples/*.fastq.gz`  
**Output:** one SAM file per library

---

## 1. Method

Paired-end reads were aligned to the TAIR10 genome with HISAT2. The `--dta` option reports alignments tailored for transcript assemblers such as StringTie (longer anchors are required for new splice sites), which is the recommended setting before StringTie assembly.

Libraries:

* **WT (0% plastic control):** `R1-1_WT`, `R2-1_WT`, `R3-1_WT`
* **Cond (5% plastic, soil pre-incubated for 120 days):** `R1-9_Cond`, `R2-9_Cond`, `R3-9_Cond`

---

## 2. Commands

One command per library (10 threads), for example:

```bash
hisat2 -p 10 --dta -x indexes/arabidopsis_index \
  -1 samples/R1-1_S0_L009_R1_001.fastq.gz \
  -2 samples/R1-1_S0_L009_R2_001.fastq.gz \
  -S R1-1_WT.sam
```

The other five libraries were run the same way with their own FASTQ files (`R2-1_S252_L004`, `R3-1_S261_L004`, `R1-9_S251_L004`, `R2-9_S0_L009`, `R3-9_S0_L009`).

---

## 3. Alignment summary

| Sample | Group | Read pairs | Concordant 0 times | Concordant, unique | Concordant, multiple | Discordant (% of non-concordant pairs) | Overall alignment |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| R1-1_WT | WT | 29,791,679 | 6,046,511 (20.30%) | 22,759,510 (76.40%) | 985,658 (3.31%) | 338,280 (5.59%) | 88.56% |
| R2-1_WT | WT | 30,293,327 | 3,896,031 (12.86%) | 25,675,061 (84.75%) | 722,235 (2.38%) | 481,753 (12.37%) | 93.35% |
| R3-1_WT | WT | 29,526,747 | 5,252,367 (17.79%) | 23,620,211 (80.00%) | 654,169 (2.22%) | 623,511 (11.87%) | 90.88% |
| R1-9_Cond | Cond | 29,061,578 | 3,758,225 (12.93%) | 24,757,144 (85.19%) | 546,209 (1.88%) | 466,708 (12.42%) | 93.78% |
| R2-9_Cond | Cond | 29,817,819 | 4,965,890 (16.65%) | 24,251,405 (81.33%) | 600,524 (2.01%) | 604,108 (12.17%) | 91.92% |
| R3-9_Cond | Cond | 29,899,345 | 3,999,998 (13.38%) | 25,397,944 (84.94%) | 501,403 (1.68%) | 563,225 (14.08%) | 93.72% |

---

## 4. Observations

* All six libraries have 29.1–30.3 million read pairs.
* Five libraries have overall alignment rates of 90.88–93.78% and unique concordant rates of 80.00–85.19%.
* `R1-1_WT` has the lowest overall alignment rate (88.56%) and unique concordant rate (76.40%), and the highest fraction of pairs that did not align concordantly (20.30%; others 12.86–17.79%). This was noted for the QC steps that follow.
* File names place `R1-1_WT`, `R2-9_Cond` and `R3-9_Cond` in one group (`S0_L009`) and the other three libraries in another (`L004`).

---

## 5. HISAT2 output

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

Next step: sorting with samtools and transcript assembly with StringTie (`modules/03_assembly_and_merge.sh`).
