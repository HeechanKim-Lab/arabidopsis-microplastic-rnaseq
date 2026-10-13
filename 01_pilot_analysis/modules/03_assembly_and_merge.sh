#!/usr/bin/env bash
# ==============================================================================
# Module 03: sorting, StringTie assembly, merge and gffcompare
#            (logs 2026-03-17 and 2026-03-18)
#
# Input : results/sam/<sample>.sam (module 02)
#         data/genome/TAIR10_GFF3_genes.gtf (module 01)
# Output: results/bam/<sample>.bam
#         results/gtf/<sample>.gtf, results/gtf/mergelist_at.txt, results/gtf/at_merged.gtf
#         results/gtf/at_compare.* (gffcompare)
# ==============================================================================

set -euo pipefail

THREADS="${THREADS:-10}"
SORT_THREADS="${SORT_THREADS:-6}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOTATION_GTF="${BASE_DIR}/data/genome/TAIR10_GFF3_genes.gtf"
SAM_DIR="${BASE_DIR}/results/sam"
BAM_DIR="${BASE_DIR}/results/bam"
GTF_DIR="${BASE_DIR}/results/gtf"
mkdir -p "${BAM_DIR}" "${GTF_DIR}"

SAMPLES="R1-1_WT R2-1_WT R3-1_WT R1-9_Cond R2-9_Cond R3-9_Cond"

# 1. Coordinate sorting
for SAMPLE in ${SAMPLES}; do
    BAM="${BAM_DIR}/${SAMPLE}.bam"
    if [[ ! -f "${BAM}" ]]; then
        echo "[03] Sorting ${SAMPLE}"
        samtools sort -@ "${SORT_THREADS}" -o "${BAM}" "${SAM_DIR}/${SAMPLE}.sam"
    fi
done

# 2. Reference-guided assembly of each library
for SAMPLE in ${SAMPLES}; do
    GTF="${GTF_DIR}/${SAMPLE}.gtf"
    if [[ ! -f "${GTF}" ]]; then
        echo "[03] Assembling ${SAMPLE}"
        stringtie -p "${THREADS}" -G "${ANNOTATION_GTF}" -o "${GTF}" -l "${SAMPLE}" "${BAM_DIR}/${SAMPLE}.bam"
    fi
done

# 3. Merge the six assemblies (same file pattern as the original 'ls R*-*_*.gtf')
MERGELIST="${GTF_DIR}/mergelist_at.txt"
MERGED_GTF="${GTF_DIR}/at_merged.gtf"
printf '%s\n' "${GTF_DIR}"/R*-*_*.gtf > "${MERGELIST}"

if [[ ! -f "${MERGED_GTF}" ]]; then
    echo "[03] Merging assemblies"
    stringtie --merge -p "${THREADS}" -G "${ANNOTATION_GTF}" -o "${MERGED_GTF}" "${MERGELIST}"
fi

# 4. Comparison of the merged transcripts with TAIR10
if [[ ! -f "${GTF_DIR}/at_compare.stats" ]]; then
    echo "[03] Running gffcompare"
    gffcompare -r "${ANNOTATION_GTF}" -G -o "${GTF_DIR}/at_compare" "${MERGED_GTF}"
fi
cat "${GTF_DIR}/at_compare.stats"
