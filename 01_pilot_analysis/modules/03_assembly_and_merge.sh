#!/usr/bin/env bash
# ==============================================================================
# Module: 03_assembly_and_merge.sh
# Description: BAM sorting, StringTie transfrag assembly, merge, and GFFCompare
# ==============================================================================

set -euo pipefail

THREADS=10
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOTATION_GTF="${BASE_DIR}/data/genome/TAIR10_GFF3_genes.gtf"
SAM_DIR="${BASE_DIR}/results/sam"
BAM_DIR="${BASE_DIR}/results/bam"
GTF_DIR="${BASE_DIR}/results/gtf"

mkdir -p "${BAM_DIR}" "${GTF_DIR}"

SAMPLES=("R1-1_WT" "R2-1_WT" "R3-1_WT" "R1-9_Cond" "R2-9_Cond" "R3-9_Cond")

# 1. BAM Sorting & Indexing
for SAMPLE in "${SAMPLES[@]}"; do
    SAM_IN="${SAM_DIR}/${SAMPLE}.sam"
    BAM_OUT="${BAM_DIR}/${SAMPLE}.bam"
    if [[ ! -f "${BAM_OUT}" ]]; then
        echo ">>> [03] Sorting ${SAMPLE}..."
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}" "${SAM_IN}"
        samtools index "${BAM_OUT}"
    fi
done

# 2. StringTie Transcript Assembly
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_IN="${BAM_DIR}/${SAMPLE}.bam"
    GTF_OUT="${GTF_DIR}/${SAMPLE}.gtf"
    if [[ ! -f "${GTF_OUT}" ]]; then
        echo ">>> [03] Assembling ${SAMPLE}..."
        stringtie -p "${THREADS}" \
            -G "${ANNOTATION_GTF}" \
            -o "${GTF_OUT}" \
            -l "${SAMPLE}" \
            "${BAM_IN}"
    fi
done

# 3. Merge Transcripts
MERGELIST="${GTF_DIR}/mergelist_at.txt"
MERGED_GTF="${GTF_DIR}/at_merged.gtf"

find "${GTF_DIR}" -name "R*-*.gtf" | sort > "${MERGELIST}"

if [[ ! -f "${MERGED_GTF}" ]]; then
    echo ">>> [03] Merging assembly models..."
    stringtie --merge -p "${THREADS}" \
        -G "${ANNOTATION_GTF}" \
        -o "${MERGED_GTF}" \
        "${MERGELIST}"
fi

# 4. GFFCompare Validation
COMPARE_PREFIX="${GTF_DIR}/at_compare"
if [[ ! -f "${COMPARE_PREFIX}.stats" ]]; then
    echo ">>> [03] Running GFFCompare evaluation..."
    gffcompare -r "${ANNOTATION_GTF}" -G -o "${COMPARE_PREFIX}" "${MERGED_GTF}"
fi