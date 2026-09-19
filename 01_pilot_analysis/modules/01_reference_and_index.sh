#!/usr/bin/env bash
# ==============================================================================
# Module: 01_reference_and_index.sh
# Description: Reference genome header harmonization and HISAT2 HGFM indexing
# ==============================================================================

set -euo pipefail

THREADS=10
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME_DIR="${BASE_DIR}/data/genome"
INDEX_DIR="${BASE_DIR}/data/indexes"

RAW_FASTA="${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
STANDARDIZED_FASTA="${GENOME_DIR}/Arabidopsis_thaliana_Chr.fa"
RAW_GFF="${GENOME_DIR}/TAIR10_GFF3_genes.gff"
ANNOTATION_GTF="${GENOME_DIR}/TAIR10_GFF3_genes.gtf"
INDEX_PREFIX="${INDEX_DIR}/arabidopsis_index"

mkdir -p "${INDEX_DIR}"

echo ">>> [01] Converting GFF3 to GTF..."
if [[ ! -f "${ANNOTATION_GTF}" ]]; then
    gffread "${RAW_GFF}" -T -o "${ANNOTATION_GTF}"
fi

echo ">>> [01] Standardizing chromosome headers..."
if [[ ! -f "${STANDARDIZED_FASTA}" ]]; then
    sed 's/^>\([1-5]\)/>Chr\1/' "${RAW_FASTA}" > "${STANDARDIZED_FASTA}"
    sed -i '' 's/^>Mt/>ChrM/' "${STANDARDIZED_FASTA}"
    sed -i '' 's/^>Pt/>ChrC/' "${STANDARDIZED_FASTA}"
fi

echo ">>> [01] Constructing HISAT2 index..."
if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    hisat2-build -p "${THREADS}" "${STANDARDIZED_FASTA}" "${INDEX_PREFIX}"
fi

echo ">>> [01] Indexing complete: ${INDEX_PREFIX}"