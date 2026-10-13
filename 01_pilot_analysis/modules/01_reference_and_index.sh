#!/usr/bin/env bash
# ==============================================================================
# Module 01: TAIR10 reference files and HISAT2 index (log 2026-03-10)
#
# Input : data/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa (Ensembl Plants)
#         data/genome/TAIR10_GFF3_genes.gff
# Output: data/genome/TAIR10_GFF3_genes.gtf
#         data/genome/Arabidopsis_thaliana_Chr.fa (chromosome names as in the GFF)
#         data/indexes/arabidopsis_index.*.ht2
# ==============================================================================

set -euo pipefail

THREADS="${THREADS:-6}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME_DIR="${BASE_DIR}/data/genome"
INDEX_DIR="${BASE_DIR}/data/indexes"

RAW_FASTA="${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
CHR_FASTA="${GENOME_DIR}/Arabidopsis_thaliana_Chr.fa"
RAW_GFF="${GENOME_DIR}/TAIR10_GFF3_genes.gff"
ANNOTATION_GTF="${GENOME_DIR}/TAIR10_GFF3_genes.gtf"
INDEX_PREFIX="${INDEX_DIR}/arabidopsis_index"

for f in "${RAW_FASTA}" "${RAW_GFF}"; do
    [[ -f "${f}" ]] || { echo "Missing input: ${f}" >&2; exit 1; }
done
mkdir -p "${INDEX_DIR}"

# 1. GFF3 -> GTF (StringTie and prepDE.py3 read the gene_id and transcript_id attributes)
if [[ ! -f "${ANNOTATION_GTF}" ]]; then
    echo "[01] Converting GFF3 to GTF"
    gffread "${RAW_GFF}" -T -o "${ANNOTATION_GTF}"
fi

# 2. Rename chromosomes to match the GFF: 1-5 -> Chr1-Chr5, Mt -> ChrM, Pt -> ChrC
if [[ ! -f "${CHR_FASTA}" ]]; then
    echo "[01] Renaming chromosomes"
    sed -e 's/^>\([1-5]\)/>Chr\1/' -e 's/^>Mt/>ChrM/' -e 's/^>Pt/>ChrC/' "${RAW_FASTA}" > "${CHR_FASTA}"
fi
grep -m 5 ">" "${CHR_FASTA}"

# 3. HISAT2 index
if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    echo "[01] Building the HISAT2 index"
    hisat2-build -p "${THREADS}" "${CHR_FASTA}" "${INDEX_PREFIX}"
fi

echo "[01] Done: ${INDEX_PREFIX}"
