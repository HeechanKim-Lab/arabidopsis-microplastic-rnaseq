#!/usr/bin/env bash
# ==============================================================================
# Script: run_upstream_tuxedo.sh
# Project: arabidopsis-microplastic-rnaseq (01_pilot_analysis)
# Description: End-to-end upstream pipeline (HISAT2 -> StringTie -> prepDE)
# ==============================================================================

set -euo pipefail

# --- Configuration & Hardware Tuning ---
THREADS=10
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${BASE_DIR}/data"
GENOME_DIR="${DATA_DIR}/genome"
INDEX_DIR="${DATA_DIR}/indexes"
SAMPLES_DIR="${DATA_DIR}/samples"
RESULTS_DIR="${BASE_DIR}/results"
SAM_DIR="${RESULTS_DIR}/sam"
BAM_DIR="${RESULTS_DIR}/bam"
GTF_DIR="${RESULTS_DIR}/gtf"
BALLGOWN_DIR="${RESULTS_DIR}/ballgown"
MATRICES_DIR="${RESULTS_DIR}/matrices"

# Reference file paths
RAW_FASTA="${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
STANDARDIZED_FASTA="${GENOME_DIR}/Arabidopsis_thaliana_Chr.fa"
RAW_GFF="${GENOME_DIR}/TAIR10_GFF3_genes.gff"
ANNOTATION_GTF="${GENOME_DIR}/TAIR10_GFF3_genes.gtf"
INDEX_PREFIX="${INDEX_DIR}/arabidopsis_index"

# Sample Manifest (Pilot 6-library cohort)
SAMPLES=(
    "R1-1_WT"
    "R2-1_WT"
    "R3-1_WT"
    "R1-9_Cond"
    "R2-9_Cond"
    "R3-9_Cond"
)

# Explicit FASTQ mappings
declare -A FASTQ_R1
declare -A FASTQ_R2

FASTQ_R1["R1-1_WT"]="${SAMPLES_DIR}/R1-1_S0_L009_R1_001.fastq.gz"
FASTQ_R2["R1-1_WT"]="${SAMPLES_DIR}/R1-1_S0_L009_R2_001.fastq.gz"

FASTQ_R1["R2-1_WT"]="${SAMPLES_DIR}/R2-1_S252_L004_R1_001.fastq.gz"
FASTQ_R2["R2-1_WT"]="${SAMPLES_DIR}/R2-1_S252_L004_R2_001.fastq.gz"

FASTQ_R1["R3-1_WT"]="${SAMPLES_DIR}/R3-1_S261_L004_R1_001.fastq.gz"
FASTQ_R2["R3-1_WT"]="${SAMPLES_DIR}/R3-1_S261_L004_R2_001.fastq.gz"

FASTQ_R1["R1-9_Cond"]="${SAMPLES_DIR}/R1-9_S251_L004_R1_001.fastq.gz"
FASTQ_R2["R1-9_Cond"]="${SAMPLES_DIR}/R1-9_S251_L004_R2_001.fastq.gz"

FASTQ_R1["R2-9_Cond"]="${SAMPLES_DIR}/R2-9_S0_L009_R1_001.fastq.gz"
FASTQ_R2["R2-9_Cond"]="${SAMPLES_DIR}/R2-9_S0_L009_R2_001.fastq.gz"

FASTQ_R1["R3-9_Cond"]="${SAMPLES_DIR}/R3-9_S0_L009_R1_001.fastq.gz"
FASTQ_R2["R3-9_Cond"]="${SAMPLES_DIR}/R3-9_S0_L009_R2_001.fastq.gz"

echo "=== [STAGE 0] Environment Initialization ==="
mkdir -p "${GENOME_DIR}" "${INDEX_DIR}" "${SAM_DIR}" "${BAM_DIR}" "${GTF_DIR}" "${BALLGOWN_DIR}" "${MATRICES_DIR}"

# Verify sample FASTQ integrity
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ ! -f "${FASTQ_R1[$SAMPLE]}" ]] || [[ ! -f "${FASTQ_R2[$SAMPLE]}" ]]; then
        echo "Error: FASTQ files for ${SAMPLE} missing in ${SAMPLES_DIR}." >&2
        exit 1
    fi
    gzip -tv "${FASTQ_R1[$SAMPLE]}"
    gzip -tv "${FASTQ_R2[$SAMPLE]}"
done

echo "=== [STAGE 1] Reference Standardization and HISAT2 Indexing ==="
# 1. Convert GFF3 to GTF if absent
if [[ ! -f "${ANNOTATION_GTF}" ]]; then
    gffread "${RAW_GFF}" -T -o "${ANNOTATION_GTF}"
fi

# 2. Harmonize FASTA chromosome headers (>1 -> >Chr1, >Mt -> >ChrM, >Pt -> >ChrC)
if [[ ! -f "${STANDARDIZED_FASTA}" ]]; then
    sed 's/^>\([1-5]\)/>Chr\1/' "${RAW_FASTA}" > "${STANDARDIZED_FASTA}"
    sed -i '' 's/^>Mt/>ChrM/' "${STANDARDIZED_FASTA}"
    sed -i '' 's/^>Pt/>ChrC/' "${STANDARDIZED_FASTA}"
fi

# 3. Build HISAT2 Index
if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    hisat2-build -p "${THREADS}" "${STANDARDIZED_FASTA}" "${INDEX_PREFIX}"
fi

echo "=== [STAGE 2] Spliced Alignment via HISAT2 ==="
for SAMPLE in "${SAMPLES[@]}"; do
    SAM_OUT="${SAM_DIR}/${SAMPLE}.sam"
    if [[ ! -f "${SAM_OUT}" ]]; then
        echo "Aligning ${SAMPLE}..."
        hisat2 -p "${THREADS}" --dta \
            -x "${INDEX_PREFIX}" \
            -1 "${FASTQ_R1[$SAMPLE]}" \
            -2 "${FASTQ_R2[$SAMPLE]}" \
            -S "${SAM_OUT}"
    else
        echo "Alignment for ${SAMPLE} already exists. Skipping."
    fi
done

echo "=== [STAGE 3] Coordinate Sorting & BAM Conversion ==="
for SAMPLE in "${SAMPLES[@]}"; do
    SAM_IN="${SAM_DIR}/${SAMPLE}.sam"
    BAM_OUT="${BAM_DIR}/${SAMPLE}.bam"
    if [[ ! -f "${BAM_OUT}" ]]; then
        echo "Sorting ${SAMPLE}..."
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}" "${SAM_IN}"
        samtools index "${BAM_OUT}"
    fi
done

echo "=== [STAGE 4] Transcript Assembly via StringTie ==="
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_IN="${BAM_DIR}/${SAMPLE}.bam"
    GTF_OUT="${GTF_DIR}/${SAMPLE}.gtf"
    if [[ ! -f "${GTF_OUT}" ]]; then
        echo "Assembling transcripts for ${SAMPLE}..."
        stringtie -p "${THREADS}" \
            -G "${ANNOTATION_GTF}" \
            -o "${GTF_OUT}" \
            -l "${SAMPLE}" \
            "${BAM_IN}"
    fi
done

echo "=== [STAGE 5] Transcriptome Merge & Evaluation ==="
MERGELIST="${GTF_DIR}/mergelist_at.txt"
MERGED_GTF="${GTF_DIR}/at_merged.gtf"
COMPARE_PREFIX="${GTF_DIR}/at_compare"

find "${GTF_DIR}" -name "R*-*.gtf" | sort > "${MERGELIST}"

if [[ ! -f "${MERGED_GTF}" ]]; then
    stringtie --merge -p "${THREADS}" \
        -G "${ANNOTATION_GTF}" \
        -o "${MERGED_GTF}" \
        "${MERGELIST}"
fi

if [[ ! -f "${COMPARE_PREFIX}.stats" ]]; then
    gffcompare -r "${ANNOTATION_GTF}" -G -o "${COMPARE_PREFIX}" "${MERGED_GTF}"
fi

echo "=== [STAGE 6] Ballgown Abundance Estimation ==="
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_BALLGOWN_DIR="${BALLGOWN_DIR}/${SAMPLE}"
    BAM_IN="${BAM_DIR}/${SAMPLE}.bam"
    mkdir -p "${SAMPLE_BALLGOWN_DIR}"
    
    if [[ ! -f "${SAMPLE_BALLGOWN_DIR}/t_data.ctab" ]]; then
        echo "Quantifying abundances for ${SAMPLE}..."
        stringtie -e -B -p "${THREADS}" \
            -G "${MERGED_GTF}" \
            -o "${SAMPLE_BALLGOWN_DIR}/${SAMPLE}.gtf" \
            "${BAM_IN}"
    fi
done

echo "=== [STAGE 7] Count Matrix Extraction & Metadata Compilation ==="
cd "${MATRICES_DIR}"

if [[ ! -f "prepDE.py3" ]]; then
    curl -O http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
fi

python3 prepDE.py3 -i "${BALLGOWN_DIR}"

printf "ids,condition\nR1-1_WT,WT\nR2-1_WT,WT\nR3-1_WT,WT\nR1-9_Cond,Cond\nR2-9_Cond,Cond\nR3-9_Cond,Cond\n" > phenodata.csv

echo "=== [STAGE 8] Count Distribution Check ==="
python3 -c "import pandas as pd; df = pd.read_csv('gene_count_matrix.csv', index_col=0); print(df.describe())"

echo "=== Upstream Tuxedo Pipeline Completed Successfully ==="