#!/usr/bin/env bash
# ==============================================================================
# Module: 02_hisat2_alignment.sh
# Description: Spliced paired-end read mapping across the 6 pilot libraries
# ==============================================================================

set -euo pipefail

THREADS=10
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDEX_PREFIX="${BASE_DIR}/data/indexes/arabidopsis_index"
SAMPLES_DIR="${BASE_DIR}/data/samples"
SAM_DIR="${BASE_DIR}/results/sam"

mkdir -p "${SAM_DIR}"

declare -A FASTQ_R1 FASTQ_R2
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

for SAMPLE in "R1-1_WT" "R2-1_WT" "R3-1_WT" "R1-9_Cond" "R2-9_Cond" "R3-9_Cond"; do
    SAM_OUT="${SAM_DIR}/${SAMPLE}.sam"
    if [[ ! -f "${SAM_OUT}" ]]; then
        echo ">>> [02] Aligning ${SAMPLE}..."
        hisat2 -p "${THREADS}" --dta \
            -x "${INDEX_PREFIX}" \
            -1 "${FASTQ_R1[$SAMPLE]}" \
            -2 "${FASTQ_R2[$SAMPLE]}" \
            -S "${SAM_OUT}"
    else
        echo ">>> [02] Alignment for ${SAMPLE} exists. Skipping."
    fi
done