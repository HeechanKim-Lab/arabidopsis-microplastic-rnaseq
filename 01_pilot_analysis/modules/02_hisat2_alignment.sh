#!/usr/bin/env bash
# ==============================================================================
# Module 02: paired-end alignment with HISAT2 (log 2026-03-13)
#
# Input : data/indexes/arabidopsis_index.*.ht2 (module 01)
#         data/samples/<FASTQ prefix>_R1_001.fastq.gz and _R2_001.fastq.gz
# Output: results/sam/<sample>.sam
# ==============================================================================

set -euo pipefail

THREADS="${THREADS:-10}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDEX_PREFIX="${BASE_DIR}/data/indexes/arabidopsis_index"
SAMPLES_DIR="${BASE_DIR}/data/samples"
SAM_DIR="${BASE_DIR}/results/sam"
mkdir -p "${SAM_DIR}"

# Sample name : FASTQ file prefix
SAMPLES="R1-1_WT:R1-1_S0_L009
R2-1_WT:R2-1_S252_L004
R3-1_WT:R3-1_S261_L004
R1-9_Cond:R1-9_S251_L004
R2-9_Cond:R2-9_S0_L009
R3-9_Cond:R3-9_S0_L009"

for ENTRY in ${SAMPLES}; do
    SAMPLE="${ENTRY%%:*}"
    PREFIX="${ENTRY#*:}"
    R1="${SAMPLES_DIR}/${PREFIX}_R1_001.fastq.gz"
    R2="${SAMPLES_DIR}/${PREFIX}_R2_001.fastq.gz"
    SAM_OUT="${SAM_DIR}/${SAMPLE}.sam"

    if [[ -f "${SAM_OUT}" ]]; then
        echo "[02] ${SAMPLE}: ${SAM_OUT} exists, skipped"
        continue
    fi
    [[ -f "${R1}" && -f "${R2}" ]] || { echo "Missing FASTQ files for ${SAMPLE}: ${R1}, ${R2}" >&2; exit 1; }

    echo "[02] Aligning ${SAMPLE}"
    hisat2 -p "${THREADS}" --dta -x "${INDEX_PREFIX}" -1 "${R1}" -2 "${R2}" -S "${SAM_OUT}"
done
