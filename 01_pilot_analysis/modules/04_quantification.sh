#!/usr/bin/env bash
# ==============================================================================
# Module 04: StringTie quantification and count matrices (log 2026-03-23)
#
# Input : results/bam/<sample>.bam, results/gtf/at_merged.gtf (module 03)
# Output: results/ballgown/<sample>/ (Ballgown tables, including t_data.ctab)
#         results/matrices/gene_count_matrix.csv, transcript_count_matrix.csv, phenodata.csv
# ==============================================================================

set -euo pipefail

THREADS="${THREADS:-10}"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MERGED_GTF="${BASE_DIR}/results/gtf/at_merged.gtf"
BAM_DIR="${BASE_DIR}/results/bam"
BALLGOWN_DIR="${BASE_DIR}/results/ballgown"
MATRICES_DIR="${BASE_DIR}/results/matrices"
mkdir -p "${BALLGOWN_DIR}" "${MATRICES_DIR}"

SAMPLES="R1-1_WT R2-1_WT R3-1_WT R1-9_Cond R2-9_Cond R3-9_Cond"

# 1. Quantify every library against the merged transcript set
for SAMPLE in ${SAMPLES}; do
    OUT_DIR="${BALLGOWN_DIR}/${SAMPLE}"
    if [[ ! -f "${OUT_DIR}/t_data.ctab" ]]; then
        echo "[04] Quantifying ${SAMPLE}"
        mkdir -p "${OUT_DIR}"
        stringtie -e -B -p "${THREADS}" -G "${MERGED_GTF}" -o "${OUT_DIR}/${SAMPLE}.gtf" "${BAM_DIR}/${SAMPLE}.bam"
    fi
done

# 2. Gene and transcript count matrices (prepDE.py3 with its default settings, as in the original run)
cd "${MATRICES_DIR}"
if [[ ! -f prepDE.py3 ]]; then
    curl -sS -O http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
fi
python3 prepDE.py3 -i "${BALLGOWN_DIR}"

# 3. Sample table for edgeR
printf "ids,condition\nR1-1_WT,WT\nR2-1_WT,WT\nR3-1_WT,WT\nR1-9_Cond,Cond\nR2-9_Cond,Cond\nR3-9_Cond,Cond\n" > phenodata.csv

# 4. Count distribution per library
python3 -c "import pandas as pd; df = pd.read_csv('gene_count_matrix.csv', index_col=0); print(df.describe())"
