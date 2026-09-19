#!/usr/bin/env bash
# ==============================================================================
# Module: 04_quantification.sh
# Description: StringTie Ballgown estimation and raw count extraction (prepDE)
# ==============================================================================

set -euo pipefail

THREADS=10
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MERGED_GTF="${BASE_DIR}/results/gtf/at_merged.gtf"
BAM_DIR="${BASE_DIR}/results/bam"
BALLGOWN_DIR="${BASE_DIR}/results/ballgown"
MATRICES_DIR="${BASE_DIR}/results/matrices"

mkdir -p "${BALLGOWN_DIR}" "${MATRICES_DIR}"

SAMPLES=("R1-1_WT" "R2-1_WT" "R3-1_WT" "R1-9_Cond" "R2-9_Cond" "R3-9_Cond")

# 1. StringTie Quantification
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_BG="${BALLGOWN_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_BG}"
    if [[ ! -f "${SAMPLE_BG}/t_data.ctab" ]]; then
        echo ">>> [04] Quantifying abundances for ${SAMPLE}..."
        stringtie -e -B -p "${THREADS}" \
            -G "${MERGED_GTF}" \
            -o "${SAMPLE_BG}/${SAMPLE}.gtf" \
            "${BAM_DIR}/${SAMPLE}.bam"
    fi
done

# 2. Count Matrix Extraction
cd "${MATRICES_DIR}"
if [[ ! -f "prepDE.py3" ]]; then
    curl -s -O http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
fi

python3 prepDE.py3 -i "${BALLGOWN_DIR}"

# 3. Phenotype Metadata Generation
printf "ids,condition\nR1-1_WT,WT\nR2-1_WT,WT\nR3-1_WT,WT\nR1-9_Cond,Cond\nR2-9_Cond,Cond\nR3-9_Cond,Cond\n" > phenodata.csv

# 4. Quartile Distribution Audit
python3 -c "import pandas as pd; df = pd.read_csv('gene_count_matrix.csv', index_col=0); print(df.describe())"