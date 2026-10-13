#!/usr/bin/env bash
# ==============================================================================
# Pilot upstream processing: runs modules 01-04 in order.
#
# Inputs (not in the repository):
#   data/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa, data/genome/TAIR10_GFF3_genes.gff
#   data/samples/*_R1_001.fastq.gz, *_R2_001.fastq.gz (file names in modules/02_hisat2_alignment.sh)
# Outputs: results/ (SAM, BAM, GTF, Ballgown tables, count matrices)
# Each module skips the steps whose output already exists.
# Threads can be set with THREADS=<n> (default 10; 6 for hisat2-build and samtools sort).
# ==============================================================================

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for MODULE in 01_reference_and_index.sh 02_hisat2_alignment.sh 03_assembly_and_merge.sh 04_quantification.sh; do
    echo "=== ${MODULE} ==="
    bash "${BASE_DIR}/modules/${MODULE}"
done

echo "=== Upstream steps finished; next: Rscript pipeline/run_downstream_pilot.R ==="
