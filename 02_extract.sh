#!/bin/bash
#SBATCH --job-name=schic_extract
#SBATCH --output=logs/extract_%j.out
#SBATCH --error=logs/extract_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
##SBATCH --account=def-yourpi
# Step 2: find the barcode junction and write genomic segments.
# Submit from the repository root:  mkdir -p logs && sbatch slurm/02_extract.sh
set -euo pipefail
REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${CONFIG:-${REPO_DIR}/config.sh}"
setup_env
mkdir -p "${OUT_DIR}"

python3 "${REPO_DIR}/scripts/extract_schic_reads.py" \
    --r1 "${R1}" --r2 "${R2}" \
    --bcr1 "${BCR1}" --bcr2 "${BCR2}" --bcr3 "${BCR3}" \
    --min-arm-len "${MIN_ARM_LEN}" \
    --out-fastq "${OUT_DIR}/${SAMPLE}_segments.fastq.gz" \
    --stats-out "${OUT_DIR}/${SAMPLE}_extract_stats.tsv"
