#!/bin/bash
#SBATCH --job-name=schic_diagnose
#SBATCH --output=logs/diagnose_%j.out
#SBATCH --error=logs/diagnose_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
##SBATCH --account=def-yourpi
# Step 3: barcode QC on a sample of raw reads (R1 and R2).
set -euo pipefail
REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${CONFIG:-${REPO_DIR}/config.sh}"
setup_env
mkdir -p "${OUT_DIR}"

python3 "${REPO_DIR}/scripts/diagnose_barcodes.py" \
    --r1 "${R1}" --r2 "${R2}" \
    --bcr1 "${BCR1}" --bcr2 "${BCR2}" --bcr3 "${BCR3}" \
    --n-reads "${DIAG_READS}" \
    --json-out "${OUT_DIR}/${SAMPLE}_diagnose.json" \
    > "${OUT_DIR}/${SAMPLE}_diagnose.txt"
