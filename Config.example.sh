#!/bin/bash
# Copy to config.sh and edit:  cp config.example.sh config.sh
# config.sh is git-ignored, so your cluster paths stay off GitHub.

SAMPLE="ScHiC1"
DATA_DIR="/path/to/scHiC/data"
OUT_DIR="/path/to/scHiC/analysis"

R1="${DATA_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${DATA_DIR}/${SAMPLE}_R2.fastq.gz"

BCR1="/path/to/scHiC/barcode/bcr1.txt"
BCR2="/path/to/scHiC/barcode/bcr2.txt"
BCR3="/path/to/scHiC/barcode/bcr3.txt"

GENOME_FA="/path/to/genome.fasta"   # bwa index is built next to it if missing

THREADS=16
MIN_ARM_LEN=20      # shortest genomic segment kept after removing the linker
MIN_MAPQ=10         # mapping-quality filter
DIAG_READS=2000000  # reads sampled per mate by diagnose_barcodes.py

# Software for the SLURM jobs (edit for your cluster)
setup_env() {
    module load StdEnv/2023 python/3.11 bwa/0.7.17 samtools/1.18
}
