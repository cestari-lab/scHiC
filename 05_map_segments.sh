#!/bin/bash
#SBATCH --job-name=schic_map
#SBATCH --output=logs/map_%j.out
#SBATCH --error=logs/map_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
##SBATCH --account=def-yourpi
# Step 5: map segments single-end with bwa, keep unique alignments,
# and build per-cell contacts.
set -euo pipefail
REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${CONFIG:-${REPO_DIR}/config.sh}"
setup_env
cd "${OUT_DIR}"

SEG="${SAMPLE}_segments.fastq.gz"
FILTERED="${SAMPLE}_segments.filtered.bam"   # bwa order, used for contacts
SORTED="${SAMPLE}_segments.sorted.bam"       # coordinate order, for IGV / QC

[ -f "${GENOME_FA}.bwt" ] || bwa index "${GENOME_FA}"
[ -f "${GENOME_FA}.fai" ] || samtools faidx "${GENOME_FA}"

# -C copies the CB:Z:cell tag from the FASTQ header into the alignment.
# -F 0x904 drops unmapped, secondary and supplementary alignments.
bwa mem -t "${THREADS}" -C "${GENOME_FA}" "${SEG}" \
    | samtools view -h -b -q "${MIN_MAPQ}" -F 0x904 -o "${FILTERED}" -

samtools sort -@ "${THREADS}" -o "${SORTED}" "${FILTERED}"
samtools index "${SORTED}"
samtools flagstat "${SORTED}" > "${SAMPLE}_segments.flagstat.txt"

samtools view -h "${FILTERED}" \
    | python3 "${REPO_DIR}/scripts/segments_to_pairs.py" --sam - \
        --min-mapq "${MIN_MAPQ}" \
        --out-pairs "${SAMPLE}.pairs.gz" \
        --per-cell-dir per_cell_pairs \
        --stats-out "${SAMPLE}_contacts_per_cell.tsv"
