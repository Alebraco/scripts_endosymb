#!/bin/bash
#SBATCH --job-name=atb_asm_submit
#SBATCH --output=atb_asm_submit.out
#SBATCH --error=atb_asm_submit.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# Wrapper script: reads the assembly batch index to determine how many batches
# exist, then launches run_atb_worker.sh as a SLURM array job.
#
# Usage:
#   sbatch data/submit_atb_assemblies.sh
#
# Must be run after download_atb_bacteria.py has generated the batch index.
# Can be submitted at the same time as submit_atb_bakta.sh

OUTDIR="atb_bacteria"
SCRIPT_DIR="$(dirname "$0")"
INDEX="${OUTDIR}/asm_batch_index.tsv"

if [[ ! -f "$INDEX" ]]; then
    echo "ERROR: Batch index not found: $INDEX"
    echo "Run download_atb_bacteria.py first to generate the manifest and batch lists."
    exit 1
fi

N=$(wc -l < "$INDEX")

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: $INDEX is empty."
    exit 1
fi

mkdir -p "${OUTDIR}/logs"

echo "Found $N assembly batches — submitting array job 1-${N}."
sbatch --array=1-${N}%10 \
       --export=ALL,MODE=assemblies \
       --output="${OUTDIR}/logs/asm_%A_%a.out" \
       --error="${OUTDIR}/logs/asm_%A_%a.err" \
       "${SCRIPT_DIR}/run_atb_worker.sh"
