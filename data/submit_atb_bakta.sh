#!/bin/bash
#SBATCH --job-name=atb_bkt_submit
#SBATCH --output=atb_bkt_submit.out
#SBATCH --error=atb_bkt_submit.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# Wrapper script: reads the bakta batch index to determine how many batches
# exist, then launches run_atb_bakta.sh as a SLURM array job.
#
# One array task = one bakta tar.xz batch = N annotation files to extract.
# Tasks are throttled to 10 running at a time (%10) to avoid hammering
# the ATB download servers.
#
# Usage:
#   sbatch data/submit_atb_bakta.sh
#
# Must be run after download_atb_bacteria.py has generated the batch index.
# Can be submitted at the same time as submit_atb_assemblies.sh — the two
# pipelines are independent.
# Edit OUTDIR below to match the --outdir you passed to that script.

# Root output directory — must match what was passed to download_atb_bacteria.py
OUTDIR="atb_bacteria"

# Directory containing this script; used to locate run_atb_bakta.sh
SCRIPT_DIR="$(dirname "$0")"

INDEX="${OUTDIR}/bkt_batch_index.tsv"

if [[ ! -f "$INDEX" ]]; then
    echo "ERROR: Batch index not found: $INDEX"
    echo "Run download_atb_bacteria.py first to generate the manifest and batch lists."
    exit 1
fi

# Count lines = number of bakta batches = upper bound for the array
N=$(wc -l < "$INDEX")

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: $INDEX is empty."
    exit 1
fi

mkdir -p "${OUTDIR}/logs"

echo "Found $N bakta batches — submitting array job 1-${N}."
sbatch --array=1-${N}%10 \
       --output="${OUTDIR}/logs/bkt_%A_%a.out" \
       --error="${OUTDIR}/logs/bkt_%A_%a.err" \
       "${SCRIPT_DIR}/run_atb_bakta.sh"
