#!/bin/bash
#SBATCH --job-name=atb_unpack_submit
#SBATCH --output=atb_unpack_submit.out
#SBATCH --error=atb_unpack_submit.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# Wrapper: builds a species list from atb_bacteria/annotations/, then submits
# a SLURM array where each task unpacks one species' Bakta JSONs into the
# feature-pipeline layout
#
# Run after submit_atb_assemblies.sh and submit_atb_bakta.sh have both finished.

OUTDIR="atb_bacteria"
SCRIPT_DIR="$(dirname "$0")"
ANN_DIR="${OUTDIR}/annotations"
SPECIES_LIST="${OUTDIR}/unpack_species_list.txt"

if [[ ! -d "$ANN_DIR" ]]; then
    echo "ERROR: annotations directory not found: $ANN_DIR"
    exit 1
fi

mkdir -p "${OUTDIR}/logs"

find "$ANN_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$SPECIES_LIST"
N=$(wc -l < "$SPECIES_LIST")

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: no species subdirectories under $ANN_DIR"
    exit 1
fi

echo "Found $N species — submitting array job 1-${N}."
sbatch --array=1-${N}%20 \
       --job-name=atb_unpack \
       --time=02:00:00 \
       --mem=4G \
       --cpus-per-task=1 \
       --output="${OUTDIR}/logs/unpack_%A_%a.out" \
       --error="${OUTDIR}/logs/unpack_%A_%a.err" \
       --export=ALL,SPECIES_LIST="${SPECIES_LIST}",OUTDIR="${OUTDIR}" \
       --wrap "source ~/.bashrc; \
               conda activate endosymb; \
               SPECIES=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"${SPECIES_LIST}\"); \
               python ${SCRIPT_DIR}/unpack_bakta_json.py \
                   --annotations-dir ${OUTDIR}/annotations \
                   --genomes-dir ${OUTDIR}/genomes \
                   --proteins-dir ${OUTDIR}/proteins \
                   --species \"\${SPECIES}\""
