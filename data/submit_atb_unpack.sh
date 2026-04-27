#!/bin/bash
#BSUB -J atb_unpack_submit
#BSUB -o atb_unpack_submit.out
#BSUB -e atb_unpack_submit.err
#BSUB -W 00:05
#BSUB -R "rusage[mem=2GB] span[hosts=1]"
#BSUB -n 1

# Wrapper: builds a species list from atb_bacteria/annotations/, then submits
# an LSF array where each task unpacks one species' Bakta JSONs into the
# feature-pipeline layout
#
# Run after submit_atb_assemblies.sh and submit_atb_bakta.sh have both finished.

OUTDIR="/rsstu/users/l/ljbobay/recombination/asoneto/atb_bacteria"
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
bsub -J "atb_unpack[1-${N}]%20" \
     -W 02:00 \
     -M 4000 \
     -n 1 \
     -o "${OUTDIR}/logs/unpack_%J_%I.out" \
     -e "${OUTDIR}/logs/unpack_%J_%I.err" \
     -env "SPECIES_LIST=${SPECIES_LIST},OUTDIR=${OUTDIR}" \
     "source ~/.bashrc; \
      conda activate endosymb; \
      SPECIES=\$(sed -n \"\${LSB_JOBINDEX}p\" \"${SPECIES_LIST}\"); \
      python ${SCRIPT_DIR}/unpack_bakta_json.py \
          --annotations-dir ${OUTDIR}/annotations \
          --genomes-dir ${OUTDIR}/genomes \
          --proteins-dir ${OUTDIR}/proteins \
          --species \"\${SPECIES}\""
