#!/bin/bash
#BSUB -J atb_unpack_submit
#BSUB -o atb_unpack_submit.out
#BSUB -e atb_unpack_submit.err
#BSUB -W 24:00
#BSUB -R "rusage[mem=2GB] span[hosts=1]"
#BSUB -n 1

# Wrapper: builds a species list from atb_bacteria/annotations/, then submits
# LSF array chunks where each task unpacks one species' Bakta JSONs into the
# feature-pipeline layout. Chunks are submitted one at a time; the script
# waits when pending jobs are at or above MAX_PENDING to avoid hitting the
# per-user pending job threshold.
#
# Run after submit_atb_assemblies.sh and submit_atb_bakta.sh have both finished.

OUTDIR="/rsstu/users/l/ljbobay/recombination/asoneto/atb_bacteria"
SCRIPT_DIR="$(dirname "$0")"
ANN_DIR="${OUTDIR}/annotations"
SPECIES_LIST="${OUTDIR}/unpack_species_list.txt"

CHUNK_SIZE=50   # tasks per submitted sub-array
MAX_PENDING=100 # wait if this many PEND jobs are already queued
CONCURRENT=20   # max running tasks within each sub-array (%N)

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

echo "Found $N species: submitting in chunks of ${CHUNK_SIZE} (max pending: ${MAX_PENDING})."

start=1
while [[ $start -le $N ]]; do
    end=$((start + CHUNK_SIZE - 1))
    [[ $end -gt $N ]] && end=$N

    # Wait until the pending job count drops below the threshold
    while true; do
        pending=$(bjobs -noheader 2>/dev/null | awk '$3=="PEND"' | wc -l)
        if [[ $pending -lt $MAX_PENDING ]]; then
            break
        fi
        echo "$(date): ${pending} pending jobs: waiting 60s before submitting chunk ${start}-${end}."
        sleep 60
    done

    echo "Submitting chunk ${start}-${end}..."
    bsub -J "atb_unpack[${start}-${end}]%${CONCURRENT}" \
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

    start=$((end + 1))
done

echo "All chunks submitted."
