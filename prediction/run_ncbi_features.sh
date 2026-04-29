#!/bin/bash

#BSUB -J ncbi_feature_wrapper
#BSUB -W 24:00
#BSUB -q bobay
#BSUB -e ncbi_feature_wrapper.err
#BSUB -o ncbi_feature_wrapper.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

mkdir -p prediction_logs/logs_IS

# CONFIGURE
BASE_PATH="${1}"
INPUT_DIRS=("$BASE_PATH/proteins")
echo "Using BASE_PATH=$BASE_PATH"
TRANS_DB="files/IS_db_fixed.dmnd"


# COUNT FILES
NUM_FILES=$(find "${INPUT_DIRS[@]}" -type f -name "*.faa" | wc -l)
NUM_JOBS=$(( (NUM_FILES + 59) / 60 ))
echo "Found $NUM_FILES files, corresponding to $NUM_JOBS array jobs"

if [[ "$NUM_FILES" -eq 0 ]]; then
    echo "No .faa files found in ${INPUT_DIRS[*]}. Exiting."
    exit 1
fi

# Chunked submission to keep pending jobs under MAX_PENDING.
# Cluster has a 1000 pending-job hard cap; we stay well below at 500.
CHUNK_SIZE=500   # array tasks per submitted sub-array
MAX_PENDING=500  # wait if this many PEND jobs are already queued
CONCURRENT=30    # max running tasks within each sub-array (%N)

JOB_IDS=()
start=1
while [[ $start -le $NUM_JOBS ]]; do
    end=$((start + CHUNK_SIZE - 1))
    [[ $end -gt $NUM_JOBS ]] && end=$NUM_JOBS

    # Wait until pending job count drops below the threshold.
    while true; do
        pending=$(bjobs -noheader 2>/dev/null | awk '$3=="PEND"' | wc -l)
        if [[ $pending -lt $MAX_PENDING ]]; then
            break
        fi
        echo "$(date): ${pending} pending jobs: waiting 60s before submitting chunk ${start}-${end}."
        sleep 60
    done

    echo "Submitting transposase chunk ${start}-${end}..."
    JOB_ID=$(bsub \
            -J "transposase[${start}-${end}]%${CONCURRENT}" \
            -n 8 \
            -R "rusage[mem=16GB] span[hosts=1]" \
            -q bobay \
            -W 24:00 \
            -e prediction_logs/logs_IS/%J_%I.err \
            -o prediction_logs/logs_IS/%J_%I.out \
            -env "INPUT_DIRS='${INPUT_DIRS[*]}',TRANS_DB='$TRANS_DB'" \
            bash scripts_endosymb/analysis/run_transposase.sh | grep -oP '(?<=<)\d+(?=>)')
    echo "  -> submitted job $JOB_ID"
    JOB_IDS+=("$JOB_ID")
    start=$((end + 1))
done

echo "All chunks submitted (${#JOB_IDS[@]} chunks): ${JOB_IDS[*]}"

# Wait for all chunks to finish (no RUN/PEND in any of them).
while true; do
    any_active=0
    for jid in "${JOB_IDS[@]}"; do
        if bjobs "$jid" 2>&1 | grep -qE 'RUN|PEND'; then
            any_active=1
            break
        fi
    done
    [[ $any_active -eq 0 ]] && break
    sleep 60
done

# Check whether any chunk exited with failure.
for jid in "${JOB_IDS[@]}"; do
    if bjobs -d "$jid" 2>/dev/null | grep -q 'EXIT'; then
        echo "Job $jid failed. Check logs in prediction_logs/logs_IS/. Stopping."
        exit 1
    fi
done

echo "All jobs done, running Python script."

python -u scripts_endosymb/prediction/collect_features.py --path "$BASE_PATH" --infer
echo "Feature collection complete."
