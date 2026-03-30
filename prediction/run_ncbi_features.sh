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
BASE_PATH="ncbi_bacteria"
INPUT_DIRS=("$BASE_PATH/proteins")
TRANS_DB="files/IS_db_fixed.dmnd"


# COUNT FILES
NUM_FILES=$(find "${INPUT_DIRS[@]}" -type f -name "*.faa" | wc -l)
NUM_JOBS=$(( (NUM_FILES + 59) / 60 ))
echo "Found $NUM_FILES files, corresponding to $NUM_JOBS array jobs"

if [[ "$NUM_FILES" -eq 0 ]]; then
    echo "No .faa files found in ${INPUT_DIRS[*]}. Exiting."
    exit 1
fi

JOB_ID=$(bsub \
        -J "transposase[1-$NUM_JOBS]%30" \
        -n 8 \
        -R "rusage[mem=16GB] span[hosts=1]" \
        -q bobay \
        -W 24:00 \
        -e prediction_logs/logs_IS/%J_%I.err \
        -o prediction_logs/logs_IS/%J_%I.out \
        -env "INPUT_DIRS='${INPUT_DIRS[*]}',TRANS_DB='$TRANS_DB'" \
        bash scripts_endosymb/analysis/run_transposase.sh | grep -oP '(?<=<)\d+(?=>)')

echo "Submitted transposase array job $JOB_ID ($NUM_JOBS jobs)"

while true; do
    job=$(bjobs "$JOB_ID" 2>&1)
    if echo "$job" | grep -qE 'RUN|PEND'; then
        sleep 60
        continue
    fi
    break
done

if bjobs -d "$JOB_ID" 2>/dev/null | grep -q 'EXIT'; then
    echo "Job $JOB_ID failed. Check logs in prediction_logs/logs_IS/. Stopping."
    exit 1
fi

echo "All jobs done, running Python script."

python -u scripts_endosymb/prediction/collect_features.py --path "$BASE_PATH" --infer
echo "Feature collection complete."
