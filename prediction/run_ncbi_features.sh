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

MISSING_FILES=$(find "${INPUT_DIRS[@]}" -type f -name "*.faa" | while read -r FAA; do
    OUTFILE=$(expected_outfile "$FAA")
    [[ -f "$OUTFILE" ]] || echo "$FAA"
done | wc -l)

if [[ "$MISSING_FILES" -eq 0 ]]; then
    echo "All transposase TSV outputs already exist. Skipping transposase array submission."
else
    echo "$MISSING_FILES files are missing transposase outputs; submitting array jobs."

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

    echo "Submitted job array with ID $JOB_ID"

    while true; do
        job=$(bjobs "$JOB_ID" 2>&1)
        if echo "$job" | grep -qE 'RUN|PEND'; then
            sleep 60
            continue
        fi

        if echo "$job" | grep -q 'is not found' || [ -z "$job" ]; then
            echo "Job $JOB_ID not found. Assuming all jobs completed."
            break
        fi
        break
    done

    if bjobs -d "$JOB_ID" 2>/dev/null | grep -q 'EXIT'; then
        echo "Job $JOB_ID failed. Check logs in prediction_logs/logs_IS/. Stopping."
        exit 1
    fi
fi

echo "All jobs done, running Python script."

python -u scripts_endosymb/prediction/collect_features.py --path "$BASE_PATH" --infer
echo "Feature collection complete."
