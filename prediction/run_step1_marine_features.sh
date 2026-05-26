#!/bin/bash

#BSUB -J marine_step1
#BSUB -W 36:00
#BSUB -q bobay
#BSUB -e marine_step1.err
#BSUB -o marine_step1.out
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

BASE="marine_free_livers"
OUTDIR="$BASE/bakta_results"

# Bakta array
SUBMIT_OUT=$(bash scripts_endosymb/annotation/submit_bakta.sh "$BASE" "$OUTDIR")
echo "$SUBMIT_OUT"
BAKTA_JOB=$(echo "$SUBMIT_OUT" | grep -oP '(?<=<)\d+(?=>)' | head -n1)
if [ -z "$BAKTA_JOB" ]; then
    echo "Could not parse bakta array job ID from submit_bakta.sh output."
    exit 1
fi
echo "Bakta array submitted as job $BAKTA_JOB; waiting..."

# Wait for bakta array to finish
while true; do
    job=$(bjobs "$BAKTA_JOB" 2>&1)
    if echo "$job" | grep -qE 'RUN|PEND'; then
        sleep 60
        continue
    fi
    if echo "$job" | grep -q 'is not found' || [ -z "$job" ]; then
        echo "Bakta job $BAKTA_JOB not found; assuming complete."
        break
    fi
    break
done
if bjobs -d "$BAKTA_JOB" 2>/dev/null | grep -q 'EXIT'; then
    echo "Bakta array job $BAKTA_JOB exited with failures. Check bakta_logs/."
    exit 1
fi
echo "Bakta done."

# Organize bakta outputs
bash scripts_endosymb/annotation/organize_bakta_outputs.sh "$BASE"

# Collect features
bash scripts_endosymb/prediction/run_ncbi_features.sh "$BASE"

# Make CSV
python -u scripts_endosymb/prediction/make_marine_csv.py --path "$BASE"

echo "Step 1 complete. Review files/marine_free_livers.csv before submitting Step 2."
