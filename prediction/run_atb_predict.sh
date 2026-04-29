#!/bin/bash
#BSUB -J atb_predict
#BSUB -W 48:00
#BSUB -q bobay
#BSUB -e atb_predict.err
#BSUB -o atb_predict.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 1

# End-to-end: run the feature pipeline on atb_bacteria proteins, then apply
# the pre-trained RF model and emit high-confidence endosymbiont calls.
#
# Assumes:
#   - atb_bacteria/proteins/{Species}/*.faa already populated by
#     data/submit_atb_unpack.sh + unpack_bakta_json.py
#   - files/rf_model.joblib already produced by run_models.sh

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

BASE_PATH="atb_bacteria"

# Stage 1: transposase BLAST + feature collection (chunked submission).
bash scripts_endosymb/prediction/run_ncbi_features.sh "$BASE_PATH"

# Stage 2: apply trained RF model to the resulting features.
python -u scripts_endosymb/prediction/predict_atb.py \
    --features-csv "${BASE_PATH}/feature_files/combined_features.csv" \
    --model-path files/rf_model.joblib

echo "ATB prediction pipeline complete."
