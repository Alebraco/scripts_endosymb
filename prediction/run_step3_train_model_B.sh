#!/bin/bash
# Step 3 + 4 wrapper: train Model B and evaluate on holdout_locked.csv.
# Writes models/model_B.joblib and results/model_B_holdout.json.
# Submit from login node: bsub < scripts_endosymb/prediction/run_step3_train_model_B.sh

#BSUB -J train_model_B
#BSUB -W 2:00
#BSUB -q bobay
#BSUB -e train_model_B.err
#BSUB -o train_model_B.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python -u scripts_endosymb/prediction/train_model_B.py
