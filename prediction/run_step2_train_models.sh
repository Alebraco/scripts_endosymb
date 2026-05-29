#!/bin/bash
# Step 2 wrapper: build training_set_B.csv, train Model B, re-train Model A.
# Runs all three in one job

#BSUB -J train_models
#BSUB -W 2:00
#BSUB -q bobay
#BSUB -e train_models.err
#BSUB -o train_models.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python -u scripts_endosymb/prediction/build_training_set_B.py
python -u scripts_endosymb/prediction/train_model_B.py
python -u scripts_endosymb/prediction/train_model_A.py

echo "Step 2 complete: files/model_{A,B}_holdout_metrics.json"
