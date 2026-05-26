#!/bin/bash
# Step 3b wrapper: re-train Model A so its metrics file includes OOB error.
# Forest is identical to the original (same random_state=28); only OOB is added.
# Submit from login node: bsub < scripts_endosymb/prediction/run_step3b_train_model_A.sh

#BSUB -J train_model_A
#BSUB -W 2:00
#BSUB -q bobay
#BSUB -e train_model_A.err
#BSUB -o train_model_A.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python -u scripts_endosymb/prediction/train_model_A.py
