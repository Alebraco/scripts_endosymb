#!/bin/bash
# Step 2 wrapper: build files/training_set_B.csv.
# Submit from login node: bsub < scripts_endosymb/prediction/run_step2_build_training_set_B.sh

#BSUB -J build_train_B
#BSUB -W 0:30
#BSUB -q bobay
#BSUB -e build_train_B.err
#BSUB -o build_train_B.out
#BSUB -R "rusage[mem=4GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python -u scripts_endosymb/prediction/build_training_set_B.py
