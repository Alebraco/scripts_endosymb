#!/bin/bash
#BSUB -J models
#BSUB -W 1:00
#BSUB -q bobay
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

BASE_PATH="endosymb+relatives"
python3 -u scripts_endosymb/prediction/bayesian_data_prep.py --base_path $BASE_PATH