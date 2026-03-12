#!/bin/bash
#BSUB -J models
#BSUB -W 24:00
#BSUB -q bobay
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "rusage[mem=16GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python3 -u scripts_endosymb/prediction/plot_correlation.py