#!/bin/bash
#BSUB -J theta_corr
#BSUB -W 1:00
#BSUB -q bobay
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python3 -u scripts_endosymb/analysis/theta_correlations.py --force