#!/bin/bash
#BSUB -J pca_atb
#BSUB -W 1:00
#BSUB -q bobay
#BSUB -e pca_atb.err
#BSUB -o pca_atb.out
#BSUB -R "rusage[mem=2GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python -u scripts_endosymb/analysis/visualize_atb_pca.py