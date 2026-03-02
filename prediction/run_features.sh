#!/bin/bash
#BSUB -n 1
#BSUB -q bobay
#BSUB -W 100:00
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -J collect_features
#BSUB -o %J.out
#BSUB -e %J.err

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation
python -u scripts_endosymb/prediction/collect_features.py