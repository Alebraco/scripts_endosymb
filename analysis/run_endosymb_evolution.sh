#!/bin/bash
#BSUB -J evolution
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -W 24:00
#BSUB -q bobay

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation/
python3 -u scripts_endosymb/analysis/endosymb_evolution.py