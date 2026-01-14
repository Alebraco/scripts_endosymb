#!/bin/bash
#BSUB -J igs
#BSUB -n 1
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -W 10:00
source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation/
python3 scripts_endosymb/analysis/igs_lengths.py
python3 scripts_endosymb/analysis/igs_plots.py