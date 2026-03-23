#!/bin/bash
#BSUB -J polymorphisms_distribution
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4GB]"
#BSUB -W 24:00
#BSUB -q bobay

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

python3 -u scripts_endosymb/analysis/codon_polymorphisms.py

