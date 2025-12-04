#!/bin/bash
#BSUB -n 1
#BSUB -e %J.err
#BSUB -o %J.out 
#BSUB -W 96:00
#BSUB -R "rusage[mem=50GB]"

export PYTHONUNBUFFERED=1
python3 usearch.py
