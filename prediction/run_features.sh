#!/bin/bash
#BSUB -n 1
#BSUB -W 100:00
#BSUB -R "rusage[mem=4GB] span[hosts=1]"
#BSUB -J collect_features
#BSUB -o %J.out
#BSUB -e %J.err

python collect_features.py