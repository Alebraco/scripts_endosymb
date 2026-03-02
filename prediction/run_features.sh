#!/bin/bash
#BSUB -n 1
#BSUB -W 100:00
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -J collect_features
#BSUB -o %J.out
#BSUB -e %J.err

python scripts_endosymb/prediction/collect_features.py