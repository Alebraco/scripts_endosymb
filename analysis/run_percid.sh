#!/bin/bash
#BSUB -J all_dist
#BSUB -o all_dist.out
#BSUB -e all_dist.err
#BSUB -n 1
#BSUB -R "rusage[mem=4GB]"
#BSUB -W 01:00

python3 percid_distribution.py all all

python3 percid_distribution.py all mean
