#!/bin/bash

#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH --job-name hostft
#SBATCH --mem 50G
#SBATCH -e %j.err
#SBATCH -o %j.out

export PYTHONUNBUFFERED=1
python3 host_feature.py
