#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20GB
#SBATCH -e %j.err
#SBATCH -o %j.log


export PYTHONUNBUFFERED=1

python3 gc_correlation.py
