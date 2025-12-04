#!/bin/bash
#SBATCH --job-name=blast_ids
#SBATCH --output=blast_ids.out
#SBATCH --error=blast_ids.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

module load blast+
export PYTHONUNBUFFERED=1
python3 blast_ids.py
