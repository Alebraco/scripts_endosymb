#!/bin/bash

#SBATCH -N 1
#SBATCH --job-name=id_processing
#SBATCH --mem=100GB
#SBATCH -p standard
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH -t 2-00:00:00

export PYTHONUNBUFFERED=1
echo "Started running gb2assm.py"
python3 gb2assm.py
echo "----------------------------"
echo "Completed running gb2assm.py"

echo "Started running cand_genomes.sh"
bash cand_genomes.sh
echo "----------------------------"
echo "Completed running cand_genomes.sh"

echo "Started running unzip_candidates.sh"
bash unzip_candidates.sh
echo "----------------------------"
echo "Completed running unzip_candidates.sh"

