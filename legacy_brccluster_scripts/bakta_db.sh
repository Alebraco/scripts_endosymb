#!/bin/bash

#SBATCH -t 2-00:00:00
#SBATCH --mem=50GB
#SBATCH -e bakta.err
#SBATCH -o bakta.out
#SBATCH -p standard

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
bakta_db download --output bakta_db/ --type full
