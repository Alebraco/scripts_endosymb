#!/bin/bash 
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH --mem 128G
#SBATCH -e %j.e
#SBATCH -o %j.o
#SBATCH -p standard
#SBATCH -n 20

muscle -align 16S_sequences.fna -output precise_alignment.fasta
