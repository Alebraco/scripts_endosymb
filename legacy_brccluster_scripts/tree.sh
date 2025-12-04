#!/bin/bash 
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH --mem 128G
#SBATCH -e tree.e
#SBATCH -o tree.o
#SBATCH -p standard
#SBATCH -c 32

raxmlHPC-PTHREADS -T 32 -s precise_alignment.fasta -n 16S_precise_tree -m GTRGAMMA -p 12345 -f a -x 12345 -# 100 
