#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p standard
#SBATCH -t 5:00:00
#SBATCH --mem=64G
#SBATCH -e %j.e
#SBATCH -o %j.o

mkdir -p 16S
for dir in */; do
        sp=${dir%/}
	mkdir -p 16S/$sp
        for genome in $dir*.fna; do
		genome_name=$(basename "$genome")
                barrnap $genome -t 16 -k bac -q -o "16S/$sp/${genome_name%.fna}_16S.fna"
        done
done
