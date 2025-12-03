#!/bin/bash 
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH --mem 128G
#SBATCH -e %j.e
#SBATCH -o %j.o
#SBATCH -p standard
#SBATCH -c 16

species_dir="species_16S/*"
out_dir="species_alignments/"

export out_dir

find $species_dir -name '*_16S.fna' | parallel --jobs 16 '
	sp_file={}
	sp=$(basename ${sp_file})
	outfile="${out_dir}algn_${sp%_16S.fna}.fasta"
	muscle -align $sp_file -output $outfile
'	
