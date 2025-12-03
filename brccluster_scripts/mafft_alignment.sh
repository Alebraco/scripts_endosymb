#!/bin/bash 
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH --mem 128G
#SBATCH -e mafft_%j.e
#SBATCH -o mafft_%j.o
#SBATCH -p standard
#SBATCH -c 16

species_dir="species_16S/*"
out_dir="mafft_alignments/"
mkdir -p $out_dir
export out_dir

find $species_dir -name '*.fna' | parallel --jobs 16 '
	sp_file={}
	if [ $(grep -c ">" $sp_file) -gt 1 ]; then
		sp=$(basename ${sp_file})
		outfile="${out_dir}algn_${sp%.fna}.fasta"
        mafft --auto --adjustdirection "$sp_file" > "$outfile"
	fi
'	
