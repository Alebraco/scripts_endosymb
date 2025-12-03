#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH --mem 128G
#SBATCH -e %j.e
#SBATCH -o %j.o
#SBATCH -p standard
#SBATCH -c 32

parent="$HOME/endosymb"
alignment_dir="$parent/mafft_alignments/"

for algn in $alignment_dir/*.fasta; do   
        if [[ $(grep -c '>' $algn) -gt 2 ]]; then
                filename=$(basename $algn .fasta)
                species=${filename#algn_}
                outdir="$parent/species_trees/${species}"
                mkdir -p $outdir
                raxmlHPC-PTHREADS -T 32 -s "$algn" -n "$species" -w $outdir -m GTRGAMMA -p 12345 -f a -x 12345 -# 100
	fi
done
