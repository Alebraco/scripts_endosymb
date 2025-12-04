#!/bin/bash
rnadir='16S_coords'
complete_path='16S_seqs'
mkdir -p $rnadir $complete_path
for species in endosymb_genomes/Hodgkinia_*; do 
	sp_name=$(basename $species)
	sp_dir="$rnadir/$sp_name"
	mkdir -p $sp_dir
	for gff in "$species"/*.gff; do
		[ -e "$gff" ] || continue
		gff="$species/$gff"
		
		filename=$(basename "$gff" .gff)
		outfile="$sp_dir"/${filename}.bed

		grep '16S ribosomal RNA' $gff | awk '$3 == "rRNA" {print $1 "\t" $4-1 "\t" $5}' > $outfile
		
		[ -s $outfile ] || continue
		genome_file=endosymb_genomes/$sp_name/${filename}.fna
		seqkit subseq --bed $outfile $genome_file | seqkit seq -m 1200 | seqkit replace -p ':\.\s$' -r '' >> $complete_path/$sp_name.fna
	done
done