#!/bin/bash

for file in $(find related_species -type f -name "*accns.txt"); do
	sp=${file#related_species/}
	cleanacc=${sp%.txt}

	outdir="candidate_genomes/$(dirname "$cleanacc")"
	outfile="$outdir/$(basename "$cleanacc")_summary.tsv"
	mkdir -p "$outdir"

	# Extract accessions (first column only)
	awk '{print $1}' $file > temp_accns.txt

	# Obtain summaries for the accessions
	datasets summary genome accession $(cat temp_accns.txt) --as-json-lines --assembly-source 'RefSeq' \
	 | dataformat tsv genome --fields accession,assmstats-gc-percent,assmstats-total-sequence-len \
	 > temp_summary.tsv

	# Join the summaries with the original file to preserve ID%
	join -t $'\t' -1 1 -2 1 <(sort $file) <(sort temp_summary.tsv) > "$outfile"
	
	# Add header to the summary file
	echo -e "Accession\tID%\tGC%\tGenome Size" > $outfile.tmp
	sort -k2,2nr "$outfile" >> $outfile.tmp
	mv "$outfile.tmp" "$outfile"

	# Download the genomes
	datasets download genome accession --inputfile temp_accns.txt --include genome \
	--assembly-source 'RefSeq' --filename "$outdir/$(basename "$cleanacc").zip" 

	# Clean up the temporary files
	rm temp_accns.txt temp_summary.tsv
done
