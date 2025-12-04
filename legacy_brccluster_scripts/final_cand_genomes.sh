#!/bin/bash
#SBATCH -e %j.e
#SBATCH -o %j.o
#SBATCH -t 2-00:00:00
#SBATCH -p standard
#SBATCH --mem=10GB

#MIN_GC=45
MIN_SIZE=2500000

for file in $(find related_species/ -type f -name "*accns.txt"); do
	sp=${file#related_species/}
	cleanacc=${sp%.txt}

	outdir="candidate_genomes/$(dirname "$cleanacc")"
	outfile="$outdir/$(basename "$cleanacc")_summary.tsv"
	mkdir -p "$outdir"

    echo "Processing $file"
	
    while read -r acc _; do
        echo -e "\t>Processing $acc"
	    datasets summary genome accession $acc --as-json-lines --assembly-source 'RefSeq' \
	        	| dataformat tsv genome --fields accession,assmstats-gc-percent,assmstats-total-sequence-len \
                | tail -n +2 >> temp_summary.tsv
        sleep 1
    done < $file
	  
	
	awk -F '\t' -v g=$MIN_GC -v s=$MIN_SIZE \
	'NR>1 && ($3+0)>=s {print $1}' temp_summary.tsv > filtered_accns.txt 
	
	header="Accession\tGC%\tGenome Size\tFree-living"
    awk -F '\t' -v g=$MIN_GC -v s=$MIN_SIZE 'NR>1 {flag = ($3+0>=s ? "yes" : "no"); print $1 "\t" $2 "\t" $3 "\t" flag}' temp_summary.tsv | sort -k2,2nr > "$outfile.tmp"
    echo -e "$header" | cat - "$outfile.tmp" > "$outfile"
    rm "$outfile.tmp"

	if [[ -s filtered_accns.txt ]]; then
		# Download the genomes
		datasets download genome accession --inputfile filtered_accns.txt --include genome --filename "$outdir/$(basename "$cleanacc").zip" --assembly-source 'RefSeq' 
        echo "Genomes downloaded successfully for $file"
	else
		echo "No genomes passed filters for $file"
	fi
	# Clean up the temporary files
	rm temp_summary.tsv filtered_accns.txt
done
