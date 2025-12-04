#!/bin/bash 

groups=(
'endosymb_only'
'endosymb+relatives'
'relatives_only'
)

for group in ${groups[@]}; do
    echo "Processing group: $group"

    mkdir -p "${group}/gbff_files"
    location="$group/proteins/"

    find $location -name "*.faa" | while read -r faa_file; do
        genome_id=$(basename $faa_file .faa)
	spname=$(basename $(dirname $faa_file))
	outdir="${group}/gbff_files/$spname"
	mkdir -p $outdir
	echo "Processing $spname - $genome_id"

        datasets download genome accession $genome_id --include gbff\
        --filename ${outdir}/${genome_id}.zip

	unzip -p "${outdir}/${genome_id}.zip" "ncbi_dataset/data/*/*.gbff" > "${outdir}/${genome_id}.gbff"
    done
done
