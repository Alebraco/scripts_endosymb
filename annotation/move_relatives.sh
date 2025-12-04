#!/bin/bash
#BSUB -n 1
#BSUB -W 2:00

annotations="bakta_results"
out_dir="relatives_only"
mkdir -p $out_dir


for sp in $annotations/*; do
	echo "Processing $sp"
	sp_name=$(basename $sp)
	sp_outdir="$out_dir/$sp_name/"
	mkdir -p "$sp_outdir"/{genomes,proteins,gbff_files,cds}
	find $sp -name "*.gbff" -type f -exec cp {} "$sp_outdir"/gbff_files/ \;
	find $sp -name "*.fna" -type f -exec cp {} "$sp_outdir"/genomes/ \;
	find $sp -name "*.ffn" -type f -exec sh -c 'cp "$0" "${1}/cds/$(basename "$0" .ffn).fna"' {} "$sp_outdir" \;
	find $sp -name "*.faa" -type f -exec cp {} "$sp_outdir"/proteins/ \;
done