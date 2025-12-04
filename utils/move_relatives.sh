#!/bin/bash
#BSUB -n 1
#BSUB -W 2:00

annotations="bakta_results"
out_dir="relatives"
mkdir -p $out_dir

for sp in $annotations/*; do
	echo "Processing $sp"
	sp_name=$(basename $sp)
	sp_outdir="$out_dir/$sp_name/"
	mkdir -p $sp_outdir
	find $sp -name "*genomic.gbff" -type f -exec cp {} $sp_outdir \;
done
