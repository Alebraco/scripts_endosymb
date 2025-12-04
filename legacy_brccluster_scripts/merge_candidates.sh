#!/bin/bash

genomes="candidate_genomes"
out_dir="merged_candidates"
mkdir -p $out_dir

for sp in $genomes/*; do
	echo "Processing $sp"
	sp_name=$(basename $sp)
	sp_outdir="$out_dir/$sp_name/"
	mkdir -p $sp_outdir
	find $sp -name "*fna" -type f -exec cp {} $sp_outdir \;
done
