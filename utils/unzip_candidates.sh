#!/bin/bash

for file in candidate_genomes/*/*.zip; do
	outdir=${file%.zip}
	mkdir -p $outdir
	unzip -j $file "ncbi_dataset/data/*/*.fna" -d $outdir  
done


