#!/bin/bash

for species in bakta_results/*; do
	sp=$(basename $species)
	for dir in $species/*; do
		dirname=$(basename $dir)
		outdir="cruncher_input/$sp/$dirname/"
		mkdir -p $outdir
		for genome in $dir/*; do
			cp $genome/*genomic.faa $outdir
		done
	done
done
