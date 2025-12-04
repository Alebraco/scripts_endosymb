#!/bin/bash
group="relatives_only"
parent_dir="$group/proteins/"
species=($(find $parent_dir -maxdepth 1 -mindepth 1 -type d))

#BSUB -q bobay
#BSUB -n 1
#BSUB -W 10:00
#BSUB -J cruncher[1-52]
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out


source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

core_cruncher="CoreCruncher/corecruncher_master.py"
input_dir="${species[$(( LSB_JOBINDEX -  1))]}"

if [[ $(ls $input_dir | wc -l) -lt 2 ]]; then
	echo "Only one available genome for $(basename $input_dir)"
else
	output_dir="$group/core/$(basename $input_dir)"
	
	mkdir -p "$output_dir"
	
	if [[ "$group" == "endosymb+relatives" ]]; then
	# IF ENDOSYMBIONTS IN THE DATASET, DEFINE RANDOM ENDOSYMB AS PIVOT
		pivot=$(find $input_dir -maxdepth 1 -type f | grep -v '_genomic' | head -n 1)
		pivot_name=$(basename $pivot)
		python3 "$core_cruncher" -in "$input_dir" -out "$output_dir" -ext .faa -score 70 -ref "$pivot_name"
	else
	# RUN CORECRUNCHER
		python3 "$core_cruncher" -in "$input_dir" -out "$output_dir" -ext .faa -score 70

	fi
fi
