#!/bin/bash

#BSUB -n 1
#BSUB -W 96:00
#BSUB -J cruncher[1-48]
#BSUB -R "rusage[mem=25GB] span[hosts=1]"
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out


source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/bakta

core_cruncher="CoreCruncher/corecruncher_master.py"
species=($(find endosymb_proteins/* -maxdepth 1 -type d))
input_dir="${species[$(( LSB_JOBINDEX -  1))]}"

output_dir="endosymb_core/$(basename $input_dir)"

if [ $(ls $input_dir | wc -l) -gt 1 ]; then
	mkdir -p $output_dir
	python3 "$core_cruncher" -in "$input_dir" -out "$output_dir" -ext .faa -score 70
fi
