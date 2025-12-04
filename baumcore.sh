#!/bin/bash

#BSUB -n 1
#BSUB -W 5:00
#BSUB -J cruncher
#BSUB -R rusage[mem=25GB] span[hosts=1]
#BSUB -e %J.err
#BSUB -o %J.out

core_cruncher="CoreCruncher/corecruncher_master.py"
#group=($(find Hodgkinia_core_input/ -mindepth 1 -type d -name '*sim'))

input_dir="endosymb_proteins/Baumannia"

output_dir="Baumannia_core_output/$(basename $input_dir)"
mkdir -p "$output_dir"
python3 "$core_cruncher" -in "$input_dir" -out "$output_dir" -ext .faa -score 70
