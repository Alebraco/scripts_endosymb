#!/bin/bash

#BSUB -n 1
#BSUB -W 5:00
#BSUB -J cruncher[1-3]
#BSUB -R rusage[mem=25GB] span[hosts=1]
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out


#source ~/.bashrc
#conda activate /usr/local/usrapps/metastrain/asoneto/bakta

core_cruncher="CoreCruncher/corecruncher_master.py"
#group=($(find Hodgkinia_core_input/ -mindepth 1 -type d -name '*sim'))
group=(
"Hodgkinia_core_input/both_highsim"
"Hodgkinia_core_input/2_core_highsim"
"Hodgkinia_core_input/both_medsim"
)

input_dir="${group[$(( LSB_JOBINDEX -  1))]}"

output_dir="Hodgkinia_core_output/$(basename $input_dir)"
mkdir -p "$output_dir"
python3 "$core_cruncher" -in "$input_dir" -out "$output_dir" -ext .faa -score 70
