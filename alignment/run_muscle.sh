#!/bin/bash
#BSUB -q bobay
#BSUB -n 1
#BSUB -J muscle[1-3]
#BSUB -W 96:00
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out
#BSUB -R "rusage[mem=20GB] span[hosts=1]"

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/tree/


job_number=$LSB_JOBINDEX


echo "This is job number: $job_number"

if [ $job_number -eq 1 ]; then
    echo 'Processing endosymbionts only:'
    python3 -u muscle.py endosymb_only

elif [ $job_number -eq 2 ]; then
    echo 'Processing endosymbionts and relatives:'
    python3 -u muscle.py endosymb+relatives

elif [ $job_number -eq 3 ]; then
    echo 'Processing relatives only:'
    python3 -u muscle.py relatives_only

else
    echo "Error: Unknown job number $job_number"
    exit 1
fi
