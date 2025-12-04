#!/bin/bash

#BSUB -n 1
#BSUB -W 96:00
#BSUB -J backtranslate[1-3]
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -q bobay

# Get the current job number (1, 2, or 3)
job_number=$LSB_JOBINDEX

echo "This is backtranslate job number: $job_number"

if [ $job_number -eq 1 ]; then
    echo 'Backtranslating alignments for endosymb_only/:'
    python3 -u backtranslate.py endosymb_only

elif [ $job_number -eq 2 ]; then
    echo 'Backtranslating alignments for endosymb+relatives/:'
    python3 -u backtranslate.py endosymb+relatives

elif [ $job_number -eq 3 ]; then
    echo 'Backtranslating alignments for relatives_only/:'
    python3 -u backtranslate.py relatives_only

else
    echo "Error: Unknown job number $job_number"
    exit 1
fi
