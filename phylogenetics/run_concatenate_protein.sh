#!/bin/bash
#BSUB -J concatenate[1-3]
#BSUB -o %J_%I.out
#BSUB -e %J_%I.err
#BSUB -n 1
#BSUB -R "rusage[mem=4GB]"
#BSUB -W 02:00

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

# Get the current job number (1, 2, or 3)
job_number=$LSB_JOBINDEX

echo "This is concatenate job number: $job_number"

if [ $job_number -eq 1 ]; then
    echo "Processing endosymb+relatives concatenates"
    python3 concatenate.py endosymb+relatives protein

elif [ $job_number -eq 2 ]; then
    echo "Processing endosymb_only concatenates"
    python3 concatenate.py endosymb_only protein

elif [ $job_number -eq 3 ]; then
    echo "Processing relatives_only concatenates"
    python3 concatenate.py relatives_only protein

else
    echo "Error: Unknown job number $job_number"
    exit 1
fi
