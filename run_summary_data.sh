#!/bin/bash
#BSUB -J summary
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4GB]"
#BSUB -W 02:00

python3 summary_data.py endosymb+relatives size sp

python3 summary_data.py endosymb_only size sp

python3 summary_data.py relatives_only size sp
