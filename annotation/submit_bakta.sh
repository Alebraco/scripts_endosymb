#!/bin/bash
#BSUB -J bakta_wrapper
#BSUB -W 100:00
#BSUB -q bobay
#BSUB -e bakta_wrapper.err
#BSUB -o bakta_wrapper.out
#BSUB -R "rusage[mem=24GB] span[hosts=1]"
#BSUB -n 4

INPUT_DIR="ncbi_bacteria"
OUTDIR="bakta_results"

NGENOMES=$(find "$INPUT_DIR" -type f -name "*.fna" | wc -l)
if [ "$NGENOMES" -eq 0 ]; then
    echo "No .fna files found in $INPUT_DIR. Aborting."
    exit 1
fi

NJOBS=$(( (NGENOMES + 9) / 10 ))

echo "Found $NGENOMES genomes in $INPUT_DIR — submitting $NJOBS array jobs."
bsub -J "bakta[1-${NJOBS}]%30" \
     -n 8 \
     -R "rusage[mem=20GB] span[hosts=1]" \
     -e "bakta_logs/%J_%I.err" \
     -o "bakta_logs/%J_%I.out" \
     -W 50:00 \
     -env "INPUT_DIR='$INPUT_DIR',OUTDIR='$OUTDIR'" \
     bash scripts_endosymb/annotation/bakta_annotation.sh
