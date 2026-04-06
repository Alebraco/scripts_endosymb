#!/bin/bash
#BSUB -J bakta_wrapper
#BSUB -W 0:30
#BSUB -q bobay
#BSUB -e bakta_wrapper.err
#BSUB -o bakta_wrapper.out
#BSUB -R "rusage[mem=1GB] span[hosts=1]"
#BSUB -n 1

WORKDIR="/rs1/researchers/l/ljbobay/asoneto/endosymb"
INPUT_DIR="ncbi_bacteria"
OUTDIR="bakta_results"
GENOME_LIST="$WORKDIR/bakta_genome_list.txt"

find "$WORKDIR/$INPUT_DIR" -type f -name "*.fna" | sort > "$GENOME_LIST"
NGENOMES=$(wc -l < "$GENOME_LIST")

if [ "$NGENOMES" -eq 0 ]; then
    echo "No .fna files found in $WORKDIR/$INPUT_DIR. Aborting."
    exit 1
fi

NJOBS=$(( (NGENOMES + 99) / 100 ))

mkdir -p "$WORKDIR/bakta_logs"

echo "Found $NGENOMES genomes — submitting $NJOBS array jobs."
bsub -J "bakta[1-${NJOBS}]%30" \
     -n 8 \
     -R "rusage[mem=20GB] span[hosts=1]" \
     -e "$WORKDIR/bakta_logs/%J_%I.err" \
     -o "$WORKDIR/bakta_logs/%J_%I.out" \
     -W 200:00 \
     -env "INPUT_DIR='$WORKDIR/$INPUT_DIR',OUTDIR='$OUTDIR',GENOME_LIST='$GENOME_LIST',WORKDIR='$WORKDIR'" \
     bash scripts_endosymb/annotation/bakta_annotation.sh
