#!/bin/bash
#BSUB -J bakta_wrapper
#BSUB -W 0:30
#BSUB -q bobay
#BSUB -e bakta_wrapper.err
#BSUB -o bakta_wrapper.out
#BSUB -R "rusage[mem=1GB] span[hosts=1]"
#BSUB -n 1

WORKDIR="/rs1/researchers/l/ljbobay/asoneto/endosymb"
INPUT_DIR="${1:-ncbi_bacteria}"
OUTDIR="${2:-bakta_results}"
BATCH_SIZE="${3:-100}"
GENOME_LIST="$WORKDIR/bakta_genome_list_${INPUT_DIR//\//_}.txt"

if [ ! -f "$GENOME_LIST" ]; then
    find "$WORKDIR/$INPUT_DIR" -type f -name "*.fna" | sort > "$GENOME_LIST"
fi
NGENOMES=$(wc -l < "$GENOME_LIST")

if [ "$NGENOMES" -eq 0 ]; then
    echo "No .fna files found in $WORKDIR/$INPUT_DIR. Aborting."
    exit 1
fi

NJOBS=$(( (NGENOMES + BATCH_SIZE - 1) / BATCH_SIZE ))

mkdir -p "$WORKDIR/bakta_logs"

echo "Found $NGENOMES genomes — submitting $NJOBS array jobs (batch size $BATCH_SIZE)."
bsub -J "bakta[1-${NJOBS}]%30" \
     -n 8 \
     -R "rusage[mem=20GB] span[hosts=1]" \
     -e "$WORKDIR/bakta_logs/%J_%I.err" \
     -o "$WORKDIR/bakta_logs/%J_%I.out" \
     -W 200:00 \
     -q bobay \
     -env "INPUT_DIR='$WORKDIR/$INPUT_DIR',OUTDIR='$OUTDIR',BATCH_SIZE='$BATCH_SIZE',GENOME_LIST='$GENOME_LIST',WORKDIR='$WORKDIR'" \
     bash scripts_endosymb/annotation/bakta_annotation.sh
