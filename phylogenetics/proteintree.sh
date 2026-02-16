#!/bin/bash

#BSUB -n 4
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -J "iqtree[1-131]%5"
#BSUB -W 96:00
#BSUB -e tree_logs/%J_%I.err
#BSUB -o tree_logs/%J_%I.out
#BSUB -q bobay

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/tree

ROOT="/rs1/researchers/l/ljbobay/asoneto/endosymb"

# 1. Look for the new protein concatenates (.faa)
FILES=($ROOT/endosymb+relatives/protein_concatenates/*.faa \
    $ROOT/endosymb_only/protein_concatenates/*.faa \
    $ROOT/relatives_only/protein_concatenates/*.faa)
    
CONCATENATE=${FILES[$((LSB_JOBINDEX - 1))]}

# 2. Output to a new directory
TREE_DIR=protein_tree_results
GROUP=$(basename $(dirname $(dirname $CONCATENATE)))

SPECIES=$(basename $CONCATENATE .faa)
SPECIES=${SPECIES#concatenate_}

echo "Processing Job Index: $LSB_JOBINDEX"
echo "Target File: $CONCATENATE"

if [[ ! -f "$CONCATENATE" ]]; then
    echo "ERROR: The file $CONCATENATE cannot be found by the script."
    exit 1
fi

OUTDIR="$ROOT/$GROUP/$TREE_DIR/$SPECIES"
TREEFILE="$OUTDIR/${SPECIES}.treefile"

if [[ -f "$TREEFILE" ]]; then
    echo "Skipping $SPECIES: $TREEFILE already exists."
    exit 0
fi

mkdir -p $OUTDIR/TEMP_FILES

TEMP_FASTA="$OUTDIR/TEMP_FILES/${SPECIES}_tmp.faa"

trimal -in "$CONCATENATE" -out "$TEMP_FASTA" -fasta -seqoverlap 50 -resoverlap 0.5

SEQS=$(grep -c '>' "$TEMP_FASTA")

if [[ $SEQS -gt 3 ]]; then
    iqtree -s $TEMP_FASTA -m MFP -bb 1000 -T 4 -pre "$OUTDIR/$SPECIES"
elif [[ $SEQS -le 2 ]]; then
    echo "Skipping $SPECIES, only $SEQS sequences to compare."
else
    iqtree -s $TEMP_FASTA -m MFP -T 4 -pre "$OUTDIR/$SPECIES"
fi