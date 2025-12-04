#!/bin/bash

#BSUB -n 4
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -J iqtree
#BSUB -W 192:00
#BSUB -e %J.err
#BSUB -o %J.out

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/core

# FILES=(protein_concatenates/*)
# CONCATENATE=${FILES[$((LSB_JOBINDEX - 1))]}
# TREE_DIR=tree_results
# SPECIES=$(basename $CONCATENATE .fasta)
# SPECIES=${SPECIES#concatenate_}


# OUTDIR="$TREE_DIR/$SPECIES"
# mkdir -p $OUTDIR

# iqtree -s $CONCATENATE -m MFP -bb 1000 -T AUTO -pre "$OUTDIR/$SPECIES"

FILE="protein_concatenates/concatenate_Wolbachia_endosymbiont.fasta"
TREE_DIR=tree_results
SPECIES=$(basename $FILE .fasta)
SPECIES=${SPECIES#concatenate_}

OUTDIR="$TREE_DIR/$SPECIES"
mkdir -p $OUTDIR

iqtree -s $FILE -m MFP -T 4 -pre "$OUTDIR/$SPECIES"
