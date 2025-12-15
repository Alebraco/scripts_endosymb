#!/bin/bash

#BSUB -n 4
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -J iqtree[1-131]
#BSUB -W 96:00
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out
#BSUB -q bobay

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/tree

FILES=(endosymb+relatives/dna_concatenates/*.fasta \
	endosymb_only/dna_concatenates/*.fasta \
	relatives_only/dna_concatenates/*.fasta)
CONCATENATE=${FILES[$((LSB_JOBINDEX - 1))]}
TREE_DIR=dna_tree_results
GROUP=$(basename $(dirname $(dirname $CONCATENATE)))
SPECIES=$(basename $CONCATENATE .fasta)
SPECIES=${SPECIES#concatenate_}

OUTDIR="$GROUP/$TREE_DIR/$SPECIES"
mkdir -p $OUTDIR

TEMP_FASTA="$OUTDIR/${SPECIES}_tmp.fasta"
seqkit seq -g 0.5 "$CONCATENATE" > "$TEMP_FASTA"

SEQS=$(grep -c '>' "$TEMP_FASTA")

if [[ $SEQS -gt 3 ]]; then
	iqtree -s $TEMP_FASTA -m GTR+G -bb 1000 -T 4 -pre "$OUTDIR/$SPECIES"
elif [[ $SEQS -le 2 ]]; then
	echo "Skipping $SPECIES, only $SEQS sequences to compare."
else
	iqtree -s $TEMP_FASTA -m GTR+G -T 4 -pre "$OUTDIR/$SPECIES"
fi

rm "$TEMP_FASTA"
