#!/bin/bash

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation/

BAKTA_DB="/rs1/researchers/l/ljbobay/asoneto/endosymb/bakta_db/db/"
INPUT_DIR=($INPUT_DIR)
OUTDIR=($OUTDIR)
OUTPATH="$INPUT_DIR/$OUTDIR"

GENOMES=($(find $INPUT_DIR -type f -name "*.fna" | sort))

START_INDEX=$((($LSB_JOBINDEX - 1) * 10))
END_INDEX=$(($START_INDEX + 9))

for ((i=START_INDEX; i<=END_INDEX && i<${#GENOMES[@]}; i++)); do
    GENOME=${GENOMES[$i]}
    SUBPATH=${GENOME#$INPUT_DIR/}
    SUBPATH=${SUBPATH%.fna}
    ANNOTPATH=$OUTPATH/$SUBPATH/

    PREFIX=$(basename $GENOME .fna)

    echo "Batch $LSB_JOBINDEX: Annotating $PREFIX"
    bakta --db $BAKTA_DB \
        --output $ANNOTPATH \
        --prefix $PREFIX \
        --threads 8 \
        $GENOME
done

conda deactivate
