#!/bin/bash

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation/

BAKTA_DB="/rs1/researchers/l/ljbobay/asoneto/endosymb/archived/bakta_db/db/"
all_files=($(cat "$GENOME_LIST"))
# BATCH_SIZE controls how many genomes each array task processes (default 100)
BATCH_SIZE="${BATCH_SIZE:-100}"
start_index=$(( (LSB_JOBINDEX - 1) * BATCH_SIZE ))
end_index=$(( start_index + BATCH_SIZE - 1 ))

for ((i=start_index; i<=end_index && i<${#all_files[@]}; i++)); do
    genome=${all_files[$i]}
    prefix=$(basename $genome .fna)
    subpath=${genome#$INPUT_DIR/}
    subpath=${subpath%.fna}

    scratch_out="/share/metastrain/asoneto/endosymb/$OUTDIR/$subpath"
    dest="$WORKDIR/$OUTDIR/$subpath"
    mkdir -p $scratch_out
    mkdir -p $dest

    echo "Batch $LSB_JOBINDEX: Annotating $prefix"
    # --keep-contig-headers preserves the original FASTA contig IDs so they
    # match between the output .fna and .gff
    
    bakta --db $BAKTA_DB \
        --output $scratch_out \
        --prefix $prefix \
        --threads 8 \
        --keep-contig-headers \
        --force \
        $genome

    cp -a $scratch_out/. $dest/
done

conda deactivate
