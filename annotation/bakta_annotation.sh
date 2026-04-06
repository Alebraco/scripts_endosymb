#!/bin/bash

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation/

BAKTA_DB="/rs1/researchers/l/ljbobay/asoneto/endosymb/bakta_db/db/"
all_files=($(cat "$GENOME_LIST"))
start_index=$((($LSB_JOBINDEX - 1) * 100))
end_index=$(($start_index + 99))

for ((i=start_index; i<=end_index && i<${#all_files[@]}; i++)); do
    genome=${all_files[$i]}
    prefix=$(basename $genome .fna)
    subpath=${genome#$INPUT_DIR/}
    subpath=${subpath%.fna}

    scratch_out="/share/metastrain/asoneto/endosymb/bakta_results/$subpath"
    dest="$WORKDIR/$OUTDIR/$subpath"
    mkdir -p $scratch_out
    mkdir -p $dest

    echo "Batch $LSB_JOBINDEX: Annotating $prefix"
    bakta --db $BAKTA_DB \
        --output $scratch_out \
        --prefix $prefix \
        --threads 8 \
        --force \
        $genome

    cp -a $scratch_out/. $dest/
done

conda deactivate
