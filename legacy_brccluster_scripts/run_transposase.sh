#!/bin/bash
#BSUB -J IS[1-111]
#BSUB -W 20:00
#BSUB -q bobay
#BSUB -e logs_diamond/%J_%I.err
#BSUB -o logs_diamond/%J_%I.out
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/transposase

TRANS_DB="files/IS_db.dmnd"
INPUT_DIRS=('endosymb_only/faa_files' 'relatives_only/faa_files')

START_IND=$((($LSB_JOBINDEX - 1) * 60))
END_IND=$(($START_IND + 59))
ALL_FILES=($(find ${INPUT_DIRS[@]} -type f -name "*.faa" | sort))

for ((i=START_IND; i<=END_IND && i<${#ALL_FILES[@]}; i++)); do
    FAA_FILE=${ALL_FILES[$i]}
    ACCESSION=$(basename $FAA_FILE .faa)
    SP_NAME=$(basename $(dirname $FAA_FILE))
    GROUP=$(basename $(dirname $(dirname $(dirname $FAA_FILE))))

    TRANS_DIR="${GROUP}/transposase"
    OUTFILE="${TRANS_DIR}/${SP_NAME}/${ACCESSION}.tsv"
    mkdir -p $(dirname $OUTFILE)

    diamond blastp \
        --db $TRANS_DB \
        --query $FAA_FILE \
        --out $OUTFILE \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --threads 8 \
        --evalue 1e-5 \
        --max-target-seqs 1
done