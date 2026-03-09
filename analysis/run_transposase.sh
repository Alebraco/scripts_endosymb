#!/bin/bash

#BSUB -J IS[1-111]%10
#BSUB -W 20:00
#BSUB -q bobay
#BSUB -e logs_IS/%J_%I.err
#BSUB -o logs_IS/%J_%I.out
#BSUB -R "rusage[mem=14GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/transposase

TRANS_DB=($TRANS_DB)
INPUT_DIRS=($INPUT_DIRS)

START_IND=$((($LSB_JOBINDEX - 1) * 60))
END_IND=$(($START_IND + 59))
ALL_FILES=($(find ${INPUT_DIRS[@]} -type f -name "*.faa" | sort))

for ((i=START_IND; i<=END_IND && i<${#ALL_FILES[@]}; i++)); do
    FAA_FILE=${ALL_FILES[$i]}
    ACCESSION=$(basename "$FAA_FILE" .faa)

    IN_DIR=""
    for cand in "${INPUT_DIRS[@]}"; do
        cand="${cand%/}"
        [[ "$FAA_FILE" == "$cand/"* ]] && IN_DIR="$cand" && break
    done

    if [[ -z "$IN_DIR" ]]; then
        OUTFILE="transposase/Unknown/${ACCESSION}.tsv"
    else
        # Get relative path from input directory to FAA file
        REL_PATH="${FAA_FILE#"$IN_DIR"/}"
        SP_NAME=$(dirname "$REL_PATH")
        
        [[ "$SP_NAME" == "." ]] && SP_NAME="Unknown"

        OUTPUT_ROOT="$(dirname "$IN_DIR")/transposase"
        OUTFILE="${OUTPUT_ROOT}/${SP_NAME}/${ACCESSION}.tsv"
    fi

    mkdir -p "$(dirname "$OUTFILE")"

    if [[ -f "$OUTFILE" ]]; then
        echo "Skipping existing transposase file: $OUTFILE"
        continue
    fi

    echo "Processing $FAA_FILE -> $OUTFILE"

    diamond blastp \
        --db $TRANS_DB \
        --query $FAA_FILE \
        --out $OUTFILE \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --threads 8 \
        --evalue 1e-5 \
        --max-target-seqs 1
done
