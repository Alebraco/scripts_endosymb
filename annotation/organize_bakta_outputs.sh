#!/bin/bash

BASE="${1:-ncbi_bacteria}"
BAKTA_DIR="$BASE/bakta_results/genomes"
GENOMES_DIR="$BASE/genomes"
PROTEINS_DIR="$BASE/proteins"

if [ ! -d "$BAKTA_DIR" ]; then
    echo "Bakta results directory not found: $BAKTA_DIR"
    exit 1
fi

n_gff=0
n_faa=0

for sp_dir in "$BAKTA_DIR"/*/; do
    sp_name=$(basename "$sp_dir")
    mkdir -p "$PROTEINS_DIR/$sp_name"

    for accession_dir in "$sp_dir"/*/; do
        for gff_file in "$accession_dir"/*.gff3; do
            [ -f "$gff_file" ] || continue
            dest="$GENOMES_DIR/$sp_name/$(basename "${gff_file%.gff3}.gff")"
            cp "$gff_file" "$dest"
            n_gff=$((n_gff + 1))
        done

        # Copy .faa into proteins/
        for faa_file in "$accession_dir"/*.faa; do
            [ -f "$faa_file" ] || continue
            cp "$faa_file" "$PROTEINS_DIR/$sp_name/"
            n_faa=$((n_faa + 1))
        done
    done
done

echo "Organized $n_gff GFF files into $GENOMES_DIR/"
echo "Organized $n_faa FAA files into $PROTEINS_DIR/"
