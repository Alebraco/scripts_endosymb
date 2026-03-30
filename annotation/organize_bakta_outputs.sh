#!/bin/bash
# Organize Bakta annotation outputs into the directory layout expected by
# the feature collection pipeline.
#
# Bakta outputs:   {BASE}/bakta_results/genomes/{species}/{accession}/{accession}.*
# Expected layout: {BASE}/genomes/{species}/{accession}.gff   (co-located with .fna)
#                  {BASE}/proteins/{species}/{accession}.faa
#
# Usage: bash organize_bakta_outputs.sh [base_dir]
#   Defaults to ncbi_bacteria if no argument given.

BASE="ncbi_bacteria"
BAKTA_DIR="$BASE/bakta_results/genomes"
GENOMES_DIR="$BASE/genomes"
PROTEINS_DIR="$BASE/proteins"

if [ ! -d "$BAKTA_DIR" ]; then
    echo "Bakta results directory not found: $BAKTA_DIR"
    exit 1
fi

n_gff=0
n_faa=0

for sp_dir in "$BAKTA_DIR"/*/; do                  # iterate over each species subdirectory
    sp_name=$(basename "$sp_dir")                  # extract species name from the directory path
    mkdir -p "$PROTEINS_DIR/$sp_name"              # create matching species subdir under proteins/ (genome dir already exists)

    for accession_dir in "$sp_dir"/*/; do          # iterate over each accession subdirectory within the species
        for gff_file in "$accession_dir"/*.gff3; do
            [ -f "$gff_file" ] || continue             # skip if no match
            dest="$GENOMES_DIR/$sp_name/$(basename "${gff_file%.gff3}.gff")"  
            cp "$gff_file" "$dest"                     # copy the annotation file
            n_gff=$((n_gff + 1))                       # increment counter
        done

        # Copy .faa into proteins/
        for faa_file in "$accession_dir"/*.faa; do     
            [ -f "$faa_file" ] || continue             # skip if no match
            cp "$faa_file" "$PROTEINS_DIR/$sp_name/"   # copy the protein file
            n_faa=$((n_faa + 1))                       # increment counter
        done
    done
done

echo "Copied $n_gff GFF files into $GENOMES_DIR/"   # total GFF files moved
echo "Copied $n_faa FAA files into $PROTEINS_DIR/"  # total FAA files moved
