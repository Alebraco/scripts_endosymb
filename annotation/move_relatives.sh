#!/bin/bash
#BSUB -n 1
#BSUB -W 2:00
#BSUB -J move_relatives
#BSUB -o move_relatives.out
#BSUB -e move_relatives.err

annotations="bakta_results"
out_dir="relatives_only"

# Loop through each species directory
for sp in "$annotations"/*; do
    sp_name=$(basename "$sp")
    echo "Processing $sp_name"

    # Output Structure: relatives_only / file_type / species_name /
    target_gbff="$out_dir/gbff_files/$sp_name"
    target_genome="$out_dir/genomes/$sp_name"
    target_cds="$out_dir/cds/$sp_name"
    target_prot="$out_dir/proteins/$sp_name"

    mkdir -p "$target_gbff" "$target_genome" "$target_cds" "$target_prot"

    find "$sp" -name "*.gbff" -type f -exec cp {} "$target_gbff"/ \;
    find "$sp" -name "*.fna" -type f -exec cp {} "$target_genome"/ \;
    find "$sp" -name "*.faa" -type f -exec cp {} "$target_prot"/ \;
    find "$sp" -name "*.ffn" -type f -exec sh -c 'cp "$0" "${1}/$(basename "$0" .ffn).fna"' {} "$target_cds" \;
done