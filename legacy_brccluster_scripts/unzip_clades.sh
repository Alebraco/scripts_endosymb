#!/bin/bash

# for file in clade_genomes/*/*.zip; do
# 	unzip -j $file "ncbi_dataset/data/*/*.faa" -d $(dirname $file)  
# 	mv $(dirname $file)
# done

for file in endosymb_cds/*/*.zip; do
    parent_dir=$(dirname "$file")
    unzip "$file" "ncbi_dataset/data/*" -d "$parent_dir"
    for extracted_file in "$parent_dir"/ncbi_dataset/data/*/*; do
        if [ -f "$extracted_file" ]; then
            subdir=$(basename "$(dirname "$extracted_file")")
            new_name="${subdir}.fna"
            mv "$extracted_file" "$parent_dir/$new_name"
        fi
    done
done
