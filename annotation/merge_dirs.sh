#!/bin/bash
#BSUB -n 1
#BSUB -W 2:00
#BSUB -J merge_dirs
#BSUB -o merge_dirs.out
#BSUB -e merge_dirs.err

types=("proteins" "genomes" "cds" "gbff_files")

for dir_type in "${types[@]}"; do

    echo "--- Merging $dir_type ---"

    # Define origin and destination based on the directory type
    origin=endosymb_only/$dir_type
    dest=endosymb+relatives/$dir_type

    # Check if directories exist before proceeding
    if [ ! -d "$origin" ] || [ ! -d "$dest" ]; then
        echo "Directories for $dir_type do not exist."
        continue
    fi

    find "$dest" -type d -maxdepth 1 -mindepth 1 -printf "%P\n" | while read -r subdir; do
        if [ -d "$origin/$subdir" ] && [ -d "$dest/$subdir" ]; then
            echo "merging $origin/$subdir/ with $dest/$subdir"
            cp "$origin/$subdir"/* "$dest/$subdir"/
        fi
    done
done