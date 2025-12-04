#!/bin/bash
origin=endosymb_only/gbff_files
dest=endosymb+relatives/gbff_files

find $dest -type d -maxdepth 1 -printf "%P\n" | while read -r subdir; do
    if [ -d $origin/$subdir ] && [ -d $dest/$subdir ]; then
        echo "merging $origin/$subdir/ with $dest/$subdir"
	cp $origin/$subdir/* $dest/$subdir/
    fi
done
