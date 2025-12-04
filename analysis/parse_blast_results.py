#!/usr/bin/env python3

import os
import subprocess

# exclude = ['endosymbiont', 'Wolbachia', 'Buchnera', 'Blochmannia']
exclude = ['uncultured', 'unclassified', 'unidentified']
related_species = 'related_species/'
os.makedirs(related_species, exist_ok=True)
blast_results = 'blast_results/'

for sp in os.listdir(blast_results): # for each species
    sp_path = os.path.join(blast_results, sp)
    related_sp_dir = os.path.join(related_species, sp)
    os.makedirs(related_sp_dir, exist_ok=True)
    
    for file in os.listdir(sp_path): # for each 16S sequence
        file_path = os.path.join(sp_path, file)
        write_path = os.path.join(related_species, sp, file)
        with open(write_path, 'w') as out:
            with open(file_path, 'r') as f:
                for line in f:
                    elements = line.strip().split('\t')
                    sseqid = elements[1]
                    pident = elements[2]
                    stitle = elements[6]
                    # add 'not' hereunder if excluding
                    if not any(keyword.lower() in stitle.lower() for keyword in exclude):
                        out.write(f'{sseqid}\t{pident}\n')

        cleaned_path = write_path + '.cleaned'
        subprocess.run(f"awk -F'\t' '{{split($1,a,\".\"); print a[1] \"\t\" $2}}' {write_path} | sort -u > {cleaned_path}", shell=True)
        os.replace(cleaned_path, write_path)

                

    
