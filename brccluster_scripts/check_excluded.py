#!/usr/bin/env python3

import os
import subprocess

include = ['endosymbiont', 'Wolbachia', 'Buchnera', 'Blochmannia']
write_path = 'excluded_titles.txt'
blast_results = 'blast_results/'

with open(write_path, 'w') as out:
    for sp in os.listdir(blast_results): # for each species
        sp_path = os.path.join(blast_results, sp)    
        for file in os.listdir(sp_path): # for each 16S sequence
            file_path = os.path.join(sp_path, file)
            out.write(f'\n>{file} in {sp}\n')
            with open(file_path, 'r') as f:
                for line in f:
                    elements = line.strip().split('\t')
                    sseqid = elements[1]
                    pident = elements[2]
                    stitle = elements[6]
                    # add 'not' hereunder if excluding
                    if any(keyword.lower() in stitle.lower() for keyword in include):
                        reduced = stitle.split(';')[-1]
                        out.write(f'{pident}\t{reduced}\n')
                

    