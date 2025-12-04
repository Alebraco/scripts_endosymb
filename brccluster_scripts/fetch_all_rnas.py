#!/usr/bin/env python3

import os
species_dirs = '16S/'
with open('16S_sequences.fna','w') as output:
    for sp in os.listdir(species_dirs):
        sp_path = os.path.join(species_dirs, sp)
        # variable sp contains the species name / name of directory
        # variable sp_path contains the path to the directory
        for file in os.listdir(sp_path): 
            # variable file contains the name of the file (accession + _16S.fna)
            accession = file.replace('_16S.fna', '')
            with open(os.path.join(sp_path, file), 'r') as f:
                counter = 0 
                # counter for number of 16S rRNA sequences

                for line in f:
                    if line.startswith('>16S_rRNA'):
                        sequence = next(f).strip()
                        counter += 1
                        name = f'{sp}_{accession}_{counter}'
                        output.write(f'>{name}\n{sequence}\n')


    
