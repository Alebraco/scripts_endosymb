#!/usr/bin/env python3

import os
species_dirs = '16S/'
os.makedirs('species_16S', exist_ok=True)

outdir = 'species_16S/'
for sp in os.listdir(species_dirs):
    sp_path = os.path.join(species_dirs, sp)
    with open(os.path.join(outdir, f'{sp}_16S.fna'), 'w') as output:
        for file in os.listdir(sp_path):
            accession = file.replace('_16S.fna', '')
            with open(os.path.join(sp_path, file), 'r') as f:
                counter = 0
                for line in f:
                    if line.startswith('>16S_rRNA'):
                        sequence = next(f).strip()
                        counter += 1
                        name = f'{sp}_{accession}_{counter}'
                        output.write(f'>{name}\n{sequence}\n')
