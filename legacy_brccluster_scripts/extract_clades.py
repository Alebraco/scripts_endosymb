#!/usr/bin/env python3

import os
from Bio import SeqIO

dna_location = 'species_16S/'
outdir = 'clades/'
os.makedirs(outdir, exist_ok=True)

species_dict={}
with open('main_clades_ids.txt', 'r') as f:
    for line in f:
        if line.startswith('>'):
            sp_fna = line.strip().split('>')[1]
            species_dict[sp_fna] = []
        elif line != '\n':
            accession = line.strip()
            species_dict[sp_fna].append(accession)

for sp, accn in species_dict.items():
    if os.path.exists(dna_location + sp):
        path_to_fna = os.path.join(dna_location, sp)
        with open(path_to_fna, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                if record.id in accn:
                    sp_name = sp.split('.fna')[0]
                    clade_dir = os.path.join(outdir, sp_name)
                    os.makedirs(clade_dir, exist_ok=True)
                    output_file = os.path.join(clade_dir, f'{record.id}.fna')
                    with open(output_file, 'w') as out:
                        SeqIO.write(record, out, 'fasta')
                        print(f'Extracted {record.id} from {path_to_fna} to {output_file}')