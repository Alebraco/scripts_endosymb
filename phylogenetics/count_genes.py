#!/usr/bin/env python3

import os
import pandas as pd
import sys

data = []
if len(sys.argv) != 2:
    print('Usage: script.py <group>')
    sys.exit(1)
main_dir = os.path.join(sys.argv[1],'core')
for sp in os.listdir(main_dir):
    sp_path = os.path.join(main_dir, sp)
    for genome in os.listdir(sp_path):
        genome_path = os.path.join(sp_path, genome)
        file = os.path.join(genome_path, 'summary.txt')
        with open(file, 'r') as f:
            first = f.readline()
            f.readline()
            third = f.readline()
        genes = first.split(':')[1].strip().split(' ')[0]
        genomes = third.split(':')[1].strip()
        name = genome.split('accns')[0]
        data.append({
            'Species': sp,
            'Reference': name,
            'Core_Genes': genes,
            'Genomes': genomes
        })
df = pd.DataFrame(data)
df.to_csv('core_gene_summary.csv', index=False)
