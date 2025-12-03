#!/usr/bin/env python3 
import os
from Bio import Phylo
import pandas as pd

df = pd.read_csv('metadata_gcsize.tsv', sep = '\t')

tree_dir = 'dna_tree_results/'
out_dir = 'labeled_dna_trees/'

for sp in os.listdir(tree_dir):
    print(f'Processing {sp}')
    tree_path = os.path.join(tree_dir, sp, f'{sp}.treefile')
    if os.path.isfile(tree_path):
        tree = Phylo.read(tree_path, 'newick')

        for leaf in tree.get_terminals():
            if leaf.name in df['Genome'].values:
                row = df[df['Genome'] == leaf.name].iloc[0]
                size, gc = row[['Size', 'GC%']]
                leaf.name = f'{leaf.name}/Size_{size}/GC_{gc}'

        out_file = os.path.join(out_dir, f'labeled_{sp}.treefile')
        Phylo.write(tree, out_file, 'newick')
    else:
        print(f'No tree file for {sp}')
        continue

