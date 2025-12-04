#!/usr/bin/env python3 
import matplotlib.pyplot as plt
from Bio import Phylo
import os

input_dir = 'labeled_dna_trees/'
out_dir = 'labeled_dna_trees/plots/'
os.makedirs(out_dir, exist_ok=True)
for file in os.listdir(input_dir):
    if file.endswith('.treefile'):
        print(f'Processing {file}')
        file_path = os.path.join(input_dir, file)
        tree = Phylo.read(file_path, 'newick')
        n = len(tree.get_terminals())
        print(f'{n} terminals')
        height = min(n, 50)
        fig, ax = plt.subplots(figsize=(20,height))

        Phylo.draw(tree, axes=ax)
        out_file = os.path.join(out_dir, file)
        plt.savefig(f'{out_file}.pdf', dpi=300, bbox_inches="tight")
        plt.close(fig)
