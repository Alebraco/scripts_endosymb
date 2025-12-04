#!/usr/bin/env python3

from Bio import Phylo
import os
tree_dir = 'tree_files/'
for tree in os.listdir(tree_dir):
    tree_path = os.path.join(tree_dir, tree)
    print(f'Processing {tree_path}')
    tree_info = Phylo.read(tree_path, 'newick')

    for clade in tree_info.find_clades():
        if clade.branch_length:
            if float(clade.branch_length) > 0.5:
                if clade.is_terminal():
                    print(f'>Accession: {clade.name}:{clade.branch_length}')
                else:
                    accessions = [leaf.name for leaf in clade.get_terminals()]
                    print(f'>Group: {accessions}:{clade.branch_length}')

                print('\n')

    print('-'*40)

