#!/usr/bin/env python3
from delta_matrix import delta_matrix
from distance_matrix import distance_matrix
from metadata_gcsize import genome_gcsize
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import subprocess
import pickle
import sys
import os
import json

#The objective of this script is to visualize the distribution of the ratios (ΔGC/Size / Patristic Distance)
#This will be done for each species, each pair of accessions will have a ratio of their own.
#Also, this will be carried out for each group (relatives_only, endosymb_only, endosymb+relatives)
#For the last group, we will only compare endosymbionts against relatives. The in-group comparisons are not needed.
#We can have a group-level distribution and a species-level distribution of ratios
#The latter will be the most useful, since it is more specific.

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <matrix (gc_genome/size)>')
    sys.exit(1)
matrix_type = sys.argv[1]

def load_or_compute(filename, compute_function, *args, **kwargs):
    """Load data from file if exists, otherwise compute and save."""
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    else:
        data = compute_function(*args, **kwargs)
        with open(filename, 'w') as f:
            json.dump(data, f)
        return data

def load_or_compute_pickle(filename, compute_function, *args, **kwargs):
    if os.path.isfile(filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)
    else:
        data = compute_function(*args, **kwargs)
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
        return data
    
titles = {
    'distance': 'Patristic Distance',
    'size': 'ΔGenome Size',
    'gc_genome': 'Genome ΔGC%',
}
group_names = {
    'endosymb+relatives': 'Endosymbionts and Free-Living Relatives',
    'relatives_only': 'Free-Living Relatives Only',
    'endosymb_only': 'Endosymbionts Only'
}
all_ratios = []
metric = titles[matrix_type]


for group in group_names.keys():
    print(f'>Processing {group}:')
    gc_json = f'files/gc_codon_data_{group}.json'
    genome_json = f'files/gcsize_genome_{group}.json'
    distance_pkl = f'files/distances_{group}.pkl'

    genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
    genome_matrices = delta_matrix(genome_dataset, matrix_type)
    distance_matrices = load_or_compute_pickle(distance_pkl, distance_matrix, group)
    print('All data has been loaded.')

    tree_loc = os.path.join(group, 'dna_tree_results')

    for species in distance_matrices.keys(): #Use distance matrices keys as they include fewer species than genome matrices
        sp_name = ' '.join(species.split('_endosymbiont')[0].split('_'))
        sp_distmatrix = distance_matrices[species]
        sp_genmatrix = genome_matrices[species]
        i,j = np.triu_indices_from(sp_distmatrix, k=1)
        ids = sp_distmatrix.index.tolist()
        sp_ratios = []
        #If distance matrix is available, the species must have a tree file
        sp_tree_log = os.path.join(tree_loc, species, f'{species}.log')
        #Retrieve the number of nucleotides used in the alignment

        totalnt_out = subprocess.run([
            'grep', '-m 1', '-oP', r'Alignment has \d+ sequences with \K\d+', sp_tree_log],
            capture_output=True,
            text=True
            )
        totalnt = int(totalnt_out.stdout)
        print(totalnt)

        for index1, index2 in zip(i, j):
            id1 = ids[index1]
            id2 = ids[index2]
            if group == 'endosymb+relatives':
                # Skip comparisons within the same group (e.g., genomic vs genomic)
                if ('genomic' in id1) == ('genomic' in id2):
                    continue
            distance = sp_distmatrix.loc[id1, id2]
            if distance == 0:
                continue
            gmetric = sp_genmatrix.loc[id1, id2]
            
            substitutions = distance * totalnt
            ratio = gmetric / substitutions

            sp_ratios.append(ratio)
            all_ratios.append({'species': sp_name, 'group': group_names[group], 'ratio': ratio, 'pair': f'{id1}vs{id2}'})
        
        print(f'{sp_name} ({group}): {len(sp_ratios)} ratios')

final_df = pd.DataFrame(all_ratios)
print(f'Total number of ratios: {len(final_df)}')

outdir = os.path.join('ratio_plots', matrix_type)
os.makedirs(outdir, exist_ok=True)

group_colors = {
'Endosymbionts and Free-Living Relatives' : '#1f77b4',
'Endosymbionts Only' : '#ff7f0e',
'Free-Living Relatives Only' : '#2ca02c',
}

sp_median_df = final_df.groupby(['group','species'])['ratio'].median().reset_index()
sp_median_df['log_median_ratio'] = np.log10(sp_median_df['ratio'])

plt.figure(figsize=(10, 6))
sns.kdeplot(
    data=sp_median_df,
    x='log_median_ratio',
    # x='ratio',
    hue='group',
    palette=group_colors,
    fill=True,
    alpha=0.5,
    hue_order=list(group_colors.keys())
    )

plt.title(f'Distribution Median Rates ({matrix_type})')
plt.xlabel(f'Log Median Ratio ({metric} / Substitutions)')
plt.ylabel('Density (Frequency of Species)')
plt.tight_layout()
plot_filename = os.path.join(outdir, f'{matrix_type}_ratio.pdf')
plt.savefig(os.path.join('ratio_plots', f'all_{matrix_type}_ratio.pdf'))
plt.close()

# for species in final_df['species'].unique():
#     sp_df = final_df[final_df['species'] == species]
    
#     if len(sp_df['group'].unique()) == 3:
#         plt.figure(figsize=(8, 6))
#         sns.violinplot(
#             data = sp_df,
#             x = 'group',
#             y = 'ratio',
#             hue = 'group',
#             palette = group_colors,
#             hue_order = list(group_colors.keys())
#         )

#         plt.title(f'Rate of Change Distribution for {species}')
#         plt.xlabel('Group')
#         plt.ylabel(f'log-Ratio ({metric} / Distance)')
#         plt.yscale('log')

#         plt.tight_layout()
#         sp_safe = species.replace(' ', '_')
#         plot_filename = os.path.join(outdir, f"{sp_safe}_ratio.pdf")
#         plt.savefig(plot_filename)
#         plt.close()
