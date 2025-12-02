#!/usr/bin/env python3
import pickle
import sys
import os
import json

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt

from delta_matrix import delta_matrix
from distance_matrix import distance_matrix
from GC_pipeline.metadata_gcsize import genome_gcsize

# Calculates (Delta / Max Size of Pair) / Patristic Distance
# Percentage change in metric per unit of evolutionary divergence

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <matrix (gc_genome/size)>')
    sys.exit(1)
matrix_type = sys.argv[1]

def load_or_compute(filename, compute_function, *args, **kwargs):
    if os.path.isfile(filename):
        with open(filename, 'r') as f: return json.load(f)
    else:
        data = compute_function(*args, **kwargs)
        with open(filename, 'w') as f: json.dump(data, f)
        return data

def load_or_compute_pickle(filename, compute_function, *args, **kwargs):
    if os.path.isfile(filename):
        with open(filename, 'rb') as f: return pickle.load(f)
    else:
        data = compute_function(*args, **kwargs)
        with open(filename, 'wb') as f: pickle.dump(data, f)
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
group_colors = {
    'Endosymbionts and Free-Living Relatives' : '#1f77b4',
    'Endosymbionts Only' : '#ff7f0e',
    'Free-Living Relatives Only' : '#2ca02c'
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

    for species in distance_matrices.keys(): 
        sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
        distance_data = distance_matrices[species]
        delta_data = genome_matrices[species]
        genome_data = genome_dataset.get(species)

        if genome_data is None:
            continue

        i,j = np.triu_indices_from(distance_data, k=1)
        ids = distance_data.index.tolist()
        sp_ratios = []
        
        for index1, index2 in zip(i, j):
            id1 = ids[index1]
            id2 = ids[index2]
            
            if group == 'endosymb+relatives':
                if ('genomic' in id1) == ('genomic' in id2):
                    continue
            
            distance = distance_data.loc[id1, id2]
            if distance > 0:
                delta_val = delta_data.loc[id1, id2]

                if matrix_type == 'size':
                    val1 = genome_data[id1][matrix_type]
                    val2 = genome_data[id2][matrix_type]
                    max_size = max(val1, val2)
                    if max_size == 0: 
                        continue

                    percent_diff = delta_val / max_size
                    ratio = percent_diff / distance
                else:
                    ratio = delta_val / distance

                sp_ratios.append(ratio)
                all_ratios.append({
                    'species': sp_name, 
                    'group': group_names[group], 
                    'ratio': ratio
                })

        print(f'{sp_name} ({group}): {len(sp_ratios)} ratios')

final_df = pd.DataFrame(all_ratios)
print(f'Total number of ratios: {len(final_df)}')

outdir = os.path.join('ratio_plots', matrix_type)
os.makedirs(outdir, exist_ok=True)

sp_median_df = final_df.groupby(['group','species'])['ratio'].median().reset_index()

# Filter zeros for log scale
sp_median_df = sp_median_df[sp_median_df['ratio'] > 0]

plt.figure(figsize=(10, 6))

sns.histplot(
    data=sp_median_df,
    x='ratio',
    hue='group',
    palette=group_colors,
    element="step",
    stat="proportion",
    common_norm=False,
    fill=True,
    alpha=0.3,
    hue_order=list(group_colors.keys()),
    log_scale=True
    )

plt.title(f'Rate of Evolution: {titles[matrix_type]}')
if matrix_type == 'size':
    plt.xlabel(f'Log Ratio: Fractional Change (Δ Size / Max Size) per substitution)')
elif matrix_type == 'gc_genome':
    plt.xlabel(f'Log Ratio: {titles[matrix_type]} / Substitution')
plt.ylabel('Proportion')
plt.tight_layout()
plt.savefig(os.path.join('ratio_plots', f'all_{matrix_type}_ratio_hist.pdf'))
plt.close()