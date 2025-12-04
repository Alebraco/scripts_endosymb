#!/usr/bin/env python3
import os
import sys
import json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from gc_metadata_size import genome_gcsize
from gc_delta_matrix import delta_matrix
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from gc_utils import (
    titles,
    group_names,
    load_or_compute,
    genome_gcsize_json_path,
)

if len(sys.argv) != 4:
    print(f'Usage: {sys.argv[0]} <group> <matrix_type> <sp/all>')
    sys.exit(1)

group = sys.argv[1]
mat_type = sys.argv[2] 
crit = sys.argv[3]

allowed = {'size', 'gc_genome'}
if mat_type not in allowed:
    print(f'Error: matrix_type must be one of {allowed}.')
    sys.exit(1)

genome_json = genome_gcsize_json_path(group)

genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
print('All data has been loaded')

matrix = delta_matrix(genome_dataset, mat_type)
print(f'Created {mat_type} matrix.')

if crit == 'sp':
    all_data = []
    for species in matrix.keys():
        sp_name = ' '.join(species.split('_endosymbiont')[0].split('_'))
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i,j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i,j):
            if group == 'endosymb+relatives':
                if ('genomic' in ids[index1]) == ('genomic' in ids[index2]):
                    continue
            value = sp_matrix.iloc[index1, index2]
            sp_values.append(value)
        
        for value in sp_values:
            all_data.append({'species':sp_name, 'value': value})
        
    df = pd.DataFrame(all_data)
    df_summary = df.groupby('species')['value'].agg(['std', 'var']).reset_index()
    df_summary.to_csv(f'variation_{group}_{mat_type}.csv', index = False)

    plt.figure(figsize=(12,8))
    sns.boxplot(data = df, y = 'species', x ='value', fliersize=3)
    plt.title(f'{titles[mat_type]}\nAcross {group_names[group]}')
    plt.xlabel(f'{titles[mat_type]}')
    if mat_type == 'size':
        plt.gca().xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb')
            )
    plt.tight_layout()
    plt.savefig(f'sp_{mat_type}_{group}.pdf')
    plt.close()

elif crit == 'all':
    all_values = []
    for species in matrix.keys():
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i,j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i,j):
            if group == 'endosymb+relatives':
                if ('genomic' in ids[index1]) == ('genomic' in ids[index2]):
                    continue
            value = sp_matrix.iloc[index1, index2]
            sp_values.append(value)
        mean_sp = np.mean(sp_values)
        all_values.append(mean_sp)
    
    plt.figure()
    sns.boxplot(data = all_values)
    plt.title(f'{titles[mat_type]}\n{group_names[group]}', fontsize=14)
    plt.xlabel(f'{titles[mat_type]}')
    if mat_type == 'size':
        plt.gca().xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb')
            )
    plt.savefig(f'all_{mat_type}_{group}.pdf')
    plt.close()