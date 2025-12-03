#!/usr/bin/env python3
from gc_codon_dict import gc_codon_dict
from delta_matrix import delta_matrix
from distance_matrix import distance_matrix
from metadata_gcsize import genome_gcsize
from matrix_correlation import matrix_correlation
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
import sys
import pickle

if len(sys.argv) != 5:
    print('Usage: main.py <group> <matrix1_type> <matrix2_type> <mean/matrix>')
    sys.exit(1)
if '/' in sys.argv[1]:
    print('Defined group must be a string, not a folder.')
    sys.exit(1)

group = sys.argv[1]
matrix1_type = sys.argv[2]
matrix2_type = sys.argv[3]
mode = sys.argv[4]

# Update these when adding new species
gc_json = f'files/gc_codon_data_{group}.json'
genome_json = f'files/gcsize_genome_{group}.json'
distance_pkl = f'files/distances_{group}.pkl'

titles = {
    'distance': 'Patristic Distance',
    'size': 'ΔGenome Size',
    'gc_genome': 'Genome ΔGC%',
    'gc_third': '4D Site ΔGC%', 
    'gc_all': 'Core ΔGC%'
}
group_names = {
    'endosymb+relatives': 'Endosymbionts and Free-Living Relatives',
    'relatives_only': 'Free-Living Relatives Only',
    'endosymb_only': 'Endosymbionts Only'
}

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


def align_matrices(mat1, mat2):
    """Align two symmetric matrices to common genomes"""
    common_genomes = sorted(set(mat1.index) & set(mat2.index))
    
    if len(common_genomes) < 2:
        return None, None
    
    mat1_aligned = mat1.loc[common_genomes, common_genomes]
    mat2_aligned = mat2.loc[common_genomes, common_genomes]
    
    return mat1_aligned, mat2_aligned

def distance_scatterplot(x_df, y_df, species, correlation, p_value):
    # Flatten dataframes
    i, j = np.triu_indices_from(x_df, k=1)

    x = x_df.values[i,j]
    y = y_df.values[i,j]

        
    scatter = plt.scatter(x, y, c=y, alpha=0.6) #cmap='inferno'
    # plt.colorbar(scatter, label='ΔGC value')

    sns.regplot(x=x, y=y, scatter=False, color="red")
    plt.xlabel(titles[matrix1_type] if matrix1_type != 'distance' else titles[matrix2_type])
    plt.ylabel(titles[matrix2_type] if matrix2_type != 'distance' else titles[matrix1_type])
    plt.title(f"{species}: {titles[matrix1_type]} vs {titles[matrix2_type]}")

    plt.text(
        0.05, 0.95,
        f'Correlation: {round(correlation, 2)}\nP-value: {round(p_value, 4)}\nN: {len(x)}',
        ha = 'left', va = 'top', transform=plt.gca().transAxes
    )

#gc_dataset = load_or_compute(gc_json, gc_codon_dict, group)
genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
print('All data has been loaded.')

matrices = {}
for type in [matrix1_type, matrix2_type]:
    if type == 'distance':
        matrices['distance'] = load_or_compute_pickle(distance_pkl, distance_matrix, group)
    elif type in ['size', 'gc_genome']:
        matrices[type] = delta_matrix(genome_dataset, type)
    else:
        matrices[type] = delta_matrix(gc_dataset, type)
print(f'Created {matrix1_type} and {matrix2_type} matrices.')

n_species = len(matrices[matrix1_type])
print(n_species)
ncols = 6
nrows = n_species // ncols + (n_species % ncols > 0)
c = 1
results = {}

if mode == 'matrix':
    plt.figure(figsize=(5*6, 4*7))
elif mode == 'mean':
    plt.figure(figsize=(8,8))
for species in matrices[matrix1_type].keys(): 
    print(f'Processing {species}')
    if species in matrices[matrix2_type]:
        mat1 = matrices[matrix1_type][species]
        mat2 = matrices[matrix2_type][species]

        # Align matrices so we have the same dimensions / order
        mat1_aligned, mat2_aligned = align_matrices(mat1, mat2)
        if mat1_aligned is None:
            print(f'Less than 2 common genomes for {species}')
            continue
        

        if mode == 'mean':
            i, j = np.triu_indices_from(mat1_aligned, k=1)
            x = mat1_aligned.values[i,j]
            y = mat2_aligned.values[i,j]
            x_mean = np.mean(x)
            y_mean = np.mean(y)
            results[species] = {matrix1_type: x_mean, matrix2_type: y_mean}

        else:
            result = matrix_correlation(mat1_aligned, mat2_aligned)
            if result is None:
                print(f'\t>Less than 2 genomes for {species}')
                continue
            correlation, p_value, n = result
            results[species] = correlation
            plt.subplot(nrows, ncols, c)
            distance_scatterplot(mat1, mat2, species, correlation, p_value)
            c += 1
if mode == 'mean':
    plt.rcParams.update({'font.size': 11})
    df = pd.DataFrame.from_dict(results, orient = 'index')
    df.to_csv(f"{mode}_{group}_{matrix1_type}vs{matrix2_type}.csv")
    sns.regplot(data = df, x = matrix1_type, y = matrix2_type, ci = None)
    plt.xlabel(titles[matrix1_type])
    plt.ylabel(titles[matrix2_type])
    plt.title(f"{group_names[group]}: {titles[matrix1_type]} vs {titles[matrix2_type]}")    
    corr, p_value = spearmanr(df[matrix1_type], df[matrix2_type])
    plt.text(
    0.05, 0.95,
    f'Correlation: {round(corr, 2)}\n$R^2$: {round(corr**2, 2)}\np-value: {round(p_value, 4)}',
    ha = 'left', va = 'top', transform=plt.gca().transAxes
    )
if matrix1_type == 'size':
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
elif matrix2_type == 'size':
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))

plt.tight_layout()
plt.savefig(f"{mode}_{group}_{matrix1_type}vs{matrix2_type}.pdf")
plt.close()
