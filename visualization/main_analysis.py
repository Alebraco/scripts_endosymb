#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from delta_matrix import delta_matrix
from distance_matrix import distance_matrix
from gcsize_dict import genome_gcsize
from matrix_correlation import matrix_correlation
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from utils import (
    titles,
    group_names,
    load_or_compute,
    load_or_compute_pickle,
    genome_gcsize_json_path,
    files_dir,
)

if len(sys.argv) != 5:
    print('Usage: main_analysis.py <group> <matrix1_type> <matrix2_type> <mean/matrix>')
    sys.exit(1)
if '/' in sys.argv[1]:
    print('Defined group must be a string, not a folder.')
    sys.exit(1)

group = sys.argv[1]
matrix1_type = sys.argv[2]
matrix2_type = sys.argv[3]
mode = sys.argv[4]

allowed = {"distance", "size", "gc_genome"}


if matrix1_type not in allowed or matrix2_type not in allowed:
    print(f"Error: Supported matrix types are: {', '.join(allowed)}")
    sys.exit(1)

genome_json = genome_gcsize_json_path(group)
distance_pkl = os.path.join(files_dir, f'distances_{group}.pkl')

def align_matrices(mat1, mat2):
    '''Align two symmetric matrices to common genomes'''
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
    sns.regplot(x=x, y=y, scatter=False, color="red")
    plt.xlabel(titles[matrix1_type] if matrix1_type != 'distance' else titles[matrix2_type])
    plt.ylabel(titles[matrix2_type] if matrix2_type != 'distance' else titles[matrix1_type])
    plt.title(f"{species}: {titles[matrix1_type]} vs {titles[matrix2_type]}")

    plt.text(
        0.05, 0.95,
        f'Correlation: {round(correlation, 2)}\nP-value: {round(p_value, 4)}\nN: {len(x)}',
        ha = 'left', va = 'top', transform=plt.gca().transAxes
    )

genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
print('All data has been loaded.')

matrices = {}
for t in [matrix1_type, matrix2_type]:
    if t == 'distance':
        matrices['distance'] = load_or_compute_pickle(distance_pkl, distance_matrix, group)
    elif t in ['size', 'gc_genome']:
        matrices[t] = delta_matrix(genome_dataset, t)
    else:
        print(f'Unsupported matrix type: {t}')

print(f'Created {matrix1_type} and {matrix2_type} matrices.')

n_species = len(matrices[matrix1_type])
print(n_species)
ncols = 6
nrows = n_species // ncols + (n_species % ncols > 0)
c = 1
results = {}

if mode == 'matrix':
    plt.figure(figsize=(5 * 6, 4 * 7))
elif mode == 'mean':
    plt.figure(figsize=(8, 8))
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