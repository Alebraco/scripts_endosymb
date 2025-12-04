#!/usr/bin/env python3
"""
Script to compute and visualize summary statistics for delta matrices
(e.g., genome size differences, GC content differences) across species
or groups of endosymbionts and free-living relatives.

Usage:
    ./summary_data.py <group> <matrix_type> <sp/all>

Arguments:
    <group>       : The group to analyze (e.g., 'endosymb+relatives', 'relatives_only', 'endosymb_only').
    <matrix_type> : The type of delta matrix to compute (e.g., 'size', 'gc_genome', 'gc_third', 'gc_all', 'distance').
    <sp/all>      : 'sp' for species-specific analysis, 'all' for combined analysis.

Outputs:
    - Species-specific or combined boxplots summarizing the delta matrix values.
    - CSV file summarizing variation (standard deviation and variance) for each species.
    - PDF plots saved to the current directory.
"""

from delta_matrix import delta_matrix
from metadata_gcsize import genome_gcsize
from gc_codon_dict import gc_codon_dict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import os
import sys
import json

if len(sys.argv) != 4:
    print(f'Usage: {sys.argv[0]} <group> <matrix_type> <sp/all>')
    sys.exit(1)

group = sys.argv[1]
type = sys.argv[2] 
crit = sys.argv[3]

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

gc_json = f'gc_codon_data_{group}.json'
genome_json = f'gcsize_genome_{group}.json'

def load_or_compute(filename, compute_function, *args, **kwargs):
    """
    Load data from a file if it exists; otherwise, compute the data and save it to the file.

    Args:
        filename (str): Path to the file to load or save.
        compute_function (function): Function to compute the data if the file does not exist.
        *args: Positional arguments to pass to the compute function.
        **kwargs: Keyword arguments to pass to the compute function.

    Returns:
        dict: The loaded or computed data.
    """
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    else:
        data = compute_function(*args, **kwargs)
        with open(filename, 'w') as f:
            json.dump(data, f)
        return data

# Load or compute datasets
genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
print('All data has been loaded')

# Generate delta matrix based on the specified type
if type in ['size', 'gc_genome']:
    matrix = delta_matrix(genome_dataset, type)
else:
    matrix = delta_matrix(gc_dataset, type)
print(f'Created {type} matrix.')

if crit == 'sp':
    """
    Species-specific analysis:
    - Computes delta matrix values for each species.
    - Summarizes variation (standard deviation and variance) for each species.
    - Generates boxplots for each species and saves them as a PDF.
    """
    all_data = []
    for species in matrix.keys():
        sp_name = ' '.join(species.split('_endosymbiont')[0].split('_'))
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i, j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i, j):
            if group == 'endosymb+relatives':
                # Skip comparisons within the same group (e.g., genomic vs genomic)
                if ('genomic' in ids[index1]) == ('genomic' in ids[index2]):
                    continue
            value = sp_matrix.iloc[index1, index2]
            sp_values.append(value)
        
        for value in sp_values:
            all_data.append({'species': sp_name, 'value': value})
    
    # Save summary statistics to a CSV file
    df = pd.DataFrame(all_data)
    df_summary = df.groupby('species')['value'].agg(['std', 'var']).reset_index()
    df_summary.to_csv(f'variation_{group}_{type}.csv', index=False)

    # Generate and save boxplot
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=df, y='species', x='value', fliersize=3)
    plt.title(f'{titles[type]}\nAcross {group_names[group]}')
    plt.xlabel(f'{titles[type]}')
    if type == 'size':
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.tight_layout()
    plt.savefig(f'sp_{type}_{group}.pdf')
    plt.close()

elif crit == 'all':
    """
    Combined analysis:
    - Aggregates delta matrix values across all species.
    - Computes mean values for each species if required.
    - Generates a combined boxplot and saves it as a PDF.
    """
    all_values = []
    for species in matrix.keys():
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i, j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i, j):
            if group == 'endosymb+relatives':
                # Skip comparisons within the same group (e.g., genomic vs genomic)
                if ('genomic' in ids[index1]) == ('genomic' in ids[index2]):
                    continue
            value = sp_matrix.iloc[index1, index2]
            sp_values.append(value)
        mean_sp = np.mean(sp_values)
        all_values.append(mean_sp)
    
    # Generate and save combined boxplot
    plt.figure()
    sns.boxplot(data=all_values)
    plt.title(f'{titles[type]}\n{group_names[group]}', fontsize=14)
    plt.xlabel(f'{titles[type]}')
    if type == 'size':
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.savefig(f'all_{type}_{group}.pdf')
    plt.close()