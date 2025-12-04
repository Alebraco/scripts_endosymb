#!/usr/bin/env python3
"""
Script to compute and visualize summary statistics for delta matrices
(e.g., genome size differences, GC content differences) across multiple groups
of endosymbionts and free-living relatives.

Usage:
    ./summary_across_groups.py <matrix_type> [mean/all]

Arguments:
    <matrix_type> : The type of delta matrix to compute (e.g., 'size', 'gc_genome', 'gc_third', 'gc_all', 'distance').
    [mean/all]    : Optional. Use 'mean' to aggregate species-level means for the combined distribution.
                    Defaults to 'all' (all pairwise distances included).

Outputs:
    - CSV file summarizing statistics (min, mean, median, max, std, count) for each group.
    - Boxplot visualizing the distribution of delta matrix values across groups.
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

if len(sys.argv) < 2:
    print(f'Usage: {sys.argv[0]} <matrix_type> <all(default)/mean>')
    sys.exit(1)

type = sys.argv[1]
mean = False
if len(sys.argv) == 3:
    if sys.argv[2] == 'mean':
        mean = True

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

# Initialize data storage
all_data = []
groups = ['endosymb+relatives', 'endosymb_only', 'relatives_only']
group_key = {
    'endosymb+relatives': 'Endosymbionts and Free-Living',
    'endosymb_only': 'Endosymbionts Only',
    'relatives_only': 'Free-Living Only'
}

# Process each group
for group in groups:
    gc_json = f'gc_codon_data_{group}.json'
    genome_json = f'gcsize_genome_{group}.json'

    # Load or compute datasets
    genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
    print('All data has been loaded')

    # Generate delta matrix based on the specified type
    if type in ['size', 'gc_genome']:
        matrix = delta_matrix(genome_dataset, type)
    else:
        matrix = delta_matrix(gc_dataset, type)
    print(f'Created {type} matrix.')

    # Extract values from the delta matrix
    all_values = []
    for species in matrix.keys():
        sp_name = ' '.join(species.split('_endosymbiont')[0].split('_'))
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i, j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i, j):
            # Skip comparisons within the same group for 'endosymb+relatives'
            if group == 'endosymb+relatives':
                if ('genomic' in ids[index1]) == ('genomic' in ids[index2]):
                    continue
            value = sp_matrix.iloc[index1, index2]
            sp_values.append(value)
        if mean:
            mean_sp = np.mean(sp_values)
            all_data.append({'group': group_names[group], 'value': mean_sp, 'species': sp_name})
        else:
            for value in sp_values:
                all_data.append({'group': group_names[group], 'value': value, 'species': sp_name})

# Create a DataFrame from the collected data
df = pd.DataFrame(all_data)

# Print and save summary statistics
print('--Summary Statistics--')
summary_table = df.groupby('group')['value'].agg(['min', 'mean', 'median', 'max', 'std', 'count']).round(4)
print(summary_table)
if mean:
    summary_table.to_csv(f'mean{type}_groups_stats.csv')
else:
    summary_table.to_csv(f'{type}_groups_stats.csv')

# Generate and save boxplot
plt.figure(figsize=(12, 8))
sns.boxplot(data=df, x='group', y='value', fliersize=2)
plt.title(f'{titles[type]}\nAcross Groups')
plt.xlabel(f'{titles[type]}')
if type == 'size':
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))

plt.tight_layout()
if mean:
    plt.savefig(f'all_{type}_boxplot-mean.pdf')
else:
    plt.savefig(f'all_{type}_boxplot.pdf')
plt.close()