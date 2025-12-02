#!/usr/bin/env python3
from delta_matrix import delta_matrix
from metadata_gcsize import genome_gcsize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys

from gc_utils import (
    titles,
    group_names,
    load_or_compute,
    genome_gcsize_json_path,
)

if len(sys.argv) < 2:
    print(f'Usage: {sys.argv[0]} <matrix_type> <all/mean>')
    sys.exit(1)

mat_type = sys.argv[1]
mean = sys.argv[2].lower() == 'mean'

allowed = {'size', 'gc_genome'}
if mat_type not in allowed:
    print(f'Error: matrix_type must be one of {allowed}')
    sys.exit(1)

all_data = []
for group in group_names.keys():
    genome_json = genome_gcsize_json_path(group)

    genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
    print('All data has been loaded')

    matrix = delta_matrix(genome_dataset, mat_type)
    print(f'Created {mat_type} matrix.')

    for species in matrix.keys():
        sp_name = ' '.join(species.split('_endosymbiont')[0].split('_'))
        print(f'>Processing {species}:')
        sp_matrix = matrix[species]
        i,j = np.triu_indices_from(sp_matrix, k=1)
        ids = sp_matrix.index.tolist()

        sp_values = []

        for index1, index2 in zip(i,j):
            genome1 = ids[index1]
            genome2 = ids[index2]
            if group == 'endosymb+relatives':
                if ('genomic' in genome1) == ('genomic' in genome2):
                    continue
            value = sp_matrix.iloc[index1, index2]
            if mean:
                sp_values.append(value)
            else:
                all_data.append({
                    'group':group_names[group], 
                    'value':value, 
                    'species':sp_name,
                    'genome1':genome1,
                    'genome2':genome2})
        if mean:
            mean_sp = np.mean(sp_values)
            all_data.append({'group':group_names[group], 'value':mean_sp, 'species':sp_name})


df = pd.DataFrame(all_data)

print('--Summary Statistics--')
summary_table = df.groupby('group')['value'].agg(['min', 'mean', 'median', 'max', 'std', 'count']).round(4)
print(summary_table)
if mean:
    summary_table.to_csv(f'mean{mat_type}_groups_stats.csv')

    #Boxplots, only produced when using mean values
    plt.figure(figsize=(12,8))
    sns.boxplot(data = df, x = 'group', y = 'value', fliersize=2)
    plt.title(f'{titles[mat_type]}\nAcross Groups')
    plt.ylabel(f'{titles[mat_type]}')
    if mat_type == 'size':
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.tight_layout()
    if mean:
        plt.savefig(f'all_{mat_type}_boxplot-mean.pdf')
    else:
        plt.savefig(f'all_{mat_type}_boxplot.pdf')
    plt.close()

else:
    summary_table.to_csv(f'{mat_type}_groups_stats.csv')
    order_list = [
    'Endosymbionts and Free-Living Relatives',
    'Endosymbionts Only',
    'Free-Living Relatives Only']

    #Violin plots to visualize distribution of absolute values
    plt.rcParams.update({'font.size': 15})

    plt.figure(figsize=(12,8))
    sns.violinplot(data = df, x = 'group', y = 'value', inner = 'box', cut = 0, order = order_list)
    plt.title(f'{titles[mat_type]} - Distribution Across Groups')
    plt.xlabel('Group')
    plt.ylabel(f'{titles[mat_type]}')
    if mat_type == 'size':
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.tight_layout()
    plt.savefig(f'all_{mat_type}_violin.pdf')
    plt.close()

    #Produce boxplots of all species together, different color coding for groups
    df = df[df['group'].isin(['Endosymbionts Only', 'Free-Living Relatives Only'])]

    group_colors = {
    # 'Endosymbionts and Free-Living Relatives' : '#1f77b4',
    'Endosymbionts Only' : '#ff7f0e',
    'Free-Living Relatives Only' : '#2ca02c',
    }
    
    plt.figure(figsize=(12,24))
    sns.boxplot(
        data = df, 
        x = 'value', 
        y = 'species', 
        hue  = 'group', 
        palette=group_colors,
        hue_order = list(group_colors.keys()),
        fliersize=2, 
        dodge = True, 
        gap = 0.1
    )

    plt.title(f'{titles[mat_type]}\nAcross Groups')
    plt.xlabel(f'{titles[mat_type]}')
    if mat_type == 'size':
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.legend(fontsize='small', loc='upper right')
    plt.tight_layout()
    plt.savefig(f'all_{mat_type}_sp_group.pdf')
    plt.close()

    # Look for outliers in the endosymbionts group
    endosymb_data = df[df['group'] == 'Endosymbionts Only']
    Q1 = endosymb_data['value'].quantile(0.25)
    Q3 = endosymb_data['value'].quantile(0.75)
    IQR = Q3 - Q1
    outlier_threshold = Q3 + 1.5 * IQR

    outliers = endosymb_data[endosymb_data['value'] > outlier_threshold]
    outliers.to_csv(
        f'endosymb_outlier_pairs.csv', index=False)