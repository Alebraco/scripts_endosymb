#!/usr/bin/env python3
from metadata_gcsize import genome_gcsize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys

from gc_utils import (
    abs_titles as titles,
    group_names,
    load_or_compute,
    genome_gcsize_json_path,
)

if len(sys.argv) < 3:
    print(f'Usage: {sys.argv[0]} <type> <mean/all>')
    sys.exit(1)

mat_type = sys.argv[1]
mean = sys.argv[2] == 'mean'
allowed = {'size', 'gc_genome'}
if mat_type not in allowed:
    print(f"Error: Supported types are: {', '.join(allowed)}")
    sys.exit(1)

all_data = []

for group in group_names.keys():
    genome_json = genome_gcsize_json_path(group)

    genome_dataset = load_or_compute(genome_json, genome_gcsize, group)

    data_dict = {}
    for sp, genomes in genome_dataset.items():
        sp_name = ' '.join(sp.split('_endosymbiont')[0].split('_'))
        data_dict = {genome:data.get(mat_type, 0) for genome,data in genomes.items()}
        if not mean:
            for value in data_dict.values():
                all_data.append({'group': group_names[group], 'value': value, 'species': sp_name})
        else:
            value = np.mean(list(data_dict.values()))
        all_data.append({'group': group_names[group], 'value': value, 'species': sp_name})

df = pd.DataFrame(all_data)
group_colors = {
'Endosymbionts Only' : '#ff7f0e',
'Free-Living Relatives Only' : '#2ca02c',
}

if mean is False:
#Plot all species, when mean = False
    plt.rcParams.update({'font.size': 15})

    df = df[df['group'].isin(['Endosymbionts Only', 'Free-Living Relatives Only'])]
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
    plt.savefig(f'all_abs{mat_type}_sp_group.pdf')
    plt.close()
else:
#When mean = True, plot mean genome size of endosymb_only vs. endosymb+relatives
    plt.rcParams.update({'font.size': 15})

    pivot_df = df.pivot_table(index='species', columns='group', values='value').reset_index()
    pivot_df.to_csv(f'mean_abs{mat_type}.csv', index = False)
    plt.figure(figsize=(10,8))
    sns.scatterplot(
        data = pivot_df, 
        x = 'Endosymbionts Only', 
        y = 'Free-Living Relatives Only'
    )

    x_max, x_min = pivot_df['Endosymbionts Only'].max(), pivot_df['Endosymbionts Only'].min()
    y_max, y_min = pivot_df['Free-Living Relatives Only'].max(), pivot_df['Free-Living Relatives Only'].min()
    common_min = min(x_min, y_min)
    common_max = max(x_max, y_max)
    padding = (common_max - common_min) * 0.1
    common_max, common_min = common_max + padding, common_min - padding
    plt.xlim(common_min, common_max)
    plt.ylim(common_min, common_max)

    plt.axline((0, 0), slope=1, color='red', linestyle='--', linewidth=1, label='y=x')
    plt.title(f'Mean {titles[mat_type]}\nEndoymbionts vs Free-Living Relatives')
    plt.xlabel('Endosymbionts')
    plt.ylabel('Free-Living Relatives')
    if mat_type == 'size':
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))
    plt.tight_layout()
    plt.savefig(f'scatter_mean_abs{mat_type}.pdf')
    plt.close()