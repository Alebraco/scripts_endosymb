#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from delta_matrix import delta_matrix
from gcsize_dict import genome_gcsize
from utils import files_dir, genome_gcsize_json_path, load_or_compute, group_names

plot_dir = os.path.join('plots', 'igs_lengths')
os.makedirs(plot_dir, exist_ok=True)
group_colors = {
    'Endosymbionts' : '#FC8D62',
    'Free-Living Relatives' : '#66C2A5',
}

#Box Plots with all data points
df_boxplot = pd.read_csv(os.path.join(files_dir, 'all_IGS_data.csv'))
df_boxplot = df_boxplot[df_boxplot['Group'].isin(['Endosymbionts','Free-Living Relatives'])]
df_boxplot_median = df_boxplot.groupby(['Group', 'Species', 'File'])['IGS_Size'].median().reset_index()

plt.figure(figsize=(12,12))
sns.boxplot(data = df_boxplot_median, x = 'IGS_Size', y = 'Species', hue = 'Group',
            palette=group_colors, hue_order = list(group_colors.keys()),
            fliersize=2, dodge = True, gap = 0.1)
plt.xlabel('Species')
plt.ylabel('IGS Size (bp)')
plt.title('Intergenic Space Size Distribution by Species and Group')
plt.legend(title = 'Group')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_total_boxplot.pdf'))
plt.close()


summary_df = pd.read_csv(os.path.join(files_dir, 'meanIGS.csv'))
summary_df['Group'] = summary_df['Group'].replace({
    'endosymb_only': 'Endosymbionts',
    'relatives_only': 'Free-Living Relatives',
    'Endosymbionts Only': 'Endosymbionts',
    'Free-Living Relatives Only': 'Free-Living Relatives'
})
summary_df = summary_df[summary_df['Group'].isin(['Endosymbionts', 'Free-Living Relatives'])]

threshold = summary_df['mean_mean_IGS'].quantile(0.90)
summary_df = summary_df[summary_df['mean_mean_IGS'] <= threshold]

pivot_df = summary_df.pivot_table(index='Species', columns='Group', values='mean_mean_IGS').reset_index()
outlier_threshold = pivot_df['Endosymbionts'].quantile(0.95)
pivot_outliers = pivot_df[pivot_df['Endosymbionts'] > outlier_threshold]

pivot_df = pivot_df[pivot_df['Endosymbionts'] <= outlier_threshold]

plt.figure(figsize=(12,12))
sns.boxplot(data = summary_df, x = 'Group', y = 'mean_mean_IGS', hue='Group', palette=group_colors, hue_order=list(group_colors.keys()), fliersize=0)
sns.stripplot(data = summary_df, x = 'Group', y = 'mean_mean_IGS', color='black', alpha=0.7)
plt.title('Intergenic Space (IGS) Size by Group', fontsize = 24, fontweight='bold')
plt.xticks(fontsize=18)
plt.xlabel('Group', fontsize = 20, fontweight='bold')
plt.yticks(fontsize=18)
plt.ylabel('IGS Size (bp)', fontsize = 20, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_mean_boxplot.pdf'))
plt.close()

plt.figure(figsize=(12,12))
sns.scatterplot(data = pivot_df, x = 'Endosymbionts', y = 'Free-Living Relatives')
x_max, x_min = pivot_df['Endosymbionts'].max(), pivot_df['Endosymbionts'].min()
y_max, y_min = pivot_df['Free-Living Relatives'].max(), pivot_df['Free-Living Relatives'].min()
common_min = min(x_min, y_min)
common_max = max(x_max, y_max)
padding = (common_max - common_min) * 0.1
common_max, common_min = common_max + padding, common_min - padding
plt.xlim(common_min, common_max)

plt.axline((0, 0), slope=1, color='red', linestyle='--', linewidth=1, label='y=x')

plt.xlabel('Endosymbionts', fontsize = 20, fontweight='bold')
plt.ylabel('Free-Living Relatives', fontsize = 20, fontweight='bold')
plt.margins(x=0.1)
plt.title('Intergenic Space Size (bp)\nEndosymbionts vs Free-Living Relatives', fontsize = 24, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_scatterplot.pdf'))
plt.close()

# Correlation of intergenic size and Delta GC%
genome_json_path = genome_gcsize_json_path('endosymb+relatives')
genome_data = load_or_compute(genome_json_path, genome_gcsize, 'endosymb+relatives')
delta_df = delta_matrix(genome_data, mat_type='gc_genome')
transposase = pd.read_csv(os.path.join(files_dir, 'median_transposases_species.csv'))

data_list = []
for species, matrix in delta_df.items():
    sp_name = species.replace('_endosymbiont','').replace('_', ' ')
    if sp_name not in summary_df['Species'].values:
        print(f'Skipping {species} as no IGS data is available.')
        continue
    print(f'Processing {species} for IGS vs Delta GC% correlation.')

    igs_size = summary_df[summary_df['Species'] == sp_name]['mean_mean_IGS'].values[0]
    
    ids = matrix.index.tolist()
    i,j = np.triu_indices_from(matrix, k=1)    
    dist_list = []
    for index1, index2 in zip(i, j):
        id1, id2 = ids[index1], ids[index2]
        if ('_genomic' in id1) == ('_genomic' in id2):
            continue
        dist_list.append(matrix.loc[id1, id2])
    if not dist_list:
        print(f'No valid distance pairs for {species}, skipping.')
        continue
        
    if species not in genome_data:
        print(f'No genome GC% or size data for {species}, skipping.')
        continue
    else:
        gc_values = [genome_data[species].get(genome_id, {}).get('gc_genome') for genome_id in genome_data[species] if not '_genomic' in genome_id and isinstance(genome_data[species], dict)]
        size_values = [genome_data[species].get(genome_id, {}).get('size') for genome_id in genome_data[species] if not '_genomic' in genome_id and isinstance(genome_data[species], dict)]
        gc_value = np.median(gc_values) if gc_values else np.nan
        size = np.median(size_values) if size_values else np.nan
    
    if sp_name not in transposase['Species'].values:
        print(f'No transposase data for {sp_name}, skipping.')
        continue
    else:
        transposases = transposase[transposase['Species'] == sp_name]['Median_Transposases'].values[0]

    delta_gc_value = np.median(dist_list)

    data_list.append({'species': sp_name, 'mean_IGS': igs_size, 'DeltaGC': delta_gc_value, 'genome_gc': gc_value, 'genome_size': size, 'transposases': transposases})

corr_df = pd.DataFrame(data_list)

def scatterplot(corr_df, x_col):
    labels = {'DeltaGC': 'Delta GC%', 'genome_gc': 'Genome GC%', 'genome_size': 'Genome Size (bp)', 'transposases': 'Median Number of Transposases'}

    plt.figure(figsize=(6,6))
    sns.scatterplot(data=corr_df, x=x_col, y='mean_IGS', alpha=0.6)
    plt.xlabel(labels[x_col])
    plt.ylabel('Mean IGS Size (bp)')
    plt.title(f'IGS Size and {labels[x_col]} in Endosymbionts')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'IGS_vs_{x_col}_endosymb.pdf'))
    plt.close()

for xcol in corr_df.columns:
    if xcol in ['DeltaGC', 'genome_gc', 'genome_size', 'transposases']:
        scatterplot(corr_df, xcol)
    else:
        continue






