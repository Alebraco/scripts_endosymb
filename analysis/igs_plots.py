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
from utils import files_dir, genome_gcsize_json_path, load_or_compute

plot_dir = os.path.join('plots', 'igs_lengths')
os.makedirs(plot_dir, exist_ok=True)
group_colors = {
# 'Endosymbionts and Free-Living Relatives' : '#1f77b4',
'Endosymbionts Only' : '#FC8D62FF',
'Free-Living Relatives Only' : '#66C2A5FF',
}

#Box Plots with all data points
df_boxplot = pd.read_csv(os.path.join(files_dir, 'all_IGS_data.csv'))
df_boxplot = df_boxplot[df_boxplot['group'].isin(['Endosymbionts Only','Free-Living Relatives Only'])]
df_boxplot_median = df_boxplot.groupby(['group', 'species', 'file'])['IGS_Size'].median().reset_index()

plt.figure(figsize=(12,12))
sns.boxplot(data = df_boxplot_median, x = 'IGS_Size', y = 'species', hue = 'group',
            palette=group_colors, hue_order = list(group_colors.keys()),
            fliersize=2, dodge = True, gap = 0.1)
plt.xlabel('Species')
plt.ylabel('IGS Size (bp)')
plt.title('Intergenic Space Size Distribution by Species and Group')
plt.legend(title = 'Group')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_total_boxplot.pdf'))
plt.close()


#Scatterplot with mean IGS only
summary_df = pd.read_csv(os.path.join(files_dir, 'meanIGS.csv'))
summary_df = summary_df[summary_df['group'].isin(['Endosymbionts Only','Free-Living Relatives Only'])]
pivot_df = summary_df.pivot_table(index='species', columns='group', values='mean_mean_IGS').reset_index()
outliers = pivot_df[pivot_df['Endosymbionts Only'] > pivot_df['Endosymbionts Only'].quantile(0.90)]

plt.figure(figsize=(8,8))
sns.boxplot(data = summary_df, x = 'group', y = 'mean_mean_IGS', hue='group', palette=group_colors, hue_order=list(group_colors.keys()), fliersize=0)
sns.stripplot(data = summary_df, x = 'group', y = 'mean_mean_IGS', color='black', alpha=0.7)
plt.title('Intergenic Space (IGS) Size by Group')
plt.xlabel('Group')
plt.ylabel('IGS Size (bp)')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_mean_boxplot.pdf'))
plt.close()

plt.figure(figsize=(8,8))
sns.scatterplot(data = pivot_df, x = 'Endosymbionts Only', y = 'Free-Living Relatives Only')
x_max, x_min = pivot_df['Endosymbionts Only'].max(), pivot_df['Endosymbionts Only'].min()
y_max, y_min = pivot_df['Free-Living Relatives Only'].max(), pivot_df['Free-Living Relatives Only'].min()
common_min = min(x_min, y_min)
common_max = max(x_max, y_max)
padding = (common_max - common_min) * 0.1
common_max, common_min = common_max + padding, common_min - padding
plt.xlim(common_min, common_max)

plt.axline((0, 0), slope=1, color='red', linestyle='--', linewidth=1, label='y=x')

# Label outliers
for idx, row in outliers.iterrows():
    plt.annotate(row['species'].replace('endosymbiont', '').replace('_', ' '),
                xy=(row['Endosymbionts Only'], row['Free-Living Relatives Only']),
                xytext=(5, 5), textcoords='offset points',
                fontsize=8, alpha=0.7)

plt.xlabel('Endosymbionts Only')
plt.ylabel('Free-Living Relatives Only')
plt.margins(x=0.1)
plt.title('Intergenic Space Size (bp)\nEndosymbionts vs Free-Living Relatives')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'IGS_scatterplot.pdf'))
plt.close()

# Correlation of intergenic size and Delta GC%
genome_json = genome_gcsize_json_path('endosymb_only')
delta_df = delta_matrix(load_or_compute(genome_json, genome_gcsize, 'endosymb_only'), mat_type='gc_genome', save_to_file=True)

print(f"Total species in delta_df: {len(delta_df)}")
print(f"Total species in summary_df: {len(summary_df['species'].unique())}")
print(f"Species in delta_df: {list(delta_df.keys())}")
print(f"Species in summary_df: {summary_df['species'].unique()}")

size, delta_gc, species_list = [], [], []
for species, matrix in delta_df.items():
    if species not in summary_df['species'].values:
        print(f'Skipping {species} as no IGS data is available.')
        continue

    igs_size = summary_df[summary_df['species'] == species]['mean_mean_IGS'].values[0]

    i,j = np.triu_indices_from(matrix, k=1)
    delta_gc_value = np.median(matrix.values[i,j])

    size.append(igs_size)
    delta_gc.append(delta_gc_value)
    species_list.append(species)

corr_df = pd.DataFrame({'species': species_list, 'mean_IGS': size, 'delta_gc': delta_gc})

plt.figure(figsize=(8,8))
sns.scatterplot(data=corr_df, x='delta_gc', y='mean_IGS', alpha=0.6)
plt.xlabel('Delta GC%')
plt.ylabel('Mean IGS Size (bp)')
plt.title(f'IGS Size and Delta GC% in Endosymbionts')
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, f'IGS_vs_DeltaGC_endosymb.pdf'))
plt.close()


