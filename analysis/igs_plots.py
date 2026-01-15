#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from utils import files_dir

group_colors = {
# 'Endosymbionts and Free-Living Relatives' : '#1f77b4',
'Endosymbionts Only' : '#ff7f0e',
'Free-Living Relatives Only' : '#2ca02c',
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
plt.savefig('IGS_total_boxplot.pdf')
plt.close()


#Scatterplot with mean IGS only
summary_df = pd.read_csv(os.path.join(files_dir, 'medianIGS.csv'))
summary_df = summary_df[summary_df['group'].isin(['Endosymbionts Only','Free-Living Relatives Only'])]
pivot_df = summary_df.pivot_table(index='species', columns='group', values='mean_median_IGS').reset_index()
outliers = pivot_df[pivot_df['Endosymbionts Only'] > pivot_df['Endosymbionts Only'].quantile(0.80)]

plt.figure(figsize=(8,8))
sns.boxplot(data = summary_df, x = 'group', y = 'mean_median_IGS', hue='group', palette=group_colors, hue_order = list(group_colors.keys()))
sns.stripplot(data = summary_df, x = 'group', y = 'mean_median_IGS', color='black', alpha=0.7)
plt.title('Intergenic Space (IGS) Size by Group')
plt.xlabel('Group')
plt.ylabel('IGS Size (bp)')
plt.tight_layout()
plt.savefig('IGS_mean_boxplot.pdf')
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
plt.title('Mean Intergenic Space Size (bp)\nEndosymbionts vs Free-Living Relatives')
plt.tight_layout()
plt.savefig('IGS_scatterplot.pdf')
plt.close()
