#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from utils import files_dir

plot_dir = os.path.join('plots', 'transposase')
os.makedirs(plot_dir, exist_ok=True)

coverage_threshold = 0.8
identity_threshold = 30.0
columns = ['qseqid', 
           'sseqid', 
           'pident', 
           'length', 
           'mismatch', 
           'gapopen', 
           'qstart', 
           'qend', 
           'sstart', 
           'send', 
           'evalue', 
           'bitscore', 
           'qlen', 
           'slen']
use_columns = ['sseqid', 'pident', 'sstart', 'send', 'slen']

def processing_transposase():
    results = []

    for group in ['endosymb_only', 'relatives_only']:
        group_path = os.path.join(group, 'transposase')

        for species in os.listdir(group_path):
            sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
            species_path = os.path.join(group_path, species)

            tsv_files = [file for file in os.listdir(species_path) if file.endswith('.tsv')]
            total_genomes = len(tsv_files)

            for file in tsv_files:
                file_path = os.path.join(species_path, file)

                complete, partial, total, families = 0,0,0,0
                fam_list = None

                if os.path.getsize(file_path) > 0:
                    df = pd.read_csv(file_path, sep='\t', names=columns, usecols=use_columns)

                    # Filter hits
                    df = df[df['pident'] >= identity_threshold]
                    df = df[df['sseqid'].str.contains('Transposase', case=False)]
                    df = df[~df['sseqid'].str.contains('Accessory', case=False)]

                    if not df.empty:

                        df['coverage'] = (df['send'] - df['sstart'] + 1) / df['slen']
                        complete = len(df[df['coverage'] >= coverage_threshold])
                        partial = len(df[df['coverage'] < coverage_threshold])
                        total = complete + partial

                        # Extract IS families
                        fam_list = df['sseqid'].str.split('_').str[0].str.split('//').str[1].tolist()
                        families = len(fam_list)

                results.append({
                    'Group': group,
                    'Species': sp_name,
                    'File': file.replace('.tsv',''),
                    'Total_Genomes': total_genomes,
                    'Complete_Transposases': complete,
                    'Partial_Transposases': partial,
                    'Total_Transposases': total,
                    'IS_Families': fam_list,
                    'Families_Count': families
                })

    if results:
        return pd.DataFrame(results)
    else:
        return None

def abundance_plot(df_master):
    palette = {
        'endosymb_only': 'lightcoral',
        'relatives_only': 'lightblue'
    }
    plt.figure(figsize=(12, 6))

    ax = sns.boxplot(x='Group', y='Total_Transposases', hue='Group', data=df_master, 
                        palette=palette, dodge=False, showfliers=False)

    sns.stripplot(x='Group', y='Total_Transposases', hue='Group',
                data=df_master, palette=['black','grey'], 
                alpha=0.7, jitter=True)

    plt.title("Number of Transposases per Genome by Group")
    plt.ylabel("Total Transposases per Genome")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'transposase_abundance.pdf'))
    plt.close()

# def heatmap_is_families(df_master):
#     # Expand copies into rows
#     df_exploded = df_master.explode('IS_Families')
    
#     # Count occurrences per species (IS Abundance)
#     matrix = pd.crosstab(index=[df_exploded['Group'], df_exploded['Species']], 
#                          columns=df_exploded['IS_Families'])
    
#     full_index = pd.MultiIndex.from_frame(df_master[['Group', 'Species']].drop_duplicates())
#     matrix = matrix.reindex(full_index, fill_value=0)
    
#     # Normalize by number of genomes to get average copies per genome
#     genome_counts = df_master.groupby(['Group', 'Species'])['Total_Genomes'].max()
#     matrix = matrix.div(genome_counts, axis=0)

#     matrix['Total_Abundance'] = matrix.sum(axis=1)

#     # Sort by group and total abundance    
#     matrix = matrix.sort_values(by=['Group', 'Total_Abundance'], ascending=[True, False])

#     plot_data = matrix.drop(['Total_Abundance'], axis=1)

#     # Select top 20 families
#     top_families = plot_data.sum().sort_values(ascending=False).head(20).index
#     plot_data = plot_data[top_families]

#     plot_data.index = [idx[1] for idx in plot_data.index]

#     plt.figure(figsize=(14, 10))
#     sns.heatmap(plot_data, cmap="Reds", linewidths=0.5, linecolor='whitesmoke',
#                 cbar_kws={'label': 'Average Copies per Genome'})
    
#     plt.title("Transposase Family Abundance: Endosymbionts vs Relatives")
#     plt.xlabel("IS Family (Top 20)")
#     plt.ylabel("Species")
#     plt.tight_layout()
#     plt.savefig(os.path.join(plot_dir, 'family_heatmap.pdf'))
#     plt.close()

import numpy as np
import matplotlib.lines as mlines

def heatmap_is_families(df_master):
    # 1. Expand and Crosstab
    df_exploded = df_master.explode('IS_Families')
    matrix = pd.crosstab(index=[df_exploded['Group'], df_exploded['Species']], 
                         columns=df_exploded['IS_Families'])
    
    # 2. Safety Net (Reindex)
    full_index = pd.MultiIndex.from_frame(df_master[['Group', 'Species']].drop_duplicates())
    matrix = matrix.reindex(full_index, fill_value=0)
    
    # 3. Normalize
    genome_counts = df_master.groupby(['Group', 'Species'])['Total_Genomes'].max()
    matrix = matrix.div(genome_counts, axis=0)

    # 4. Sorting
    matrix['Total_Abundance'] = matrix.sum(axis=1)
    matrix = matrix.sort_values(by=['Group', 'Total_Abundance'], ascending=[True, False])

    plot_data = matrix.drop(['Total_Abundance'], axis=1)

    # 5. Filter Top 40 (Increased from 20 to capture more rare families)
    top_families = plot_data.sum().sort_values(ascending=False).head(40).index
    plot_data = plot_data[top_families]
    
    # 6. LOG TRANSFORM (The Fix for the "Empty" look)
    # We use log1p (log(1+x)) to handle zeros gracefully
    plot_data_log = np.log1p(plot_data)

    # Clean Index
    plot_data_log.index = [idx[1] for idx in plot_data_log.index]
    n_endo = sum(matrix.index.get_level_values('Group') == 'endosymb_only')
    n_total = len(plot_data_log)

    # 7. Plot
    plt.figure(figsize=(16, 10)) # Wider for 40 families
    ax = sns.heatmap(plot_data_log, cmap="Reds", linewidths=0.5, linecolor='whitesmoke',
                     cbar_kws={'label': 'Log(Avg Copies + 1)'}) # Updated label
    
    # --- BRACKETS ---
    trans = ax.get_yaxis_transform()
    def draw_bracket(ax, y_start, y_end, label, x_pos=1.01):
        ax.add_line(mlines.Line2D([x_pos, x_pos], [y_start, y_end], transform=trans, color='black', clip_on=False))
        ax.add_line(mlines.Line2D([x_pos, x_pos-0.01], [y_start, y_start], transform=trans, color='black', clip_on=False))
        ax.add_line(mlines.Line2D([x_pos, x_pos-0.01], [y_end, y_end], transform=trans, color='black', clip_on=False))
        ax.text(x_pos + 0.02, (y_start + y_end) / 2, label, transform=trans, va='center', ha='left', fontsize=12)

    draw_bracket(ax, 0, n_endo, "Endosymbionts")
    draw_bracket(ax, n_endo, n_total, "Relatives")

    plt.title("Transposase Family Abundance (Log Scale)")
    plt.xlabel("IS Family (Top 40)")
    plt.ylabel("Species")
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    plt.savefig(os.path.join(files_dir, 'family_heatmap_log.pdf'))
    plt.close()

if __name__ == "__main__":
    df_master = processing_transposase()
    if df_master is not None:
        print(df_master.head())
        df_master.to_csv(os.path.join(files_dir, 'transposase_summary.csv'), index=False)
        abundance_plot(df_master)
        heatmap_is_families(df_master)

