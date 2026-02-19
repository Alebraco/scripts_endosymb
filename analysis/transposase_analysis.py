#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from utils import files_dir, group_names

plot_dir = os.path.join('plots', 'transposase')
os.makedirs(plot_dir, exist_ok=True)

coverage_threshold = 0.8
identity_threshold = 25.0
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

def processing_transposase(path, group_label = None, auto_classify = False):
    """
    base_path: The directory containing the 'transposase' folder.
               e.g. 'path/to/group' or 'path/to/general_folder'
    group_label: Optional string to label the 'Group' column. 
                 If None, defaults to 'Ungrouped'.
    """
    results = []

    target_dir = os.path.join(path, 'transposase')

    if not os.path.exists(target_dir):
        print(f"Directory not found: {target_dir}")
        return None
    
    default_group = group_label if group_label else 'Ungrouped'

    for species in os.listdir(target_dir):
        sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
        species_path = os.path.join(target_dir, species)

        if not os.path.isdir(species_path):
            continue

        tsv_files = [file for file in os.listdir(species_path) if file.endswith('.tsv')]
        total_genomes = len(tsv_files)

        for file in tsv_files:
            file_path = os.path.join(species_path, file)

            complete, partial, total, families = 0,0,0,0
            fam_list = None

            if auto_classify:
                if '_genomic' in file:
                    current_group = 'relatives_only'
                else:
                    current_group = 'endosymb_only'
            else:
                current_group = default_group

            if os.path.getsize(file_path) > 0:
                df = pd.read_csv(file_path, sep='\t', names=columns, usecols=use_columns)

                df = df[df['pident'] >= identity_threshold]
                df = df[df['sseqid'].str.contains('Transposase', case=False)]
                df = df[~df['sseqid'].str.contains('Accessory', case=False)]

                if not df.empty:
                    df['coverage'] = (df['send'] - df['sstart'] + 1) / df['slen']
                    complete = len(df[df['coverage'] >= coverage_threshold])
                    partial = len(df[df['coverage'] < coverage_threshold])
                    total = complete + partial

                    fam_list = df['sseqid'].str.split('_').str[0].str.split('//').str[1].tolist()
                    families = len(set(fam_list)) if fam_list else 0

            results.append({
                'Group': current_group,
                'Species': sp_name,
                'File': file.replace('.tsv',''),
                'Total_Genomes': total_genomes,
                'Complete_Transposases': complete,
                'Partial_Transposases': partial,
                'Total_Transposases': total,
                'IS_Families': fam_list,
                'Unique_Families': families
            })

    if results:
        return pd.DataFrame(results)
    else:
        return None

def transposase_plot(df_master, metric, title, label, filename):
    endo_order = sorted(df_master[df_master['Group'] == 'endosymb_only']['Species'].unique())
    rel_order = sorted(df_master[df_master['Group'] == 'relatives_only']['Species'].unique())

    colors = sns.color_palette('Set2', n_colors=2)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 10), sharex=True)

    sns.stripplot(y='Species', x=metric,
                  data=df_master[df_master['Group'] == 'endosymb_only'],
                  alpha=0.7, jitter=True, ax=axes[0], color=colors[0],
                  order=endo_order)
    
    axes[0].set_title('Endosymbionts')
    axes[0].tick_params(axis='x', rotation=90) 
    axes[0].set_ylabel('')
    axes[0].set_xlabel(label)
    axes[0].grid(axis='x', linestyle='--', alpha=0.5)

    sns.stripplot(y='Species', x=metric,
                  data=df_master[df_master['Group'] == 'relatives_only'],
                  alpha=0.7, jitter=True, ax=axes[1], color=colors[1],
                  order=rel_order)
   
    axes[1].set_title('Relatives')
    axes[1].grid(axis='x', linestyle='--', alpha=0.5)
    axes[1].set_xlabel(label)
    axes[1].set_ylabel('') 

    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, filename))
    plt.close()

def transposase_group_count(df_master):
    plt.figure(figsize=(12,12))

    metric = 'Transposases_per_Gene' if 'Transposases_per_Gene' in df_master.columns else 'Total_Transposases'
    

    sns.boxplot(x='Group', y=metric, hue='Group', data=df_master,
                    palette='Set2', showfliers=False, width=0.5, boxprops={'alpha': 0.4},
                    order=['relatives_only', 'endosymb_only'], hue_order=['relatives_only', 'endosymb_only'])

    sns.stripplot(x='Group', y=metric, hue='Group',
                  data=df_master, alpha=0.7, jitter=True, palette='Set2', 
                  order=['relatives_only', 'endosymb_only'], hue_order=['relatives_only', 'endosymb_only'])
    
    if metric == 'Transposases_per_Gene':
        plt.title('Normalized Transposase Abundance by Group', fontsize=24, fontweight='bold')
        plt.ylabel('Transposases per Gene', fontsize=22, fontweight='bold')
    else:
        plt.title('Total Transposase Count by Group', fontsize=24, fontweight='bold')
        plt.ylabel('Number of Transposases', fontsize=22, fontweight='bold')
    plt.xlabel('Group', fontsize=22, fontweight='bold')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.grid(axis='y', linestyle='--', alpha=0.5)

    ax = plt.gca()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Endosymbionts', 'Relatives'])

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'transposase_group_count.pdf'))
    plt.close()

def transposase_completeness(df_master, group_name):

    df = df_master[df_master['Group'] == group_name]
    df_agg = df.groupby('Species')[['Complete_Transposases', 'Partial_Transposases']].mean()
    df_agg['Total'] = df_agg['Complete_Transposases'] + df_agg['Partial_Transposases']
    df_sorted = df_agg.sort_values('Total', ascending=True)

    ax = df_sorted[['Complete_Transposases', 'Partial_Transposases']].rename(
        columns={
            'Complete_Transposases': 'Complete', 
            'Partial_Transposases': 'Partial'
        }).plot(
        kind='barh', 
        stacked=True, 
        figsize=(8, 14),
        color=['black', 'grey'],
        edgecolor='white'
    )
    
    plt.title(f'Transposase Integrity: {group_names[group_name].replace(" Only", "")}', fontsize=14)
    plt.xlabel('Number of Transposases')
    plt.ylabel('')
    plt.legend(title='Integrity', loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    
    plt.savefig(os.path.join(plot_dir, f'transposase_completeness_{group_name}.pdf'))
    plt.close()

def transposase_completeness_perc(df_master):
    df = df_master.groupby(['Group'])[['Complete_Transposases', 'Partial_Transposases']].sum().reset_index()
    df['Total'] = df['Complete_Transposases'] + df['Partial_Transposases']
    df['Complete_Percentage'] = (df['Complete_Transposases'] / df['Total']) * 100
    df['Partial_Percentage'] = 100 - df['Complete_Percentage']
    df_melted = df.melt(id_vars='Group', value_vars=['Complete_Percentage', 'Partial_Percentage'],
                        var_name='Integrity', value_name='Percentage')
    df_melted['Integrity'] = df_melted['Integrity'].replace({'Complete_Percentage': 'Complete', 'Partial_Percentage': 'Partial'})

    plt.figure(figsize=(6,6))
    sns.barplot(x='Group', y='Percentage', hue='Integrity', data=df_melted,
                palette=['black', 'grey'], edgecolor='white')
    
    plt.title('Transposase Completeness Percentage by Group')
    plt.ylabel('Percentage (%)')
    plt.xlabel('Group')
    plt.ylim(0, 100)
    plt.legend(title='Integrity')

    ax = plt.gca()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Endosymbionts', 'Relatives'])

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'transposase_completeness_perc.pdf'))
    plt.close()

if __name__ == "__main__":
    dfs = []

    for group in ['endosymb_only', 'relatives_only']:
        df = processing_transposase(group, group_label=group)
        if df is not None:
            dfs.append(df)
    if dfs:
        df_master = pd.concat(dfs, ignore_index=True)

        gene_counts_path = os.path.join(files_dir, 'gene_counts.csv')
        if os.path.exists(gene_counts_path):
            gene_counts_df = pd.read_csv(gene_counts_path)
            df_master = pd.merge(df_master, gene_counts_df, on=['Group', 'Species', 'File'], how='left')

            df_master['Transposases_per_Gene'] = df_master['Total_Transposases'] / df_master['Gene_Count']
            print('Calculated Transposases per Gene and added to DataFrame.')
        else:
            print(f'Warning: gene_counts.csv not found at {gene_counts_path}.')

        df_master.to_csv(os.path.join(files_dir, 'transposase_summary.csv'), index=False)

        # Compute median number of transposases per species
        median_transposases = df_master.groupby('Species')['Total_Transposases'].median().reset_index().rename(columns={'Total_Transposases': 'Median_Transposases'})
        median_transposases.to_csv(os.path.join(files_dir, 'median_transposases_species.csv'), index=False)

        parameters = [
            ('Total_Transposases', 'Abundance of Transposases per Genome', 'Total Transposases', 'transposase_abundance.pdf'),
            ('Unique_Families', 'Diversity of Transposases (Unique Families)', 'Unique IS Families', 'transposase_diversity.pdf')
        ]

        for metric, title, label, filename in parameters:
            transposase_plot(df_master, metric, title, label, filename)

        transposase_group_count(df_master)

        for group in ['endosymb_only', 'relatives_only']:
            transposase_completeness(df_master, group)

        transposase_completeness_perc(df_master)
    else:
        print("No transposase data found to process.")

