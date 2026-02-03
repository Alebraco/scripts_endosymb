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
                        families = len(set(fam_list)) if fam_list else 0

                results.append({
                    'Group': group,
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
    plt.figure(figsize=(6,6))

    sns.stripplot(x='Group', y='Total_Transposases', hue='Group',
                  data=df_master, alpha=0.7, jitter=True, palette='Set2', 
                  order=['endosymb_only', 'relatives_only'])
    
    plt.title('Total Transposase Count by Group')
    plt.ylabel('Number of Transposases')
    plt.xlabel('Group')
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'transposase_group_count.pdf'))
    plt.close()

if __name__ == "__main__":
    df_master = processing_transposase()
    if df_master is not None:
        print(df_master.head())
        df_master.to_csv(os.path.join(files_dir, 'transposase_summary.csv'), index=False)

        parameters = [
            ('Total_Transposases', 'Abundance of Transposases per Genome', 'Total Transposases', 'transposase_abundance.pdf'),
            ('Unique_Families', 'Diversity of Transposases (Unique Families)', 'Unique IS Families', 'transposase_diversity.pdf')
        ]
        for metric, title, label, filename in parameters:
            transposase_plot(df_master, metric, title, label, filename)
        transposase_group_count(df_master)

