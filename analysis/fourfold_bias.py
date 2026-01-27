#!/usr/bin/env python3
import os
import statistics
from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from utils import files_dir
from collections import Counter

csv_path = os.path.join(files_dir,'fourfold_gc_content.csv')
detailed_csv_path = os.path.join(files_dir,'fourfold_gc_detailed.csv')
plot_dir = os.path.join('plots', 'fourfold_bias')
os.makedirs(plot_dir, exist_ok=True)

def fourfold_gc_content():
    '''
    Compares global and fourfold amino acid-specific GC4 across species between two groups.
    '''
    fourfold_families = {
        'GC': 'Ala', 'GG': 'Gly', 'CC': 'Pro', 'AC': 'Thr', 
        'GT': 'Val', 'CG': 'Arg4', 'CT': 'Leu4', 'TC': 'Ser4'
    }
    groups = {'relatives_only': 'Relatives', 'endosymb_only': 'Endosymbionts'}
    all_data = []
    detailed_data = []

    for group, group_name in groups.items():
        group_path = os.path.join(group, 'dna_concatenates')
        for species in os.listdir(group_path):
            temp_data = []
            sp_name = species.replace('concatenate_', '').replace('_', ' ').replace('.fasta', '')
            species_path = os.path.join(group_path, species)
            counts = {aa: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for aa in fourfold_families.values()}
            for record in SeqIO.parse(species_path, 'fasta'):
                seq = str(record.seq).upper()
                for i in range(0, len(seq) - 2, 3):
                    codon = seq[i:i+3]
                    prefix = codon[0:2]
                    third_base = codon[2]

                    if prefix in fourfold_families.keys() and third_base in 'ACGT':
                        aa = fourfold_families[prefix]
                        counts[aa][third_base] += 1
            total_gc4 = 0
            total_count = 0
            overall_gc4 = 0
            stdev_gc4 = 0.0
            gc_values = []

            for aa, nuc_counts in counts.items():
                gc_count = nuc_counts['G'] + nuc_counts['C']
                sum_count = sum(nuc_counts.values())

                if sum_count > 0:
                    gc4_percent = round((gc_count / sum_count) * 100, 2)
                    gc_values.append(gc4_percent)

                    temp_data.append({
                        'group': group_name, 
                        'species': sp_name,
                        'amino_acid': aa,
                        'gc4_percent': gc4_percent,
                        'count': sum_count
                    })

                    total_gc4 += gc_count
                    total_count += sum_count
                    
            if total_count > 0:
                overall_gc4 = round((total_gc4 / total_count) * 100, 2)

                for item in temp_data:
                    item['overall_gc4'] = overall_gc4
                    item['bias'] = abs(round(item['gc4_percent'] - overall_gc4, 2))
                    detailed_data.append(item)

                stdev_gc4 = round(statistics.stdev(gc_values), 2) if len(gc_values) > 1 else 0.0

                all_data.append({
                    'group': group_name, 
                    'species': sp_name, 
                    'overall_gc4': overall_gc4, 
                    'stdev_gc4': stdev_gc4
                })

    df = pd.DataFrame(all_data)
    df.to_csv(csv_path, index=False)

    detailed_df = pd.DataFrame(detailed_data)
    detailed_df.to_csv(detailed_csv_path, index=False)

    print(f'Saved fourfold GC content data to {csv_path}')
    print(f'Saved detailed fourfold GC content data to {detailed_csv_path}')

def stat_analysis():
    '''
    Performs statistical analysis on the fourfold GC content data.
    '''
    df = pd.read_csv(csv_path)
    relatives = df[df['group'] == 'Relatives']['stdev_gc4']
    endosymb = df[df['group'] == 'Endosymbionts']['stdev_gc4']

    t_stat, p_val = stats.ttest_ind(relatives, endosymb, equal_var=False)

    print(f'T-test results: t-statistic = {round(t_stat, 4)}, p-value = {p_val:.4e}')
    if p_val < 0.05:
        print('Significant difference in GC4 variation between groups')
    else:
        print('No significant difference in GC4 variation between groups')
    return t_stat, p_val

def plot_results(df, detailed_df, p_val):
    
    #Boxplot of GC4 variation
    plt.figure(figsize=(6, 6))
    sns.boxplot(x='group', y='stdev_gc4', data=df, palette='Set2', hue = 'group')
    sns.stripplot(x='group', y='stdev_gc4', data=df, color='black', alpha=0.5)
    
    if p_val < 0.001:
        sig_label = '***'
    elif p_val < 0.01:
        sig_label = '**'
    elif p_val < 0.05:
        sig_label = '*'
    else:
        sig_label = 'ns'
    x1, x2 = 0, 1
    y_max = df['stdev_gc4'].max()
    y, h = y_max + (y_max * 0.05), y_max * 0.03

    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='grey')
    plt.text((x1+x2)*.5, y+h + (y_max * 0.02), f"{sig_label}\n(p = {p_val:.2e})", 
             ha='center', va='bottom', color='grey', fontsize=10)
    
    plt.ylim(top=y_max * 1.20)
    
    plt.title('Codon-specific GC4 Variation Between Groups')
    plt.xlabel('Group')
    plt.ylabel('Standard Deviation of GC4 (%)')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gc4_boxplot.pdf'))
    plt.close()

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='overall_gc4', 
        y='stdev_gc4', 
        hue='group', 
        hue_order=['Relatives','Endosymbionts'],
        data=df, 
        palette='Set2', 
        alpha=0.8
        )
    plt.title('Overall GC4 vs GC4 Variation')
    plt.xlabel('Overall GC4 (%)')
    plt.ylabel('Standard Deviation of GC4 (%)')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gc4_scatter.pdf'))
    plt.close()

    #Amino acid specific boxplots
    plt.figure(figsize=(10, 6))
    detailed_df = pd.read_csv(detailed_csv_path)
    sns.boxplot(
        x='amino_acid', 
        y='gc4_percent', 
        hue='group',
        data=detailed_df, 
        hue_order=['Relatives','Endosymbionts'],
        palette='Set2',
        fliersize=0
        )
    sns.stripplot(
        x='amino_acid', 
        y='gc4_percent', 
        hue='group', 
        dodge=True,
        hue_order=['Relatives','Endosymbionts'],
        data=detailed_df, 
        palette=['black','grey'],
        alpha=0.5,
        legend=False
        )
    plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('GC4 Percentage by Amino Acid and Group')
    plt.xlabel('Amino Acid')
    plt.ylabel('GC4 Percentage (%)')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gc4_amino_acid_boxplot.pdf'))
    plt.close()

def amino_acid_bias(detailed_df):
    '''
    Determines the most/least biased amino acid for each species and group.
    '''
    idx_max = detailed_df.groupby(['group','species'])['bias'].idxmax()
    max_bias = detailed_df.loc[idx_max]

    idx_min = detailed_df.groupby(['group','species'])['bias'].idxmin()
    min_bias = detailed_df.loc[idx_min]

    print('Most Biased Amino Acids by Group:')
    print(max_bias.groupby('group')['amino_acid'].value_counts())

    print('\nLeast Biased Amino Acids by Group:')
    print(min_bias.groupby('group')['amino_acid'].value_counts())

    fig, axes = plt.subplots(1,2, figsize=(14,6), sharey=True)

    sns.countplot(
        x='amino_acid', 
        hue='group', 
        data=max_bias, 
        palette='Set2',
        hue_order=['Relatives','Endosymbionts'],
        order=detailed_df['amino_acid'].unique(),
        ax=axes[0],
        legend=False
        )
    for container in axes[0].containers:
        axes[0].bar_label(container, padding=2) 

    axes[0].set_title('Most Biased Amino Acids by Group')
    axes[0].set_xlabel('Amino Acid')
    axes[0].set_ylabel('Number of Species')

    sns.countplot(
        x='amino_acid', 
        hue='group', 
        data=min_bias, 
        palette='Set2',
        hue_order=['Relatives','Endosymbionts'],
        order=detailed_df['amino_acid'].unique(),
        ax=axes[1]
        )
    axes[1].legend(title='Group', loc='upper right')
    for container in axes[1].containers:
        axes[1].bar_label(container, padding=2) 

    axes[1].set_title('Least Biased Amino Acids by Group')
    axes[1].set_xlabel('Amino Acid')
    axes[1].set_ylabel('Number of Species')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gc4_amino_acid_bias.pdf'))
    plt.close()

def counts_matrix(detailed_df):
    '''
    Creates a counts matrix of amino acid biases per species and group.
    '''
    matrix = detailed_df.pivot_table(
        index=['group', 'species'], 
        columns='amino_acid', 
        values='count', 
        fill_value=0
    )
    matrix['Total'] = matrix.sum(axis=1)
    matrix = matrix.astype(int)

    output_csv = os.path.join(files_dir, 'fourfold_counts_matrix.csv')
    matrix.to_csv(output_csv)
    print(f'Saved fourfold counts matrix to {output_csv}')

if __name__ == '__main__':
    if not os.path.exists(csv_path) or not os.path.exists(detailed_csv_path):
        print('Calculating fourfold GC content...')
        fourfold_gc_content()
    else:
        print('Fourfold GC content data already exists. Skipping calculation.')

    df = pd.read_csv(csv_path)
    detailed_df = pd.read_csv(detailed_csv_path)
    counts_matrix(detailed_df)

    t_stat, p_val = stat_analysis()
    plot_results(df, detailed_df, p_val)
    amino_acid_bias(detailed_df)