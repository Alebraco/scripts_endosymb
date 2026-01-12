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


csv_path = 'fourfold_gc_content.csv'

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

    for group, group_name in groups.items():
        group_path = os.path.join(group, 'dna_concatenates')
        for species in os.listdir(group_path):
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
            print(f'Results for {sp_name} in {group} ---')

            for aa, nuc_counts in counts.items():
                gc_count = nuc_counts['G'] + nuc_counts['C']
                sum_count = sum(nuc_counts.values())

                if sum_count > 0:
                    gc4_percent = round((gc_count / sum_count) * 100, 2)
                    gc_values.append(gc4_percent)
                    print(f'{aa} GC4: {gc4_percent:.2f}% (n={sum_count})')

                    total_gc4 += nuc_counts['G'] + nuc_counts['C']
                    total_count += sum_count
            if total_count > 0:
                overall_gc4 = round((total_gc4 / total_count) * 100, 2)
                stdev_gc4 = round(statistics.stdev(gc_values), 2) if len(gc_values) > 1 else 0.0
                print(f'Overall GC4: {overall_gc4:.2f}%')
                print(f'Variance (SD) between families: {stdev_gc4:.2f}')

                all_data.append({'group': group_name, 'species': sp_name, 'overall_gc4': overall_gc4, 'stdev_gc4': stdev_gc4})

    df = pd.DataFrame(all_data)
    df.to_csv(csv_path, index=False)
    print(f'Saved fourfold GC content data to {csv_path}')

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

def plot_results(df, p_val):
    plt.figure(figsize=(8, 6))
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
    plt.savefig('gc4_boxplot.pdf')
    plt.close()

    linear = sns.lmplot(x='overall_gc4', y='stdev_gc4', hue='group', data=df, palette='Set2', 
               aspect=1.2, scatter_kws={'alpha':0.6})
    linear.figure.suptitle('Overall GC4 vs GC4 Variation', y=1.02)
    linear.set_axis_labels('Overall GC4 (%)', 'Standard Deviation of GC4 (%)')

    plt.tight_layout()
    plt.savefig('gc4_correlation.pdf')
    plt.close()

if __name__ == '__main__':
    if not os.path.exists(csv_path):
        fourfold_gc_content()

    t_stat, p_val = stat_analysis()
    df = pd.read_csv(csv_path)
    plot_results(df, p_val)