#!/usr/bin/env python3

import os
from adjustText import adjust_text
from itertools import combinations
from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from distance_matrix import distance_matrix
from utils import files_dir, load_or_compute_pickle

output_csv = 'spectrum_rates.csv'
data_dir = 'files'
csv_path = os.path.join(data_dir, output_csv)
plot_dir = os.path.join('plots', 'spectrum_plots')
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)

group_colors = {
'Endosymbiont-Relative Pairs' : '#1f77b4',
'Free-Living Control Pairs' : '#2ca02c',
'Endosymbiont Control Pairs' : '#ff7f0e'
}
mut_map = {
    'rGC->AT': 'rGC→AT',
    'rAT->GC': 'rAT→GC'
}
group_names = {
    'endosymb+relatives': 'Endosymbiont-Relative Pairs',
    'relatives_only': 'Free-Living Control Pairs',
    'endosymb_only': 'Endosymbiont Control Pairs'
}
mutation_types = list(mut_map.values())

def load_data(csv_file):
    '''
    Loads CSV, checks column names and renames if needed.
    Returns a DataFrame. 
    '''
    df = pd.read_csv(csv_file)
    df = df.rename(columns = mut_map)
    return df

def polymorphism_total(alignment):
    '''
    Calculation of polymorphic sites
    '''
    # Filter empty records
    seqs = [list(rec.seq) for rec in alignment if len(rec.seq) > 0]
    if len(seqs) < 2:
        return 0
    try:
        array = np.array(seqs)
    except ValueError:
        print('Sequences in the alignment have different lengths.')
        return 0
    
    # Determine polymorphisms in each column (site)
    A_check = np.any(array == 'A', axis = 0)
    C_check = np.any(array == 'C', axis = 0)
    G_check = np.any(array == 'G', axis = 0)
    T_check = np.any(array == 'T', axis = 0)

    # Sum different nucleotides for each column
    score = A_check.astype(int) + C_check.astype(int)  + G_check.astype(int) + T_check.astype(int)

    # If sum is 1, no polymorphism. Count how many columns had polymorphisms
    poly_count = np.sum(score > 1)
    return poly_count

def calculate_rates(group_name, species, ancestor, descendant):
    '''
        Calculates r(GC->AT) and r(AT->GC) for a pair of 
        ancestor-descendant sequences using the 
        Hershberg and Petrov (2010) formula.
    
    Args:
        group (str): Group name.
        species (str): Species name.
        ancestor (SeqIO.SeqRecord)
        descendant (SeqIO.SeqRecord)
    
    Returns:
        dict: species, ancestor ID, descendant ID, and rates.
    '''

    GC_set, AT_set, gap_set = {'G','C'}, {'A','T'}, {'N', '-'}

    GC_sites, AT_sites = 0,0
    GC_AT, AT_GC = 0,0

    for ancestor_base, descendant_base in zip(ancestor.seq, descendant.seq):
        ancestor_base = ancestor_base.upper()
        descendant_base = descendant_base.upper()

        if ancestor_base in gap_set or descendant_base in gap_set:
            continue
        
        if ancestor_base in GC_set:
            GC_sites += 1
        elif ancestor_base in AT_set:
            AT_sites += 1

        if ancestor_base != descendant_base:
            if ancestor_base in GC_set and descendant_base in AT_set:
                GC_AT += 1
            elif ancestor_base in AT_set and descendant_base in GC_set:
                AT_GC += 1

    if GC_sites != 0:
        r_GC_to_AT = GC_AT / GC_sites
    else:
        r_GC_to_AT = np.nan

    if AT_sites != 0:
        r_AT_to_GC = AT_GC / AT_sites
    else:
        r_AT_to_GC = np.nan
    
    return {
            'Group': group_name,
            'Species': species,
            'Ancestor': ancestor.id, 
            'Descendant': descendant.id,
            'rGC->AT': r_GC_to_AT,
            'rAT->GC': r_AT_to_GC,
            }

def plot_distributions(df, plot_type = 'kde'):
    '''
    Plots histograms side-by-side and summary boxplot 
    '''
    median_df = df.groupby(['Group','Species'])[['rGC→AT', 'rAT→GC']].median().reset_index()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    plt.suptitle('Comparison of Mutational Spectra by Group', fontsize=16, y=0.98)

    if plot_type != 'kde':
        # Histograms
        for ax, mut in zip(axes, mutation_types):
            sns.histplot(data = median_df,
                        x = mut, 
                        hue = 'Group', 
                        palette = group_colors,
                        hue_order = list(group_colors.keys()),
                        element = 'step',
                        fill = True,
                        alpha = 0.5,
                        ax=ax,
                        bins = 50,
                        stat = 'probability',
                        common_norm = False,
                        legend = (ax == axes[1])
            )
            ax.set_title(f'{mut} Rate Distributions')
            ax.set_xlabel(f'Rate ({mut})')        
            ax.set_xlim(0, 1)
        axes[0].set_ylabel('Probability') 
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'spectrum_hist.pdf'))
        plt.close()
        print('Saved spectrum_hist.pdf')
    else:
        # KDE plots
        for ax, mut in zip(axes, mutation_types):
            sns.kdeplot(data = median_df,
                        x = mut, 
                        hue = 'Group', 
                        palette = group_colors,
                        hue_order = list(group_colors.keys()),
                        clip=(0, 1),
                        fill = True,
                        alpha = 0.5,
                        ax=ax,
                        common_norm = False,
                        legend = (ax == axes[1])
            )
            ax.set_title(f'{mut} Rate Distributions')
            ax.set_xlabel(f'Rate ({mut})')        
            ax.set_xlim(0, 1)
        axes[0].set_ylabel('Density') 
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'spectrum_kde.pdf'))
        plt.close()
        print('Saved spectrum_kde.pdf')

    # Boxplot
    df_melted = median_df.melt(
        id_vars = ['Group', 'Species'],
        value_vars = mutation_types,
        var_name = 'Mutation Type',
        value_name = 'Rate'
    ).dropna(subset=['Rate'])

    plt.figure(figsize=(8, 6))
    sns.boxplot(data = df_melted,
                x = 'Group', 
                y = 'Rate',
                hue = 'Mutation Type', 
                palette = 'Paired',
                fliersize = 1,
                width = 0.5
                )
    plt.title('Summary of Mutation Rates by Group', fontsize=16)
    plt.xlabel('Group')
    plt.ylabel('Rate')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'spectrum_boxplot.pdf'))
    plt.close()
    print('Saved spectrum_boxplot.pdf')

def plot_species_grid(df):
    '''
    Generates a grid of histogram plots for each species in the CSV.
    Only plots if species are shared by both groups.
    Produces two files: one for GC->AT rates and one for AT->GC rates.
    '''

    # Include only species present in both groups
    group1 = set(df[df['Group'] == 'Endosymbiont-Relative Pairs']['Species'])
    group2 = set(df[df['Group'] == 'Free-Living Control Pairs']['Species'])
    common_species = group1.intersection(group2)
    df = df[df['Species'].isin(common_species)].copy()
    print(f'Plotting {len(common_species)} species shared by both groups.')

    sns.set_style('whitegrid')
    
    for mut in mutation_types:
        safe_name = mut.replace('→', '_').replace('r', '')
        output = f'grid_{safe_name}.pdf'
    
        grid = sns.FacetGrid(
            df, 
            col='Species', 
            hue='Group',
            palette = group_colors,
            hue_order = list(group_colors.keys()),
            col_wrap=7,
            sharex=True,
            sharey=False,
            height=2.5,
            aspect=1.2
        )
        
        grid.map(sns.histplot, mut, stat='density',
                common_norm=False, fill=True, alpha=0.5, 
                bins='auto', element = 'step')
        grid.add_legend(title='Group')
        grid.set_titles('{col_name}')
        grid.set_axis_labels(f'Rate ({mut})', 'Density')
        grid.set(xlim=(0, 1))
        grid.figure.suptitle(f'{mut} Mutation Rate Distributions by Species', fontsize=18, y =1.02)
        
        grid.savefig(os.path.join(plot_dir, output))
        plt.close()
        print(f'Saved {output}')

def run_diagnostic():
    print('Running diagnostic...')
    total_sites = []
    for group in group_names.keys():
        print(f'Processing {group}')
        input_dir = os.path.join(group, 'third_sites')
        sp_list = os.listdir(input_dir)
        for species in sp_list:
            sp_name = species.replace('concatenate_', '').replace('.fasta', '')
            sp_path = os.path.join(input_dir, species)
            aln = list(SeqIO.parse(sp_path, 'fasta'))
            # length = next((len(record.seq) for record in aln if len(record.seq) > 0), 0)
            clean_seqs = [str(rec.seq) for rec in aln if len(rec.seq) > 0]
            lengths = set([len(record.seq) for record in aln if len(record.seq) > 0])
            
            poly_count = polymorphism_total(aln)
            dseqnum = len(set(clean_seqs))
            identical = (dseqnum == 1)
            diffs = ''
            if poly_count == 0:
                if identical:
                    diffs = 'Identical sequences.'
                else:
                    diffs = 'Variations are gaps only.'
            else:
                diffs = dseqnum

            total_sites.append({'group': group, 
                               'species': sp_name, 
                               'polymorphic sites': poly_count,
                               'aln length': lengths,
                               'differences': diffs
                               })

    df = pd.DataFrame(total_sites)
    df.to_csv(os.path.join(data_dir, 'diagnostic.csv'), index = False)

def run_mutation_analysis():
    '''
    Main analysis logic. Parses FASTA files, calculates rates,
    saves to CSV, and generates plots.
    '''

    all_data = []
    sites_threshold = 100

    for group, group_name in group_names.items():
        print(f'Processing {group_name}')
        input_dir = os.path.join(group, 'third_sites')
        if not os.path.exists(input_dir):
            print(f'Warning: Directory {input_dir} not found. Skipping.')
            continue

        sp_list = os.listdir(input_dir)
        for species in sp_list:
            try:
                sp_name = species.split('concatenate_')[1].split('.fasta')[0].replace('_',' ').replace(' endosymbiont', '')
            except IndexError:
                sp_name = species.replace('.fasta', '')

            sp_path = os.path.join(input_dir, species)
            if not os.path.exists(sp_path):
                print(f'No alignment found for {sp_name} in {group_name} group.')
                continue
            
            aln = list(SeqIO.parse(sp_path, 'fasta'))
            print(f'>Processing {sp_name}')

            count = polymorphism_total(aln)
            if count < sites_threshold:
                print(f'  Skipping {sp_name}: only {count} polymorphic sites found (threshold: {sites_threshold}).')
                continue

            pairs = combinations(aln, 2)

            for seq1, seq2 in pairs:
                if group == 'endosymb+relatives':
                    # Take free-living and endosymbiont pairs only
                    if ('_genomic' in seq1.id) == ('_genomic' in seq2.id):
                        continue
                    if '_genomic' in seq1.id:
                        ancestor, endosymb = seq1, seq2
                    else:
                        ancestor, endosymb = seq2, seq1
                    rates = calculate_rates(group_name, sp_name, ancestor, endosymb)
                    all_data.append(rates)
                elif group == 'relatives_only':
                    if '_genomic' not in seq1.id or '_genomic' not in seq2.id:
                        continue
                    # Since there is no ancestor, call function twice
                    rates1 = calculate_rates(group_name, sp_name, seq1, seq2)
                    all_data.append(rates1)

                    rates2 = calculate_rates(group_name, sp_name, seq2, seq1)
                    all_data.append(rates2)
                elif group == 'endosymb_only':
                    if '_genomic' in seq1.id or '_genomic' in seq2.id:
                        continue
                    rates1 = calculate_rates(group_name, sp_name, seq1, seq2)
                    all_data.append(rates1)

                    rates2 = calculate_rates(group_name, sp_name, seq2, seq1)
                    all_data.append(rates2)
    final_df = pd.DataFrame(all_data)
    final_df.to_csv(csv_path, index = False)

# def rate_shift_plot(df):
#     '''
#     Regression plots comparing rates between control groups to determine mutation rate shifts by species.
#     '''
#     distance_matrices = load_or_compute_pickle(os.path.join(files_dir, 'distances_endosymb+relatives.pkl'))
#     df_pivot = df.pivot_table(index = 'Species', 
#                               columns = 'Group', 
#                               values = mutation_types, 
#                               aggfunc = 'median').dropna()
#     for mut in mutation_types:
#         relative_rates = df_pivot[(mut, 'Free-Living Control Pairs')]
#         endosymbiont_rates = df_pivot[(mut, 'Endosymbiont Control Pairs')]

#     texts = []
#     median_distance = []
#     for species_key, distance_matrix in distance_matrices.items():
#             sp_name = species_key.replace('_endosymbiont', '').replace('_', ' ')
#             dist_list = []
#             distance_data = distance_matrix
#             i,j = np.triu_indices_from(distance_data, k=1)
#             ids = distance_data.index.tolist()

#             for index1, index2 in zip(i, j):
#                 id1 = ids[index1]
#                 id2 = ids[index2]
#                 if ('_genomic' in id1) == ('_genomic' in id2):
#                     continue
#                 dist_list.append(distance_data.loc[id1, id2])

#             median_distance.append(np.median(dist_list))
#             x = relative_rates.loc[sp_name]
#             y = endosymbiont_rates.loc[sp_name]     
            
#             t = plt.text(x, y, 
#                          sp_name, 
#                          fontsize=8, 
#                          alpha=0.7)
#             texts.append(t)

#     plt.figure(figsize=(8,8))
#     scatter = plt.scatter(x = relative_rates, 
#                     y = endosymbiont_rates,
#                     c = median_distance,
#                     cmap = 'viridis',
#                     s = 50,
#                     alpha = 0.7,
#                     label = 'Species')
#     cbar = plt.colorbar(scatter)
#     cbar.set_label('Median Genetic Distance\n(Endosymbiont-Relative Pairs)', rotation=270, labelpad=15)
    
#     plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (No shift)')
#     plt.xlim(-0.05,1.05)
#     plt.ylim(-0.05,1.05)
#     adjust_text(texts, 
#                 arrowprops=dict(arrowstyle="-", color='black', lw=0.5, alpha=0.7, shrinkA = 10),
#                 )

#     plt.title(f'{mut.replace("r","")} Median Rate Shift\nEndosymbionts vs Free-Living Relatives', fontsize=16)
#     plt.legend()
#     plt.tight_layout()
#     plt.xlabel('Free-Living Relatives')
#     plt.ylabel('Endosymbionts')
#     plt.savefig(os.path.join(plot_dir, f'rate_shift_{mut.replace("→","_")}.pdf'))
#     plt.close()

def rate_shift_plot(df):
    '''
    Regression plots comparing rates between control groups to determine mutation rate shifts by species.
    '''
    distance_matrices = load_or_compute_pickle(
                        os.path.join(files_dir, 'distances_endosymb+relatives.pkl'),
                        distance_matrix,
                        'endosymb+relatives'
                        )
    
    species_distances = {}
    for species_key, matrix in distance_matrices.items():
        sp_name = species_key.replace('_endosymbiont', '').replace('_', ' ')
        
        i, j = np.triu_indices_from(matrix, k=1)
        ids = matrix.index.tolist()
        
        dist_list = []
        for index1, index2 in zip(i, j):
            id1, id2 = ids[index1], ids[index2]
            # Only compare across endosymbiont-relative pairs
            if ('_genomic' in id1) == ('_genomic' in id2):
                continue
            dist_list.append(matrix.loc[id1, id2])
        
        species_distances[sp_name] = np.median(dist_list)

    df_pivot = df.pivot_table(index='Species', 
                              columns='Group', 
                              values=mutation_types, 
                              aggfunc='median').dropna()

    for mut in mutation_types:
        plt.figure(figsize=(8, 8))
        current_data = df_pivot[mut].dropna()
        
        relative_rates = current_data[('Free-Living Control Pairs')]
        endosymbiont_rates = current_data[('Endosymbiont Control Pairs')]
        
        # We need to align distances with the species present in the pivot table
        current_distances = [species_distances.get(sp, np.nan) for sp in current_data.index]
        
        scatter = plt.scatter(x=relative_rates, 
                              y=endosymbiont_rates,
                              c=current_distances,
                              cmap='magma',
                              edgecolors='black',
                              s=60,
                              alpha=0.8,
                              zorder=2,
                              label='Species')
        
        texts = []
        for sp_name in current_data.index:
            x_val = relative_rates.loc[sp_name]
            y_val = endosymbiont_rates.loc[sp_name]
            t = plt.text(x_val, y_val, sp_name, fontsize=8, alpha=0.8, zorder = 3)
            texts.append(t)
            
        cbar = plt.colorbar(scatter)
        cbar.set_label('Median Genetic Distance\n(Endosymbiont-Relative Pairs)', rotation=270, labelpad=25)
        
        plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (No shift)')
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)
        
        adjust_text(texts, 
            only_move={'points':'y', 'texts':'xy'}, 
            force_points=0.3,
            force_text=0.5, 
            expand_points=(1.8, 1.8),
            expand_text=(1.2, 1.2),
            add_objects=texts,
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5, alpha=0.3, zorder = 1))

        plt.title(f'{mut.replace("r","")} Median Rate Shift\nEndosymbionts vs Free-Living Relatives', fontsize=16)
        plt.xlabel('Free-Living Relatives')
        plt.ylabel('Endosymbionts')
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(os.path.join(plot_dir, f'rate_shift_{mut.replace("→","_")}.pdf'))
        plt.close()

if __name__ == '__main__':
    if os.path.exists(csv_path):
        # run_diagnostic()  
        df = load_data(csv_path)
        plot_distributions(df)
        plot_species_grid(df)
        rate_shift_plot(df)
    else:
        print('File not found, run mutation analysis first.')
        run_mutation_analysis()
        df = load_data(csv_path)
        plot_distributions(df)
        plot_species_grid(df)
