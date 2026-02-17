#!/usr/bin/env python3

import os
from adjustText import adjust_text
from itertools import combinations
from Bio import SeqIO
import pandas as pd
import matplotlib

from gcsize_dict import genome_gcsize
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from distance_matrix import distance_matrix
from gcsize_dict import genome_gcsize
from gc_codon_dict import gc_codon_dict
from utils import files_dir, genome_gcsize_json_path, load_or_compute_pickle, load_or_compute, gc_codon_json_path

plot_dir = os.path.join('plots', 'spectrum_plots')
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(files_dir, exist_ok=True)

site_labels = {
    'first_sites': 'First Codon Positions',
    'second_sites': 'Second Codon Positions',
    'third_sites': 'Fourfold Degenerate Third Sites'
}

# short filename prefixes for each analyzed site
site_file_prefix = {
    'first_sites': 'first',
    'second_sites': 'second',
    'third_sites': 'third'
}

group_colors = {
'Endosymbiont-Relative Pairs' : '#8da0cb',
'Free-Living Control Pairs' : '#66c2a5',
'Endosymbiont Control Pairs' : '#fc8d62'
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

def plot_distributions(df, position, plot_type = 'kde'):
    '''
    Plots histograms side-by-side and summary boxplot 
    '''
    median_df = df.groupby(['Group','Species'])[['rGC→AT', 'rAT→GC']].median().reset_index()

    site_label = site_labels[position]
    fig, axes = plt.subplots(1, 2, figsize=(18, 9), sharey=True)
    plt.suptitle(f'Comparison of Mutational Spectra: {site_label}', fontsize=20, y=0.98)

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
            ax.set_title(f'{mut} Rate Distributions', fontsize=16)
            ax.set_xlabel(f'Rate ({mut})', fontsize=14)
            ax.set_xlim(0, 1)
            ax.tick_params(axis='both', labelsize=12)
        axes[0].set_ylabel('Probability', fontsize=14)
        plt.tight_layout()
        
        outdir = os.path.join(plot_dir, position)
        prefix = site_file_prefix.get(position, position)
        plt.savefig(os.path.join(outdir, f'{prefix}_spectrum_hist.pdf'))
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
            ax.set_title(f'{mut} Rate Distributions', fontsize=16)
            ax.set_xlabel(f'Rate ({mut})', fontsize=14)
            ax.set_xlim(0, 1)
            ax.tick_params(axis='both', labelsize=12)
        axes[0].set_ylabel('Density', fontsize=14)
        plt.tight_layout()
        prefix = site_file_prefix.get(position, position)
        plt.savefig(os.path.join(plot_dir, position, f'{prefix}_spectrum_kde.pdf'))
        plt.close()
        print('Saved spectrum_kde.pdf')

    # Boxplot
    df_melted = median_df.melt(
        id_vars = ['Group', 'Species'],
        value_vars = mutation_types,
        var_name = 'Mutation Type',
        value_name = 'Rate'
    ).dropna(subset=['Rate'])

    plt.figure(figsize=(12, 10))
    sns.boxplot(data = df_melted,
                x = 'Group', 
                y = 'Rate',
                hue = 'Mutation Type', 
                palette = 'Paired',
                fliersize = 1,
                width = 0.5
                )
    plt.title('Summary of Mutation Rates by Group', fontsize=20)
    plt.xlabel('Group', fontsize=16)
    plt.ylabel('Rate', fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()

    prefix = site_file_prefix.get(position, position)
    outpath = os.path.join(plot_dir, position, f'{prefix}_spectrum_boxplot.pdf')
    plt.savefig(outpath)
    plt.close()
    print(f'Saved {outpath}')

def plot_separated_group_boxplots(df, position):
    '''
    Creates side-by-side boxplots for Endosymbionts and Relatives
    X-axis is Mutation Type (GC→AT vs AT→GC) 
    '''
    site_label = site_labels[position]
    
    # 1. Filter groups and clean up names for the presentation
    keep_groups = ['Endosymbiont Control Pairs', 'Free-Living Control Pairs']
    df_filtered = df[df['Group'].isin(keep_groups)].copy()
    df_filtered['Group'] = df_filtered['Group'].replace({
        'Endosymbiont Control Pairs': 'Endosymbionts',
        'Free-Living Control Pairs': 'Free-Living Relatives'
    })
    
    # 2. Take the median per species (Crucial for smoothing out noise)
    median_df = df_filtered.groupby(['Group', 'Species'])[mutation_types].median().reset_index()

    # 3. Melt the dataframe for Seaborn
    df_melted = median_df.melt(
        id_vars = ['Group', 'Species'],
        value_vars = mutation_types,
        var_name = 'Mutation Type',
        value_name = 'Rate'
    ).dropna(subset=['Rate'])

    groups_to_plot = ['Endosymbionts', 'Free-Living Relatives']

    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)

    for ax, grp in zip(axes, groups_to_plot):
        subset = df_melted[df_melted['Group'] == grp]
        
        sns.boxplot(data = subset,
                    x = 'Mutation Type', 
                    y = 'Rate',
                    hue = 'Mutation Type',
                    palette = 'Paired', 
                    width = 0.5,
                    showfliers = False,
                    ax = ax
                    )
        
        ax.set_title(f'{grp}', fontsize=24, fontweight='bold')
        ax.set_xlabel('Mutation Type', fontsize=20, fontweight='bold')
        
        display_muts = [mut.replace('r', '') for mut in mutation_types]
        ax.set_xticks([0, 1])
        ax.set_xticklabels(display_muts, fontsize=18)
        
        if ax == axes[0]:
            ax.set_ylabel('Median Rate per Species', fontsize=20, fontweight='bold')
            ax.tick_params(axis='y', labelsize=18)
        else:
            ax.set_ylabel('')
            ax.tick_params(axis='y', length=0) 
            
        ax.set_ylim(-0.05, 1) 
        ax.grid(True, which="major", axis="y", ls="-", alpha=0.5)

    plt.suptitle(f'Mutational Spectra by Group:\n{site_label}', fontsize=26)
    plt.tight_layout()

    prefix = site_file_prefix.get(position, position)
    outpath = os.path.join(plot_dir, position, f'{prefix}_separated_spectrum_boxplot.pdf')
    
    plt.savefig(outpath)
    plt.close()
    print(f'Saved {outpath}')

def plot_species_grid(df, position):
    '''
    Generates a grid of histogram plots for each species in the CSV.
    Only plots if species are shared by both groups.
    Produces two files: one for GC->AT rates and one for AT->GC rates.
    '''

    site_label = site_labels[position]
    # Include only species present in both groups
    group1 = set(df[df['Group'] == 'Endosymbiont-Relative Pairs']['Species'])
    group2 = set(df[df['Group'] == 'Free-Living Control Pairs']['Species'])
    common_species = group1.intersection(group2)
    df = df[df['Species'].isin(common_species)].copy()
    print(f'Plotting {len(common_species)} species shared by both groups.')

    sns.set_style('whitegrid')
    
    for mut in mutation_types:
        safe_name = mut.replace('→', '_').replace('r', '')
        prefix = site_file_prefix.get(position, position)
        output = f'{prefix}_grid_{safe_name}.pdf'
    
        grid = sns.FacetGrid(
            df, 
            col='Species', 
            hue='Group',
            palette = group_colors,
            hue_order = list(group_colors.keys()),
            col_wrap=7,
            sharex=True,
            sharey=False,
            height=3.5,
            aspect=1.25
        )
        
        grid.map(sns.histplot, mut, stat='density',
                common_norm=False, fill=True, alpha=0.5, 
                bins='auto', element = 'step')
        grid.add_legend(title='Group', title_fontsize=14, fontsize=12)
        grid.set_titles('{col_name}')
        grid.set_axis_labels(f'Rate ({mut})', 'Density')
        grid.set(xlim=(0, 1))
        grid.figure.suptitle(f'{mut} Rates by Species\n{site_label}', fontsize=20, y =1.02)
        for ax in grid.axes.flatten():
            ax.tick_params(axis='both', labelsize=12)
            ax.set_xlabel(ax.get_xlabel(), fontsize=12)
            ax.set_ylabel(ax.get_ylabel(), fontsize=12)
        outpath = os.path.join(plot_dir, position, output)
        grid.savefig(outpath)
        plt.close()
        print(f'Saved {outpath}')

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
    df.to_csv(os.path.join(files_dir, 'diagnostic.csv'), index = False)

def run_mutation_analysis(position):
    '''
    Main analysis logic. Parses FASTA files, calculates rates,
    saves to CSV, and generates plots.
    Args:
        position (str): Position type to analyze ('first_sites', 'second_sites', 'third_sites').
    '''

    output_csv = f'{position}_spectrum_rates.csv'
    csv_path = os.path.join(files_dir, output_csv)
    all_data = []
    sites_threshold = 100

    for group, group_name in group_names.items():
        print(f'Processing {group_name}')
        input_dir = os.path.join(group, position)
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

def rate_shift_plot(df, position, color_by = 'distance', annotate=False):
    '''
    Plots comparing mutation rates between endosymbionts and relatives, colored by genetic distance or GC content.
    '''
    site_label = site_labels[position]
    species_metric = {}

    if color_by not in ['distance', 'gc_genome']:
        print('Invalid color_by option. Choose "distance" or "gc_genome".')
        return
    if color_by == 'distance':
        matrices = load_or_compute_pickle(
                            os.path.join(files_dir, 'distances_endosymb+relatives.pkl'),
                            distance_matrix,
                            'endosymb+relatives'
                            )
        
        for species_key, matrix in matrices.items():
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
            
            if dist_list:
                species_metric[sp_name] = np.median(dist_list)

        cmap = 'magma'
        cbar_label = 'Median Genetic Distance\n(Endosymbiont-Relative Pairs)'

    elif color_by == 'gc_genome':
        genome_dataset = load_or_compute(
                            genome_gcsize_json_path('endosymb+relatives'),
                            genome_gcsize,
                            'endosymb+relatives'
                            )

        for sp, genomes in genome_dataset.items():
            sp_name = ' '.join(sp.split('_endosymbiont')[0].split('_'))

            endosymb_gcs = [g['gc_genome'] for genome_id, g in genomes.items() if '_genomic' not in genome_id]
            relatives_gcs = [g['gc_genome'] for genome_id, g in genomes.items() if '_genomic' in genome_id]

            if endosymb_gcs and relatives_gcs:
                species_metric[sp_name] = abs(np.median(endosymb_gcs) - np.median(relatives_gcs))

        cmap = 'Greys'
        cbar_label = 'Absolute ΔGC (%)'
    

    df_pivot = df.pivot_table(index='Species', 
                              columns='Group', 
                              values=mutation_types, 
                              aggfunc='median').dropna()

    for mut in mutation_types:

        plt.figure(figsize=(12, 12))
        current_data = df_pivot[mut].dropna()
        
        relative_rates = current_data[('Free-Living Control Pairs')]
        endosymbiont_rates = current_data[('Endosymbiont Control Pairs')]
        
        current_metric = [species_metric.get(sp, np.nan) for sp in current_data.index]

        if color_by == 'gc_genome':
            limit = np.nanmax(np.abs(current_metric))
            vmax_threshold = limit
            vmin_val = 0
        elif color_by == 'distance':
            # Cap color scale at 95th percentile to avoid outlier compression
            vmax_threshold = np.nanpercentile(current_metric, 95)
            vmin_val = np.nanmin(current_metric)

        scatter = plt.scatter(x=relative_rates, 
                              y=endosymbiont_rates,
                              c=current_metric,
                              cmap=cmap,
                              edgecolors='black',
                              s=100,
                              alpha=0.85,
                              vmin=vmin_val,
                              vmax=vmax_threshold,
                              zorder=2,
                              label='Species' if annotate else None)

        if annotate:
            texts = []
            for sp_name in current_data.index:
                x_val = relative_rates.loc[sp_name]
                y_val = endosymbiont_rates.loc[sp_name]
                t = plt.text(x_val, y_val, sp_name, fontsize=10, alpha=0.9, zorder = 3)
                texts.append(t)
            adjust_text(texts, 
                only_move={'points':'y', 'texts':'xy'}, 
                force_points=0.03,
                force_text=2.0, 
                expand_points=(2.1, 2.1),
                expand_text=(1.5, 1.5),
                autoalign='y',
                add_objects=[scatter]
                )

        cbar = plt.colorbar(scatter)
        cbar.set_label(cbar_label, rotation=270, labelpad=25, fontsize=12)
        cbar.ax.tick_params(labelsize=12)
        
        plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (No shift)')
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)

        plt.title(f'{mut.replace("r","")} Rate Shift\n{site_label}', fontsize=20)
        plt.xlabel('Free-Living Relatives', fontsize=16)
        plt.ylabel('Endosymbionts', fontsize=16)
        plt.tick_params(axis='both', labelsize=14)
        plt.legend(fontsize=12)
        plt.tight_layout()
        mut_file = mut.replace('→', '_')
        prefix = site_file_prefix.get(position, position)
        outpath = os.path.join(plot_dir, position, f'{prefix}_{mut_file}_rate_shift_{color_by}.pdf')
        plt.savefig(outpath)
        plt.close()
        print(f'Saved {outpath}')

def rate_shift_plot_gc(df, position, annotate=False):
    site_label = site_labels[position]

    genome_json = genome_gcsize_json_path('endosymb_only')
    genome_dataset = load_or_compute(genome_json, genome_gcsize, 'endosymb_only')

    sp_gc_values = {}
    for sp, genomes in genome_dataset.items():
        sp_name = ' '.join(sp.split('_endosymbiont')[0].split('_'))

        gc_data = [data.get('gc_genome', np.nan) for data in genomes.values()]
        if gc_data:
            sp_gc_values[sp_name] = np.mean(gc_data)
        else:
            sp_gc_values[sp_name] = np.nan
    
    df_pivot = df.pivot_table(index='Species', 
                              columns='Group', 
                              values=mutation_types, 
                              aggfunc='median').dropna()
    for mut in mutation_types:
        plt.figure(figsize=(12, 12))
        current_data = df_pivot[mut].dropna()
        
        relative_rates = current_data[('Free-Living Control Pairs')]
        endosymbiont_rates = current_data[('Endosymbiont Control Pairs')]
        
        current_gcs = [sp_gc_values.get(sp, np.nan) for sp in current_data.index]

        scatter = plt.scatter(x=relative_rates, 
                              y=endosymbiont_rates,
                              c=current_gcs,
                              cmap='viridis',
                              edgecolors='black',
                              s=100,
                              alpha=0.85,
                              zorder=2,
                              label='Species' if annotate else None)
        if annotate:
            texts = []
            for sp_name in current_data.index:
                x_val = relative_rates.loc[sp_name]
                y_val = endosymbiont_rates.loc[sp_name]
                t = plt.text(x_val, y_val, sp_name, fontsize=10, alpha=0.9, zorder = 3)
                texts.append(t)
            adjust_text(texts, 
            only_move={'points':'y', 'texts':'xy'}, 
            force_points=0.03,
            force_text=2.0, 
            expand_points=(2.1, 2.1),
            expand_text=(1.5, 1.5),
            autoalign='y',
            add_objects=[scatter]
            )
            
        cbar = plt.colorbar(scatter)
        cbar.set_label('Endosymbiont GC Content (%)', rotation=270, labelpad=25, fontsize=12)
        cbar.ax.tick_params(labelsize=12)
        
        plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (No shift)')
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)

        plt.title(f'{mut.replace("r","")} Rate Shift\n{site_label}', fontsize=20)
        plt.xlabel('Free-Living Relatives', fontsize=16)
        plt.ylabel('Endosymbionts', fontsize=16)
        plt.tick_params(axis='both', labelsize=14)
        plt.legend(fontsize=12)
        plt.tight_layout()

        mut_file = mut.replace('→', '_')
        prefix = site_file_prefix.get(position, position)
        outpath = os.path.join(plot_dir, position, f'{prefix}_gcabs_rate_shift_{mut_file}.pdf')
        plt.savefig(outpath)
        plt.close()
        print(f'Saved {outpath}')

def gc_equilibrium(df):
    '''
    GC equilibrium plot comparing observed GC content to predicted equilibrium GC based on mutation rates.
    '''

    # Use the mapped mutation column names (they contain the unicode arrow)
    r_gc_at, r_at_gc = mutation_types[0], mutation_types[1]
    median_rates = df.groupby(['Group','Species'])[[r_gc_at, r_at_gc]].median().reset_index()
    median_rates = median_rates[median_rates['Group'] != 'Endosymbiont-Relative Pairs'].copy()
    # Calculate predicted equilibrium GC content
    median_rates['GC_eq'] = (median_rates[r_at_gc] / (median_rates[r_at_gc] + median_rates[r_gc_at])) * 100
    
    # Map group names to match fourfold_gc
    group_map = {
        'Endosymbiont Control Pairs': 'Endosymbionts',
        'Free-Living Control Pairs': 'Relatives'
    }
    median_rates['Group'] = median_rates['Group'].map(group_map)
    
    # Merge with observed GC content
    fourfold_gc = pd.read_csv(os.path.join(files_dir, 'fourfold_gc_content.csv'))
    fourfold_gc = fourfold_gc.rename(columns={'group': 'Group', 'species': 'Species', 'overall_gc4': 'GC_obs'})

    merged_df = pd.merge(median_rates, fourfold_gc, on=['Species', 'Group'], how='inner')
    
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=merged_df,
                    x='GC_eq',
                    y='GC_obs',
                    hue='Group',
                    palette = {'Endosymbionts': '#FC8D62FF', 'Relatives': '#66C2A5FF'},
                    alpha = 0.8, 
                    s = 100,
                    edgecolor='black'
                    )
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='Observed = Equilibrium')
    plt.gca().set_aspect('equal')
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.title('GC Equilibrium vs Observed GC4', fontsize=16)
    plt.ylabel('Observed GC4 (Fourfold Sites)')
    plt.xlabel('Predicted Equilibrium GC ($GC_{eq}$)')

    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.legend()

    prefix = site_file_prefix.get('third_sites', 'third')
    outpath = os.path.join(plot_dir, 'third_sites', f'{prefix}_gc_equilibrium_plot.pdf')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
    print(f'Saved {outpath}')


def combined_sites_boxplot(strip=False):
    """
    Combined faceted boxplots comparing mutation rates across sites
    for the two control groups: Endosymbiont Control Pairs and Free-Living Control Pairs.

    Produces a faceted plot (one facet per mutation type) as well as individual plots for each mutation type.
    """

    positions = ['first_sites', 'second_sites', 'third_sites']
    dfs = []
    for site in positions:
        csv_name = f'{site}_spectrum_rates.csv'
        csv_path = os.path.join(files_dir, csv_name)
        if not os.path.exists(csv_path):
            print(f'Warning: {csv_path} not found, skipping {site} in combined plot.')
            continue
        temp = load_data(csv_path)
        if temp.empty:
            continue
        temp['Site'] = site_file_prefix[site]
        dfs.append(temp)

    if not dfs:
        print('No data available for combined sites plot.')
        return

    combined = pd.concat(dfs, ignore_index=True)

    keep_groups = ['Endosymbiont Control Pairs', 'Free-Living Control Pairs']
    combined = combined[combined['Group'].isin(keep_groups)].copy()
    combined['Group'] = combined['Group'].replace({
        'Endosymbiont Control Pairs': 'Endosymbionts',
        'Free-Living Control Pairs': 'Free-Living Relatives'
    })

    med = combined.groupby(['Site', 'Group', 'Species'])[mutation_types].median().reset_index()

    med_melt = med.melt(id_vars=['Site', 'Group', 'Species'], value_vars=mutation_types,
                        var_name='Mutation Type', value_name='Rate')

    site_order = ['first', 'second', 'third']
    display_labels = ['First', 'Second', 'Third (4-fold)']

    # Create faceted boxplot
    fig, axes = plt.subplots(1, 2, figsize=(18, 12), sharey=True)

    for ax, mut in zip(axes, mutation_types):
        sns.boxplot(data=med_melt[med_melt['Mutation Type'] == mut],
                    x='Site', 
                    y='Rate', 
                    hue='Group', 
                    palette = {'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
                    order=site_order,
                    ax=ax,
                    width=0.5,
                    showfliers=False
                    )
        if strip:
            sns.stripplot(data=med_melt[med_melt['Mutation Type'] == mut],
                        x='Site', 
                        y='Rate', 
                        hue='Group', 
                        palette = {'Endosymbionts': 'gray', 'Free-Living Relatives': 'gray'},
                        order=site_order,
                        ax=ax,
                        dodge=True,
                        alpha=0.4,
                        legend=False
                        )
        ax.set_title(f'{mut.replace("r","")} Rates by Site', fontsize=18, fontweight='bold')
        ax.set_xlabel('')
        ax.set_xticks([0,1,2])
        ax.set_xticklabels(display_labels, fontsize=18)
        
        ax.set_ylabel('Mutation Rate', fontsize=22, fontweight='bold')
        ax.tick_params(axis='y', labelsize=18)
        ax.set_ylim(-0.05, 1)

        ax.grid(True, which="major", axis="y", ls="-", alpha=0.5)
        ax.grid(True, which="minor", axis="y", ls=":", alpha=0.2)

        if ax == axes[1]:
            ax.legend(title='Group', fontsize=16, title_fontsize=18)
        else:
            ax.get_legend().remove()
    fig.suptitle('Comparison of Mutation Rates Across Codon Sites', fontsize=26, fontweight='bold')
    fig.supxlabel('Codon Site', fontsize=22, fontweight='bold')

    plt.tight_layout()

    outpath = os.path.join(plot_dir, f'combined_sites_boxplot.pdf')
    plt.savefig(outpath)
    plt.close()
    print(f'Saved {outpath}')

    for mut in mutation_types:
        plt.figure(figsize=(14, 12))
        
        ax = sns.boxplot(data=med_melt[med_melt['Mutation Type'] == mut],
                    x='Site', 
                    y='Rate', 
                    hue='Group', 
                    palette = {'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
                    order=site_order,
                    width=0.5,
                    showfliers=False
                    )
        if strip:
            sns.stripplot(data=med_melt[med_melt['Mutation Type'] == mut],
                        x='Site', 
                        y='Rate', 
                        hue='Group', 
                        palette = {'Endosymbionts': 'gray', 'Free-Living Relatives': 'gray'},
                        order=site_order,
                        dodge=True,
                        alpha=0.4,
                        legend=False,
                        ax = ax
                        )
        plt.title(f'{mut.replace("r","")} Rates by Site', fontsize=30, fontweight='bold')
        plt.xlabel('Codon Site', fontsize=26, fontweight='bold')
        plt.xticks([0,1,2], display_labels, fontsize=22)

        plt.ylabel('Mutation Rate', fontsize=26, fontweight='bold')
        plt.yticks(fontsize=22)
        plt.ylim(-0.05, 1)

        plt.grid(True, which="major", axis="y", ls="-", alpha=0.5)
        plt.grid(True, which="minor", axis="y", ls=":", alpha=0.2)

        plt.legend(title='Group', fontsize = 18, title_fontsize = 20)
        mut_file = mut.replace('→', '_').replace('r', '')
        outpath_single = os.path.join(plot_dir, f'{mut_file}_sites_boxplot.pdf')

        plt.tight_layout()

        plt.savefig(outpath_single)
        plt.close()
        print(f'Saved {outpath_single}')

    
if __name__ == '__main__':
    for position in ['first_sites', 'second_sites', 'third_sites']:

        output_csv = f'{position}_spectrum_rates.csv'
        csv_path = os.path.join(files_dir, output_csv)

        if not os.path.exists(csv_path):
            print('File not found, run mutation analysis first.')
            run_mutation_analysis(position)
        
        if not os.path.exists(csv_path):
            print(f'Error: {csv_path} could not be generated. Skipping {position}.')
            continue
        
        print(f'Loading data from {csv_path}')
        df = load_data(csv_path)

        if df.empty:
            print(f'No data found in {csv_path}, skipping {position}.')
            continue

        outdir = os.path.join(plot_dir, position)
        os.makedirs(outdir, exist_ok=True)

        plot_distributions(df, position)
        plot_species_grid(df, position)
        rate_shift_plot(df, position, color_by='distance')
        rate_shift_plot(df, position, color_by='gc_genome')
        rate_shift_plot_gc(df, position)

        if position == 'third_sites':
            plot_separated_group_boxplots(df, position)
            gc_equilibrium(df)
            combined_sites_boxplot()