from itertools import combinations
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from utils import group_names, fourfold_families
import os

def codon_polymorphisms(alignment_dir):
    all_data = []
    gap_set = {'-', 'N'}

    for folder, label in group_names.items():
        if folder in ['relatives_only', 'endosymb_only']: 
            dna_dir = os.path.join(folder, alignment_dir)
            if not os.path.isdir(dna_dir):
                continue

            for filename in os.listdir(dna_dir):
                if not filename.endswith('.fasta'):
                    continue
                file_path = os.path.join(dna_dir, filename)
                sp_name = filename.replace('_endosymbiont.fasta', '').replace('concatenate_', '').replace('_', ' ')
                sequences = list(SeqIO.parse(file_path, 'fasta'))

                # Iterate over all pairs of sequences in the alignment
                for rec1, rec2 in combinations(sequences, 2):
                    s1, s2 = str(rec1.seq).upper(), str(rec2.seq).upper()

                    # Identify the number of matches between sequences
                    matches = sum(1 for a, b in zip(s1, s2) if a == b and a != '-')
                    # Identify the number of valid sites in both sequences
                    total_sites = sum(1 for a, b in zip(s1, s2) if a != '-' and b != '-')
                    if total_sites == 0:
                        continue
                    # Calculate the percent identity value
                    pid = matches / total_sites
                    # Filter too similar and too distant sequences
                    if not (0.7 < pid < 0.95):
                        continue

                    # Initialize variables for all positions
                    first_total = first_diff = 0
                    second_total = second_diff = 0
                    fourfold_total = fourfold_diff = 0

                    # Iterate over each codon in the sequence
                    for i in range(0, len(s1) - 2, 3):
                        # Capture the current codon in both sequences
                        c1, c2 = s1[i:i+3], s2[i:i+3]
                        # Filter incomplete codons
                        if len(c1) < 3 or len(c2) < 3:
                            continue

                        # First position
                        if c1[0] not in gap_set and c2[0] not in gap_set:
                            first_total += 1
                            if c1[0] != c2[0]:
                                first_diff += 1

                        # Second position
                        if c1[1] not in gap_set and c2[1] not in gap_set:
                            second_total += 1
                            if c1[1] != c2[1]:
                                second_diff += 1

                        # Fourfold-degenerate third position only
                        # Both codons must have a fourfold-degenerate prefix
                        prefix1, prefix2 = c1[0:2], c2[0:2]
                        if (prefix1 in fourfold_families and prefix2 in fourfold_families
                                and c1[2] not in gap_set and c2[2] not in gap_set):
                            fourfold_total += 1
                            if c1[2] != c2[2]:
                                fourfold_diff += 1

                    all_data.append({
                        'group': label,
                        'species': sp_name,
                        'id1': rec1.id,
                        'id2': rec2.id,
                        'pid': round(pid, 4),
                        'poly_first':    round(first_diff / first_total, 4) if first_total else float('nan'),
                        'poly_second':   round(second_diff / second_total, 4) if second_total else float('nan'),
                        'poly_fourfold': round(fourfold_diff / fourfold_total, 4) if fourfold_total else float('nan'),
                    })

    df = pd.DataFrame(all_data)
    out = os.path.join('files', 'codon_polymorphisms.csv')
    os.makedirs('files', exist_ok=True)
    df.to_csv(out, index=False)
    print(f'Saved {out}')
    return df


def plot_codon_polymorphisms(df):
    '''
    Boxplot comparing polymorphism frequencies across codon positions and groups.
    '''

    plot_dir = os.path.join('plots', 'codon_polymorphisms')
    os.makedirs(plot_dir, exist_ok=True)

    med = df.groupby(['group', 'species'])[['poly_first', 'poly_second', 'poly_fourfold']].median().reset_index()
    melted = med.melt(
        id_vars=['group', 'species'],
        value_vars=['poly_first', 'poly_second', 'poly_fourfold'],
        var_name='position', value_name='poly_freq'
    )
    position_labels = {
        'poly_first': 'First',
        'poly_second': 'Second',
        'poly_fourfold': 'Fourfold (3rd)'
    }
    melted['position'] = melted['position'].map(position_labels)

    plt.figure(figsize=(12, 9))
    sns.boxplot(
        data=melted, x='position', y='poly_freq', hue='group',
        order=['First', 'Second', 'Fourfold (3rd)'],
        palette = {'Endosymbionts': '#FC8D62FF', 'Relatives': '#66C2A5FF'},
        width=0.5, showfliers=False
    )
    plt.title('Polymorphism Frequency by Codon Position and Group', fontsize=20, fontweight='bold')
    plt.xlabel('Codon Position', fontsize=16, fontweight='bold')
    plt.ylabel('Polymorphism Frequency', fontsize=16, fontweight='bold')
    plt.legend(title='Group', fontsize=14, title_fontsize=14)
    plt.tight_layout()
    outpath = os.path.join(plot_dir, 'codon_polymorphisms_boxplot.pdf')
    plt.savefig(outpath)
    plt.close()
    print(f'Saved {outpath}')


def plot_codon_polymorphisms_scatter(df):
    '''
    Scatter plot comparing polymorphism frequencies
    Plots endosymbionts vs. free-living relatives for each codon position.
    '''

    plot_dir = os.path.join('plots', 'codon_polymorphisms')
    os.makedirs(plot_dir, exist_ok=True)

    g_rel  = "Free-Living Relatives Only"
    g_endo = "Endosymbionts Only"

    positions = ['poly_first', 'poly_second', 'poly_fourfold']
    position_labels = {
        'poly_first':    'First',
        'poly_second':   'Second',
        'poly_fourfold': 'Fourfold (3rd)',
    }

    df_pivot = df.pivot_table(index='species', columns='group', values=positions, aggfunc='median').dropna()

    out = os.path.join('files', 'codon_polymorphisms_scatter.csv')
    df_pivot.to_csv(out)
    print(f'Saved {out}')

    for pos in positions:
        current_data = df_pivot[pos].dropna()

        relative_rates = current_data[g_rel]
        endosymbiont_rates = current_data[g_endo]

        plt.figure(figsize=(8, 8))
        plt.scatter(x=relative_rates, y=endosymbiont_rates,
                    edgecolors='black', s=80, alpha=0.85)
        
        plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (No change)')
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)

        plt.xlabel('Free-Living Relatives', fontsize=14, fontweight='bold')
        plt.ylabel('Endosymbionts', fontsize=14, fontweight='bold')
        plt.title(f'Polymorphism Frequency\n{position_labels[pos]} Position',
                  fontsize=16, fontweight='bold')
        plt.tight_layout()
        outpath = os.path.join(plot_dir, f'scatter_{pos}.pdf')
        plt.savefig(outpath, bbox_inches='tight')
        plt.close()
        print(f'Saved {outpath}')


if __name__ == '__main__':
    csv_path = os.path.join('files', 'codon_polymorphisms.csv')
    if not os.path.exists(csv_path):
        alignment_dir = 'dna_concatenates'
        df = codon_polymorphisms(alignment_dir)
    else:
        df = pd.read_csv(csv_path)
    plot_codon_polymorphisms(df)
    plot_codon_polymorphisms_scatter(df)
