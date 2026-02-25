#!/usr/bin/env python3
import argparse
import os

from distance_matrix import distance_matrix
from gcsize_dict import genome_gcsize
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import (
    load_or_compute,
    load_or_compute_pickle,
    genome_gcsize_json_path,
    files_dir
)

plot_dir = os.path.join('plots', 'endosymb_evolution')
os.makedirs(plot_dir, exist_ok=True)
# Size - GC - Evolutionary Distance Plot
def sge_data(mode='mean', seq_type='protein'):

    group = 'endosymb+relatives'
    genome_json = genome_gcsize_json_path(group)
    suffix = f'_{seq_type}' if seq_type != 'dna' else ''
    distance_pkl = os.path.join(files_dir, f'distances_{group}{suffix}.pkl')

    genome_dataset = load_or_compute(genome_json, genome_gcsize, group)
    distance_matrices = load_or_compute_pickle(distance_pkl, distance_matrix, group, mode=seq_type)
    print('All data has been loaded.')

    all_data = []


    for species, matrix in distance_matrices.items():
        if species not in genome_dataset:
            continue

        sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
        
        genome_data = genome_dataset.get(species)
        ids = matrix.index.tolist()

        endosymb_ids = [genome_id for genome_id in ids if '_genomic' not in genome_id]

        i, j = np.triu_indices_from(matrix, k=1)        
        dist_list = []
        for index1, index2 in zip(i, j):
            id1, id2 = ids[index1], ids[index2]
            # Only compare across endosymbiont-relative pairs
            if ('_genomic' in id1) == ('_genomic' in id2):
                continue
            dist_list.append(matrix.loc[id1, id2])
            
        median_distance = np.median(dist_list)

        if mode == 'mean':
            gc = np.median([genome_data[genome_id]['gc_genome'] for genome_id in endosymb_ids])
            size = np.median([genome_data[genome_id]['size'] for genome_id in endosymb_ids])
            
            all_data.append({
                'species': sp_name,
                'gc': gc,
                'size': size,
                'distance': median_distance,
                'n_endosymbionts': len(endosymb_ids)
            })

        else:
            for genome_id in endosymb_ids:
                all_data.append({
                    'species': sp_name,
                    'gc': genome_data[genome_id]['gc_genome'],
                    'size': genome_data[genome_id]['size'],
                    'distance': median_distance
                })

    df = pd.DataFrame(all_data)
    suffix = f'_{seq_type}' if seq_type != 'dna' else ''
    df.to_csv(os.path.join(files_dir, f'endosymb_evolution_data{suffix}.csv'), index=False)
    return df

def seg_plot(df, seq_type='protein'):
    suffix = f'_{seq_type}' if seq_type != 'dna' else ''
    plt.figure(figsize=(16,12))
    scatter = plt.scatter(
        x = df['distance'],         
        y = df['size'],  
        c = df['gc'],
        cmap = 'viridis', 
        s = 50, 
        edgecolors = 'black',
        linewidths = 0.5
        )
    cbar = plt.colorbar(scatter)
    cbar.set_label('GC Content (%)', fontsize=18, fontweight='bold')
    cbar.ax.tick_params(labelsize=14)

    plt.xlabel('Median Evolutionary Distance', fontsize=22, fontweight='bold')
    plt.ylabel('Genome Size (Mb)' , fontsize=22, fontweight='bold')
    plt.tick_params(axis='both', labelsize=14)
    plt.title('The Evolution of Endosymbionts\nGenome Size, GC Content, and Evolutionary Distance', fontsize=24, fontweight='bold')
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'endosymb_evolution{suffix}.pdf'))
    plt.close()

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--seq-type', default='protein', choices=['protein', 'dna'],
                    help='Distance matrix type to use (default: protein)')
    args = ap.parse_args()

    suffix = f'_{args.seq_type}' if args.seq_type != 'dna' else ''
    csv_path = os.path.join(files_dir, f'endosymb_evolution_data{suffix}.csv')
    if os.path.exists(csv_path):
        print('Loading existing data.')
        df = pd.read_csv(csv_path)
    else:
        print('Computing data.')
        df = sge_data(mode='mean', seq_type=args.seq_type)
    seg_plot(df, seq_type=args.seq_type)
    print('Plot saved.')


