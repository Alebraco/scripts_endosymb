#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from transposase_analysis import processing_transposase
from sequence_features import collect_codon_stats
from igs_lengths import collect_gff_stats
from utils import files_dir

def main():
    parser = argparse.ArgumentParser(description='Collect and merge genome features.')
    parser.add_argument('--path', required=True, help='Path to process')
    parser.add_argument('--infer', action='store_true', help='Inference mode: do not auto-classify by filename, default to Ungrouped')
    args = parser.parse_args()

    path = args.path
    if not os.path.exists(path):
        raise FileNotFoundError(f'Input path does not exist: {path}')

    print('Starting feature collection process.')

    if args.infer:
        classify_flag = False
        print('Running in inference mode (auto classification disabled; Group defaults to Ungrouped).')
    else:
        classify_flag = True
        print('Running in training mode (auto classification enabled).')

    transposase_df = processing_transposase(path, auto_classify=classify_flag)
    if transposase_df is None or transposase_df.empty:
        print('No transposase records found. Exiting.')
        return
    print('Done with transposase analysis.')

    igs_data, gene_data = collect_gff_stats(path, auto_classify=classify_flag)
    if not igs_data or not gene_data:
        print('No IGS or gene records found. Exiting.')
        return

    igs_df = pd.DataFrame(igs_data).groupby(['File'])['IGS_Size'].mean().reset_index()
    igs_df = igs_df.rename(columns={'IGS_Size': 'Mean_IGS_Size'})
    print('Done with IGS analysis.')

    gene_df = pd.DataFrame(gene_data)
    gene_count_df = gene_df.groupby(['File']).size().reset_index(name='Gene_Count')
    gene_length_df = gene_df.groupby(['File'])['Gene_Length'].mean().reset_index(name='mean_gene_length')
    gene_df = pd.merge(gene_count_df, gene_length_df, on='File')

    print('Done with gene length and count analysis.')

    codon_stats_df = collect_codon_stats(path, auto_classify=classify_flag)
    print('Done with GC metrics.')

    merged_df = pd.merge(transposase_df, igs_df, on='File', how='inner')
    merged_df = pd.merge(merged_df, gene_df, on='File', how='inner')
    merged_df = pd.merge(merged_df, codon_stats_df, on='File', how='inner')

    merged_df = merged_df[merged_df['Gene_Count'] > 0]
    merged_df['Transposase_Per_Gene'] = merged_df['Total_Transposases'] / merged_df['Gene_Count']
    merged_df['Delta_GC2_4'] = merged_df['GC2'] - merged_df['GC4']


    merged_df.to_csv(os.path.join(path, 'combined_features.csv'), index=False)
    print('All features saved to combined_features.csv')


if __name__ == '__main__':
    main()