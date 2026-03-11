#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from transposase_analysis import processing_transposase
from sequence_features import collect_codon_stats
from igs_lengths import collect_gff_stats

def main():
    parser = argparse.ArgumentParser(description='Collect and merge genome features.')
    parser.add_argument('--path', required=True, help='Path to process')
    parser.add_argument('--infer', action='store_true', help='Inference mode: do not auto-classify by filename, default to Ungrouped')
    parser.add_argument('--force', action='store_true', help='Recompute features even if output file already exists')
    args = parser.parse_args()

    path = args.path
    if not os.path.exists(path):
        raise FileNotFoundError(f'Input path does not exist: {path}')

    feature_dir = os.path.join(path, 'feature_files')
    output_csv = os.path.join(feature_dir, 'combined_features.csv')

    if os.path.exists(output_csv) and not args.force:
        print(f'Found existing feature file: {output_csv}')
        print('Skipping analysis. Use --force to recompute.')
        return

    print('Starting feature collection process.')

    if args.infer:
        classify_flag = False
        print('Running in inference mode (auto classification disabled; Group defaults to Ungrouped).')
    else:
        classify_flag = True
        print('Running in training mode (auto classification enabled).')

    transposase_df = processing_transposase(path, auto_classify=classify_flag, get_fam_list = False)
    if transposase_df is None or transposase_df.empty:
        print('No transposase records found. Exiting.')
        return
    transposase_df = transposase_df.rename(columns={'Total': 'Total_Transposases'})
    print('Done with transposase analysis.')

    igs_data, gene_data = collect_gff_stats(path, auto_classify=classify_flag)
    if not igs_data or not gene_data:
        print('No IGS or gene records found. Exiting.')
        return

    igs_df = pd.DataFrame(igs_data).groupby(['Group', 'Species', 'File'])['IGS_Size'].mean().reset_index()
    igs_df = igs_df.rename(columns={'IGS_Size': 'Mean_IGS_Size'})
    print('Done with IGS analysis.')

    gene_df = pd.DataFrame(gene_data)
    gene_count_df = gene_df.groupby(['Group', 'Species', 'File']).size().reset_index(name='Gene_Count')
    gene_length_df = gene_df.groupby(['Group', 'Species', 'File'])['Gene_Length'].mean().reset_index(name='mean_gene_length')
    gene_df = pd.merge(gene_count_df, gene_length_df, on=['Group', 'Species', 'File'])

    print('Done with gene length and count analysis.')

    codon_stats_df = collect_codon_stats(path, auto_classify=classify_flag)
    print('Done with GC metrics.')

    merged_df = pd.merge(transposase_df, igs_df, on=['Group', 'Species', 'File'], how='inner')
    merged_df = pd.merge(merged_df, gene_df, on=['Group', 'Species', 'File'], how='inner')
    merged_df = pd.merge(merged_df, codon_stats_df, on=['Group', 'Species', 'File'], how='inner')

    merged_df = merged_df[merged_df['Gene_Count'] > 0]
    merged_df['Transposase_Per_Gene'] = merged_df['Total_Transposases'] / merged_df['Gene_Count']
    merged_df['Delta_GC2_4'] = merged_df['GC2'] - merged_df['GC4']

    keep_cols = [
        'Group', 
        'Species', 
        'File',   
        'GC4',
        'Delta_GC2_4', 
        'AV_Bias', 
        'Rest_Bias',
        'Transposase_Per_Gene', 
        'Mean_IGS_Size'
    ]
    merged_df = merged_df[keep_cols]

    transposase_files = set(transposase_df['File'].unique())
    final_files = set(merged_df['File'].unique())
    dropped_files = transposase_files - final_files
    if dropped_files:
        print(f'Warning: {len(dropped_files)} files were dropped due to missing data.')
    else:
        print('All files successfully merged.')

    os.makedirs(feature_dir, exist_ok=True)
    merged_df.to_csv(output_csv, index=False)
    print('All features saved to combined_features.csv')

if __name__ == '__main__':
    main()