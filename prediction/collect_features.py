#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

from diamond_job import run_diamond
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from gcsize_dict import genome_gcsize
from transposase_analysis import processing_transposase
from igs_lengths import collect_gff_stats
from utils import files_dir, load_or_compute, genome_gcsize_json_path

def main():
    parser = argparse.ArgumentParser(description='Collect and merge genome features.')
    parser.add_argument('--path', help='Path to process')
    args = parser.parse_args()

    path = args.path

    print('Starting feature collection process.')
    frames = []
    gcsize_data = load_or_compute(genome_gcsize_json_path(path), genome_gcsize, path)
    for species, accessions in gcsize_data.items():
        for accn, metadata in accessions.items():
            frames.append({
                'Species': species.replace('_endosymbiont', '').replace('_', ' '),
                'File': accn,
                'GC_Content': metadata['gc_genome'],
                'Genome_Size': metadata['size']
            })
    gcsize_df = pd.DataFrame(frames)
    print('Done with GC content and genome size analysis.')
    
    transposase_df = processing_transposase(path, auto_classify=True)
    print('Done with transposase analysis.')

    igs_data, gene_data = collect_gff_stats(path, auto_classify=True)
    igs_df = pd.DataFrame(igs_data).groupby(['File'])['IGS_Size'].mean().reset_index()
    igs_df = igs_df.rename(columns={'IGS_Size': 'Mean_IGS_Size'})
    print('Done with IGS analysis.')

    merged_df = pd.merge(gcsize_df, transposase_df, on=['File', 'Species'], how='inner')
    merged_df = pd.merge(merged_df, igs_df, on='File', how='inner')

    merged_df.to_csv(os.path.join(files_dir, 'combined_features.csv'), index=False)
    print('All features saved to combined_features.csv')


if __name__ == '__main__':
    main()