import os
import re
import subprocess
import sys
import pandas as pd

from diamond_job import run_diamond
from analysis.gcsize_dict import genome_gcsize
from analysis.transposase_analysis import processing_transposase
from analysis.igs_lengths import collect_gff_stats
from analysis.utils import files_dir

if __name__ == '__main__':
    # Clustering with labeled groups (use auto_classify to determine group based on filename)
    path = 'endosymb+relatives'

    frames = []
    gcsize_data = genome_gcsize(path)
    for species, accessions in gcsize_data.items():
        for accn, metadata in accessions.items():
            frames.append({
                'Species': species.replace('_endosymbiont', '').replace('_', ' '),
                'File': accn,
                'GC_Content': metadata['gc_genome'],
                'Genome_Size': metadata['size']
            })
    gcsize_df = pd.DataFrame(frames)

    # Structure must be path/proteins/species/*.faa
    run_diamond(path, threads=8, wait=True, max_parallel=20)
    print('Done with DIAMOND searches of transposases.')
    
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