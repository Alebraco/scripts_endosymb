import os
import subprocess
import sys
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from gcsize_dict import genome_gcsize
from transposase_analysis import processing_transposase
from igs_lengths import collect_gff_stats
from utils import files_dir

# Structure must be path/proteins/species/*.faa

def run_diamond(path):
    # The path must be the directory containing the 'proteins' folder
    TRANS_DB = os.path.join(files_dir, 'IS_db_fixed.dmnd')

    for species in os.listdir(path):
        species_path = os.path.join(path, 'proteins', species)
        if not os.path.isdir(species_path):
            continue

        out_dir = os.path.join(path, 'transposase', species)
        os.makedirs(out_dir, exist_ok=True)

        for file in os.listdir(species_path):
            if file.endswith('.faa'):
                accession = file.split('.faa')[0]
                faa_file = os.path.join(species_path, file)

                outfile = os.path.join(out_dir, f'{accession}.tsv')

                cmd = [
                    "diamond", "blastp",
                    "--db", TRANS_DB,
                    "--query", faa_file,
                    "--out", outfile,
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", 
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
                    "--threads", "8",
                    "--evalue", "1e-5",
                    "--max-target-seqs", "1"
                ]
    
                subprocess.run(cmd, check=True)



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

    run_diamond(path)
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






