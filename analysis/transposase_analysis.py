#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from utils import files_dir

coverage_threshold = 0.8
identity_threshold = 30.0
columns = ['qseqid', 
           'sseqid', 
           'pident', 
           'length', 
           'mismatch', 
           'gapopen', 
           'qstart', 
           'qend', 
           'sstart', 
           'send', 
           'evalue', 
           'bitscore', 
           'qlen', 
           'slen']
use_columns = ['sseqid', 'pident', 'sstart', 'send', 'slen']

def processing_transposase():
    results = []

    for group in ['endosymb_only', 'relatives_only']:
        group_path = os.path.join(group, 'transposase')

        for species in os.listdir(group_path):
            sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
            species_path = os.path.join(group_path, species)

            tsv_files = [file for file in os.listdir(species_path) if file.endswith('.tsv')]
            total_genomes = len(tsv_files)

            for file in tsv_files:
                file_path = os.path.join(species_path, file)

                complete, partial, total, families = 0,0,0,0
                fam_list = None

                if os.path.getsize(file_path) > 0:
                    df = pd.read_csv(file_path, sep='\t', names=columns, usecols=use_columns)

                    # Filter hits
                    df = df[df['pident'] >= identity_threshold]
                    df = df[df['sseqid'].str.contains('Transposase', case=False)]
                    df = df[~df['sseqid'].str.contains('Accessory', case=False)]

                    if not df.empty:

                        df['coverage'] = (df['send'] - df['sstart'] + 1) / df['slen']
                        complete = len(df[df['coverage'] >= coverage_threshold])
                        partial = len(df[df['coverage'] < coverage_threshold])
                        total = complete + partial

                        # Extract IS families
                        fam_list = df['sseqid'].str.split('_').str[0].unique().tolist()
                        families = len(fam_list)

                results.append({
                    'Group': group,
                    'Species': sp_name,
                    'File': file.replace('.tsv',''),
                    'Total_Genomes': total_genomes,
                    'Complete_Transposases': complete,
                    'Partial_Transposases': partial,
                    'Total_Transposases': total,
                    'IS_Families': fam_list,
                    'Families_Count': families
                })

    if results:
        return pd.DataFrame(results)
    else:
        return None
    
if __name__ == "__main__":
    df_master = processing_transposase()
    if df_master is not None:
        print(df_master.head())
        df_master.to_csv(os.path.join(files_dir, 'transposase_summary.csv'), index=False)
