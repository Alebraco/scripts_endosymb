#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np

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

            for file in os.listdir(species_path):
                file_path = os.path.join(species_path, file)

                if os.path.getsize(file_path) > 0:
                    df = pd.read_csv(file_path, sep='\t', names=columns, usecols=use_columns)
                    df['Group'] = group
                    df['Species'] = sp_name
                    df['Accession'] = file.replace('.tsv','')
                    df['Total Genomes'] = len(os.listdir(species_path))

                    df = df[df['pident'] >= identity_threshold]

                    df = df[df['sseqid'].str.contains('Transposase', case=False)]
                    df = df[~df['sseqid'].str.contains('Accessory', case=False)]

                    df['coverage'] = (df['send'] - df['sstart'] + 1) / df['slen']

                    df['status'] = np.where(df['coverage'] >= coverage_threshold, 'complete', 'partial')

                    df['IS_Family'] = df['sseqid'].str.split('_').str[0]

                    results.append(df)

    if results:
        final_df = pd.concat(results, ignore_index=True)
        return final_df
    else:
        return None
    
if __name__ == "__main__":
    df_master = processing_transposase()
    # if df_master is not None:
