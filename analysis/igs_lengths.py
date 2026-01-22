#!/usr/bin/env python3

import os
import gffutils
import pandas as pd
from utils import files_dir

def gff_processing(file_path):
    db = gffutils.create_db(file_path, dbfn=':memory:', force = True, 
                            keep_order = True, merge_strategy = 'merge',
                            sort_attribute_values = True)
    
    IGS_sizes = []
    gene_lengths = []
    IGS_start = None

    for feature in db.all_features(order_by = ('seqid', 'start', 'end')):
        start = feature.start
        end = feature.end

        # Skip region features (they span entire sequences)
        if feature.featuretype == 'region':
            continue

        # Record gene lengths
        if feature.featuretype == 'gene':
            gene_length = end - start + 1
            gene_lengths.append(gene_length)

        # Calculate IGS sizes
        if IGS_start:
            IGS_size = start - IGS_start - 1
            if IGS_size > 0:
                IGS_sizes.append(IGS_size)
        IGS_start = end

    return IGS_sizes, gene_lengths

all_igs_data = []
all_gene_data = []

group_names = {'endosymb_only': 'Endosymbionts Only', 
               'endosymb+relatives': 'Endosymbionts and Free-Living Relatives', 
               'relatives_only': 'Free-Living Relatives Only'}
print('Processing GFF files.')

for group in group_names.keys():
    gff_folder = os.path.join(group, 'genomes')
    species = os.listdir(gff_folder)
    for sp in species:
        sp_path  = os.path.join(gff_folder, sp)
        spname = sp.replace('_endosymbiont','').replace('_', ' ')
        for file in os.listdir(sp_path):
            if file.endswith('gff'):
                file_path = os.path.join(sp_path, file)
                IGS_sizes, gene_lengths = gff_processing(file_path)
                for size in IGS_sizes:
                    all_igs_data.append({'group': group_names[group],
                                     'species': sp, 'IGS_Size': size, 
                                     'file': file})
                for length in gene_lengths:
                    all_gene_data.append({'group': group_names[group],
                                     'species': sp, 'Gene_Length': length,
                                     'file': file})

if all_igs_data:                    
    #Show overall median IGS per species, for each group
    df = pd.DataFrame(all_igs_data)
    summary_df = df.groupby(['group', 'species', 'file'])['IGS_Size'].median().reset_index()
    summary_df = summary_df.groupby(['group', 'species']).agg({
        'IGS_Size': 'mean',
        'file': 'count'
    }).rename(columns = {'IGS_Size': 'mean_median_IGS', 'file': 'num_genomes'}).reset_index()

    summary_df.to_csv(os.path.join(files_dir, 'medianIGS.csv'), index = False)
    df.to_csv(os.path.join(files_dir, 'all_IGS_data.csv'), index = False)

if all_gene_data:
    #Show average Gene Length per species, for each group
    df_gene = pd.DataFrame(all_gene_data)
    summary_gene_df = df_gene.groupby(['group', 'species']).agg({
        'Gene_Length': 'mean',
        'file': 'count'
    }).rename(columns = {'Gene_Length': 'mean_gene_length', 'file': 'num_genomes'}).reset_index()

    summary_gene_df.to_csv(os.path.join(files_dir, 'mean_gene_lengths.csv'), index = False)