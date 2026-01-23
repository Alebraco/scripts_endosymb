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
    cds_lengths = []
    IGS_start = None

    target_features = ['gene', 'CDS']

    for feature in db.all_features(order_by = ('seqid', 'start', 'end')):
        if feature.featuretype not in target_features:
            continue
        
        # Determine pseudogene status
        pseudo = False
        if 'pseudo' in feature.attributes:
            if feature.attributes['pseudo'][0].lower() == 'true':
                pseudo = True
        
        start = feature.start
        end = feature.end
        length = end - start + 1

        # Only consider true genes
        if not pseudo:
            # Record CDS lengths
            if feature.featuretype == 'CDS':
                cds_lengths.append(length)
            # Record gene lengths
            if feature.featuretype == 'gene':
                gene_lengths.append(length)

        # Calculate IGS sizes
        if IGS_start:
            IGS_size = start - IGS_start - 1
            if IGS_size > 0:
                IGS_sizes.append(IGS_size)

        if IGS_start is None or end > IGS_start:
            IGS_start = end

    final_lengths = gene_lengths if len(gene_lengths) > 0 else cds_lengths

    return IGS_sizes, final_lengths

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
                                     'species': spname, 'IGS_Size': size, 
                                     'file': file})
                for length in gene_lengths:
                    all_gene_data.append({'group': group_names[group],
                                     'species': spname, 'Gene_Length': length,
                                     'file': file})

if all_igs_data:                    
    df = pd.DataFrame(all_igs_data)
    # Median per file
    summary_df = df.groupby(['group', 'species', 'file'])['IGS_Size'].median().reset_index()
    # Mean of Medians per Species
    summary_df = summary_df.groupby(['group', 'species']).agg({
        'IGS_Size': 'mean',
        'file': 'count'
    }).rename(columns = {'IGS_Size': 'mean_median_IGS', 'file': 'num_genomes'}).reset_index()

    summary_df.to_csv(os.path.join(files_dir, 'medianIGS.csv'), index = False)
    df.to_csv(os.path.join(files_dir, 'all_IGS_data.csv'), index = False)

if all_gene_data:
    df_gene = pd.DataFrame(all_gene_data)
    # Mean per file
    summary_gene_df = df_gene.groupby(['group', 'species', 'file']).agg({
        'Gene_Length': 'mean'}).reset_index()
    # Mean of Mean per Species
    summary_gene_df = summary_gene_df.groupby(['group', 'species']).agg({
        'Gene_Length': 'mean',
        'file': 'count'
    }).rename(columns = {'Gene_Length': 'mean_gene_length', 'file': 'num_genomes'}).reset_index()

    summary_gene_df.to_csv(os.path.join(files_dir, 'mean_gene_lengths.csv'), index = False)