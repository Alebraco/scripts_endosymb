#!/usr/bin/env python3

import os
import gffutils
import pandas as pd
from utils import files_dir

def gff_processing(file_path):
    db = gffutils.create_db(file_path, dbfn=':memory:', force = True, keep_order = True, merge_strategy = 'merge', sort_attribute_values = True)
    IGS_sizes = []
    IGS_start = None
    feature_count = 0

    for feature in db.all_features(order_by = ('seqid', 'start', 'end')):
        start = feature.start
        end = feature.end
        feature_count += 1

        if feature.featuretype == 'region':
            continue

        if IGS_start:
            IGS_size = start - IGS_start - 1
            if IGS_size > 0:
                IGS_sizes.append(IGS_size)
        IGS_start = end

    return IGS_sizes

all_data = []
group_names = {'endosymb_only': 'Endosymbionts Only', 
               'endosymb+relatives': 'Endosymbionts and Free-Living Relatives', 
               'relatives_only': 'Free-Living Relatives Only'}

for group in group_names.keys():
    gff_folder = os.path.join(group, 'genomes')
    species = os.listdir(gff_folder)
    for sp in species:
        IGS_lengths = []
        sp_path  = os.path.join(gff_folder, sp)
        spname = sp.replace('_endosymbiont','').replace('_', ' ')
        for file in os.listdir(sp_path):
            if file.endswith('gff'):
                file_path = os.path.join(sp_path, file)
                IGS_sizes = gff_processing(file_path)
                for size in IGS_sizes:
                    all_data.append({'group': group_names[group],
                                     'species': sp, 'IGS_Size': size, 
                                     'file': file})
#Show overall median IGS per species, for each group
df = pd.DataFrame(all_data)
summary_df = df.groupby(['group', 'species', 'file'])['IGS_Size'].median().reset_index()
summary_df = summary_df.groupby(['group', 'species']).agg({
    'IGS_Size': 'mean',
    'file': 'count'
}).rename(columns = {'IGS_Size': 'mean_median_IGS', 'file': 'num_genomes'}).reset_index()

summary_df.to_csv(os.path.join(files_dir, 'medianIGS.csv'), index = False)
df.to_csv(os.path.join(files_dir, 'all_IGS_data.csv'), index = False)