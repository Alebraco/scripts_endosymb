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
    current_seqid = None

    target_features = ['gene', 'CDS', 'rRNA']

    for feature in db.all_features(order_by = ('seqid', 'start', 'end')):
        if feature.featuretype not in target_features:
            continue
        
        if current_seqid != feature.seqid:
            current_seqid = feature.seqid
            IGS_start = None
            
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

            # Calculate IGS sizes if not pseudo gene
            if IGS_start:
                IGS_size = start - IGS_start - 1
                if IGS_size > 0:
                    IGS_sizes.append(IGS_size)

            if IGS_start is None or end > IGS_start:
                IGS_start = end

    final_lengths = gene_lengths if len(gene_lengths) > 0 else cds_lengths

    return IGS_sizes, final_lengths

def collect_gff_stats(path, group_label = None, auto_classify = False):
    '''
    Process a directory of genomes in GFF format
    Extracts IGS sizes and gene lengths.
    '''
    all_igs_data = []
    all_gene_data = []

    default_group = group_label if group_label else 'Ungrouped'

    target_dir = os.path.join(path, 'genomes')
    if not os.path.exists(target_dir):
        print(f"Directory not found: {target_dir}")
        return [], []
    print(f"Processing GFF files in: {target_dir}")

    for sp in os.listdir(target_dir):
        sp_path = os.path.join(target_dir, sp)
        
        if not os.path.isdir(sp_path):
            continue
        
        spname = sp.replace('_endosymbiont','').replace('_', ' ')

        for file in os.listdir(sp_path):
            if not file.lower().endswith('.gff'):
                continue
            
            filename = file.replace('.gff', '')
            file_path = os.path.join(sp_path, file)

            IGS_sizes, gene_lengths = gff_processing(file_path)

            if auto_classify:
                if '_genomic' in file:
                    current_group = 'relatives_only'
                else:
                    current_group = 'endosymb_only'
            else:
                current_group = default_group

            for size in IGS_sizes:
                all_igs_data.append({
                    'Group': current_group,
                    'Species': spname, 
                    'IGS_Size': size, 
                    'File': filename
                })
            for length in gene_lengths:
                all_gene_data.append({
                    'Group': current_group,
                    'Species': spname, 
                    'Gene_Length': length,
                    'File': filename
                })
    return all_igs_data, all_gene_data

def save_summary_stats(all_igs_data, all_gene_data, prefix=''):
    '''
    Saves summary statistics for IGS sizes and gene lengths to CSV files.
    '''

    if all_igs_data:                    
        df = pd.DataFrame(all_igs_data)

        summary_df = df.groupby(['Group', 'Species', 'File'])['IGS_Size'].mean().reset_index()

        summary_df = summary_df.groupby(['Group', 'Species']).agg({
            'IGS_Size': 'mean',
            'File': 'count'
        }).rename(columns = {'IGS_Size': 'mean_mean_IGS', 'File': 'num_genomes'}).reset_index()

        outfile = f'{prefix}meanIGS.csv' if prefix else 'meanIGS.csv'
        summary_df.to_csv(os.path.join(files_dir, outfile), index = False)

        raw_name = f'{prefix}all_IGS_data.csv' if prefix else 'all_IGS_data.csv'
        df.to_csv(os.path.join(files_dir, raw_name), index = False)

    if all_gene_data:
        df_gene = pd.DataFrame(all_gene_data)

        gene_counts = df_gene.groupby(['Group', 'Species', 'File']).size().reset_index(name='Gene_Count')
        count_name = f'{prefix}gene_counts.csv' if prefix else 'gene_counts.csv'
        gene_counts.to_csv(os.path.join(files_dir, count_name), index = False)

        summary_gene_df = df_gene.groupby(['Group', 'Species', 'File']).agg({
            'Gene_Length': 'mean'}).reset_index()

        summary_gene_df = summary_gene_df.groupby(['Group', 'Species']).agg({
            'Gene_Length': 'mean',
            'File': 'count'
        }).rename(columns = {'Gene_Length': 'mean_gene_length', 'File': 'num_genomes'}).reset_index()

        outname = f'{prefix}mean_gene_lengths.csv' if prefix else 'mean_gene_lengths.csv'
        summary_gene_df.to_csv(os.path.join(files_dir, outname), index = False)

if __name__ == "__main__":

    master_igs = []
    master_gene = []

    csv_names = ['meanIGS.csv', 'all_IGS_data.csv', 'mean_gene_lengths.csv']
    files = all(os.path.exists(os.path.join(files_dir, f)) for f in csv_names)

    if not files:
        print('IGS and Gene Length CSV files not found. Generating data...')
        for group in ['endosymb_only', 'relatives_only']:
            if os.path.exists(group):
                print(f"Processing group: {group}")
                igs_data, gene_data = collect_gff_stats(group, group_label=group)
                master_igs.extend(igs_data)
                master_gene.extend(gene_data)
        save_summary_stats(master_igs, master_gene)
        print('Done.')

    else:
        print('IGS and Gene Length CSV files found.')