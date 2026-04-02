#!/usr/bin/env python3

import os
import pandas as pd
from joblib import Parallel, delayed
from utils import files_dir

TARGET_TYPES = {'gene', 'CDS', 'rRNA'}

def gff_features(file_path):

    seqs = {}

    with open(file_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in TARGET_TYPES:
                continue

            # Parse attributes column for pseudo=true
            pseudo = False
            for attr in parts[8].split(';'):
                if attr.strip().lower().startswith('pseudo='):
                    pseudo = attr.split('=', 1)[1].strip().lower() == 'true'
                    break
            if pseudo:
                continue

            seqid  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]

            if seqid not in seqs:
                seqs[seqid] = []
            seqs[seqid].append((start, end, ftype, strand))

    IGS_sizes   = []
    gene_lengths = []
    cds_lengths  = []
    gene_coords  = []
    cds_coords   = []

    for seqid, features in seqs.items():
        features.sort(key=lambda x: (x[0], x[1]))  # sort by start, then end

        igs_frontier = None  # rightmost end position seen so far on this contig
        for start, end, ftype, strand in features:
            length = end - start + 1

            if ftype == 'CDS':
                cds_lengths.append(length)
                cds_coords.append((seqid, start, end, strand))
            elif ftype == 'gene':
                gene_lengths.append(length)
                gene_coords.append((seqid, start, end, strand))

            if igs_frontier is not None:
                IGS_sizes.append(max(start - igs_frontier - 1, 0))

            if igs_frontier is None or end > igs_frontier:
                igs_frontier = end

    final_lengths = gene_lengths if gene_lengths else cds_lengths
    final_coords  = gene_coords  if gene_coords  else cds_coords

    return IGS_sizes, final_lengths, final_coords


def _process_single_gff(file_path, sp_name, auto_classify, default_group):

    filename = os.path.basename(file_path).replace('.gff', '')

    IGS_sizes, gene_lengths, _ = gff_features(file_path)

    if auto_classify:
        current_group = 'relatives_only' if '_genomic' in filename else 'endosymb_only'
    else:
        current_group = default_group

    igs_rows = [
        {'Group': current_group, 'Species': sp_name, 'IGS_Size': size, 'File': filename}
        for size in IGS_sizes
    ]
    gene_rows = [
        {'Group': current_group, 'Species': sp_name, 'Gene_Length': length, 'File': filename}
        for length in gene_lengths
    ]
    return igs_rows, gene_rows


def collect_gff_stats(path, group_label=None, auto_classify=False, n_jobs=-1):

    default_group = group_label if group_label else 'Unknown'

    target_dir = os.path.join(path, 'genomes')
    if not os.path.exists(target_dir):
        print(f"Directory not found: {target_dir}")
        return [], []
    print(f"Processing GFF files in: {target_dir}")

    gff_tasks = []

    species_dirs = [
        entry for entry in os.listdir(target_dir)
        if os.path.isdir(os.path.join(target_dir, entry))
    ]
    for sp in species_dirs:
        sp_path = os.path.join(target_dir, sp)
        sp_name = sp.replace('_endosymbiont', '').replace('_', ' ')
        for fname in os.listdir(sp_path):
            if fname.lower().endswith('.gff'):
                gff_tasks.append((os.path.join(sp_path, fname), sp_name))

    flat_gff_files = [f for f in os.listdir(target_dir) if f.lower().endswith('.gff')]
    for fname in flat_gff_files:
        gff_tasks.append((os.path.join(target_dir, fname), 'Unknown'))

    if not gff_tasks:
        return [], []

    print(f"  Found {len(gff_tasks)} GFF files to process.")

    results = Parallel(n_jobs=n_jobs, backend='loky')(
        delayed(_process_single_gff)(fp, sp, auto_classify, default_group)
        for fp, sp in gff_tasks
    )

    all_igs_data  = [row for igs_rows, _         in results for row in igs_rows]
    all_gene_data = [row for _,        gene_rows  in results for row in gene_rows]

    return all_igs_data, all_gene_data


def save_summary_stats(all_igs_data, all_gene_data, prefix=''):
    '''
    Saves summary statistics for IGS sizes and gene lengths to CSV files.
    '''

    if all_igs_data:
        df = pd.DataFrame(all_igs_data)

        genome_stats = df.groupby(['Group', 'Species', 'File'])['IGS_Size'].agg(
            mean_IGS='mean', std_IGS='std'
        ).reset_index()
        std_name = f'{prefix}genome_std_IGS.csv' if prefix else 'genome_std_IGS.csv'
        genome_stats.to_csv(os.path.join(files_dir, std_name), index=False)

        summary_df = genome_stats.groupby(['Group', 'Species']).agg(
            mean_mean_IGS=('mean_IGS', 'mean'),
            num_genomes=('File', 'count')
        ).reset_index()

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

    csv_names = ['meanIGS.csv', 'all_IGS_data.csv', 'genome_std_IGS.csv', 'mean_gene_lengths.csv', 'gene_counts.csv']
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
