#!/usr/bin/env python3

import os
import pandas as pd

from igs_lengths import gff_features
from Bio import SeqIO
from Bio.Seq import Seq
from gc_calculate import calculate_gc_content
from fourfold_bias import compute_fourfold_bias
from utils import fourfold_codons
from joblib import Parallel, delayed


def parse_fna(file_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(file_path, 'fasta')}

def compute_gc_metrics(coordinates, seq_dict):

    all_second = []
    all_third = []
    gene_seqs = []

    for seqid, start, end, strand in coordinates:
        if seqid not in seq_dict:
            continue
        gene = seq_dict[seqid][start - 1 : end]
        if strand == '-':
            gene = str(Seq(gene).reverse_complement())

        gene_seqs.append(gene)

        for i in range(0, len(gene) - 2, 3):
            codon = gene[i:i+3]
            if len(codon) == 3:
                all_second.append(codon[1])
                if codon in fourfold_codons:
                    all_third.append(codon[2])

    aa_gc4_dict = compute_fourfold_bias(''.join(gene_seqs))
    gc2 = calculate_gc_content(''.join(all_second))
    gc4 = calculate_gc_content(''.join(all_third))

    return gc2, gc4, aa_gc4_dict


def _process_single_genome(gff_path, sp_name, auto_classify, default_group):

    filename = os.path.basename(gff_path)
    clean_filename = filename.replace('.gff', '')
    fna_path = gff_path.replace('.gff', '.fna')

    if not os.path.exists(fna_path):
        print(f'Warning: .fna file not found for {gff_path}.')
        return None

    if auto_classify:
        current_group = 'relatives_only' if '_genomic' in filename else 'endosymb_only'
    else:
        current_group = default_group

    print(f'Processing {gff_path}')
    _, _, coordinates = gff_features(gff_path)
    seq_dict = parse_fna(fna_path)
    gc2, gc4, aa_gc4_dict = compute_gc_metrics(coordinates, seq_dict)

    AV_bias   = sum(aa_gc4_dict[aa]['bias'] for aa in ['Val', 'Ala'] if aa in aa_gc4_dict) / 2
    rest_bias = sum(aa_gc4_dict[aa]['bias'] for aa in aa_gc4_dict if aa not in ['Val', 'Ala']) / 6

    return {
        'Group':   current_group,
        'Species': sp_name,
        'File':    clean_filename,
        'GC2':     gc2,
        'GC4':     gc4,
        'AV_Bias': AV_bias,
        'Rest_Bias': rest_bias,
    }


def collect_codon_stats(path, group_label=None, auto_classify=False, n_jobs=-1):

    default_group = group_label if group_label else 'Unknown'
    target_dir = os.path.join(path, 'genomes')

    if not os.path.exists(target_dir):
        print(f'Warning: Genomes directory not found at {target_dir}. Skipping codon stats collection.')
        target_dir = path

    tasks = []
    for roots, dirs, files in os.walk(target_dir):
        for filename in files:
            if filename.endswith('.gff'):
                rel_path = os.path.relpath(roots, target_dir)
                sp_name = 'Unknown' if rel_path == '.' else rel_path.replace('_endosymbiont', '').replace('_', ' ')
                tasks.append((os.path.join(roots, filename), sp_name))

    if not tasks:
        return pd.DataFrame()

    print(f"  Found {len(tasks)} GFF files for codon stats — running with n_jobs={n_jobs}")

    results = Parallel(n_jobs=n_jobs, backend='loky')(
        delayed(_process_single_genome)(gff_path, sp_name, auto_classify, default_group)
        for gff_path, sp_name in tasks
    )

    data = [r for r in results if r is not None]
    return pd.DataFrame(data)
