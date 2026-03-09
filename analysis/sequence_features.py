#!/usr/bin/env python3

import pandas as pd

from igs_lengths import gff_features
from Bio import SeqIO
from Bio.Seq import Seq
from gc_calculate import calculate_gc_content
from fourfold_bias import compute_fourfold_bias
from utils import fourfold_codons
import os

def parse_fna(file_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(file_path, 'fasta')}

def compute_gc_metrics(coordinates, seq_dict):
    # Works on one file at a time, one seq_dict entry
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

def collect_codon_stats(path, group_label = None, auto_classify = False):
    
    data = []
    for roots, dirs, files in os.walk(path):
        for filename in files:
            if filename.endswith('.gff'):

                if auto_classify:
                    current_group = 'relatives_only' if '_genomic' in filename else 'endosymb_only'
                else:
                    current_group = group_label if group_label else 'Ungrouped'

                file_path = os.path.join(roots, filename)
                _, _, coordinates = gff_features(file_path)
                fna_path = file_path.replace('.gff', '.fna')
                clean_filename = filename.replace('.gff', '')
                if not os.path.exists(fna_path):
                    print(f'Warning: .fna file not found for {file_path}.')
                    continue
                else:
                    print(f'Processing {file_path} and corresponding {fna_path}')
                    seq_dict = parse_fna(fna_path)
                    gc2, gc4, aa_gc4_dict = compute_gc_metrics(coordinates, seq_dict)

                    AV_bias = sum(aa_gc4_dict[aa]['bias'] for aa in ['Val', 'Ala'] if aa in aa_gc4_dict) / 2
                    rest_bias = sum(aa_gc4_dict[aa]['bias'] for aa in aa_gc4_dict if aa not in ['Val', 'Ala']) / 6

                    print(f'File: {filename}, GC2: {gc2:.4f}, GC4: {gc4:.4f}')
                    data.append({
                        'Group': current_group,
                        'File': clean_filename, 
                        'GC2': gc2, 
                        'GC4': gc4, 
                        'AV_Bias': AV_bias, 
                        'Rest_Bias': rest_bias
                        })

    return pd.DataFrame(data)
                         
