#!/usr/bin/env python3

from gc_calculate import calculate_gc_content
from utils import genome_gcsize_json_path
from Bio import SeqIO
import json
import os

groups = ['endosymb_only', 'endosymb+relatives', 'relatives_only']

def fetch_gc_size(sequence):
    genome_size = len(sequence)
    gc_content = calculate_gc_content(sequence)
    return genome_size, gc_content

def genome_gcsize(group):
    '''
    Process genomes in input_dir and compute genome-wide size + GC content.

    Args:
        group (str): endosymb_only, endosymb+relatives, relatives_only.

    Returns:
        dict: Nested dictionary with species → accession → metadata.
    '''

    genome_dir = os.path.join(group, 'genomes')

    all_data = {}

    for sp in os.listdir(genome_dir):
        print(f'Processing {sp}')
        all_data[sp] = {}

        genome_path = os.path.join(genome_dir, sp)
        genome_list = [g for g in os.listdir(genome_path) if g.endswith('.fna')]

        for file in genome_list:
            accn = file.split('.fna')[0]
            file_path = os.path.join(genome_path, file)

            seqs = []
            for rec in SeqIO.parse(file_path, 'fasta'):
                seqs.append(str(rec.seq))

            
            final_seq = ''.join(seqs)
            genome_size, gc_content = fetch_gc_size(final_seq)

            all_data[sp][accn] = {
                'gc_genome': gc_content,
                'size': genome_size
            }
    out_path = genome_gcsize_json_path(group)
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(all_data, f, indent=4)

    return all_data

if __name__ == '__main__':
    for group in groups:
        genome_gcsize(group)