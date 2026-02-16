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

def genome_gcsize(path):
    '''
    Process genomes in input_dir and compute genome-wide size + GC content.

    Args:
        path (str): Path to the directory containing genome groups.

    Returns:
        dict: Nested dictionary with species → accession → metadata.
    '''

    genome_dir = os.path.join(path, 'genomes')
    group_name = os.path.basename(path)

    all_data = {}
    print()

    for sp in os.listdir(genome_dir):
        print(f'Processing {sp}')
        all_data[sp] = {}

        species_path = os.path.join(genome_dir, sp)
        genome_list = [g for g in os.listdir(species_path) if g.endswith('.fna')]

        for file in genome_list:
            accn = file.split('.fna')[0]
            file_path = os.path.join(species_path, file)

            seqs = []
            try:
                for rec in SeqIO.parse(file_path, 'fasta'):
                    seqs.append(str(rec.seq))
                
                final_seq = ''.join(seqs)
                if final_seq:
                    genome_size, gc_content = fetch_gc_size(final_seq)
                    all_data[sp][accn] = {
                        'gc_genome': gc_content,
                        'size': genome_size
                    }
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue

    if all_data and group_name in groups:
        try:
            out_path = genome_gcsize_json_path(path)
            with open(out_path, 'w') as f:
                json.dump(all_data, f, indent=4)
            print(f"Data saved to {out_path}")
        except Exception as e:
            print(f"Error saving data to {out_path}: {e}")

    return all_data

if __name__ == '__main__':
    for group in groups:
        if os.path.exists(group):
            genome_gcsize(group)