#!/usr/bin/env python3

from gc_calculate import calculate_gc_content
from utils import gc_codon_json_path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import os
import json

groups = ['endosymb_only', 'endosymb+relatives', 'relatives_only']

def gc_codon_dict(group, save_to_file = True):
    '''Create a dictionary with the GC content
    Produces fourfold degenerate and overall GC content 
    Computes GC for each core gene alignment concatenate
    args:
        group (str): endosymb_only, endosymb+relatives, relatives_only
        save_to_file (boolean): Saves dictionary as json file
    '''
    fourfold_codons = {
        'GCT', 'GCC', 'GCA', 'GCG',   # Ala
        'CGT', 'CGC', 'CGA', 'CGG',   # Arg
        'GGT', 'GGC', 'GGA', 'GGG',   # Gly
        'CCT', 'CCC', 'CCA', 'CCG',   # Pro
        'ACT', 'ACC', 'ACA', 'ACG',   # Thr
        'GTT', 'GTC', 'GTA', 'GTG',   # Val
        'CTT', 'CTC', 'CTA', 'CTG',   # Leu
        'TCT', 'TCC', 'TCA', 'TCG'    # Ser
    }
    gc_data = {}

    output_dirs = {
        'first': os.path.join(group, 'first_sites'),
        'second': os.path.join(group, 'second_sites'),
        'third': os.path.join(group, 'third_sites'),
    }
    for path in output_dirs.values():
        os.makedirs(path, exist_ok=True)

    input_dir = os.path.join(group, 'dna_concatenates')

    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        species_name = file.split('concatenate_')[1].split('.fasta')[0]
        out_paths = {site: os.path.join(output_dirs[site], file) for site in output_dirs}
        
        print(f'Processing {species_name}')

        records_lists = {site: [] for site in output_dirs}
        alignment_gc = {}

        for rec in SeqIO.parse(path, 'fasta'):
            sites = {site: [] for site in output_dirs}

            # Iterate over codon positions 
            for i in range(0, len(rec.seq), 3): 
                codon = str(rec.seq[i:i+3]) # Retrieve each codon sequence
                for site in output_dirs:
                    sites['first'].append(codon[0])
                    sites['second'].append(codon[1])
                    sites['third'].append(codon[2] if codon in fourfold_codons else '-')

            for site in output_dirs:
                seq_str = ''.join(sites[site]) # Concatenate the nucleotides for each record (ID)
                new_rec = SeqRecord(Seq(seq_str), id = rec.id) # Convert to SeqRecord object
                records_lists[site].append(new_rec) # Append SeqRecord to records list

            # Calculate GC content 
            gc_third = calculate_gc_content(Seq(''.join(sites['third'])))
            gc_all = calculate_gc_content(rec.seq)

            alignment_gc[rec.id] = {
                'gc_third': gc_third,
                'gc_all' : gc_all,
            }
        gc_data[species_name] = alignment_gc

        # Write all records at once, for this species
        for site in output_dirs:
            with open(out_paths[site], 'w') as output: 
                SeqIO.write(records_lists[site], output, 'fasta')

    if save_to_file:
        out_json = gc_codon_json_path(group)
        os.makedirs(os.path.dirname(out_json) or ".", exist_ok=True)
        with open(out_json, 'w') as save:
            json.dump(gc_data, save)
    return gc_data

if __name__ == '__main__':
    for group in groups:
        gc_codon_dict(group)