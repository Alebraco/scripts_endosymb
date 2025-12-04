#!/usr/bin/env python3

from gc_calculate import calculate_gc_content
from gc_utils import gc_codon_json_path
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
    output_dir = os.path.join(group, 'third_sites')
    input_dir = os.path.join(group, 'dna_concatenates')

    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        species_name = file.split('concatenate_')[1].split('.fasta')[0]
        out_path = os.path.join(output_dir, file)
        print(f'Processing {species_name}')

        records_list = []
        alignment_gc = {}

        for rec in SeqIO.parse(path, 'fasta'):
            third_sites_degen = []

            # Iterate over codon positions 
            for i in range(0, len(rec.seq), 3): 
                codon = str(rec.seq[i:i+3]) # Retrieve each codon sequence
                if codon in fourfold_codons: # Check if codon is fourfold degenerate
                    third_sites_degen.append(codon[2]) # Append only third base
                else:
                    third_sites_degen.append('-')

            seq_str = ''.join(third_sites_degen) # Concatenate the nucleotides for each record (ID)
            new_rec = SeqRecord(Seq(seq_str), id = rec.id) # Convert to SeqRecord object
            records_list.append(new_rec) # Append SeqRecord to records list

            # Calculate GC content
            gc_third = calculate_gc_content(new_rec.seq)
            gc_all = calculate_gc_content(rec.seq)

            alignment_gc[rec.id] = {
                'gc_third': gc_third,
                'gc_all' : gc_all,
            }
        gc_data[species_name] = alignment_gc

        # Write all records at once, for this species
        with open(out_path, 'w') as output: 
            SeqIO.write(records_list, output, 'fasta')

    if save_to_file:
        out_json = gc_codon_json_path(group)
        os.makedirs(os.path.dirname(out_json) or ".", exist_ok=True)
        with open(out_json, 'w') as save:
            json.dump(gc_data, save)
    return gc_data

if __name__ == '__main__':
    for group in groups:
        gc_codon_dict(group)