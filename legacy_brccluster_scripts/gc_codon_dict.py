#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from calculate_gc import calculate_gc_content
import os
import json

def gc_codon_dict(group, save_to_file = True):
    '''Create a dictionary with the GC content
    Produces fourfold degenerate and overall GC content 
    Computes GC for each core gene alignment concatenate
    args:
        group (str): endosymb_only, endosymb+relatives, relatives_only
        save_to_file (boolean): Saves dictionary as json file
    '''
    fourfold_codons = {
        "GCT", "GCC", "GCA", "GCG",   # Ala
        "CGT", "CGC", "CGA", "CGG",   # Arg
        "GGT", "GGC", "GGA", "GGG",   # Gly
        "CCT", "CCC", "CCA", "CCG",   # Pro
        "ACT", "ACC", "ACA", "ACG",   # Thr
        "GTT", "GTC", "GTA", "GTG",   # Val
        "CTT", "CTC", "CTA", "CTG",   # Leu
        "TCT", "TCC", "TCA", "TCG"    # Ser
    }
    gc_data = {}
    output_dir = 'third_sites/'
    input_dir = os.path.join(group, 'dna_concatenates/')

    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        path = os.path.join(input_dir, file)
        species_name = file.split('concatenate_')[1].split('.fasta')[0]
        out_path = os.path.join(output_dir, file)

        records_list = []
        alignment_gc = {}
        id_list = []

        for rec in SeqIO.parse(path, 'fasta'):
            third_sites_degen = []

            for i in range(0, len(rec.seq), 3): # Iterate over codon positions 
                codon = str(rec.seq[i:i+3]) # Retrieve each codon sequence
                if codon in fourfold_codons: # Check if codon is fourfold degenerate
                    third_sites_degen.append(codon[2]) # Append only third base

            seq_str = ''.join(third_sites_degen) # Concatenate the nucleotides for each record (ID)
            new_rec = SeqRecord(Seq(seq_str), id = rec.id) # Convert to SeqRecord object
            records_list.append(new_rec) # Append SeqRecord to records list
            id_list.append(rec.id) # Append record ID to ID list

            # Calculate GC content
            gc_third = calculate_gc_content(new_rec.seq)
            gc_all = calculate_gc_content(rec.seq)

            alignment_gc[rec.id] = {
                'gc_third': gc_third,
                'gc_all' : gc_all,
            }
        gc_data[species_name] = alignment_gc

        with open(out_path, 'w') as output: # Write all records at once, for each file
            SeqIO.write(records_list, output, 'fasta')

    if save_to_file:
        with open(f'gc_codon_data_{group}.json', 'w') as save:
            json.dump(gc_data, save)
    return gc_data


