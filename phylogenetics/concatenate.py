#!/usr/bin/env python3

import os
from Bio import SeqIO, AlignIO
import sys

id_dict = {}
sp_list = []

if len(sys.argv) not in [2, 3]:
    print(f'Usage: {sys.argv[0]} <group> [protein]')
    sys.exit(1)

group = sys.argv[1]

mode = 'protein' if len(sys.argv) == 3 and sys.argv[2] == 'protein' else 'dna'

core_dir = os.path.join(group,'core_alignments') # Used to list species only
genome_dir = os.path.join(group, 'proteins') # Used to retrieve IDs only

if mode == 'protein':
    alignment_dir = os.path.join(group, 'core_alignments') # Fetch the actual sequences
    output_dir = os.path.join(group, 'protein_concatenates')
    file_extension = '.faa'
else:
    alignment_dir = os.path.join(group, 'backtranslated') # Fetch the actual sequences
    output_dir = os.path.join(group, 'dna_concatenates')
    file_extension = '.fna'

os.makedirs(output_dir, exist_ok = True)

# List species
for sp in os.listdir(core_dir):
    sp_list.append(sp)

# List genomes for each species
for sp in sp_list:
    genome_path = os.path.join(genome_dir, sp)
    genome_files = os.listdir(genome_path)
    genomes = list(set(map(lambda x : x.split('.faa')[0], genome_files)))
    id_dict[sp] = genomes

# Concatenate alignments for each species
for sp in os.listdir(alignment_dir):
    concatenate = {genome : '' for genome in id_dict[sp]}
    print(f'Processing {sp}, {mode} mode')
    sp_path = os.path.join(alignment_dir, sp)
    alignments = []
    # Read all alignments for the species
    for file in os.listdir(sp_path):
        aln_path = os.path.join(sp_path, file)
        aln = AlignIO.read(aln_path, 'fasta')
        alignments.append(aln)
    
    # Concatenate sequences for each genome
    for aln in alignments:
        aln_length = aln.get_alignment_length()
        aln_ids = {}
        for rec in aln:
            clean_id = rec.id.split('.faa&')[0]
            aln_ids[clean_id] = str(rec.seq)

        # Add sequences to the concatenate dictionary
        for genome in concatenate:
            if genome in aln_ids:
                concatenate[genome] += aln_ids[genome]
            else:
                concatenate[genome] += '-' * aln_length
    outpath = os.path.join(output_dir, f'concatenate_{sp}{file_extension}')
    with open(outpath, 'w') as out:
        for ids, seq in concatenate.items():
            out.write(f'>{ids}\n{seq}\n')

