#!/usr/bin/env python3

import os
from Bio import SeqIO, AlignIO
import sys

id_dict = {}
sp_list = []

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <group>')
    sys.exit()

group = sys.argv[1]
core_dir = os.path.join(group,'core_alignments') # Used to list species only
genome_dir = os.path.join(group, 'proteins') # Used to retrieve IDs only
alignment_dir = os.path.join(group, 'backtranslated') # Fetch the actual sequences

for sp in os.listdir(core_dir):
    sp_list.append(sp)

for sp in sp_list:
    genome_path = os.path.join(genome_dir, sp)
    genome_files = os.listdir(genome_path)
    genomes = list(set(map(lambda x : x.split('.faa')[0], genome_files)))
    id_dict[sp] = genomes



for sp in os.listdir(alignment_dir):
    concatenate = {genome : '' for genome in id_dict[sp]}
    print(f'Processing {sp}')
    sp_path = os.path.join(alignment_dir, sp)
    alignments = []
    for file in os.listdir(sp_path):
        aln_path = os.path.join(sp_path, file)
        aln = AlignIO.read(aln_path, 'fasta')
        alignments.append(aln)
    
    for aln in alignments:
        aln_length = aln.get_alignment_length()
        aln_ids = {}
        for rec in aln:
            clean_id = rec.id.split('.faa&')[0]
            aln_ids[clean_id] = str(rec.seq)
            
        for genome in concatenate:
            if genome in aln_ids:
                concatenate[genome] += aln_ids[genome]
            else:
                concatenate[genome] += '-' * aln_length

    with open(f'concatenate_{sp}.fasta', 'w') as output:
        for ids, seq in concatenate.items():
            output.write(f'>{ids}\n{seq}\n')
