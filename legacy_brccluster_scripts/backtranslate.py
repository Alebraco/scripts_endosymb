#!/usr/bin/env python3

import os
import re
import sys

if len(sys.argv) != 2:
    print(f'Usage {sys.argv[0]} <group>')
    sys.exit(1)
group = sys.argv[1]

cds_dir = os.path.join(group, 'cds')
alignment_dir = os.path.join(group, 'core_alignments')
out_dir = os.path.join(group, 'backtranslated')
species = os.listdir(alignment_dir)
nt_dict = {}

for sp in species:
    sp_path = os.path.join(cds_dir, sp)
    genomes = os.listdir(sp_path)
    seq_dict = {}
    for genome in genomes:
        genome_path = os.path.join(sp_path, genome)
        with open(genome_path, 'r') as f:
            protein_id = None
            temp = ''
            for line in f:
                if line.startswith('>'):
                    if protein_id and temp:
                        seq_dict[protein_id] = temp

                    if 'protein_id' in line:
                        protein_id = re.search(r'\[protein_id=([^\]]+)\]', line).group(1)
                        
                    else:
                        if 'genomic' in genome:
                            protein_id = line.split()[0].lstrip('>')
                        else:
                            protein_id = None
                    temp = ''
                elif protein_id:
                    temp += line.strip().upper()
            if protein_id and temp:
                seq_dict[protein_id] = temp
    nt_dict[sp] = seq_dict
    
    print(f"CDS IDs loaded for {sp}: {len(seq_dict)}")
    if len(seq_dict) > 0:
        print("Sample CDS IDs:", list(seq_dict.keys())[:3])

    aln_path = os.path.join(alignment_dir, sp)
    alignments = os.listdir(aln_path)
    backtranslate = os.path.join(out_dir, sp)
    os.makedirs(backtranslate, exist_ok=True)

    for fam in alignments:
        fam_path = os.path.join(aln_path, fam)
        fam_name = fam.split('_aln')[0]
        print(f"\n=== Processing {fam} ===") 

        aa_seq = ''
        header = ''
        id = None

        outfile = os.path.join(backtranslate, f'{fam_name}.fasta')

        with open(fam_path, 'r') as f:
            with open(outfile, 'w') as out:
                for line in f:
                    if line.startswith('>'):
                        if id:
                            if id in nt_dict[sp]:
                                if len(nt_dict[sp][id]) >= 3 * (len(aa_seq) - aa_seq.count('-')):
                                    nt_aln = ''
                                    u = 0
                                    for i in range(len(aa_seq)):
                                        if aa_seq[i] == '-':
                                            nt_aln += '---'
                                            
                                        else:
                                            nt_aln += nt_dict[sp][id][u*3:u*3+3]
                                            u += 1
                                    out.write(header)
                                    out.write(nt_aln + '\n')
                                    
                                else:
                                    with open(f'missing_cds_{sp}.log', 'a') as log:
                                        log.write(f'\nInsufficient CDS length for {id} in {sp}, {fam}\n')
                                        log.write(f'{len(nt_dict[sp][id])} >= {3 * (len(aa_seq) - aa_seq.count("-"))}')
                            else:
                                with open(f'missing_cds_{sp}.log', 'a') as log:
                                    log.write(f'\nMissing CDS for {id} in {sp}, {fam}\n')

                        id = line.split('.faa&')[1].strip()
                        header = line
                        aa_seq = ''

                        print(f"Looking for ID: '{id}'")
                        print(f"Is it in seq_dict? {id in seq_dict}")
                        if id not in seq_dict:
                            print("Available similar IDs:", [k for k in seq_dict.keys() if id in k][:3])
                        

                    elif line.strip():
                        aa_seq += line.strip().upper()
                    
                # last sequence
                if id and header and aa_seq:
                    nt_aln = ''
                    if id in nt_dict[sp]:
                        if len(nt_dict[sp][id]) >= 3 * (len(aa_seq) - aa_seq.count('-')):
                            u = 0
                            for i in range(len(aa_seq)):
                                if aa_seq[i] == '-':
                                    nt_aln += '-' * 3
                                else:
                                    nt_aln += nt_dict[sp][id][u*3:u*3+3]
                                    u += 1
                            out.write(header)
                            out.write(nt_aln + '\n')

                        else:
                            with open(f'missing_cds_{sp}.log', 'a') as log:
                                log.write(f'\nInsufficient CDS length for {id} in {sp}, {fam}\n')
                                log.write(f'{len(nt_dict[sp][id])} >= {3 * (len(aa_seq) - aa_seq.count("-"))}')
                    else:
                        with open(f'missing_cds_{sp}.log', 'a') as log:
                            log.write(f'\nMissing nucleotide sequence for {id} in {sp}, {fam}\n') 