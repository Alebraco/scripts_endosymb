#!/usr/bin/env python3

import os
import re

cds_dir = 'endosymb_cds'
alignment_dir = 'core_alignments'
out_dir = 'backtranslated'

species = os.listdir(alignment_dir)
nt_dict = {}

for sp in species: # Iterate over species in the alignment directory
    sp_path = os.path.join(cds_dir, sp) # Obtain path to CDS for the current species 
    genomes = os.listdir(sp_path) # List CDS sequences in path 
    seq_dict = {} # Initialize a sequence dictionary 
    for genome in genomes: # Iterate over the CDS sequences of the species
        genome_path = os.path.join(sp_path, genome) # Obtain path to a particular CDS file
        with open(genome_path, 'r') as f: # Read the CDS file
            protein_id = None # Initialize protein ID variable 
            temp = '' 
            for line in f: # Iterate over each line of the file
                if line.startswith('>'): 
                    if protein_id and temp: # Store previous sequence
                        seq_dict[protein_id] = temp 
                    if 'protein_id' in line: # Obtain protein_id from CDS
                        protein_id = re.search(r'\[protein_id=([^\]]+)\]', line).group(1)
                        temp = '' # Reset temp to remove previous sequence 
                    else: # If no ID is present, skip
                        protein_id = None
                        temp = ''
                elif protein_id: # If currently processing a protein_id from CDS, store sequence in temp
                    temp += line.strip().upper() 
            if protein_id and temp: # Process last CDS in file
                seq_dict[protein_id] = temp
    nt_dict[sp] = seq_dict # Add CDS dictionary to species dictionary
    
    aln_path = os.path.join(alignment_dir, sp) # Obtain alignments path for the current species   
    alignments = os.listdir(aln_path) # List elements in the path
    backtranslate = os.path.join(out_dir, sp) # Define the output directory for the current species
    os.makedirs(backtranslate, exist_ok=True) # Create said directory

    for fam in alignments: # Iterate over the alignments available for the current species
        fam_path = os.path.join(aln_path, fam) # Obtain path of the alignment
        fam_name = fam.split('_aln')[0] # Retrive family name
        
        aa_seq = '' # Initialize amino acid sequence variable
        header = '' 
        id = None

        outfile = os.path.join(backtranslate, f'{fam_name}.fasta') # Define output file

        with open(fam_path, 'r') as f: # Read alignment file
            with open(outfile, 'w') as out:  # Start writing output file
                for line in f: 
                    if line.startswith('>'): # If new identifier starts 
                        if id: # Store previous ID (for the first header, this won't run)
                            if id in nt_dict[sp]: # If the ID is in the species dictionary for the current species
                                if len(nt_dict[sp][id]) >= 3 * (len(aa_seq) - aa_seq.count('-')): # Backtranslation
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
                                    
                                else: # If DNA sequence does not have sufficient length
                                    with open(f'missing_cds_{sp}.log', 'a') as log:
                                        log.write(f'\nInsufficient CDS length for {id} in {sp}, {fam}\n')
                                        log.write(f'{len(nt_dict[sp][id])} >= {3 * (len(aa_seq) - aa_seq.count("-"))}')
                            else: # If no match is found
                                with open(f'missing_cds_{sp}.log', 'a') as log:
                                    log.write(f'\nMissing CDS for {id} in {sp}, {fam}\n')
                        
                        id = line.split('.faa&')[1].strip() # Retrieve accession ID from the header (to match CDS and prot. alignment)
                        header = line # Store header to know if currently processing a sequence
                        aa_seq = '' # Reset amino acid sequence variable to remove previous
                        

                    elif line.strip(): # If not a header line, add sequence to variable
                        aa_seq += line.strip().upper()
                    
                # Process the last sequence
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