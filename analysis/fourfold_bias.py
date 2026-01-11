#!/usr/bin/env python3
import os
from Bio import SeqIO

fourfold_families = {
    'GC': 'Ala', 'GG': 'Gly', 'CC': 'Pro', 'AC': 'Thr', 
    'GT': 'Val', 'CG': 'Arg4', 'CT': 'Leu4', 'TC': 'Ser4'
}
groups = ['relatives_only', 'endosymb_only']

for group in groups:
    group_path = os.path.join(group, 'dna_concatenates')
    for species in os.listdir(group_path):
        species_path = os.path.join(group_path, species)
        counts = {aa: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for aa in fourfold_families.values()}
        for record in SeqIO.parse(species_path, 'fasta'):
            seq = str(record.seq).upper()
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                prefix = codon[0:2]
                third_base = codon[2]

                if prefix in fourfold_families.keys() and third_base in 'ACGT':
                    aa = fourfold_families[prefix]
                    counts[aa][third_base] += 1
        print(f'--- Results for {species} in {group} ---')
        total_gc4 = 0
        total_count = 0

        for aa, nuc_counts in counts.items():
            gc_count = nuc_counts['G'] + nuc_counts['C']
            sum_count = sum(nuc_counts.values())

            if sum_count > 0:
                gc4_percent = round((gc_count / sum_count) * 100, 2)
                print(f'{aa} GC4: {gc4_percent}% (n={sum_count})')

                total_gc4 += nuc_counts['G'] + nuc_counts['C']
                total_count += sum_count
        if total_count > 0:
            overall_gc4 = round((total_gc4 / total_count) * 100, 2)
            print(f'Overall GC4: {overall_gc4}%')
    print('----------------------------------')


        


            
        
    
