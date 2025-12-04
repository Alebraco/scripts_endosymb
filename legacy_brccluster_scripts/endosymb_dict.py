#!/usr/bin/env python3

import pandas as pd
import json

species_list = []
with open('endosymbionts_list.txt', 'r') as f:
    for line in f:
        species_list.append(line.strip())

df = pd.read_csv('genome_metadata.tsv', sep='\t', usecols=[0, 1])

genome_dict = {}
for species in species_list:
    sp_df = df[df['Organism Name'].str.contains(species, case=False)]

    accession_list = sp_df['Assembly Accession'].tolist()

    genome_dict[species] = accession_list
with open('endosymbiont_dict.json', 'w') as f:
    json.dump(genome_dict, f, indent=4)


    