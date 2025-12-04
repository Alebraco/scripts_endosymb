#!/usr/bin/env python3

import os
from Bio import Entrez
import time
import json

Entrez.email = 'alekey039@hotmail.com'
clades_path = 'files/main_clades_ids.txt'
# objective: to fetch genomes and get them into a structured folder
clades_dict = {}
clades_assm = {}
with open(clades_path, 'r') as f:
    for line in f:
        if line.startswith('>'):
            key = line.split('>')[1].replace('.fna', '').strip()
            clades_dict[key] = []
        elif line.strip():
            value = line.strip()
            accn = value.split('.')[0]
            clades_dict[key].append(accn)          
            
parent_dir = 'clades_genomes'
for species in clades_dict.keys():
    assm_ids = []
    assm_accns = []
    print(f'> {species}')
    for accn in clades_dict.get(species):
        print(f'\tProcessing {accn}')
        try:
            handle = Entrez.elink(
            dbfrom='nuccore',
            db='assembly',
            id=accn
            )
            record = Entrez.read(handle)
            handle.close()
            time.sleep(1)
            if record and len(record) > 0 and 'LinkSetDb' in record[0]:
                for linkset in record[0]['LinkSetDb']:
                    for link in linkset['Link']:
                        assm_ids.append(link['Id'])            

        except Exception as e:
            print(f'\t >Failed to process {species}: {e}')
            continue

    handle2 = Entrez.esummary(db = 'assembly', id = assm_ids)
    summary = Entrez.read(handle2)
    handle2.close()
    time.sleep(1)
    for doc in summary['DocumentSummarySet']['DocumentSummary']:
        assm_accns.append(doc.get('AssemblyAccession'))
        
    clades_assm[species] = assm_accns

with open('clade_dict.json', 'w') as f:
    json.dump(clades_assm, f)