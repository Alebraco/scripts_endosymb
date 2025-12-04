#!/usr/bin/env python3

import os
from Bio import Entrez, SeqIO
Entrez.email = "asoneto@ncsu.edu"

def get_host(file):
    with open(file, 'r') as f:
        ids = [line.strip() for line in f]

    write_file = f'{file.split('.tsv')[0]}_host.tsv'
    with open(write_file, 'w') as f:
        batch_size = 200
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]
            handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "genbank")
            
            for record in records:
                host = None
                for feature in record.features:
                    if feature.type == "source":
                        host = feature.qualifiers.get('host', [None])[0]
                        break
                f.write(f"{record.id}\t{host}\n")

related_sp = 'related_species/'
for sp in os.listdir(related_sp):
    sp_path = os.path.join(related_sp, sp)
    for file in os.listdir(sp_path):
        if 'host' not in file and file.endswith('.tsv'):
            file_path = os.path.join(sp_path, file)
            get_host(file_path)
            print(f"Processed {file_path}")