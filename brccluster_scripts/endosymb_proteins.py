#!/usr/bin/env python3

import json
import os
import subprocess

with open('files/clade_dict.json', 'r') as f:
    clades_assm = json.load(f)

output_dir = 'test'
for sp, accns in clades_assm.items():
    try:
        sp_path = os.path.join(output_dir, sp)
        os.makedirs(sp_path, exist_ok = True)
        cmd = f'datasets download genome accession {",".join(accns)} --include protein --filename "{sp_path}/{sp}.zip"'
        output = subprocess.run(cmd, shell=True, check=True)
        print(f'Downloaded genomes for {sp}')
    except Exception as e:
        print(f'Error downloading genomes for {sp}: {e}')
