#!/usr/bin/env python3

import json
import os
import subprocess

output_dir = 'endosymb_cds'
prot_dir = 'endosymb_proteins/'
for sp in os.listdir(prot_dir):
    sp_path = os.path.join(prot_dir, sp)
    file_list = []
    for file in os.listdir(sp_path):
        gcf = file.split('.faa')[0]
        file_list.append(gcf)
    out_path = os.path.join(output_dir, sp)
    os.makedirs(out_path, exist_ok = True)
    try:
        cmd = f'datasets download genome accession {",".join(file_list)} --include cds --filename "{out_path}/{sp}.zip"'
        output = subprocess.run(cmd, shell = True, check = True)
        print(f'Downloaded genomes for {sp}')
    except: 
        print(f'Error downloading genomes for {sp}')