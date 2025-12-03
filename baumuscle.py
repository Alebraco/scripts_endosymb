#!/usr/bin/env python3

import os
import subprocess

out_dir = 'core_alignments/'
input_dir = 'endosymb_core/'
dirlist = ['Baumannia']

for sp in dirlist:
    core_path = os.path.join(input_dir, sp, 'core')
    print(f'Processing {sp}:')
    for file in os.listdir(core_path):
        in_file = os.path.join(core_path, file)
        file_name = file.split('.faa')[0]
        out_path = os.path.join(out_dir, sp)
        os.makedirs(out_path, exist_ok=True)
        out_file = os.path.join(out_path, f'{file_name}_aln.faa')


        cmd = ['muscle','-align', in_file, '-output', out_file]
        print(f'> Running alignment on {file}')
        subprocess.run(cmd, check=True)
