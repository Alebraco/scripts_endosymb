#!/usr/bin/env python3

import os
from itertools import combinations
import subprocess

genome_dir = 'clade_genomes/'
outdir = 'usearch_results/'
for sp in os.listdir(genome_dir):
    sp_path = os.path.join(genome_dir, sp)
    sp_output_dir = os.path.join(outdir, sp)
    genomes = os.listdir(sp_path)
    if len(genomes) >= 2:
        print(f'Processing {sp}')
        os.makedirs(sp_output_dir, exist_ok = True)
        for comb in combinations(genomes, 2):
            genome1_path = os.path.join(sp_path, comb[0])
            genome2_path = os.path.join(sp_path, comb[1])
            name1 = os.path.splitext(comb[0])[0]
            name2 = os.path.splitext(comb[1])[0]
            output_file = os.path.join(sp_output_dir, f"{name1}_vs_{name2}.txt")
            print(f'\n>Running USEARCH: {output_file}\n')
            cmd = [
                'usearch61',
                '-usearch_global', genome1_path,
                '-db', genome2_path,
                '-id', '0.7',
                '-maxaccepts', '1',
                '-blast6out', output_file,
                '-strand', 'both'
            ] 
            subprocess.run(cmd, check=True)
