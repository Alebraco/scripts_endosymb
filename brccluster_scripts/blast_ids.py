#!/usr/bin/env python3

import os
import subprocess

blast_dir = 'blast_results/'
clades_dir = 'clades_16S/'
silva_db = "/home6/asoneto/endosymb/silvadb/SILVA_16S"

for sp in os.listdir(clades_dir):
    sp_path = os.path.join(clades_dir, sp)
    for file in os.listdir(sp_path): # Contains accession names (e.g. NZ_CACTIE1000084.1_766-2227.fna)
        file_name = file.split('.fna')[0] # Get the file name without the extension
        file_path = os.path.join(sp_path, file) # Full path to the file

        outdir = os.path.join(blast_dir, sp)
        os.makedirs(outdir, exist_ok=True)
        print(f'Blasting {sp}: {file_name}')
        outfile = os.path.join(outdir, f'{file_name}.tsv')

        cmd = [
            'blastn',
            '-query', file_path,
            '-db', silva_db,
            '-out', outfile,
            '-outfmt', '6 qseqid sseqid pident length evalue bitscore stitle',
            '-evalue', '1e-10',
            '-num_threads', '8',
        ]
        subprocess.run(cmd, check=True)

