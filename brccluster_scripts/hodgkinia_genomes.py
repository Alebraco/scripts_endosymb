#!/usr/bin/env python3
import subprocess
import os

hodgk_dir = 'endosymb_genomes/'
dirs = ['Hodgkinia_1', 'Hodgkinia_2', 'Hodgkinia_3']

for dir in dirs:
    path = os.path.join(hodgk_dir, dir)
    accns = [x.replace('.gff', '') for x in os.listdir(path)]
    for accn in accns:
        if accn:
            name = accn.replace(' ', '_')
            filename = f'{name}.zip'
            print(filename)
            cmd = ['datasets', 'download', 'genome', 'accession']
            cmd.append(accn)
            cmd.extend(['--include', 'genome'])
            cmd.extend(['--filename', filename]) 
            subprocess.run(cmd, check=True)
            subprocess.run([
                'unzip', '-j', filename, f'ncbi_dataset/data/{accn}/*.fna', '-d', path
            ], check=True)
            subprocess.run(f'mv {path}/{name}*genomic.fna {path}/{name}.fna', shell = True)
            os.remove(filename) 
