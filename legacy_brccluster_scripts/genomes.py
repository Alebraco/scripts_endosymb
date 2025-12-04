#!/usr/bin/env python3
import time
import subprocess
import os
import json

with open('endosymbiont_dict.json', 'r') as f:
    genome_dict = json.load(f)

for species,accns in genome_dict.items():
    try:
        print(f'Processing {species}.')
        if accns:
            name = species.replace(' ', '_')
            filename = f'{name}.zip'
            cmd = ['datasets', 'download', 'genome', 'accession']
            cmd.extend(accns)
            cmd.extend(['--include', 'genome,gff3'])
            cmd.extend(['--assembly-source', "RefSeq", '--filename', filename]) 
            subprocess.run(cmd, check=True)
            for accn in accns:
                subprocess.run([
                    'unzip', '-j', filename, f'ncbi_dataset/data/{accn}/*.fna', f'ncbi_dataset/data/{accn}/*.gff', '-d', name
                ], check=True)
                subprocess.run(f'mv {name}/*genomic.fna {name}/{accn}.fna', shell = True)
                subprocess.run(f'mv {name}/genomic.gff {name}/{accn}.gff', shell = True)
            os.remove(filename)

            print(f'Finished processing {species}.')
        else:
            print(f'No genome accession found for {species}.')
            continue

    except subprocess.CalledProcessError as e:
        print(f'Failed to fetch genome information for {species}', e)
    except Exception as e:
        print(f'Failed to fetch genome information for {species}', e)
    
    time.sleep(1)
