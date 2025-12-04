#!/usr/bin/env python3

import os
import pandas as pd
import sys

data = []
if len(sys.argv) != 2:
    print('Usage: script.py <group>')
    sys.exit(1)
group = sys.argv[1]
main_dir = os.path.join(group,'core')
for sp in os.listdir(main_dir):
    sp_path = os.path.join(main_dir, sp)
    file = os.path.join(sp_path, 'summary.txt')
    if not os.path.isfile(file):
        continue
    with open(file, 'r') as f:
        first = f.readline()
        f.readline()
        third = f.readline()
    genes = first.split(':')[1].strip().split(' ')[0]
    genomes = third.split(':')[1].strip()
    data.append({
        'Species': sp,
        'Core_Genes': genes,
        'Genomes': genomes
    })
df = pd.DataFrame(data)
df.to_csv(f'core_summary_{group}.csv', index=False)
