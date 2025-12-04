#!/usr/bin/env python3 

import os
import pandas as pd

output = 'Hodgkinia_core_output/'

for sp in os.listdir(output):
    sp_path = os.path.join(output, sp, 'CC')
    pair_data = {}
    data = []
    for pair in os.listdir(sp_path):
        genome1 = pair.split('.faa-')[0] + '.faa'
        genome2 = pair.split('.faa-')[1]
        filepath = os.path.join(sp_path, pair)
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if not lines:
            continue

        hit_count = len(lines)
        pids = [float(line.strip().split('\t')[2]) for line in lines]
        total_pid = sum(pids)
        min_id = min(pids) 

        key = tuple(sorted((genome1, genome2)))

        if key not in pair_data:
            pair_data[key] = {
                'hits': hit_count,
                'total_pid': total_pid,
                'min_id': min_id
            }
        else:
            pair_data[key]['hits'] += hit_count
            pair_data[key]['total_pid'] += total_pid
            pair_data[key]['min_id'] = min(pair_data[key]['min_id'], min_id)

    allgenomes = [genome for pair in pair_data.keys() for genome in pair]
    genome_counts = pd.Series(allgenomes).value_counts()
    pivot = genome_counts.index[0]

    for (genome1, genome2), stats in pair_data.items():
        avg_pid = round(stats['total_pid'] / stats['hits'], 1)
        nonpivot = genome2 if genome1 == pivot else genome1
        data.append({
            'genome': nonpivot,
            'hits': stats['hits'],
            'avg_pid': avg_pid,
            'min_id': stats['min_id'],
        })

    df = pd.DataFrame(data).sort_values(['hits', 'avg_pid'], ascending = [False, False])

    df.to_csv(f'{sp}.tsv', sep = '\t', index=False)



