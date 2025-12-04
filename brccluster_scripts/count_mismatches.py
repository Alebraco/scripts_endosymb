#!usr/bin/env python3

from Bio import AlignIO
import os
from itertools import combinations
import sys

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} <mafft/muscle alignments>')
    sys.exit(1)

output = f'{sys.argv[1]}_mm_counts.txt'
alg_dir = f'{sys.argv[1]}_alignments'

def pairwise_mismatches(alignment):
    pairs = list(combinations(alignment, 2))
    results = []
    for seq1, seq2 in pairs:
        mismatches = sum(base1 !=  base2 for base1,base2 in zip(seq1.seq, seq2.seq))
        results.append(f'{seq1.id}\t{seq2.id}\t{mismatches}')
    return results

with open(output, 'w') as f:
    for alg in os.listdir(alg_dir):
        if alg.endswith('.fasta'):
            path = os.path.join(alg_dir, alg)
            algobj = AlignIO.read(path, 'fasta')
            f.write(f'\n>{alg.split(".fasta")[0]}')
            pairwise_results = pairwise_mismatches(algobj)
            for line in pairwise_results:
                f.write(f'\n{line}\n')
