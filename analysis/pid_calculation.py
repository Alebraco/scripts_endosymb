#!/usr/bin/env python3 
from itertools import combinations

def pairwise_pid(alignment):
    '''Receives alignment object and calculates percent identity'''
    pairs = list(combinations(alignment, 2))
    results = []
    gaps = False
    for seq1, seq2 in pairs:
        if ('genomic' in seq1.id) == ('genomic' in seq2.id):
            continue

        seq1_str = seq1.seq
        seq2_str = seq2.seq

        if gaps:
            mismatches = sum(1 for a,b in zip(seq1_str, seq2_str) if a != b)
            length = len(seq1_str)
        else:
            mismatches = 0
            length = 0 
            for a,b in zip(seq1_str, seq2_str):
                if a != '-' and b != '-':
                    length += 1
                    if a != b:
                        mismatches += 1

        pid = mismatches / length if length > 0 else 0
        results.append(pid)
    return results
