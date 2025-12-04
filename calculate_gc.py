#!/usr/bin/env python3

def calculate_gc_content(seq):
    '''Calculate GC content
    args: 
        seq (str): nucleotide sequence
    '''
    gc_count = seq.count('G') + seq.count('C')
    total_length = seq.count('G') + seq.count('C') + seq.count('A') + seq.count('T')
    if total_length > 0:
        return round((gc_count/total_length)*100, 2)
    else:
        return 0
