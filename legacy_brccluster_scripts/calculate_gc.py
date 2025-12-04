#!/usr/bin/env python3

def calculate_gc_content(seq):
    '''Calculate GC content
    args: 
        seq (str): nucleotide sequence
    '''
    gc_count = seq.count('G') + seq.count('C')
    if len(seq) > 0:
        return round((gc_count/len(seq))*100, 2)
    else:
        return 0