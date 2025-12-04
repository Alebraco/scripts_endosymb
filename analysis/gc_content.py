#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    file = sys.argv[1]
    print(f'Processing {file}')
else: 
    print('No file included')

with open(file, 'r') as f:
    content = f.read()

gc_content = []
for seq in content.strip().split('>')[1:]:
    lines = seq.split('\n')
    id = lines[0]
    sequence = ''.join(lines[1:]).replace('\n', '')
    gc_count = sequence.count('G') + sequence.count('C')
    seq_length = len(sequence) - sequence.count('-')
    gc_perc = (gc_count / seq_length) * 100
    gc_content.append((id, gc_perc))

gc_list = [perc for _,perc in gc_content]
average_gc = sum(perc for _,perc in gc_content) / len(gc_content)
min_gc = min(gc_list)
max_gc = max(gc_list)
min_index, max_index = gc_list.index(min_gc), gc_list.index(max_gc)
min_id, max_id = gc_content[min_index], gc_content[max_index]

#print(min_id, max_id)
#print('\nGC Content:')
#print(gc_content)

print(f'Minimum GC%: {min_gc:.2f}')
print(f'Average GC%: {average_gc:.2f}')
print(f'Maximum GC%: {max_gc:.2f}')
