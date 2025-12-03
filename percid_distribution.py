
#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from pid_calculation import pairwise_pid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import math
group = 'endosymb+relatives/core_alignments'

mean = False
if len(sys.argv) not in (2,3):
    print(f'Usage: {sys.argv[0]} <sp/all> optional(mean/all)')
    sys.exit(1)
elif len(sys.argv) == 3:
    agg = sys.argv[2]
    if sys.argv[1] == 'sp':
        print('Warning: Last parameter ignored for sp mode.')
    elif agg == 'mean':
        mean = True
    else:
        mean = False
crit = sys.argv[1]

if crit == 'sp':
    n = len(os.listdir(group))
    ncols = int(math.sqrt(n))
    nrows = n // ncols + (n % ncols > 0)

    plt.figure(figsize = (5 * ncols, 4 * nrows))
    plt.suptitle('Pairwise Distance Distributions Across Endosymbionts and Free-Living Relatives', fontsize=16, y=0.98)
    i = 1
    for sp in os.listdir(group):
        pid_sp = []
        print(f'Processing {sp}:')
        sp_path = os.path.join(group, sp)
        aln_files = os.listdir(sp_path)
        print(f'{len(aln_files)} alignments found.')
        
        for aln_file in aln_files:
            print(f'>Processing {aln_file}')
            aln_path = os.path.join(sp_path, aln_file)
            aln = SeqIO.parse(aln_path, 'fasta')
            
            pid_aln = pairwise_pid(aln)
            pid_sp.extend(pid_aln)

        n_bins = min(int(math.sqrt(len(pid_sp))), 30) 

        plt.subplot(nrows, ncols, i)
        sns.histplot(x = pid_sp, bins = n_bins, stat = 'probability') 
        #plt.hist(pid_sp, bins = n_bins, density = True)
        plt.xlabel('Distance (% ID)')
        plt.ylabel('Probability')
        plt.title(sp)
        
        plt.text(0.95, 0.95, f'bins: {n_bins}\nn:{len(pid_sp)}', transform=plt.gca().transAxes,
                ha='right', va='top')
        i += 1

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('sp_distributions.pdf')
    plt.close()
elif crit == 'all':
    pid_all = []
    sp_info = []
    for sp in os.listdir(group):
        pid_sp = []
        print(f'Processing {sp}:')
        sp_path = os.path.join(group, sp)
        aln_files = os.listdir(sp_path)
        
        for aln_file in aln_files:
            print(f'>Processing {aln_file}')
            aln_path = os.path.join(sp_path, aln_file)
            aln = SeqIO.parse(aln_path, 'fasta')
            
            pid_aln = pairwise_pid(aln)
            pid_sp.extend(pid_aln)
        if mean:
            mean_pid = np.mean(pid_sp)
            pid_all.append(mean_pid)

            sp_info.append({'species': sp, 'mean pid': mean_pid})
            df = pd.DataFrame(sp_info)
            df.sort_values('mean pid')
            df.to_csv('pid_per_sp.csv', index = False)
        else:
            pid_all.extend(pid_sp)


    plt.figure()
    sns.histplot(x = pid_all, stat = 'probability') 
    plt.xlabel('Distance (ID)')
    plt.ylabel('Probability')
    plt.title('Pairwise Distance Distributions\nEndosymbionts and Free-Living Relatives')
    plt.savefig(f'combined_distribution_{sys.argv[2]}.pdf')
    plt.close()



