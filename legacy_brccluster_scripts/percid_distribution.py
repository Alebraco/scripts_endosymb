#!/usr/bin/env python3
"""
Script to calculate and visualize pairwise percentage identity (PID) distributions
for alignments of endosymbionts and free-living relatives.

Usage:
    ./percid_distribution.py <sp/all> [mean/all]

Arguments:
    <sp/all> : 'sp' to generate species-specific distributions, 'all' for combined distribution.
    [mean/all] : Optional. Use 'mean' to aggregate species-level means for the combined distribution.
                 Defaults to 'all' (all pairwise distances included).

Outputs:
    - Species-specific PID distribution plots (if 'sp' mode is selected).
    - Combined PID distribution plot (if 'all' mode is selected).
    - Plots are saved as PDF files.
"""

import os
import sys
from Bio import SeqIO
from pid_calculation import pairwise_pid
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math

group = 'endosymb+relatives/core_alignments'  # Directory containing alignment files

mean = False
if len(sys.argv) not in (2, 3):
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
crit = sys.argv[1]  # Mode: 'sp' for species-specific, 'all' for combined

if crit == 'sp':
    """
    Species-specific mode:
    - Generates a subplot for each species, showing the distribution of pairwise PIDs.
    - Saves the combined figure as 'sp_distributions.pdf'.
    """
    n = len(os.listdir(group))  # Number of species
    ncols = int(math.sqrt(n))  # Number of columns in the subplot grid
    nrows = n // ncols + (n % ncols > 0)  # Number of rows in the subplot grid

    plt.figure(figsize=(5 * ncols, 4 * nrows))  # Set figure size based on grid dimensions
    plt.suptitle('Pairwise Distance Distributions Across Endosymbionts and Free-Living Relatives', fontsize=16, y=0.98)
    i = 1  # Subplot index

    for sp in os.listdir(group):
        """
        Process each species:
        - Reads alignment files for the species.
        - Calculates pairwise PIDs for each alignment.
        - Plots the distribution as a histogram.
        """
        pid_sp = []  # List to store pairwise PIDs for the species
        print(f'Processing {sp}:')
        sp_path = os.path.join(group, sp)
        aln_files = os.listdir(sp_path)
        print(f'{len(aln_files)} alignments found.')
        
        for aln_file in aln_files:
            """
            Process each alignment file:
            - Reads the alignment in FASTA format.
            - Calculates pairwise PIDs using the `pairwise_pid` function.
            """
            print(f'>Processing {aln_file}')
            aln_path = os.path.join(sp_path, aln_file)
            aln = SeqIO.parse(aln_path, 'fasta')  # Parse alignment file
            
            pid_aln = pairwise_pid(aln)  # Calculate pairwise PIDs
            pid_sp.extend(pid_aln)

        n_bins = min(int(math.sqrt(len(pid_sp))), 30)  # Determine number of bins for the histogram

        plt.subplot(nrows, ncols, i)  # Create subplot for the species
        sns.histplot(x=pid_sp, bins=n_bins, stat='probability')  # Plot histogram
        plt.xlabel('Distance (% ID)')
        plt.ylabel('Probability')
        plt.title(sp)
        
        plt.text(0.95, 0.95, f'bins: {n_bins}\nn:{len(pid_sp)}', transform=plt.gca().transAxes,
                 ha='right', va='top')  # Add metadata to the plot
        i += 1

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit all subplots
    plt.savefig('sp_distributions.pdf')  # Save the figure
    plt.close()

elif crit == 'all':
    """
    Combined mode:
    - Aggregates pairwise PIDs across all species.
    - Optionally calculates species-level means if 'mean' is specified.
    - Saves the combined distribution plot as a PDF file.
    """
    pid_all = []  # List to store all pairwise PIDs
    for sp in os.listdir(group):
        """
        Process each species:
        - Reads alignment files for the species.
        - Calculates pairwise PIDs for each alignment.
        - Aggregates PIDs into a combined list.
        """
        pid_sp = []  # List to store pairwise PIDs for the species
        print(f'Processing {sp}:')
        sp_path = os.path.join(group, sp)
        aln_files = os.listdir(sp_path)
        
        for aln_file in aln_files:
            """
            Process each alignment file:
            - Reads the alignment in FASTA format.
            - Calculates pairwise PIDs using the `pairwise_pid` function.
            """
            print(f'>Processing {aln_file}')
            aln_path = os.path.join(sp_path, aln_file)
            aln = SeqIO.parse(aln_path, 'fasta')  # Parse alignment file
            
            pid_aln = pairwise_pid(aln)  # Calculate pairwise PIDs
            pid_sp.extend(pid_aln)
        if mean:
            mean_pid = np.mean(pid_sp)  # Calculate mean PID for the species
            pid_all.append(mean_pid)
        else:
            pid_all.extend(pid_sp)  # Add all PIDs to the combined list

    plt.figure()
    sns.histplot(x=pid_all, stat='probability')  # Plot combined histogram
    plt.xlabel('Distance (ID)')
    plt.ylabel('Probability')
    plt.title('Pairwise Distance Distributions\nEndosymbionts and Free-Living Relatives')
    plt.savefig(f'combined_distribution_{sys.argv[2]}.pdf')  # Save the figure
    plt.close()


