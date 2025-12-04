#!/usr/bin/env python3

from calculate_gc import calculate_gc_content
from Bio import SeqIO
import pandas as pd
import json
import os

def fetch_gc_size(sequence):
    genome_size = len(sequence)
    gc_content = calculate_gc_content(sequence)
    return genome_size, gc_content

def genome_gcsize(group):
    """
    Process genomes in input_dir and compute genome-wide size + GC content.

    Args:
        group (str): endosymb_only, endosymb+relatives, relatives_only.
        output_file (str): Path for output JSON file.

    Returns:
        dict: Nested dictionary with species → accession → metadata.
    """

    tree_dir = os.path.join(group, "dna_tree_results")
    genome_dir = os.path.join(group, "genomes")

    all_data = {}

    for sp in os.listdir(genome_dir):
        print(f"Processing {sp}")
        all_data[sp] = {}

        genome_path = os.path.join(genome_dir, sp)
        genome_list = [g for g in os.listdir(genome_path) if g.endswith(".fna")]

        for file in genome_list:
            accn = file.split(".fna")[0]
            file_path = os.path.join(genome_path, file)

            seqs = []
            for rec in SeqIO.parse(file_path, "fasta"):
                seqs.append(str(rec.seq))

            
            final_seq = "".join(seqs)
            genome_size, gc_content = fetch_gc_size(final_seq)

            all_data[sp][accn] = {
                "gc_genome": gc_content,
                "size": genome_size
            }

    with open(f'gcsize_genome_{group}.json', "w") as f:
        json.dump(all_data, f, indent=4)

    return all_data

