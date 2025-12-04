#!/usr/bin/env python3

"""
Utilities for the GC pipeline.

Contains:
- title mappings
- group name mappings
- helpers for loading or computing JSON and pickle data
- paths for GC-related JSON files
"""

import json
import os
import pickle

# Different matrix titles

files_dir = "files"

titles = {
    "distance": "Patristic Distance",
    "size": "ΔGenome Size",
    "gc_genome": "Genome ΔGC%",
    "gc_third": "4D Site ΔGC%",
    "gc_all": "Core ΔGC%",
}

abs_titles = {
    'size': 'Genome Size',
    'gc_genome': 'Genome GC%',
    'gc_third': '4D Site GC%', 
    'gc_all': 'Core GC%'
}

group_names = {
    "endosymb+relatives": "Endosymbionts and Free-Living Relatives",
    "relatives_only": "Free-Living Relatives Only",
    "endosymb_only": "Endosymbionts Only",
}

def gc_codon_json_path(group):
    """
    Path to the GC codon (GC3) JSON for a given group.
    """
    return os.path.join(files_dir, f"gc_codon_data_{group}.json")

def genome_gcsize_json_path(group):
    """
    Path to the genome size and GC JSON for a given group.
    """
    return os.path.join(files_dir, f"gcsize_genome_{group}.json")

def load_or_compute(
    filename,
    compute_function,
    *args,
    **kwargs,
):
    """
    Load JSON data from `filename` if it exists, otherwise run
    `compute_function(*args, **kwargs)`, save to `filename` (JSON),
    and return the result.
    """
    if os.path.isfile(filename):
        with open(filename, "r") as f:
            return json.load(f)
    else:
        data = compute_function(*args, **kwargs)
        os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
        with open(filename, "w") as f:
            json.dump(data, f)
        return data
    
def load_or_compute_pickle(
    filename,
    compute_function,
    *args,
    **kwargs,
):
    """
    Pickle-based version  of load_or_compute, used for
    distance matrices.
    """
    if os.path.isfile(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)
    else:
        data = compute_function(*args, **kwargs)
        os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)
        with open(filename, "wb") as f:
            pickle.dump(data, f)
        return data


