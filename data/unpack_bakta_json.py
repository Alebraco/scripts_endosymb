#!/usr/bin/env python3
"""Unpack ATB Bakta JSONs into .faa and .gff files.

Reads   {annotations_dir}/{Species}/{SAMPLE}.bakta.json
Writes  {proteins_dir}/{Species}/{SAMPLE}.faa
        {genomes_dir}/{Species}/{SAMPLE}.gff
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile


def process_sample(json_path, sample, species, proteins_dir, genomes_dir, force):
    faa_path = os.path.join(proteins_dir, species, f'{sample}.faa')
    gff_path = os.path.join(genomes_dir, species, f'{sample}.gff')

    if not force and os.path.exists(faa_path) and os.path.exists(gff_path):
        return 'existing'

    os.makedirs(os.path.dirname(faa_path), exist_ok=True)
    os.makedirs(os.path.dirname(gff_path), exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        result = subprocess.run(
            ['bakta_io', '--output', tmp, '--prefix', sample, json_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f'ERROR: bakta_io failed for {json_path}: {result.stderr}', file=sys.stderr)
            return 'error'

        tmp_faa = os.path.join(tmp, f'{sample}.faa')
        tmp_gff = os.path.join(tmp, f'{sample}.gff')
        if not os.path.exists(tmp_faa) or not os.path.exists(tmp_gff):
            print(f'ERROR: bakta_io did not produce expected files for {json_path}', file=sys.stderr)
            return 'error'

        shutil.move(tmp_faa, faa_path)
        shutil.move(tmp_gff, gff_path)

    return 'extracted'


def process_species(species, annotations_dir, genomes_dir, proteins_dir, force):
    ann_sp = os.path.join(annotations_dir, species)
    if not os.path.isdir(ann_sp):
        print(f'ERROR: species dir not found: {ann_sp}', file=sys.stderr)
        return

    counts = {'extracted': 0, 'existing': 0, 'error': 0}
    for name in sorted(os.listdir(ann_sp)):
        if name.endswith('.bakta.json.gz'):
            sample = name[:-len('.bakta.json.gz')]
        elif name.endswith('.bakta.json'):
            sample = name[:-len('.bakta.json')]
        else:
            continue
        json_path = os.path.join(ann_sp, name)
        status = process_sample(json_path, sample, species, proteins_dir, genomes_dir, force)
        counts[status] = counts.get(status, 0) + 1

    print(f'{species}: extracted={counts["extracted"]} '
          f'skipped_existing={counts["existing"]} '
          f'errors={counts["error"]}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotations-dir', required=True)
    p.add_argument('--genomes-dir', required=True)
    p.add_argument('--proteins-dir', required=True)
    p.add_argument('--species', help='Process only this species subdirectory. If omitted, process all.')
    p.add_argument('--force', action='store_true', help='Re-extract even if outputs already exist.')
    args = p.parse_args()

    if args.species:
        species_list = [args.species]
    else:
        species_list = sorted(
            d for d in os.listdir(args.annotations_dir)
            if os.path.isdir(os.path.join(args.annotations_dir, d))
        )

    for species in species_list:
        process_species(species, args.annotations_dir, args.genomes_dir, args.proteins_dir, args.force)


if __name__ == '__main__':
    main()
