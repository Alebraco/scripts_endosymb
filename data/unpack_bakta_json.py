#!/usr/bin/env python3
"""Unpack ATB Bakta JSONs into .faa and .gff files matching the feature pipeline layout.

Reads   {annotations_dir}/{Species}/{SAMPLE}.bakta.json
Writes  {proteins_dir}/{Species}/{SAMPLE}.faa
        {genomes_dir}/{Species}/{SAMPLE}.gff    (co-located with existing {SAMPLE}.fna)

Samples whose matching {SAMPLE}.fna is missing from genomes_dir are skipped, because
downstream codon-usage analysis requires the co-located assembly.
"""

import argparse
import gzip
import json
import os
import sys

GFF_TYPE_MAP = {
    'cds': 'CDS',
    'trna': 'tRNA',
    'rrna': 'rRNA',
    'tmrna': 'tmRNA',
    'ncrna': 'ncRNA',
    'ncrna-region': 'ncRNA',
    'crispr': 'repeat_region',
    'oric': 'origin_of_replication',
    'oriv': 'origin_of_replication',
    'orit': 'origin_of_transfer',
    'gap': 'gap',
}

FAA_WRAP = 60


def load_json(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return json.load(fh)


def get_aa(feature):
    for key in ('aa', 'translation', 'sequence'):
        val = feature.get(key)
        if isinstance(val, str) and val:
            return val.rstrip('*')
    return None


def feature_attrs(feat):
    parts = []
    locus = feat.get('locus') or feat.get('locus_tag') or feat.get('id')
    if locus:
        parts.append(f'ID={locus}')
        parts.append(f'locus_tag={locus}')
    gene = feat.get('gene')
    if gene:
        parts.append(f'gene={gene}')
    product = feat.get('product')
    if product:
        parts.append(f'product={escape_gff(product)}')
    if feat.get('pseudogene') or feat.get('pseudo'):
        parts.append('pseudo=true')
    return ';'.join(parts) if parts else '.'


def escape_gff(text):
    return (text.replace('%', '%25')
                .replace(';', '%3B')
                .replace('=', '%3D')
                .replace('&', '%26')
                .replace(',', '%2C')
                .replace('\t', ' ')
                .replace('\n', ' '))


def write_gff(features, out_path):
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for feat in features:
            raw_type = str(feat.get('type', '')).lower()
            gff_type = GFF_TYPE_MAP.get(raw_type)
            if gff_type is None:
                continue
            contig = feat.get('contig')
            start = feat.get('start')
            stop = feat.get('stop')
            strand = feat.get('strand', '.')
            if contig is None or start is None or stop is None:
                continue
            if strand not in ('+', '-', '.'):
                strand = '.'
            phase = '0' if gff_type == 'CDS' else '.'
            fh.write('\t'.join([
                str(contig),
                'Bakta',
                gff_type,
                str(start),
                str(stop),
                '.',
                strand,
                phase,
                feature_attrs(feat),
            ]) + '\n')


def write_faa(features, out_path):
    count = 0
    with open(out_path, 'w') as fh:
        for feat in features:
            if str(feat.get('type', '')).lower() != 'cds':
                continue
            if feat.get('pseudogene') or feat.get('pseudo'):
                continue
            aa = get_aa(feat)
            if not aa:
                continue
            locus = feat.get('locus') or feat.get('locus_tag') or feat.get('id') or f'cds_{count}'
            product = feat.get('product', '')
            header = f'>{locus} {product}'.rstrip()
            fh.write(header + '\n')
            for i in range(0, len(aa), FAA_WRAP):
                fh.write(aa[i:i + FAA_WRAP] + '\n')
            count += 1
    return count


def process_sample(json_path, fna_path, faa_path, gff_path, force):
    if not os.path.exists(fna_path):
        return 'no_fna'
    if not force and os.path.exists(faa_path) and os.path.exists(gff_path):
        return 'existing'

    try:
        data = load_json(json_path)
    except Exception as e:
        print(f'ERROR: failed to read {json_path}: {e}', file=sys.stderr)
        return 'error'

    features = data.get('features') or []
    if not features:
        print(f'WARN: no features in {json_path}', file=sys.stderr)
        return 'error'

    os.makedirs(os.path.dirname(faa_path), exist_ok=True)
    os.makedirs(os.path.dirname(gff_path), exist_ok=True)

    tmp_faa = faa_path + '.tmp'
    tmp_gff = gff_path + '.tmp'
    try:
        n_cds = write_faa(features, tmp_faa)
        write_gff(features, tmp_gff)
        if n_cds == 0:
            print(f'WARN: no CDS with translation in {json_path}', file=sys.stderr)
        os.replace(tmp_faa, faa_path)
        os.replace(tmp_gff, gff_path)
    finally:
        for p in (tmp_faa, tmp_gff):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
    return 'extracted'


def process_species(species, annotations_dir, genomes_dir, proteins_dir, force):
    ann_sp = os.path.join(annotations_dir, species)
    if not os.path.isdir(ann_sp):
        print(f'ERROR: species dir not found: {ann_sp}', file=sys.stderr)
        return

    counts = {'extracted': 0, 'existing': 0, 'no_fna': 0, 'error': 0}
    for name in sorted(os.listdir(ann_sp)):
        if not name.endswith('.bakta.json') and not name.endswith('.bakta.json.gz'):
            continue
        sample = name[:-len('.bakta.json.gz')] if name.endswith('.bakta.json.gz') else name[:-len('.bakta.json')]
        json_path = os.path.join(ann_sp, name)
        fna_path = os.path.join(genomes_dir, species, f'{sample}.fna')
        faa_path = os.path.join(proteins_dir, species, f'{sample}.faa')
        gff_path = os.path.join(genomes_dir, species, f'{sample}.gff')
        status = process_sample(json_path, fna_path, faa_path, gff_path, force)
        counts[status] = counts.get(status, 0) + 1
        if status == 'no_fna':
            print(f'WARN: skipping {species}/{sample} — no .fna in {genomes_dir}', file=sys.stderr)

    print(f'{species}: extracted={counts["extracted"]} '
          f'skipped_existing={counts["existing"]} '
          f'skipped_no_fna={counts["no_fna"]} '
          f'errors={counts["error"]}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotations-dir', required=True, help='e.g. atb_bacteria/annotations')
    p.add_argument('--genomes-dir', required=True, help='e.g. atb_bacteria/genomes (holds .fna, will receive .gff)')
    p.add_argument('--proteins-dir', required=True, help='e.g. atb_bacteria/proteins (will receive .faa)')
    p.add_argument('--species', help='Process only this species subdirectory (by name). If omitted, process all.')
    p.add_argument('--force', action='store_true', help='Re-extract even if .faa and .gff already exist.')
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
