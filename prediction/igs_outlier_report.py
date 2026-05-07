"""Diagnostic report on Mean_IGS_Size outliers in the RF training set.

Lists every genome with Mean_IGS_Size > 500 bp and looks up genome size and
contig count by parsing the assembly FASTA under annotation/genomes/. Flags
values > 1000 bp as likely artifacts (assembly fragmentation, annotation gaps,
partitioned genomes such as Hodgkinia).

"""
import argparse
import os
import sys

import pandas as pd

DEFAULT_INPUT = os.path.join(
    'endosymb+relatives', 'feature_files', 'combined_features.csv'
)
DEFAULT_OUTPATH = os.path.join('files', 'igs_outlier_report.md')
DEFAULT_GENOMES_DIR = os.path.join('endosymb+relatives', 'genomes')
FASTA_EXTS = ('.fna', '.fasta', '.fa', '.fna.gz', '.fasta.gz', '.fa.gz')

FLAG_THRESHOLD = 500
ARTIFACT_THRESHOLD = 1000


def fasta_stats(path):
    """Return (contig_count, genome_size_bp) for a FASTA, or (None, None)."""
    if path.endswith('.gz'):
        import gzip
        opener = lambda p: gzip.open(p, 'rt')
    else:
        opener = lambda p: open(p, 'r')
    contigs = 0
    size = 0
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith('>'):
                    contigs += 1
                else:
                    size += len(line.strip())
    except OSError:
        return None, None
    return contigs, size


def find_fasta(genomes_dir, file_stem):
    try:
        species_dirs = [
            os.path.join(genomes_dir, d)
            for d in os.listdir(genomes_dir)
            if os.path.isdir(os.path.join(genomes_dir, d))
        ]
    except OSError:
        species_dirs = []
    for folder in [genomes_dir] + species_dirs:
        for ext in FASTA_EXTS:
            candidate = os.path.join(folder, file_stem + ext)
            if os.path.exists(candidate):
                return candidate
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=DEFAULT_INPUT)
    parser.add_argument('--genomes-dir', default=DEFAULT_GENOMES_DIR)
    parser.add_argument('--output', default=DEFAULT_OUTPATH)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input CSV not found: {args.input}")

    df = pd.read_csv(args.input)

    s = df['Mean_IGS_Size']
    summary = (
        f"n={len(s)}  median={s.median():.1f}  p90={s.quantile(.90):.1f}  "
        f"p95={s.quantile(.95):.1f}  p99={s.quantile(.99):.1f}  max={s.max():.1f}"
    )

    flagged = (
        df[df['Mean_IGS_Size'] > FLAG_THRESHOLD]
        .sort_values('Mean_IGS_Size', ascending=False)
        .reset_index(drop=True)
    )

    rows = []
    for _, r in flagged.iterrows():
        fasta = find_fasta(args.genomes_dir, r['File'])
        if fasta is None:
            contigs, gsize = 'N/A', 'N/A'
        else:
            c, g = fasta_stats(fasta)
            contigs = 'N/A' if c is None else c
            gsize = 'N/A' if g is None else g
        flag = 'likely artifact' if r['Mean_IGS_Size'] > ARTIFACT_THRESHOLD else ''
        rows.append({
            'Species': r['Species'],
            'File': r['File'],
            'Mean_IGS_Size': f"{r['Mean_IGS_Size']:.1f}",
            'genome_size_bp': gsize,
            'contigs': contigs,
            'flag': flag,
        })

    n_flagged = len(flagged)
    n_artifact = (flagged['Mean_IGS_Size'] > ARTIFACT_THRESHOLD).sum()

    lines = []
    lines.append('# Mean_IGS_Size outlier report')
    lines.append('')
    lines.append(f"Source: `{args.input}`")
    lines.append('')
    lines.append(f"**Distribution:** {summary}")
    lines.append('')
    lines.append(
        f"**Flagged:** {n_flagged} genomes with Mean_IGS_Size > {FLAG_THRESHOLD} bp; "
        f"{n_artifact} of these > {ARTIFACT_THRESHOLD} bp marked likely artifact."
    )
    lines.append('')
    lines.append(
        '| Species | File | Mean_IGS_Size | genome_size_bp | contigs | flag |'
    )
    lines.append('| --- | --- | ---: | ---: | ---: | --- |')
    for r in rows:
        lines.append(
            f"| {r['Species']} | {r['File']} | {r['Mean_IGS_Size']} | "
            f"{r['genome_size_bp']} | {r['contigs']} | {r['flag']} |"
        )

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    with open(args.output, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    print(f"Distribution: {summary}")
    print(f"Flagged (>{FLAG_THRESHOLD} bp): {n_flagged}  "
          f"likely artifact (>{ARTIFACT_THRESHOLD} bp): {n_artifact}")
    print(f"Wrote {args.output}")


if __name__ == '__main__':
    main()
