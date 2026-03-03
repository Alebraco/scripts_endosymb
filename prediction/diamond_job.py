#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
import shlex
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir


def run_diamond(path, threads=8, max_parallel=20, wait=True, queue='bobay'):
    '''Run DIAMOND via a single LSF array job based on discovered .faa files.'''

    path = os.path.abspath(path)
    TRANS_DB = os.path.join(files_dir, 'IS_db_fixed.dmnd')
    proteins_root = os.path.join(path, 'proteins')

    if not os.path.isdir(proteins_root):
        raise FileNotFoundError(f'Proteins directory not found: {proteins_root}')

    if not os.path.isfile(TRANS_DB):
        raise FileNotFoundError(f'DIAMOND database not found: {TRANS_DB}')

    tasks = []
    for species in os.listdir(proteins_root):
        species_path = os.path.join(proteins_root, species)
        if not os.path.isdir(species_path):
            continue

        out_dir = os.path.join(path, 'transposase', species)
        os.makedirs(out_dir, exist_ok=True)

        for filename in os.listdir(species_path):
            if filename.endswith('.faa'):
                accession = filename.rsplit('.faa', 1)[0]
                faa_file = os.path.join(species_path, filename)

                outfile = os.path.join(out_dir, f'{accession}.tsv')

                if not os.path.isfile(outfile):  # skip already completed
                    tasks.append((faa_file, outfile))

    print(f'Found {len(tasks)} pending DIAMOND tasks.')
    if not tasks:
        print('No new DIAMOND jobs to run.')
        return {'submitted': 0, 'job_id': None}
    
    task_file = os.path.join(os.getcwd(), 'diamond_tasks.tsv')
    with open(task_file, 'w') as f:
        for faa_file, outfile in tasks:
            f.write(f'{faa_file}\t{outfile}\n')
    
    n_jobs = len(tasks)
    job_name = f'transposase_batch[1-{n_jobs}]%{max_parallel}'
    logs_dir = os.path.join(os.getcwd(), 'transpred_logs')
    os.makedirs(logs_dir, exist_ok=True)

    quoted_task_file = shlex.quote(task_file)
    quoted_db = shlex.quote(TRANS_DB)

    cmd = (
        f'line=$(sed -n "${{LSB_JOBINDEX}}p" {quoted_task_file}); '
        'query=$(echo "$line" | cut -f1); '
        'out=$(echo "$line" | cut -f2); '
        f'diamond blastp --db {quoted_db} --query "$query" --out "$out" '
        '--outfmt 6 qseqid sseqid pident length mismatch gapopen '
        'qstart qend sstart send evalue bitscore qlen slen '
        f'--threads {threads} --evalue 1e-5 --max-target-seqs 1 --quiet'
    )
    
    bsub_cmd = [
        'bsub', 
        '-n', str(threads),
        '-q', queue,
        '-R', f'rusage[mem=8GB] span[hosts=1]',
        '-J', job_name,
        '-e', os.path.join(logs_dir, '%J_%I.err'),
        '-o', os.path.join(logs_dir, '%J_%I.out'),
        cmd,
    ]

    result = subprocess.run(bsub_cmd, check=True, capture_output=True, text=True)

    match = re.search(r'Job <(\d+)>', result.stdout)
    if not match:
        raise ValueError(f'Failed to parse Job ID from bsub output: {result.stdout}')
    
    job_id = match.group(1)
    print(f'Submitted {n_jobs} DIAMOND jobs with Job ID {job_id}.')

    if wait:
        print('Waiting for DIAMOND jobs to finish.')
        subprocess.run(['bwait', '-w', f'done({job_id})'], check=True)

    return {'submitted': n_jobs, 'job_id': job_id}


def main():
    parser = argparse.ArgumentParser(description='Submit DIAMOND LSF array jobs for transposase searches.')
    parser.add_argument('--path', required=True, help='Group path containing proteins/ and transposase/ dirs')
    parser.add_argument('--threads', type=int, default=8, help='Threads per DIAMOND array task')
    parser.add_argument('--max-parallel', type=int, default=20, help='Maximum concurrent array tasks')
    parser.add_argument('--queue', default='bobay', help='LSF queue name')
    parser.add_argument('--no-wait', action='store_true', help='Do not wait for array completion')
    parser.add_argument('--done-file', default=None, help='Optional marker file to write at end')
    args = parser.parse_args()

    result = run_diamond(
        path=args.path,
        threads=args.threads,
        max_parallel=args.max_parallel,
        wait=not args.no_wait,
        queue=args.queue,
    )

    if args.done_file:
        with open(args.done_file, 'w') as f:
            f.write(f"submitted={result['submitted']}\n")
            f.write(f"job_id={result['job_id']}\n")


if __name__ == '__main__':
    main()