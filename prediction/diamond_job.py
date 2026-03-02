import os
import re
import subprocess
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir



def run_diamond(path, threads = 8, max_parallel = 20, wait=True):
    '''
    Run DIAMOND using LSF array jobs
    '''
    # The path must be the directory containing the 'proteins' folder
    TRANS_DB = os.path.join(files_dir, 'IS_db_fixed.dmnd')
    proteins_root = os.path.join(path, 'proteins')

    tasks = []
    for species in os.listdir(proteins_root):
        species_path = os.path.join(proteins_root, species)
        if not os.path.isdir(species_path):
            continue

        out_dir = os.path.join(path, 'transposase', species)
        os.makedirs(out_dir, exist_ok=True)

        for filename in os.listdir(species_path):
            if filename.endswith('.faa'):
                accession = filename.split('.faa')[0]
                faa_file = os.path.join(species_path, filename)

                outfile = os.path.join(out_dir, f'{accession}.tsv')

                if not os.path.isfile(outfile):  # skip already completed
                    tasks.append((faa_file, outfile))
    if not tasks:
        print('No new DIAMOND jobs to run.')
        return
    
    task_file = os.path.join(path, 'diamond_tasks.tsv')
    with open(task_file, 'w') as f:
        for faa_file, outfile in tasks:
            f.write(f'{faa_file}\t{outfile}\n')
    
    n_jobs = len(tasks)
    job_name = f'transposase_batch[1-{n_jobs}]%{max_parallel}'
    logs_dir = os.path.join(path, 'transpred_logs')
    os.makedirs(logs_dir, exist_ok=True)

    cmd = (
        f'line=$(sed -n "${{LSB_JOBINDEX}}p" {task_file}); '
        'query=$(echo "$line" | cut -f1); '
        'out=$(echo "$line" | cut -f2); '
        f'diamond blastp --db {TRANS_DB} --query "$query" --out "$out" '
        '--outfmt 6 qseqid sseqid pident length mismatch gapopen '
        'qstart qend sstart send evalue bitscore qlen slen '
        f'--threads {threads} --evalue 1e-5 --max-target-seqs 1 --quiet'
    )
    
    bsub_cmd = [
        'bsub', 
        '-n', str(threads),
        '-q', 'bobay',
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

if __name__ == '__main__':
    path = 'endosymb+relatives'
    run_diamond(path)
    print('Done with DIAMOND searches of transposases.')