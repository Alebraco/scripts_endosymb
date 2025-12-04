#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem 50G
#SBATCH -N 1
#SBATCH -e %j.e
#SBATCH -o %j.o
#SBATCH -p standard

datasets summary genome taxon 2 --as-json-lines --assembly-source 'RefSeq'| dataformat tsv genome --fields accession,organism-name,assmstats-total-sequence-len,assminfo-level,annotinfo-featcount-gene-protein-coding,assmstats-contig-n50,annotinfo-busco-complete,checkm-marker-set,checkm-completeness,checkm-contamination
