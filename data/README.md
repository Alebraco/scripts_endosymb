# Data Retrieval Scripts

This directory contains scripts for identifying, retrieving, and processing endosymbiont and related bacterial genomes from NCBI databases.

## Scripts

### genus.py
Filter bacterial genomes and identify endosymbionts based on metadata and quality criteria.
- **Input**: `genome_metadata.tsv` with genome metadata from NCBI
- **Output**: Filtered list of genomes including endosymbionts and genome size/gene ratio outliers
- **Filters**: CheckM completeness > 85%, contamination < 10%; identifies genome size outliers (30% below median)
- **Function**: Identifies endosymbionts by name and detects outliers within each genus
- **Dependencies**: pandas, numpy

### extract_16S.sh
Extract 16S ribosomal RNA sequences from annotated genomes.
- **Input**: Annotated genome GFF files with 16S rRNA annotations
- **Output**: 16S sequences in FASTA format (minimum 1200 bp)
- **Dependencies**: seqkit, awk
- **Usage**: Scans for '16S ribosomal RNA' features in GFF files

### blast_ids.py
BLAST 16S rRNA sequences against SILVA database for taxonomic identification.
- **Input**: 16S sequences from `clades_16S/` directory
- **Output**: BLAST results in TSV format with top hits
- **Database**: SILVA 16S rRNA reference database
- **Parameters**: E-value 1e-10, 8 threads
- **Dependencies**: BLAST+, subprocess

### parse_blast_results.py
Parse BLAST results and identify related species (non-endosymbionts).
- **Input**: BLAST output TSV files from `blast_results/`
- **Output**: Filtered species list with percent identity values
- **Function**: Excludes endosymbionts, Wolbachia, Buchnera, Blochmannia, and unclassified hits
- **Dependencies**: subprocess, awk

### gb2assm.py
Convert GenBank accession IDs to assembly accession IDs via NCBI Entrez.
- **Input**: GenBank accession IDs from BLAST results
- **Output**: RefSeq assembly accessions with percent identity mapping
- **Function**: Links nucleotide records to assembly database entries
- **Dependencies**: BioPython (Entrez), time (rate limiting)
- **Rate Limit**: 0.3 seconds between API calls

### final_cand_genomes.sh
Download candidate genome metadata and filter by quality criteria.
- **Input**: Assembly accession lists from `related_species/`
- **Output**: Filtered genome summary TSV files
- **Filters**: Minimum genome size (default: 2,500,000 bp)
- **Dependencies**: NCBI datasets CLI tool
- **Configuration**: SLURM batch system (2 days runtime, 10GB memory)

### unzip_candidates.sh
Unzip downloaded genome datasets and extract FASTA files.
- **Input**: Genome ZIP files from `candidate_genomes/`
- **Output**: Extracted FASTA genome files (.fna)
- **Function**: Extracts only FASTA files from NCBI dataset packages

### merge_candidates.sh
Merge candidate genomes from multiple sources into organized directories.
- **Input**: Genome FASTA files from `candidate_genomes/`
- **Output**: Merged genome collection in `merged_candidates/`
- **Function**: Consolidates genomes by species

## Typical Workflow

```bash
# 1. Filter species from metadata
python3 genus.py

# 2. Extract 16S sequences from annotated genomes
./extract_16S.sh

# 3. BLAST 16S against SILVA database
python3 blast_ids.py

# 4. Parse BLAST results to identify relatives
python3 parse_blast_results.py

# 5. Convert GenBank IDs to assembly IDs
python3 gb2assm.py

# 6. Download and filter candidate genomes
sbatch final_cand_genomes.sh

# 7. Extract genome files
./unzip_candidates.sh

# 8. Merge all candidates
./merge_candidates.sh
```

## Dependencies

- Python 3 with BioPython, pandas, numpy
- BLAST+ (blastn)
- seqkit
- NCBI datasets CLI tool
- SILVA 16S rRNA database
