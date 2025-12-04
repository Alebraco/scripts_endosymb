# Annotation Scripts

This directory contains scripts for genome annotation and feature extraction.

## Scripts

### bakta_annotation.sh
Batch genome annotation using Bakta. Configured for LSF batch system with array jobs.
- **Input**: FASTA genome files (.fna)
- **Output**: Annotated genomes in GenBank format
- **Dependencies**: Bakta, conda environment

### extract_cds.sh
Extract coding sequences (CDS) from annotated genomes.
- **Input**: Annotated genome files
- **Output**: CDS FASTA files

### endosymb_cds.py
Specialized CDS extraction for endosymbiont genomes.
- **Input**: JSON configuration and genome files
- **Output**: Endosymbiont-specific CDS sequences

### gb2assm.py
Convert GenBank files to assembly format.
- **Input**: GenBank files (.gb, .gbff)
- **Output**: Assembly format files

### cand_genomes.sh
Identify and process candidate genomes for analysis.
- **Input**: Genome directories
- **Output**: Filtered candidate genome list

### remove_genomes.py
Filter and remove genomes based on quality criteria.
- **Input**: Genome list and quality metrics
- **Output**: Filtered genome set

## Usage Example

```bash
# Annotate all genomes
./bakta_annotation.sh

# Extract CDS sequences
./extract_cds.sh

# Process for endosymbiont analysis
python3 endosymb_cds.py
```
