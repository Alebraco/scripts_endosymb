# Annotation Scripts

This directory contains scripts for genome annotation and organizing annotated genome files.

## Scripts

### bakta_annotation.sh
Batch genome annotation using Bakta. Configured for LSF batch system with array jobs.
- **Input**: FASTA genome files (.fna) from genome directories
- **Output**: Annotated genomes in GenBank format (.gbff), protein files (.faa), CDS files (.ffn), and genome files (.fna)
- **Dependencies**: Bakta, conda environment
- **Configuration**: LSF batch system with array jobs for parallel processing

### move_relatives.sh
Organize Bakta annotation results into structured directories for relative species.
- **Input**: Bakta annotation results from `bakta_results/` directory
- **Output**: Organized directory structure in `relatives_only/` with subdirectories for genomes, proteins, gbff_files, and CDS
- **Function**: Creates species-specific directories and copies annotation files to appropriate subdirectories
- **Configuration**: LSF batch system (2 hour runtime, 1 CPU)
- **Dependencies**: bash, find

### merge_dirs.sh
Merge endosymbiont-only annotations with endosymbiont+relatives dataset.
- **Input**: GenBank files from `endosymb_only/gbff_files/` and `endosymb+relatives/gbff_files/`
- **Output**: Combined dataset in `endosymb+relatives/gbff_files/` with both endosymbionts and relatives
- **Function**: Merges annotation files for species present in both directories
- **Use case**: Creating combined datasets for comparative genomics analysis

## Workflow

```bash
# 1. Annotate all genomes with Bakta
bsub < bakta_annotation.sh

# 2. Organize relative species annotations
bsub < move_relatives.sh

# 3. Merge endosymbiont and relative annotations
./merge_dirs.sh
```

## Directory Structure

After running these scripts, the output will be organized as:

```
relatives_only/
├── species_name/
│   ├── genomes/       # Genome FASTA files (.fna)
│   ├── proteins/      # Protein sequences (.faa)
│   ├── gbff_files/    # GenBank format files (.gbff)
│   └── cds/           # Coding sequences (.fna)

endosymb+relatives/
└── gbff_files/
    └── species_name/  # Combined GenBank files
```

## Dependencies

- Bakta (genome annotation tool)
- LSF or SLURM batch system (for array jobs)
- bash with standard utilities (find, cp)
