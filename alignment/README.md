# Alignment Scripts

This directory contains scripts for multiple sequence alignment and processing.

## Scripts

### backtranslate.py
Backtranslate protein alignments to nucleotide alignments using original CDS sequences.
- **Usage**: `python3 backtranslate.py <group>`
- **Input**: Protein alignments and CDS sequences
- **Output**: Nucleotide alignments preserving codon structure

### muscle.py
MUSCLE multiple sequence alignment wrapper.
- **Input**: Unaligned sequences (FASTA)
- **Output**: Aligned sequences (FASTA)
- **Dependencies**: MUSCLE aligner

### concatenate.py
Concatenate multiple aligned gene sequences into a supermatrix.
- **Usage**: `python3 concatenate.py <group>`
- **Input**: Individual gene alignments
- **Output**: Concatenated alignment file

### Workflow Scripts

#### run_backtranslate.sh
Batch workflow for backtranslation of multiple groups.
- Processes: endosymb_only, endosymb+relatives, relatives_only
- Configured for LSF batch system

#### run_concatenate.sh
Batch workflow for concatenating alignments across groups.
- Processes all groups in parallel
- Array job configuration

#### run_muscle.sh
Batch MUSCLE alignment workflow.
- Aligns multiple sequence files
- Parallel processing support

## Usage Example

```bash
# Align sequences
./run_muscle.sh

# Backtranslate protein alignments
./run_backtranslate.sh

# Concatenate alignments
./run_concatenate.sh

# Or run individually
python3 muscle.py
python3 backtranslate.py endosymb_only
python3 concatenate.py endosymb_only
```

## Workflow

1. Align protein sequences with MUSCLE
2. Backtranslate to nucleotides
3. Concatenate genes for phylogenetic analysis
