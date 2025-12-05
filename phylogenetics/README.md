# Phylogenetics Scripts

This directory contains scripts for phylogenetic analysis, multiple sequence alignment, tree building, and core genome analysis.

## Multiple Sequence Alignment

### muscle.py
Wrapper for MUSCLE multiple sequence alignment on core genes.
- **Input**: Core gene protein sequences from `core/` directory
- **Output**: Aligned sequences in FASTA format (_aln.faa)
- **Method**: MUSCLE5 super5 iterative refinement algorithm
- **Dependencies**: MUSCLE v5+
- **Usage**: `python3 muscle.py <group>`

### run_muscle.sh
Batch workflow for running MUSCLE alignments on all species groups.
- **Input**: Group name (endosymb_only, endosymb+relatives, relatives_only)
- **Output**: Aligned core gene sequences in `core_alignments/` directory
- **Function**: Executes muscle.py for specified group

## Backtranslation

### backtranslate.py
Backtranslate protein alignments to nucleotide sequences.
- **Input**: Protein alignments from `{group}/core_alignments/` and CDS from `{group}/cds/`
- **Output**: Backtranslated DNA alignments in `{group}/backtranslated/` directory
- **Function**: Maps aligned protein sequences back to original CDS nucleotide sequences
- **Usage**: `python3 backtranslate.py <group>`
- **Dependencies**: regex for parsing protein IDs

### run_backtranslate.sh
Batch workflow for backtranslating aligned protein sequences.
- **Input**: Group name
- **Output**: Nucleotide alignments preserving protein alignment gaps
- **Function**: Executes backtranslate.py for specified group

## Concatenation

### concatenate.py
Concatenate aligned core gene sequences into supermatrices.
- **Input**: Backtranslated alignments from `backtranslated/` directory
- **Output**: Concatenated alignments in `dna_concatenates/` directory
- **Function**: Merges multiple gene alignments into single sequence per genome
- **Usage**: `python3 concatenate.py <group>`
- **Dependencies**: BioPython (SeqIO, AlignIO)

### run_concatenate.sh
Batch workflow for concatenating aligned sequences.
- **Input**: Group name
- **Output**: Supermatrix FASTA files ready for phylogenetic inference
- **Function**: Executes concatenate.py for specified group

## Tree Building

### dnatree.sh
Build phylogenetic trees from DNA sequence concatenates.
- **Dependencies**: IQ-TREE
- **Input**: Concatenated DNA alignments from `dna_concatenates/`
- **Output**: Phylogenetic trees in Newick format (.treefile)
- **Model**: GTR+G with 1000 bootstrap replicates
- **Configuration**: LSF batch system array jobs

## Core Genome Analysis

### corecruncher.sh
Analyze core genome composition across taxa using CoreCruncher.
- **Input**: Protein sequences from `proteins/` directory
- **Output**: Core gene sets in `core/` directory
- **Function**: Identifies genes shared across all genomes in a species
- **Dependencies**: CoreCruncher (USEARCH-based)
- **Configuration**: LSF batch system (48 hour runtime)

### count_core.sh
Count core genes shared across genomes.
- **Input**: Core gene output from CoreCruncher
- **Output**: Core gene count statistics
- **Function**: Reports number of core genes per species

### count_genes.py
Count total genes across all genomes.
- **Input**: Annotated protein files (.faa)
- **Output**: Gene count statistics per genome
- **Function**: Parses FASTA headers and counts entries

### count_genes_sp.py
Count genes per species (aggregated by species).
- **Input**: Protein sequences organized by species
- **Output**: Species-level gene count summaries
- **Function**: Aggregates gene counts across genomes within each species

### count_hits.py
Count and summarize core gene hits from pairwise comparisons.
- **Input**: CoreCruncher pairwise comparison results from `CC/` directory
- **Output**: Hit count statistics with percent identity metrics
- **Function**: Analyzes pairwise BLAST results to identify core genes
- **Usage**: `python3 count_hits.py <group>`
- **Dependencies**: pandas

## Typical Workflow

```bash
# 1. Identify core genes
bsub < corecruncher.sh

# 2. Align core genes
./run_muscle.sh endosymb_only

# 3. Backtranslate to DNA
./run_backtranslate.sh endosymb_only

# 4. Concatenate alignments
./run_concatenate.sh endosymb_only

# 5. Build phylogenetic trees
bsub < dnatree.sh

# 6. Count core genes and hits
./count_core.sh
python3 count_hits.py endosymb_only
```

## Notes

- Scripts are configured for HPC environments with LSF batch system
- Adjust resource requirements (-R, -W) based on dataset size
- Tree building can be computationally intensive for large datasets
- Core gene identification requires USEARCH or CoreCruncher installed
