# Analysis Scripts

This directory contains scripts for genomic analysis and statistics.

## GC Content and Genome Size Analysis

### gc_calculate.py
Core utility for calculating GC content from sequences.
- **Function**: `calculate_gc_content(seq)` - Returns GC percentage
- **Use case**: Imported by other GC analysis scripts

### gcsize_dict.py
Compute genome-wide size and GC content metadata.
- **Function**: `genome_gcsize(group)` - Process genomes and compute metrics
- **Input**: Genome FASTA files from group directories (endosymb_only, endosymb+relatives, relatives_only)
- **Output**: JSON file with nested dictionary {species: {accession: {size, gc_genome}}}
- **Dependencies**: BioPython (SeqIO), gc_calculate.py

### gc_codon_dict.py
Create dictionaries of GC content for codon analysis.
- **Function**: Computes fourfold degenerate and overall GC content
- **Input**: Core gene alignment concatenates
- **Output**: JSON file with GC content data

### delta_matrix.py
Create delta matrices for pairwise differences in genomic metrics.
- **Function**: `delta_matrix(data, mat_type)` - Compute absolute pairwise differences
- **Input**: Nested dictionary {species: {genome_id: metric_value}}
- **Output**: Symmetric delta matrix (pandas DataFrame) for each species
- **Use case**: Calculate pairwise differences in GC content or genome size
- **Dependencies**: pandas, numpy

### distance_matrix.py
Create phylogenetic distance matrices from trees.
- **Function**: `distance_matrix(group)` - Extract patristic distances
- **Input**: Phylogenetic trees in Newick format from dna_tree_results/
- **Output**: Distance matrices (pandas DataFrame) for each species
- **Method**: Computes pairwise patristic distances between all terminals
- **Dependencies**: BioPython (Phylo), pandas

### matrix_correlation.py
Compute correlations between two distance matrices.
- **Function**: `matrix_correlation(x_df, y_df)` - Mantel test between matrices
- **Method**: Spearman correlation with 10,000 permutations
- **Output**: Tuple of (correlation coefficient, p-value, sample size)
- **Minimum**: Requires at least 3 genomes per species
- **Dependencies**: scikit-bio (mantel test)

### utils.py
Utility functions for the analysis pipeline.
- **Contains**: Title mappings, group names, JSON I/O helpers, file paths
- **Functions**: `load_or_compute()`, `load_or_compute_pickle()`, path helpers
- **Use case**: Imported by all analysis and visualization scripts

## Intergenic Spacer (IGS) Analysis

### igs_lengths.py
Analyze intergenic spacer lengths across genomes.
- **Input**: Genome annotations
- **Output**: IGS length distributions

### igs_plots.py
Generate visualization plots for IGS analysis.
- **Input**: IGS data CSV files
- **Output**: Box plots and scatter plots (PDF format)
- **Dependencies**: matplotlib, seaborn

### run_igs.sh
Batch workflow for IGS analysis.

## Gene Counting

### count_genes.py
Count genes across genomes.
- **Input**: Annotated genomes
- **Output**: Gene count statistics

### count_genes_sp.py
Count genes per species.
- **Input**: Species-organized genome annotations
- **Output**: Species-level gene counts

### count_core.sh
Count core genes shared across genomes.
- **Input**: Gene presence/absence data
- **Output**: Core gene statistics

## BLAST Analysis

### parse_blast_results.py
Parse BLAST output and extract relevant hits.
- **Input**: BLAST output files
- **Output**: Parsed hit tables

### count_hits.py
Count and summarize BLAST hits.
- **Input**: BLAST results
- **Output**: Hit count statistics

### usearch.py
USEARCH sequence similarity search wrapper.
- **Dependencies**: USEARCH
- **Input**: Query and target sequences
- **Output**: Similarity search results

### run_usearch.sh
Batch USEARCH workflow.

## Sequence Identity Analysis

### percid_distribution.py
Analyze percent identity distributions.
- **Input**: Alignment files
- **Output**: Identity distribution plots and statistics

### pid_calculation.py
Calculate pairwise percent identity.
- **Input**: Aligned sequences
- **Output**: Pairwise identity matrix

### run_percid.sh
Batch percent identity analysis workflow.

## Mutation Spectrum

### spectrum.py
Analyze mutation/substitution spectra.
- **Input**: Aligned sequences
- **Output**: Mutation spectrum analysis
- **Dependencies**: BioPython

### run_spectrum.sh
Batch spectrum analysis workflow.

## Specialized Analysis

### run_pseudo.sh
Analyze pseudogene features.

### run_transposase.sh
Analyze transposase elements.

### run_summary_data.sh
Generate summary statistics across analyses.

## Usage Examples

```bash
# Calculate GC content
python3 gc_content.py alignment.fasta

# Analyze IGS
python3 igs_lengths.py
python3 igs_plots.py

# Count genes
python3 count_genes.py
./count_core.sh

# Run complete workflows
./run_igs.sh
./run_percid.sh
./run_spectrum.sh
```
