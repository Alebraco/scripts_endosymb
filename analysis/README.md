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
- **Output**: JSON file with nested dictionary {species: {genome_filename: {size, gc_genome}}}
- **Note**: Keys are genome filenames without .fna extension, not necessarily accession IDs
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

## Pseudogene Analysis

### run_pseudo.sh
Analyze pseudogene features.
- **Dependencies**: PseudoFinder

## Transposase Analysis

### run_transposase.sh
Analyze transposase elements.