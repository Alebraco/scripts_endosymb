# Analysis Scripts

This directory contains scripts for genomic analysis and statistics.

## GC Content Analysis

### gc_content.py
Calculate GC content statistics from aligned sequences.
- **Usage**: `python3 gc_content.py <alignment_file>`
- **Output**: Min, average, and max GC% statistics

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

### run_igs.py
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
./run_igs.py
./run_percid.sh
./run_spectrum.sh
```
