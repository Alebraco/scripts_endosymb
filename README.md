# Endosymbiont Genomics Analysis Scripts

A collection of research scripts for endosymbiont genomics analysis, including genome annotation, sequence alignment, phylogenetic analysis, and comparative genomics.

## Overview

This repository contains bioinformatics scripts developed for analyzing endosymbiont genomes and comparing them with free-living relatives. The scripts support various analyses including:

- Genome annotation and feature extraction
- Multiple sequence alignment and backtranslation
- Phylogenetic tree construction
- GC content and codon usage analysis
- Intergenic spacer (IGS) analysis
- Comparative genomics and clustering

## Directory Structure

```
.
├── annotation/          # Genome annotation and feature extraction
│   ├── bakta_annotation.sh      # Batch genome annotation using Bakta
│   ├── extract_cds.sh           # Extract CDS sequences from genomes
│   ├── endosymb_cds.py          # Endosymbiont-specific CDS extraction
│   ├── gb2assm.py               # Convert GenBank to assembly format
│   ├── cand_genomes.sh          # Identify candidate genomes
│   └── remove_genomes.py        # Filter/remove genomes
│
├── alignment/           # Sequence alignment and processing
│   ├── backtranslate.py         # Backtranslate protein alignments to nucleotide
│   ├── muscle.py                # MUSCLE alignment wrapper
│   ├── concatenate.py           # Concatenate aligned sequences
│   ├── run_backtranslate.sh     # Batch backtranslation workflow
│   ├── run_concatenate.sh       # Batch concatenation workflow
│   └── run_muscle.sh            # Batch alignment workflow
│
├── analysis/            # Genomic analysis and statistics
│   ├── gc_content.py            # Calculate GC content
│   ├── igs_lengths.py           # Analyze intergenic spacer lengths
│   ├── igs_plots.py             # Generate IGS visualization plots
│   ├── spectrum.py              # Analyze mutation/substitution spectra
│   ├── count_genes.py           # Count genes in genomes
│   ├── count_genes_sp.py        # Count genes per species
│   ├── count_hits.py            # Count BLAST hits
│   ├── count_core.sh            # Count core genes
│   ├── parse_blast_results.py   # Parse BLAST output
│   ├── percid_distribution.py   # Percent identity distributions
│   ├── pid_calculation.py       # Calculate percent identity
│   ├── usearch.py               # USEARCH analysis wrapper
│   └── run_*.sh                 # Batch analysis workflows
│
├── phylogenetics/       # Phylogenetic analysis
│   ├── clusters.py              # Hierarchical clustering of genomes
│   ├── dnatree.sh               # DNA-based phylogenetic trees
│   ├── iqtree.sh                # IQ-TREE phylogenetic inference
│   ├── longbranch_outliers.py   # Detect long-branch artifacts
│   ├── corecruncher.sh          # Core genome analysis
│   └── new_path_precrunchr.sh   # Pre-processing for core analysis
│
├── visualization/       # Data visualization
│   └── variation_plot.py        # Plot genetic variation
│
├── workflows/           # High-level analysis workflows
│   └── endosymb_core.sh         # Core endosymbiont analysis workflow
│
├── utils/               # Utility scripts
│   ├── merge_dirs.sh            # Merge directory contents
│   ├── merge_relatives.sh       # Merge relative genome data
│   ├── move_relatives.sh        # Move relative genome files
│   └── unzip_candidates.sh      # Extract candidate genome archives
│
└── legacy_*/            # Legacy script collections (preserved for reference)
    ├── legacy_GCSize_pipeline/
    └── legacy_brccluster_scripts/
```

## Requirements

### Software Dependencies
- Python 3.7+
- Bakta (genome annotation)
- MUSCLE (multiple sequence alignment)
- IQ-TREE (phylogenetic inference)
- BLAST+ (sequence similarity search)
- USEARCH (sequence analysis)

### Python Packages
See `requirements.txt` for Python dependencies. Install with:
```bash
pip install -r requirements.txt
```

Common packages used:
- BioPython
- NumPy
- Pandas
- Matplotlib
- Seaborn
- SciPy

## Usage

### Genome Annotation
```bash
# Annotate genomes using Bakta
cd annotation/
./bakta_annotation.sh

# Extract CDS sequences
./extract_cds.sh
```

### Sequence Alignment
```bash
# Align sequences and backtranslate
cd alignment/
./run_muscle.sh
./run_backtranslate.sh
./run_concatenate.sh
```

### Phylogenetic Analysis
```bash
# Build phylogenetic trees
cd phylogenetics/
./dnatree.sh

# Cluster genomes
python3 clusters.py
```

### Comparative Analysis
```bash
# Calculate GC content
cd analysis/
python3 gc_content.py <alignment_file>

# Analyze intergenic spacers
python3 igs_lengths.py
python3 igs_plots.py
```

## Workflow Examples

### Complete Endosymbiont Analysis Pipeline
1. **Annotation**: Extract and annotate genomes
2. **Alignment**: Align core genes across genomes
3. **Concatenation**: Concatenate alignments
4. **Phylogenetics**: Build phylogenetic trees
5. **Analysis**: Calculate statistics (GC%, IGS, etc.)
6. **Visualization**: Generate plots

See `workflows/endosymb_core.sh` for an example pipeline.

## Input Data

Scripts expect input data in standard bioinformatics formats:
- Genome assemblies: FASTA (.fna, .fasta)
- Annotations: GenBank (.gb, .gbff)
- Alignments: FASTA format
- Trees: Newick format (.treefile)

## Output

Outputs vary by script but typically include:
- Annotated genomes (GenBank format)
- Aligned sequences (FASTA)
- Phylogenetic trees (Newick)
- Summary statistics (CSV)
- Plots (PDF, PNG)

## HPC Environment

Many scripts include LSF batch system headers (BSUB) for running on HPC clusters. These may need adjustment for your computing environment.

## Legacy Scripts

The `legacy_*` directories contain older versions of scripts organized differently. These are preserved for reference but the reorganized structure should be used for new analyses.

## Citation

If you use these scripts in your research, please cite:
[Add citation information when available]

## Contact

For questions or issues, please contact the repository maintainer or open an issue on GitHub.

## License

[Add license information]
