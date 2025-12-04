# Endosymbiont Genomics Analysis Scripts

Research scripts for endosymbiont genomics analysis, including endosymbiont identification, genome annotation, phylogenetic analysis, and comparative genomics.

## Directory Structure

```
.
├── data/            # Data Retrieval
│   ├── genus.py                 # Selects endosymbiont species  
│   ├── extract_16S.sh           # Scans genomes for 16S gene
│   ├── blast_ids.py             # BLAST 16S RNA against SILVA database
│   ├── parse_blast_results.py   # Parse BLAST results
│   ├── gb2assm.py               # Convert GenBank to assembly format
│   ├── final_cand_genomes.sh    # Download candidate genomes
│   ├── unzip_candidates.sh      # Unzip candidate genomes
│   └── merge_candidates.sh      # Merge candidate genomes
│
├── annotation/          # Genome annotation and feature extraction
│   ├── bakta_annotation.sh      # Batch genome annotation using Bakta
│   ├── move_relatives.sh        # Move annotated protein files
│   └── merge_dirs.sh            # Create new group (endosymb+relatives)
│
├── phylogenetics/       # Phylogenetic analysis
│   ├── corecruncher.sh         # Get core genes (modify for each group)
│   ├── muscle.py                # MUSCLE alignment wrapper
│   ├── backtranslate.py         # Backtranslate protein alignments to nucleotide
│   ├── concatenate.py           # Concatenate aligned sequences
│   ├── dnatree.sh               # DNA-based phylogenetic trees
│   ├── run_backtranslate.sh     # Batch backtranslation workflow
│   ├── run_concatenate.sh       # Batch concatenation workflow
│   ├── run_muscle.sh            # Batch alignment workflow
│   ├── count_hits.py            # Count core gene hits per group
│   ├── count_genes.py           # Count genes in genomes
│   ├── count_genes_sp.py        # Count genes per species
│   └── count_core.sh            # Count core genes
│
├── analysis/            # Genomic analysis and statistics
│   ├── gc_calculate.py            # Calculate GC content
│   ├── delta_matrix.py            # Create pairwise difference matrices
│   ├── distance_matrix.py         # Extract phylogenetic distances from trees
│   ├── gc_codon_dict.py           # Compute fourfold degenerate GC content
│   ├── gcsize_dict.py             # Compute genome size and GC content
│   ├── matrix_correlation.py      # Mantel test correlations between matrices
│   ├── utils.py                   # Utility functions and path helpers
│   ├── igs_lengths.py             # Analyze intergenic spacer lengths
│   ├── igs_plots.py               # Generate IGS visualization plots
│   ├── spectrum.py                # Analyze mutation/substitution spectra
│   ├── percid_distribution.py     # Percent identity distributions
│   ├── pid_calculation.py         # Calculate percent identity
│   ├── run_*.sh                   # Various wrappers for python scripts
│   ├── run_transposase.sh         # Transposase Workflow
│   └── run_pseudo.sh              # PseudoFinder Workflow
│
├── visualization/       # Data visualization
│   ├── absolute_scatter.py        # Scatter plots for GC/size across groups
│   ├── clockwise_plot.py          # Normalized metric changes vs. divergence
│   ├── main_analysis.py           # Comprehensive correlation analysis plots
│   ├── summary_across_groups.py   # Compare delta values across groups
│   ├── summary_data.py            # Group-specific summary plots
│   └── variation_plot.py          # Plot genetic variation
│
└── legacy_*/            # Legacy script collections (preserved for reference)
    └── legacy_brccluster_scripts/
```

## Requirements

### Software Dependencies
- Python 3
- Bakta (genome annotation)
- MUSCLE (multiple sequence alignment)
- IQ-TREE (phylogenetic inference)
- BLAST+ (sequence similarity search)
- USEARCH (sequence alignment)
- CoreCruncher (core genome)

### Python Packages
See `requirements.txt` for Python dependencies. Install with:
```bash
pip install -r requirements.txt
```

- BioPython
- NumPy
- Pandas
- Matplotlib
- Seaborn
- SciPy

## HPC 

Some scripts include LSF system headers (BSUB) for running on HPC clusters. 
