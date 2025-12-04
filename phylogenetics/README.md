# Phylogenetics Scripts

This directory contains scripts for phylogenetic analysis, tree building, and genome clustering.

## Tree Building

### dnatree.sh
Build phylogenetic trees from DNA sequence concatenates.
- **Dependencies**: IQ-TREE
- **Input**: Concatenated DNA alignments (.fasta)
- **Output**: Phylogenetic trees in Newick format
- **Model**: GTR+G with 1000 bootstrap replicates
- **Configuration**: LSF batch system array jobs

### iqtree.sh
General IQ-TREE phylogenetic inference.
- **Model**: MFP (ModelFinder Plus) for automatic model selection
- **Bootstrap**: 1000 ultrafast bootstrap replicates
- **Input**: Protein or DNA concatenates
- **Output**: Tree files, bootstrap support values

## Clustering

### clusters.py
Hierarchical clustering of genomes based on phylogenetic distances.
- **Input**: Phylogenetic tree (.treefile)
- **Output**: Cluster assignments (CSV)
- **Method**: Average linkage clustering
- **Dependencies**: BioPython, scipy, pandas
- **Default**: Distance threshold = 0.9

## Core Genome Analysis

### corecruncher.sh
Analyze core genome composition across taxa.
- **Input**: Gene presence/absence data
- **Output**: Core gene sets

### new_path_precrunchr.sh
Pre-processing pipeline for core genome analysis.
- **Input**: Raw genome data
- **Output**: Processed data ready for core analysis

## Outlier Detection

### longbranch_outliers.py
Detect and identify long-branch artifacts in phylogenetic trees.
- **Input**: Phylogenetic tree
- **Output**: List of potential long-branch outliers
- **Purpose**: Quality control for phylogenetic analyses

## Usage Examples

```bash
# Build DNA trees for multiple groups
./dnatree.sh

# Build a single tree with automatic model selection
./iqtree.sh

# Cluster genomes based on tree
python3 clusters.py

# Analyze core genome
./corecruncher.sh

# Detect outliers
python3 longbranch_outliers.py
```

## Workflow

1. **Pre-process** genomes (new_path_precrunchr.sh)
2. **Identify core genes** (corecruncher.sh)
3. **Build trees** (dnatree.sh or iqtree.sh)
4. **Cluster genomes** (clusters.py)
5. **QC check** (longbranch_outliers.py)

## Notes

- Scripts are configured for HPC environments with LSF batch system
- Adjust resource requirements (-R, -W) based on dataset size
- Tree building can be computationally intensive for large datasets
