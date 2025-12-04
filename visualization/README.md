# Visualization Scripts

This directory contains scripts for generating plots and visualizations of genomic data, particularly GC content, genome size, and genetic variation.

## Scripts

### absolute_scatter.py
Generate absolute scatter plots for GC content and genome size across groups.
- **Usage**: `python3 absolute_scatter.py <type> <mean/all>`
- **Input**: Genome metadata (GC content or size) from analysis/gcsize_dict.py
- **Output**: Scatter plots comparing endosymbionts vs. relatives (PDF format)
- **Types**: 'gc_genome' for GC content, 'size' for genome size
- **Dependencies**: matplotlib, seaborn, pandas, numpy

### clockwise_plot.py
Generate clockwise correlation plots for GC and genome size changes vs. phylogenetic distance.
- **Usage**: `python3 clockwise_plot.py <matrix_type>`
- **Input**: Delta matrices and phylogenetic distances
- **Output**: Clockwise plots showing normalized metric changes per unit divergence
- **Function**: Calculates (Delta / Max Size of Pair) / Patristic Distance
- **Types**: 'gc_genome' or 'size'
- **Dependencies**: matplotlib, seaborn, pandas, numpy

### main_analysis.py
Comprehensive GC and genome size correlation analysis with visualization.
- **Usage**: `python3 main_analysis.py <group> <matrix1_type> <matrix2_type> <mean/matrix>`
- **Input**: GC data, distance matrices, genome metadata from analysis/
- **Output**: Correlation plots, Mantel test statistics, scatter plots (PDF format)
- **Matrix Types**: 'gc_genome', 'size', 'distance'
- **Function**: Correlates genomic metrics with phylogenetic distances
- **Dependencies**: matplotlib, seaborn, pandas, numpy, scipy

### summary_across_groups.py
Generate summary plots comparing delta values across different groups.
- **Usage**: `python3 summary_across_groups.py <matrix_type> <all/mean>`
- **Input**: GC or genome size data across all groups
- **Output**: Comparative box plots or violin plots (PDF format)
- **Function**: Visualizes distribution of pairwise differences across endosymbionts and relatives
- **Types**: 'gc_genome' or 'size'
- **Dependencies**: matplotlib, seaborn, pandas, numpy

### summary_data.py
Generate summary data plots for specific groups.
- **Usage**: `python3 summary_data.py <group> <matrix_type> <sp/all>`
- **Input**: Group-specific GC or genome size data
- **Output**: Summary visualization plots (PDF format)
- **Function**: Displays delta matrix distributions within a group
- **Dependencies**: matplotlib, seaborn, pandas, numpy

### variation_plot.py
Generate plots for genetic variation analysis.
- **Input**: Variation data (CSV or similar format)
- **Output**: Visualization plots (PDF/PNG)
- **Dependencies**: matplotlib, seaborn
- **Use case**: Visualizing genetic diversity and variation patterns

## Usage Examples

```bash
# Genetic variation plot
python3 variation_plot.py

# GC content scatter plots (mean values)
python3 absolute_scatter.py gc_genome mean

# Genome size scatter plots (all genomes)
python3 absolute_scatter.py size all

# Clockwise plot for genome size changes
python3 clockwise_plot.py size

# Comprehensive GC analysis for endosymbionts
python3 main_analysis.py endosymb_only gc_genome distance mean

# Summary across all groups
python3 summary_across_groups.py gc_genome mean

# Group-specific summary
python3 summary_data.py endosymb_only size all
```

## Dependencies

- matplotlib
- seaborn
- pandas
- numpy
- scipy (for correlation analyses)

## Output Formats

Plots are typically saved as PDF files for publication quality, with some options for PNG for web use.

## Integration with Analysis Scripts

These visualization scripts import analysis functions from the `analysis/` directory:
- `gcsize_dict.py` - Genome size and GC content data
- `delta_matrix.py` - Pairwise difference matrices
- `distance_matrix.py` - Phylogenetic distance matrices
- `matrix_correlation.py` - Mantel test correlations
- `utils.py` - Utility functions and path helpers