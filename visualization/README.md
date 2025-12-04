# Visualization Scripts

This directory contains scripts for generating plots and visualizations.

## Scripts

### variation_plot.py
Generate plots for genetic variation analysis.
- **Input**: Variation data (CSV or similar format)
- **Output**: Visualization plots (PDF/PNG)
- **Dependencies**: matplotlib, seaborn
- **Use case**: Visualizing genetic diversity and variation patterns

### igs_plots.py
Generate visualization plots for intergenic spacer (IGS) analysis.
- **Input**: IGS data CSV files
- **Output**: Box plots and scatter plots (PDF format)
- **Dependencies**: matplotlib, seaborn, pandas

## GC Content and Genome Size Visualization

### gc_absolute_scatter.py
Generate absolute scatter plots for GC content and genome size.
- **Usage**: `python3 gc_absolute_scatter.py <type> <mean/all>`
- **Input**: GC content and genome size data
- **Output**: Scatter plots showing relationships

### gc_clockwise_plot.py
Generate clockwise plots for GC and genome size analysis.
- **Input**: Delta matrices and phylogenetic distances
- **Output**: Clockwise correlation plots

### gc_main_analysis.py
Main GC content analysis with comprehensive plotting.
- **Usage**: `python3 gc_main_analysis.py <group> <matrix_type> <dist_type> <data_type>`
- **Input**: GC data, distance matrices, genome metadata
- **Output**: Correlation plots and statistics

### gc_summary_across_groups.py
Generate summary plots across different groups (endosymbionts vs relatives).
- **Usage**: `python3 gc_summary_across_groups.py <matrix_type> <all/mean>`
- **Input**: GC data across groups
- **Output**: Comparative summary plots

### gc_summary_data.py
Generate summary data plots for GC content analysis.
- **Usage**: `python3 gc_summary_data.py <group> <matrix_type> <sp/all>`
- **Input**: Group-specific GC data
- **Output**: Summary visualization plots

## Usage Examples

```bash
# Genetic variation plot
python3 variation_plot.py

# IGS analysis plots
python3 igs_plots.py

# GC content scatter plots
python3 gc_absolute_scatter.py gc_content mean

# Comprehensive GC analysis
python3 gc_main_analysis.py endosymb_only gc_content patristic mean
```

## Dependencies

- matplotlib
- seaborn
- pandas
- numpy
- scipy (for some GC scripts)
- BioPython (for phylogenetic analyses)

## Output Formats

Plots are typically saved as PDF files for publication quality, with options for PNG for web use.

