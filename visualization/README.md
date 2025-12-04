# Visualization Scripts

This directory contains scripts for generating plots and visualizations.

## Scripts

### variation_plot.py
Generate plots for genetic variation analysis.
- **Input**: Variation data (CSV or similar format)
- **Output**: Visualization plots (PDF/PNG)
- **Dependencies**: matplotlib, seaborn
- **Use case**: Visualizing genetic diversity and variation patterns

## Additional Visualization

For more visualization scripts, see:
- **analysis/igs_plots.py** - Intergenic spacer visualizations
- **GCSize_pipeline/** - GC content and genome size plots (in legacy directory)

## Usage Example

```bash
python3 variation_plot.py
```

## Dependencies

- matplotlib
- seaborn
- pandas (for data handling)

## Output Formats

Plots are typically saved as PDF files for publication quality, with options for PNG for web use.
