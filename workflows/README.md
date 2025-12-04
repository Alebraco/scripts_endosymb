# Workflow Scripts

This directory contains high-level workflow scripts that orchestrate multiple analysis steps.

## Scripts

### endosymb_core.sh
Main workflow script for core endosymbiont genome analysis.
- **Input**: Genome datasets
- **Output**: Complete analysis results
- **Process**: 
  1. Genome quality control
  2. Core gene identification
  3. Alignment and concatenation
  4. Phylogenetic analysis
  5. Comparative statistics

## Usage

```bash
./endosymb_core.sh
```

## Notes

Workflow scripts typically call multiple scripts from other directories (annotation, alignment, analysis, phylogenetics) in a coordinated sequence. They may require adjustment of paths and parameters based on your specific dataset and computing environment.

## Creating Custom Workflows

To create your own workflow:
1. Start with this template structure
2. Source/call scripts from relevant directories
3. Ensure proper input/output file paths
4. Add error checking and logging
5. Document required inputs and expected outputs

## Related Workflows

See the `alignment/run_*.sh` and `analysis/run_*.sh` scripts for specialized sub-workflows.
