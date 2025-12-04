# Utility Scripts

This directory contains utility scripts for file management and data organization.

## Scripts

### merge_dirs.sh
Merge contents from multiple directories.
- **Input**: Multiple source directories
- **Output**: Merged directory structure
- **Use case**: Combining results from parallel analyses

### merge_relatives.sh
Merge endosymbiont data with free-living relative data.
- **Input**: Endosymbiont and relative genome directories
- **Output**: Combined dataset
- **Purpose**: Comparative analysis preparation

### move_relatives.sh
Move free-living relative genome files to appropriate locations.
- **Input**: Source directory with relative genomes
- **Output**: Organized directory structure
- **Purpose**: Data organization

### unzip_candidates.sh
Extract candidate genome archives.
- **Input**: Compressed genome files (.zip, .gz, .tar.gz)
- **Output**: Extracted genome files
- **Purpose**: Batch decompression

## Usage Examples

```bash
# Merge analysis results
./merge_dirs.sh

# Combine datasets
./merge_relatives.sh

# Organize genome files
./move_relatives.sh

# Extract archives
./unzip_candidates.sh
```

## Notes

These scripts help maintain organized data structures and prepare datasets for analysis workflows.
