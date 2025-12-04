# Legacy Script Directories

These directories contain older versions of scripts organized differently. They are preserved for reference and backward compatibility but should not be used for new analyses.

## Directories

### legacy_brccluster_scripts/
Collection of scripts from the BRC cluster analysis.
- Various analysis, annotation, and phylogenetics scripts
- Many duplicates of scripts now in the main organized structure
- Preserved for reference only

## Migration

If you were using scripts from these directories:

1. **Check the main directories first** - Most functionality has been reorganized into:
   - `/annotation/` - for genome annotation
   - `/alignment/` - for sequence alignment
   - `/analysis/` - for GC content, IGS, and other analyses
   - `/phylogenetics/` - for tree building and clustering
   - `/visualization/` - for plotting

2. **Update your paths** - Scripts have moved but functionality remains the same

3. **Check README files** - Each directory has documentation

## GC Content Pipeline

The GCSize_pipeline scripts have been fully integrated into the main repository structure:
- **Analysis scripts** (computation, matrices, utilities) → `/analysis/gc_*.py`
- **Visualization scripts** (plotting) → `/visualization/gc_*.py`

All scripts have been renamed with a `gc_` prefix for clarity and imports have been updated.

## Why Keep Legacy Directories?

- **Backward compatibility** - Existing workflows may reference these paths
- **Complete archive** - Preserves all historical script versions
- **Reference** - May contain experimental or specialized versions

## Recommendation

For new analyses, use the reorganized structure in the main directories. The legacy directories will be deprecated in future versions once all workflows have been migrated.

