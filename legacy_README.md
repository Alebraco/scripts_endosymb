# Legacy Script Directories

These directories contain older versions of scripts organized differently. They are preserved for reference and backward compatibility but should not be used for new analyses.

## Directories

### legacy_GCSize_pipeline/
Scripts for GC content and genome size analysis pipeline.
- Many scripts have been reorganized into `/analysis/` and `/visualization/`
- Contains specialized GC content analysis utilities

### legacy_brccluster_scripts/
Collection of scripts from the BRC cluster analysis.
- Various analysis, annotation, and phylogenetics scripts
- Many have been reorganized into appropriate directories
- Some duplicates of scripts now in the main organized structure

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

## Why Keep Legacy Directories?

- **Backward compatibility** - Existing workflows may reference these paths
- **Complete archive** - Preserves all historical script versions
- **Reference** - May contain experimental or specialized versions

## Recommendation

For new analyses, use the reorganized structure in the main directories. The legacy directories will be deprecated in future versions once all workflows have been migrated.
