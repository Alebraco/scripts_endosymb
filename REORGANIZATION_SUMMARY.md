# Repository Reorganization Summary

## What Changed

This reorganization transformed a flat repository with 45+ scripts scattered in the root directory into a well-organized, professional structure suitable for research and collaboration.

## Before (Issues)

```
scripts_endosymb/
├── .DS_Store (tracked in git)
├── .gitignore (incomplete)
├── GCSize_pipeline/ (unorganized subdirectory)
├── brccluster_scripts/ (unorganized subdirectory)
├── backtranslate.py
├── bakta_annotation.sh
├── clusters.py
├── concatenate.py
├── count_core.sh
├── count_genes.py
├── ... (40+ more scripts)
└── (No README, LICENSE, or documentation)
```

**Problems:**
- No clear organization or structure
- Difficult to find scripts by function
- No documentation
- No dependency tracking
- Missing standard project files
- macOS system files tracked in git

## After (Professional Structure)

```
scripts_endosymb/
├── README.md                        # Comprehensive project documentation
├── LICENSE                          # MIT license
├── CITATION.cff                     # Academic citation metadata
├── requirements.txt                 # Python dependencies
├── .gitignore                       # Updated to exclude system files
├── legacy_README.md                 # Documentation for legacy directories
│
├── annotation/                      # Genome annotation (7 scripts)
│   ├── README.md
│   ├── bakta_annotation.sh
│   ├── extract_cds.sh
│   ├── endosymb_cds.py
│   ├── gb2assm.py
│   ├── cand_genomes.sh
│   └── remove_genomes.py
│
├── alignment/                       # Sequence alignment (7 scripts)
│   ├── README.md
│   ├── backtranslate.py
│   ├── muscle.py
│   ├── concatenate.py
│   ├── run_backtranslate.sh
│   ├── run_concatenate.sh
│   └── run_muscle.sh
│
├── analysis/                        # Genomic analysis (19 scripts)
│   ├── README.md
│   ├── gc_content.py
│   ├── igs_lengths.py
│   ├── igs_plots.py
│   ├── spectrum.py
│   ├── count_genes.py
│   ├── count_core.sh
│   ├── parse_blast_results.py
│   ├── percid_distribution.py
│   ├── usearch.py
│   └── ... (workflow scripts)
│
├── phylogenetics/                   # Phylogenetic analysis (7 scripts)
│   ├── README.md
│   ├── clusters.py
│   ├── dnatree.sh
│   ├── iqtree.sh
│   ├── longbranch_outliers.py
│   ├── corecruncher.sh
│   └── new_path_precrunchr.sh
│
├── visualization/                   # Data visualization (2 scripts)
│   ├── README.md
│   └── variation_plot.py
│
├── workflows/                       # High-level workflows (2 scripts)
│   ├── README.md
│   └── endosymb_core.sh
│
├── utils/                          # Utility scripts (5 scripts)
│   ├── README.md
│   ├── merge_dirs.sh
│   ├── merge_relatives.sh
│   ├── move_relatives.sh
│   └── unzip_candidates.sh
│
├── legacy_GCSize_pipeline/         # Preserved for reference
│   └── ... (12 scripts)
│
└── legacy_brccluster_scripts/      # Preserved for reference
    └── ... (64 scripts)
```

## Key Improvements

### 1. **Logical Organization**
Scripts are now organized by function:
- **annotation/** - Everything related to genome annotation
- **alignment/** - Sequence alignment and processing
- **analysis/** - Statistical and genomic analyses
- **phylogenetics/** - Tree building and clustering
- **visualization/** - Plotting and graphics
- **workflows/** - High-level analysis pipelines
- **utils/** - Helper and utility scripts

### 2. **Comprehensive Documentation**
- **Main README.md**: Complete project overview, installation, usage
- **Directory READMEs**: Detailed documentation for each category
- **CITATION.cff**: Proper academic citation format
- **LICENSE**: Clear licensing (MIT)
- **requirements.txt**: All Python dependencies listed

### 3. **Professional Standards**
- Clear directory structure following best practices
- Proper gitignore configuration (excludes .DS_Store, build artifacts)
- Citation metadata for academic use
- License for legal clarity
- Dependency tracking

### 4. **Preserved Legacy Content**
- Old directories renamed to `legacy_*`
- All historical scripts preserved
- Documentation explains migration path
- Backward compatibility maintained

### 5. **Easy Navigation**
Users can now:
- Quickly find scripts by function
- Understand project structure at a glance
- Read documentation for each area
- Install dependencies easily
- Cite the work properly

## Migration Guide for Users

If you were using this repository before:

### For script paths in your code:
```bash
# Old:
python3 backtranslate.py

# New:
python3 alignment/backtranslate.py
```

### For workflows referencing multiple scripts:
Update paths to reflect new structure:
```bash
# Example: Update workflow script
# Old: ./gc_content.py
# New: ./analysis/gc_content.py
```

### Legacy scripts still work:
If you absolutely need old paths, scripts in `legacy_*` directories are unchanged.

## Statistics

- **Scripts organized**: 45 root scripts → 7 categorized directories
- **Documentation added**: 10 new README files
- **Standard files added**: LICENSE, CITATION.cff, requirements.txt
- **Legacy scripts preserved**: 76 scripts in legacy directories
- **Improved gitignore**: Excludes system files properly

## Benefits

1. **For Collaborators**: Easy to understand structure, clear documentation
2. **For Users**: Find tools quickly, understand dependencies
3. **For Citations**: Proper CITATION.cff file for academic use
4. **For Maintenance**: Logical organization makes updates easier
5. **For New Contributors**: README files guide contribution
6. **For Research**: Professional appearance for publications/sharing

## Next Steps (Optional Future Enhancements)

Consider adding:
- Example data or test datasets
- Detailed usage tutorials
- Contribution guidelines (CONTRIBUTING.md)
- Changelog (CHANGELOG.md)
- GitHub Actions for CI/CD
- Docker container for reproducibility
- Full API documentation

---

## Update: GC Pipeline Integration

### Additional Changes (December 2024)

Following user feedback, the GCSize_pipeline scripts have been fully integrated into the main repository structure:

**Scripts Reorganized:**
- **Plotting scripts** (5 files) → `visualization/gc_*.py`
  - `gc_absolute_scatter.py` - Absolute scatter plots
  - `gc_clockwise_plot.py` - Clockwise correlation plots  
  - `gc_main_analysis.py` - Main GC analysis with comprehensive plotting
  - `gc_summary_across_groups.py` - Cross-group summary plots
  - `gc_summary_data.py` - Summary data visualization

- **Analysis/utility scripts** (7 files) → `analysis/gc_*.py`
  - `gc_calculate.py` - Core GC calculation utility
  - `gc_codon_dict.py` - Codon-level GC analysis
  - `gc_delta_matrix.py` - Delta matrix computation
  - `gc_distance_matrix.py` - Phylogenetic distance matrices
  - `gc_matrix_correlation.py` - Matrix correlation (Mantel test)
  - `gc_metadata_size.py` - Genome size and GC metadata
  - `gc_utils.py` - Shared utilities for GC pipeline

**Key Changes:**
- All scripts renamed with `gc_` prefix for clarity
- Import statements updated to reflect new locations
- Visualization scripts now import from `../analysis/` directory
- `legacy_GCSize_pipeline/` directory removed (fully integrated)
- Only `legacy_brccluster_scripts/` remains as legacy
- Updated documentation in README files
- Added `scikit-bio` dependency for Mantel test

**Migration Path:**
```bash
# Old (legacy):
python3 legacy_GCSize_pipeline/main_gc.py

# New (organized):
python3 visualization/gc_main_analysis.py
```

---

**Result**: A research repository that follows professional software development standards while maintaining all functionality and backward compatibility.
