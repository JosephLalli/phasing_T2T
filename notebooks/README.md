# Notebooks

This directory contains Jupyter notebooks for data analysis, figure generation, and statistical computation.

## Directory Structure
- Main directory - Core analysis notebooks (88 files total)
- `notebooks_whole_genome/` - The same notebooks run on genome-wide results

## Core Analysis Notebooks

### Primary Publication Notebooks
- **calc_figures_for_paper.ipynb** - Calculates all values reported in manuscript text
  - Summary statistics across populations and chromosomes
  - Switch error rate calculations  
  - Imputation accuracy measurements
  - Designed to run on laptop with summary data
  
- **calc_per_variant_figures_for_paper.ipynb** - Per-variant detailed analysis
  - Requires loading full variant dataframes into memory
  - Variant-level switch error analysis
  - Regional and population-specific variant patterns
  - Requires high-memory server environment

- **make_plots.ipynb** - Primary figure generation notebook
  - Creates all main, extended, and supplementary figures
  - Outputs figures to `figures/` directory
  - Generates tables for `tables/` directory
  - Note: Does not generate Figure 6 (handled by R script)

## Data Dependencies

### Input Data Sources
- `intermediate_data/` - Test-region summary analysis files
- `imputation_statistics/` - Test-region imputation results
- `intermediate_data_whole_genome/` - Complete genome-wide analysis results
- `imputation_statistics_whole_genome/` - Complete genome-wide imputation results
- `SHAPEIT5_switch_output/` - Switch error calculations

### Output Locations
- `figures/` and `tables/` - Test run figures and tables
- `figures_whole_genome/` and `tables_whole_genome/` - Publication-scale figures and tables
- Inline notebook outputs - Interactive analysis results
