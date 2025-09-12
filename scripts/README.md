# Scripts

This directory contains analysis pipeline scripts and utility programs for the T2T phasing project.

## Contents (20 files)

### Main Pipeline Scripts
- **phase_T2T_clean.sh** - Main production phasing pipeline for single chromosomes
- **create_and_assess_haplotype_panels.sh** - Comprehensive phasing and validation workflow
- **test_phasing.sh** - Parameterized testing framework for pipeline validation
- **parameters.sh** - Centralized parameter configuration

### Quality Assessment Scripts  
- **assess_imputation.sh** - Imputation accuracy evaluation across reference panels
- **calc_genomewide_imputation_statistics_full.sh** - Genome-wide imputation statistics compilation
- **genotyping_accuracy_with_happy.sh** - Genotype concordance measurement with hap.py

### Data Processing Scripts
- **create_summary_phasing_dataframes_polars_regional.py** - Regional analysis data processing (Polars-based)
- **generate_regions_for_rare_phasing.py** - Genomic region generation for rare variant chunks
- **find_trio_singletons.sh** - Identifies singleton variants in trio families
- **create_syn_nonsyn_bins.sh** - Creates synonymous/non-synonymous variant bins

### Specialized Analysis Scripts
- **Figure_6_script.R** - R script for generating publication Figure 6
- **liftover_panel.sh** - Coordinate liftover between T2T and GRCh38 reference systems
- **liftover_indels.py** - Python utility for indel coordinate conversion
- **get_discordant_multiallelic_sites.py** - Identifies problematic multiallelic variants

### Testing and Validation
- **all_testing_commands.sh** - Comprehensive test suite execution
- **test_testing.sh** - Basic pipeline validation

## Key Script Usage

### Basic Test Run
```bash  
# Test chr22 region with 12 threads
./scripts/test_phasing.sh chr22_test 12 test_run CHM13v2.0 true true 0.05 true false false 1
```

### Full Chromosome Phasing
```bash
# Phase chromosome 22 with production settings
./scripts/phase_T2T_clean.sh chr22 24 production 0.05
```

### Comprehensive Analysis Pipeline
```bash
# Complete phasing and assessment workflow
./scripts/create_and_assess_haplotype_panels.sh chr22 16 analysis_run CHM13v2.0
```

### Summary Statistics Generation
```bash
# Generate summary dataframes for all results
python3 scripts/create_summary_phasing_dataframes_polars_regional.py
```

## Script Dependencies

### System Requirements
- **Bash** - Shell scripts require bash 4.0+
- **Python 3** - Python scripts require Python 3.7+ with polars, pandas, numpy
- **R** - Figure generation requires R with ggplot2, data.table
- **bcftools** - VCF manipulation (v1.15+)
- **GNU parallel** - Job parallelization

### External Tools
- Docker (for containerized GATK4 and hap.py)
- Static binaries in `bin/` directory
- Reference data in `resources/` directory

## Pipeline Architecture

### Modular Design
Scripts are designed for:
- **Independent execution** - Each script can run standalone
- **Pipeline integration** - Scripts can be chained for complex workflows  
- **Parameterization** - Flexible parameter configuration via command line
- **Error handling** - Robust error detection and logging

### Configuration Management
- `parameters.sh` - Centralized configuration
- Command-line parameter passing
- Environment variable support
- Default parameter fallbacks

## Development and Testing

### Code Style
- Comprehensive error checking with `set -euo pipefail`
- Detailed logging and progress reporting
- Modular functions for reusability
- Clear parameter documentation

### Validation Framework
- Test regions for rapid validation
- Parameter sweep capabilities  
- Cross-reference validation support
- Performance benchmarking integration