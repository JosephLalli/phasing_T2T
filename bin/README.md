# Binaries

This directory contains static executable binaries required for the T2T phasing pipeline.

## Contents

### SHAPEIT5 Executables
- **SHAPEIT5_phase_common_static_v1.1.1** - Phases common variants in large chunks
- **SHAPEIT5_phase_rare_static_v1.1.1** - Phases rare variants in smaller genomic chunks  
- **SHAPEIT5_ligate_static_v1.1.1** - Ligates phased chunks into chromosome-wide haplotypes
- **SHAPEIT5_switch_static_v1.1.1** - Calculates switch errors against ground truth
- **SHAPEIT5_switch_static_JLL** - Custom switch error calculator with enhanced metrics

### Imputation and Genotype Concordance Tools
- **GLIMPSE2_concordance_static** - Measures imputation accuracy and concordance
- **impute5_1.1.5_static** - Alternative imputation method for validation

## Installation and Setup

All binaries are pre-compiled static executables and should work without additional dependencies.

### Make Executables
After cloning the repository, ensure binaries have execute permissions:
```bash
chmod +x bin/SHAPEIT5_*
chmod +x bin/GLIMPSE2_*  
chmod +x bin/impute5_*
```


## Binary Sources

### SHAPEIT5 v1.1.1
- Official release from SHAPEIT5 development team
- Source: https://github.com/odelaneau/shapeit5
- Used for statistical phasing with pedigree information

### Custom SHAPEIT5_switch_static_JLL  
- Modified version with enhanced switch error calculations
- Measures both flip and true switch rates
- [Source code modifications available in project repository, shapeit5 branch](https://github.com/JosephLalli/shapeit5)

### GLIMPSE2 and impute5
- GLIMPSE2: https://github.com/odelaneau/GLIMPSE  
- impute5: https://jmarchini.org/impute5/
- Used for imputation accuracy assessment and validation

## Platform Compatibility

Static binaries are compiled for:
- **Architecture**: x86_64 Linux
- **Dependencies**: None (statically linked)
- **Container support**: Compatible with Docker environment

For other platforms, recompilation from source may be required.