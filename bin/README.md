# Binaries

This directory contains the static binaries used by the T2T phasing pipeline.

- `SHAPEIT5_phase_common_static_v1.1.1` phases common variants in large chunks.
- `SHAPEIT5_phase_rare_static_v1.1.1` phases rare variants in smaller chunks.
- `SHAPEIT5_ligate_static_v1.1.1` ligates phased chunks into chromosome-wide haplotypes.
- `SHAPEIT5_switch_static_v1.1.1` calculates switch errors against ground truth.
- `SHAPEIT5_switch_static_JLL` is the custom switch error calculator.

`GLIMPSE2_concordance` measures imputation accuracy and concordance. `impute5_v1.2.0_static` is used for imputation validation.

The binaries are included and executable in the Docker image. They are compiled
for x86_64 Linux and are statically linked.

SHAPEIT5 v1.1.1 is from [its source repository](https://github.com/odelaneau/shapeit5). `SHAPEIT5_switch_static_JLL` measures flip and true switch rates; its [source changes are in the project shapeit5 branch](https://github.com/JosephLalli/shapeit5). GLIMPSE2 is from [its source repository](https://github.com/odelaneau/GLIMPSE) and impute5 is from [its source page](https://jmarchini.org/impute5/).
