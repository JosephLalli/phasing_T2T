# Unphased Variant Calls

This directory is where to place unphased 1000 Genomes Project variant calls mapped to T2T-CHM13 and GRCh38.

## Directories

- `t2t/`: calls aligned to the T2T-CHM13v2.0 reference
- `grch38/`: calls aligned to the GRCh38 reference

## Download calls

From the repository root on your host machine, use the command in the [main README](../README.md#phase-1kgp-variation-and-calculate-summary-statistics) to download and prepare the inputs. It downloads the calls and their tabix indexes into the directories above. For the chr15 and chr22 test regions, use the `--test` command shown there; it uses `bcftools` in the image to obtain only the needed regional calls.

## Sources and manual downloads

The T2T-CHM13 calls are from the [T2T/HPRC bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/). The GRCh38 calls are from the [1000 Genomes Project archive](https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/).

Full-chromosome files are about 30--130 GB each. Individual files can be downloaded directly on your host machine; keep the `.vcf.gz` file and its matching `.tbi` index in the appropriate directory. Use the automated command above to prepare a complete set.
