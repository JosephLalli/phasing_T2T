# Unphased Variant Calls

This directory contains unphased variant call files from the 1000 Genomes Project mapped to both T2T-CHM13 and GRCh38 reference genomes.

## Directory Structure
- `t2t/` - Variant calls aligned to T2T-CHM13v2.0 reference
- `grch38/` - Variant calls aligned to GRCh38 reference

## Required Downloads

To obtain the raw vcf files used as inputs to this pipeline, please download the following files.

### T2T-CHM13v2.0 Variant Calls (~30-130GB per chromosome)
```bash
# Download T2T-CHM13 variant calls for all chromosomes
for chr in {1..22} X; do
  wget -P unphased_variant_calls/t2t https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz
  wget -P unphased_variant_calls/t2t https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi
done
```

### GRCh38 Variant Calls (~30-130GB per chromosome)
```bash
# Download unphased, annotated GRCh38 1kGP variant calls
for chr in {1..22} X; do
  wget -P unphased_variant_calls/grch38 https://42basepairs.com/download/s3/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz
  wget -P unphased_variant_calls/grch38 https://42basepairs.com/download/s3/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz.tbi
done
```