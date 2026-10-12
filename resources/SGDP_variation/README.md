# SGDP Variation Inputs

Ground-truth variant calls from the Simons Genome Diversity Project, used to
evaluate imputation accuracy. Data was published with the T2T HG002-Y paper
(Rhie et al., Nature 2023; https://www.nature.com/articles/s41586-023-06457-y).

From the repository root on your host machine, use the commands in the
[main README](../../README.md#phase-1kgp-variation-and-calculate-summary-statistics).
The downloader obtains both T2T-CHM13 and GRCh38 inputs and prepares their indexes with the tools in the Docker image.
Use `--test` for the chr15 and chr22 test regions.

## T2T-CHM13 coordinates

The full chromosome VCFs and indexes are publicly available. To download them
manually instead, run this on your host machine from the repository root:

```bash
for chr in {1..22} X; do
  wget -P resources/SGDP_variation/t2t/ \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
  wget -P resources/SGDP_variation/t2t/ \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi"
done
```

## GRCh38 coordinates

The automated download uses the full `SGDP_GRCh38_all_chromosomes.tar.gz` or
test-region `SGDP_GRCh38_test_regions.tar.gz` archive from
[Zenodo](https://zenodo.org/records/19371182), then stores the VCFs here and
creates missing `.tbi` indexes.

Originally sourced from Terra/AnVIL:
https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY/data

Expected location: `resources/SGDP_variation/grch38/`
Expected naming: `chr{N}.recalibrated.snp_indel.pass.vcf.gz` (with `.tbi` indexes).
