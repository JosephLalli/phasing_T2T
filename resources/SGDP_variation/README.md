# SGDP Variation Inputs

Ground-truth variant calls from the Simons Genome Diversity Project, used to
evaluate imputation accuracy. Data was published with the T2T HG002-Y paper
(Rhie et al., Nature 2023; https://www.nature.com/articles/s41586-023-06457-y).

## T2T-CHM13 coordinates

Publicly available. Run `scripts/utility/download_resources.sh` or download
manually:

```bash
for chr in {1..22} X; do
  wget -P resources/SGDP_variation/t2t/ \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
  wget -P resources/SGDP_variation/t2t/ \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/SGDP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi"
done
```

## GRCh38 coordinates

Download from Zenodo: [TBD]

Originally sourced from Terra/AnVIL:
https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY/data

Place (or symlink) per-chromosome VCFs at `resources/SGDP_variation/grch38/`.
Expected naming: `chr{N}.recalibrated.snp_indel.pass.vcf.gz` (with `.tbi` indexes).
