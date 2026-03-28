# Resources

This directory mixes two kinds of inputs:

1. Small, publication-tracked resources that ship with the repository.
2. Large external resources that are required at runtime but intentionally do
   not live in git.

## Tracked in git

- `GIABv3.6_bedfiles/`
- `pedigrees/`
- `recombination_maps/`
- `sample_subsets/`
- chain files, cytobands, BEDs, and other small annotations used by the paper

## Provided externally

These inputs are required for the full pipeline and for the Docker smoke test:

- CHM13 and GRCh38 reference FASTA files and indexes
- unphased 1KGP variant calls in T2T and GRCh38 coordinates
- phased GRCh38 reference panels
- SGDP truth data in T2T and GRCh38 coordinates
- HPRC and HGSVC3 pangenome VCFs and indexes

The publication branch no longer tracks machine-specific symlinks to those
files. Supply them in one of two ways:

1. Mount them into the container with `scripts/run_docker_smoke_test.sh`
   using a `docker.env` created from `docker.env.example`.
2. Download them manually to the runtime paths expected by the scripts.

## Canonical runtime paths

The shell pipeline expects these paths when running inside the repository or
inside the Docker container:

- `resources/chm13v2.0.fa.gz`
- `resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz`
- `resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz`
- `resources/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz`
- `resources/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz`
- `resources/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz`
- `resources/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz`
- `resources/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz`
- `resources/SGDP_variation/t2t`
- `resources/SGDP_variation/grch38`

## Source URLs

Representative download sources are below. Keep a frozen archive outside git
for the exact publication run if long-term bitwise replication matters.

### Reference genomes

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz
curl https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa | bgzip > resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
```

### HPRC pangenome VCFs

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
```

### HGSVC3 pangenome VCFs

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
```
