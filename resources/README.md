# Resources

This directory holds the reference files used by the analysis. Small annotations are included with the repository; download the larger inputs separately.

## Included files

- `GIABv3.6_bedfiles/`
- `pedigrees/`
- `recombination_maps/`
- `sample_subsets/`
- chain files, cytobands, BEDs, and other small annotations used by the paper

## Download inputs

Use the Docker download command in the [main README](../README.md#phase-1kgp-variation-and-calculate-summary-statistics). It runs `scripts/utility/download_resources.sh` inside the project image and downloads the reference genomes, unphased 1KGP calls, GRCh38 comparison panels, SGDP truth data, and HPRC/HGSVC3 pangenome VCFs; it also prepares and indexes the FASTAs. Existing files are skipped.

For the chr15 and chr22 test regions, use the `--test` Docker command shown there. It extracts test-region unphased calls and SGDP truth data, but downloads whole-chromosome GRCh38 comparison panels because liftover needs variants outside the test regions.

## Required file paths

The pipeline expects these paths in the cloned repository:

- `resources/chm13v2.0.fa.gz`
- `resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz`
- `resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz` (+`.tbi`)
- `resources/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz` (+`.tbi`)
- `resources/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz` (+`.tbi`)
- `resources/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz` (+`.tbi`)
- `resources/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz` (+`.tbi`)
- `resources/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz` (+`.tbi`)
- `resources/SGDP_variation/t2t/`
- `resources/SGDP_variation/grch38/`

All VCF files require their corresponding `.tbi` tabix index at the same path.

## Sources and manual downloads

Use the supported Docker download command above for a complete setup. The links below are provided for readers who need an individual input or wish to inspect its source.

To run the manual download commands, start a container from the repository root:

```bash
docker run --rm -it --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 bash
```

Run the commands below inside that container. The downloaded files are saved in the repository on your disk.

### Reference genomes

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz
curl https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa | bgzip > resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
```

### HPRC pangenome VCFs

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi
```

### HGSVC3 pangenome VCFs

```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi
```

### SGDP ground truth

T2T-CHM13 VCFs are available from the [T2T/HPRC bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/).

The GRCh38 archives are available on Zenodo: [test-region data](https://zenodo.org/records/19371182/files/SGDP_GRCh38_test_regions.tar.gz?download=1) or [whole-genome data](https://zenodo.org/records/19371182/files/SGDP_GRCh38_all_chromosomes.tar.gz?download=1). The automated downloader selects the appropriate archive and extracts it into `resources/SGDP_variation/grch38/`.

The GRCh38 data were originally sourced from [AnVIL T2T CHRY](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_T2T_CHRY/data). See [SGDP_variation/README.md](SGDP_variation/README.md) for expected file names.
