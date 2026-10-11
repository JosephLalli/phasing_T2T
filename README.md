[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17178669.svg)](https://doi.org/10.5281/zenodo.17178669)

**[Download phased haplotype panel at this location](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/)**

**[CHM13v2.0 recombination maps used by the pipeline](resources/recombination_maps/t2t_native_scaled_maps)**

The GRCh38 maps are in [resources/recombination_maps/grch38](resources/recombination_maps/grch38).


# Computationally phased T2T 1KGP panel

This repository contains CHM13v2-aligned 1000 Genomes Project (1KGP) variant data ([Rhie et al 2023](https://www.nature.com/articles/s41586-023-06457-y)) that have been computationally phased using SHAPEIT5 ([Hofmeister et al 2023](https://www.nature.com/articles/s41588-023-01415-w)). Phased panels (both unrelated 2504 member panels and full 3202 member panels) are available at the [T2T/HPRC aws bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/).

The previous code and analysis data are available at [Zenodo record 17178669](https://zenodo.org/record/17178669). Use the code and inputs described below to generate the current figures. Please cite our preprint, [Lalli et al. 2025](https://www.biorxiv.org/content/10.1101/2025.02.24.639687v1).

Read the whole-genome notebooks at [Phasing T2T notebooks](https://josephlalli.github.io/phasing_T2T/). Standalone HTML copies are included with the code archive and can be opened locally.

## Repository structure

```
├── bin/
│   └── [ GLIMPSE2 / impute5 / SHAPEIT5 static binaries ]
├── figures/ (Test-region figures)
│   ├── figure5/ (Figure 5 karyotype plots and assembled SVGs)
│   └── supplemental/
├── intermediate_data_whole_genome/ (Summary Parquet files)
├── imputation_statistics_whole_genome/ (Chromosome-level imputation results)
├── figures_whole_genome/ (Generated publication figures)
├── tables_whole_genome/ (Generated tables)
├── notebook_runs/ (Executed notebooks)
├── notebooks/
│   ├── calc_figures_for_paper.ipynb (Summary statistics)
│   ├── calc_per_variant_figures_for_paper.ipynb (Per-variant analysis)
│   └── notebooks_whole_genome/
│       ├── [ Summary and per-variant notebooks for whole-genome data ]
│       ├── make_plots.{py,ipynb} (Publication figures)
│       ├── utilities.py (Plotting and data helpers)
│       └── figure_composition.py (Combines figure panels)
├── resources/
│   ├── GIABv3.6_bedfiles/
│   │   └── [ GIAB v3.6 T2T/GRCh38 stratifications ]
│   ├── pedigrees/
│   │   └── 1kGP Pedigrees
│   ├── recombination_maps/
│   │   ├── grch38/ (HAPMAP GRCh38 recombination maps)
│   │   └── t2t_native_scaled_maps/ (T2T native recombination maps)
│   ├── sample_subsets/
│   │   └── [ Sample IDs, including pangenome relatives ]
│   ├── SGDP_variation/
│   │   └── [ Reference GRCh38/CHM13 Simons Genome Diversity Project variation ]
│   └── [ Additional reference files ]
├── scripts/
│   └── [Utility and analysis scripts]
├── tables/
│   └── [Test-region summary tables]
├── Dockerfile (Pipeline and figure dependencies)
└── README.md (This file)
```

## How to create this phased haplotype panel

### Clone and pull the image

Clone the repository and go into the phasing_T2T directory:

```bash
git clone https://github.com/JosephLalli/phasing_T2T
cd phasing_T2T
docker pull jlalli/phasing_t2t:v2.1
```

### Run the pipeline in Docker

Run on a server with Git and Docker.

Download the inputs with `scripts/utility/download_resources.sh`, then run the
pipeline from the repository root, providing the chromosome, number of
threads, output suffix, and reference genome. The input files are listed in
[resources/README.md](resources/README.md).

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 \
  bash scripts/utility/download_resources.sh
bash run_pipeline.sh chr22 12 run1 CHM13v2.0
bash run_pipeline.sh chr22 12 run1 GRCh38
```

For the chr15 and chr22 test regions, download the test inputs in Docker and
run with `chr15_test` or `chr22_test`:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 \
  bash scripts/utility/download_resources.sh --test
bash run_pipeline.sh chr15_test 12 test_run CHM13v2.0
bash run_pipeline.sh chr15_test 12 test_run GRCh38
bash run_pipeline.sh chr22_test 12 test_run CHM13v2.0
bash run_pipeline.sh chr22_test 12 test_run GRCh38
```

`run_pipeline.sh` runs `scripts/create_and_assess_haplotype_panels.sh` in Docker.
Put inputs at the paths listed in the resources README; results are saved in
the repository. The image and downloads require internet access.

### Publication figure code

Processing and phasing the full panel takes a lot of compute power. To
reproduce the figures without repeating those steps, we provide the
intermediate files on [Zenodo](https://zenodo.org/records/23289651): 12 summary
Parquet files and 176 chromosome-level files containing variant-level
imputation accuracy results.

Run these commands from an empty directory. They extract the code, download
the inputs into the repository's existing folders, and run the notebook in
Docker.

```bash
zenodo_files_url="https://zenodo.org/records/23289651/files"

curl -fL --retry 3 \
  "$zenodo_files_url/phasing_T2T_code.zip?download=1" \
  -o phasing_T2T_code.zip
unzip phasing_T2T_code.zip
cd phasing_T2T-NG-TR68126R

curl -fL --retry 3 \
  "$zenodo_files_url/NG-TR68126R_concordance_files.zip?download=1" \
  -o NG-TR68126R_concordance_files.zip
unzip NG-TR68126R_concordance_files.zip

mkdir -p intermediate_data_whole_genome
for filename in \
  MAF_performance_variants.parquet \
  ancestries.parquet \
  binned_maf_data.parquet \
  chroms.parquet \
  filter_summary_stats.parquet \
  methods.parquet \
  per_cytoband_variant_data.parquet \
  per_genome_cnv_regions.parquet \
  per_genome_cnv_regions_with_slop.parquet \
  per_sample_imputation_performance.parquet \
  per_variant_category_imputation_performance.parquet \
  samples.parquet
do
  curl -fL --retry 3 \
    "$zenodo_files_url/$filename?download=1" \
    -o "intermediate_data_whole_genome/$filename"
done

docker pull jlalli/phasing_t2t:v2.1
bash scripts/run_make_plots.sh
```

`run_make_plots.sh` checks the inputs against `SHA256SUMS.txt` and runs
`notebooks/notebooks_whole_genome/make_plots.ipynb` in Docker. The repository
is mounted at `/phasing_T2T_project`, which is the container's working
directory. Figures are saved in `figures_whole_genome/`, tables in
`tables_whole_genome/`, and the executed notebook in `notebook_runs/`, all
inside the repository on your disk. Internet access is needed to pull the
image and download inputs; the notebook runs without it.

The script uses 8 CPUs and 64 GiB of memory by default. Set `PHASING_T2T_THREADS`,
`PHASING_T2T_MEMORY`, or `PHASING_T2T_IMAGE` to change them.

Edit `notebooks/notebooks_whole_genome/make_plots.py`, then update the paired
`make_plots.ipynb` with Jupytext in the interactive Docker container below.

The Docker image includes Python 3.11 and the reported package versions from
`requirements.txt`, along with PyMuPDF, `rsvg-convert`, Arial, and the R packages
used for Figure 5.

`notebooks/notebooks_whole_genome/figure_composition.py` combines the figure
panels; these functions are also available through `utilities.py`. Generate
Figure 5 with `scripts/figure5/Figure_5_script.R` and combine its panels with
`scripts/figure5/stitch_svgs.py`. These scripts set the text sizes and overlap
between colorbar tiles.

The main and Extended Data figure filenames are:

| Figure | Output file under `figures_whole_genome/` |
| --- | --- |
| Figure 2 | `Figure 2.pdf` |
| Figure 3 | `Figure 3.pdf` |
| Figure 4 | `Figure 4.pdf` |
| Figure 5 | `figure5/figure5.svg` from the R script |
| Figure 6 | `Figure 6.pdf` |
| Extended Data Figure 1 | `Extended Data Figure 1.pdf` |
| Extended Data Figure 2 | `Extended Data Figure 2.pdf` |
| Extended Data Figure 3 | `Extended Data Figure 3.pdf` |
| Extended Data Figure 4 | `Extended Data Figure 4.pdf` |
| Extended Data Figure 5 | `Extended Data Figure 5.pdf` |
| Extended Data Figure 6 | `Extended Data Figure 6.pdf` |
| Extended Data Figure 7 | `Extended Data Figure 7.pdf` |
| Extended Data Figure 8 | `Extended Data Figure 8.pdf` |

Supplementary figure outputs use `supplemental/Supplementary Figure N.pdf`,
with the submission's numbering from 1 to 17. Some supplementary figures are
supplied separately and are not generated by this notebook.

Figure 1 is supplied separately. The notebook uses precomputed summary and
per-site inputs; it does not rerun variant processing or phasing. It has been
run with the reported Python package versions. The Figure 5 script has also
been run, but the original R and package versions were not reported.

Panel-wide headline SER is calculated in
`notebooks/notebooks_whole_genome/calc_figures_for_paper.ipynb`: counts are pooled
within each benchmark category, and the three category SERs are averaged using
panel-membership weights of 608 probands, 1,195 parents, and 1,399 non-trio
individuals.

### Run downloads and manual analysis in Docker

The commands from here through Figure 5 run inside one interactive container.
From the repository root, start it once:

```bash
docker run --rm -it --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 bash
```

Run the following commands at `/phasing_T2T_project` inside that container.

### Download inputs

To run the chr15 and chr22 test regions, first download the required input files:

```bash
bash scripts/utility/download_resources.sh --test
```

This downloads the test inputs, converts and indexes the FASTA files
(including `.gzi` files), and downloads the SGDP ground truth data.

Add `--parallel` to download the different datasets at the same time:

```bash
bash scripts/utility/download_resources.sh --test --parallel
```

To repeat the full analysis, download the larger files below.


### Primary 1kGP Variant Calls

**Unphased 1000 Genomes Project Variant Calls (T2T-CHM13v2.0)** (~30-130GB per chromosome):
```bash
# Download T2T-CHM13 variant calls for all chromosomes
for chr in {1..22} X; do
  wget -P unphased_variant_calls/t2t https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz
  wget -P unphased_variant_calls/t2t https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi
done
```

**Unphased 1000 Genomes Project Variant Calls (GRCh38)** (~30-130GB per chromosome):

These are only needed to compare how many variants are filtered by different
pre-phasing cutoffs.

```bash
# Download unphased, annotated GRCh38 1kGP variant calls.
for chr in {1..22} X; do
  wget -P unphased_variant_calls/grch38 https://42basepairs.com/download/s3/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz
  wget -P unphased_variant_calls/grch38 https://42basepairs.com/download/s3/1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz.tbi
done
```

**Phased GRCh38 1000 Genomes Variant Calls** (~30-130GB per chromosome):

```bash
# Download Byrska-Bishop (2022) phased GRCh38 1kGP variant calls.
for chr in {1..22} X; do
  wget -P phased_panels/grch38 https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  wget -P phased_panels/grch38 https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
done
```

### Reference Genomes

**T2T-CHM13 Reference Genome** (914MB):
```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz.fai
```

**GRCh38 Reference Genome** (783MB):
```bash
curl https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa | bgzip > resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
wget -P resources/ https://42basepairs.com/download/s3/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

### Reference Pangenome Variation:
```bash
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz.tbi

wget  -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
wget  -P resources/  https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz.tbi

wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi

wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi

wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz.tbi

wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz
wget -P resources/ https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_26_minigraph_cactus_hgsvc3/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz.tbi
```

You can also use `scripts/utility/download_resources.sh` to download the
reference genomes, pangenome VCFs, SGDP ground truth data, and FASTA indexes.
Run it from the repository root. Use `--test` for the chr15 and chr22 inputs,
and add `--parallel` to download the datasets at the same time.

### Run phasing scripts directly

To phase the chromosome 22 test region and compare the two panels, run:

```bash
# Phase chr22_test (chr22_test and chr15_test are recognized, along with chr1-22, X, PAR1, PAR2)
cd scripts
CHM13_suffix=test_phasing_CHM13
GRCh38_suffix=test_phasing_GRCh38
num_threads=12

./create_and_assess_haplotype_panels.sh chr22_test $num_threads $CHM13_suffix CHM13v2.0

#...after phasing the region, assessing the phasing quality, imputing SGDP data with the phased region, and assessing the quality of the imputation...#

# Evaluate the GRCh38 panel
./create_and_assess_haplotype_panels.sh chr22_test $num_threads $GRCh38_suffix GRCh38
```

The CHM13 panels are saved to
`phased_panels/phased_CHM13v2.0_panel_$CHM13_suffix`. The corresponding subsets
of the Byrska-Bishop GRCh38 panel are saved to
`phased_panels/phased_GRCh38_panel_$GRCh38_suffix`.

SHAPEIT5_switch results are saved to `SHAPEIT5_switch_output`.

## To replicate analysis in paper

See https://github.com/mccoy-lab/1kgp_chm13_maps for the recombination map code.
This repository contains the code for producing and evaluating the phased
T2T-CHM13 panel.

### Gather summary statistics

#### Download summary statistics

Use the concordance archive and 12 Parquet files described in
[Publication figure code](#publication-figure-code). The older Zenodo record
contains inputs for the previous version.

#### Calculate summary statistics

After running `create_and_assess_haplotype_panels.sh`, use the scripts below to
combine the SHAPEIT5_switch and GLIMPSE2_concordance results into dataframes.

These scripts can use over 350GB of RAM. Run on a server.

```bash
cd ..
# Build the syntenic/nonsyntenic label files used by the concordance aggregation.
./scripts/create_syn_nonsyn_bins.sh CHM13v2.0 $CHM13_suffix imputation_statistics/imputation_results_${CHM13_suffix} true
./scripts/create_syn_nonsyn_bins.sh GRCh38 $GRCh38_suffix imputation_statistics/imputation_results_${CHM13_suffix} true

# Collect imputation statistics in one place:
./scripts/calc_genomewide_imputation_statistics_full.sh $GRCh38_suffix $CHM13_suffix $num_threads true

# Collect all data into a series of summary parquet files
python3 ./scripts/analysis/create_summary_phasing_dataframes_polars_regional.py \
    --CHM13_run_suffix $CHM13_suffix \
    --GRCh38_run_suffix $GRCh38_suffix \
    --test
```

GLIMPSE2_concordance results for both references are saved to
`imputation_statistics/imputation_results_${CHM13_suffix}`.

Summary parquet files will be found in intermediate data

### Run data notebooks

The current publication notebooks are in `notebooks/notebooks_whole_genome/`:

- `calc_figures_for_paper.ipynb` calculates the values reported in the paper
  using the smaller summary Parquet files. Run on a laptop.
- `calc_per_variant_figures_for_paper.ipynb` calculates the values that require
  loading the full variant dataframes into memory. Run on a server.
- `make_plots.py` and its paired `.ipynb` generate the figures.

Run `make_plots.ipynb` in this Docker container.
The `test` profile writes to `figures/` and `tables/`; `whole_genome` writes
to `figures_whole_genome/` and `tables_whole_genome/`.

### Generate Figure 5

Run from the repository root inside the Docker container:

```bash
PHASING_T2T_RUN_PROFILE=whole_genome Rscript scripts/figure5/Figure_5_script.R
python3 scripts/figure5/stitch_svgs.py --batch figures_whole_genome/figure5
python3 scripts/figure5/stitch_svgs.py --grid figures_whole_genome/figure5 --pdf
```

Figure 5 is written to `figure5/figure5.svg` and `figure5/figure5.pdf`.
Set `PHASING_T2T_RUN_PROFILE=test` and use `figures/figure5` for test-region inputs.

Sync the paired publication notebook after editing it:

```bash
jupytext --sync notebooks/notebooks_whole_genome/make_plots.ipynb
```

When these manual commands finish, leave the container. Run the host
`run_pipeline.sh` or `scripts/run_make_plots.sh` wrappers only after exiting:

```bash
exit
```



## Switch Error Rates

Performance was measured by comparing the phased haplotypes of 39 samples shared between the Human Pangenome Reference Consortium (HPRC)'s draft human pangenome and the 1000 Genotypes dataset, or the 61 samples shared between HGSVC3 released assemblies and the 1000 Genomes Project.

|Panel|Source of ground truth assembly|Trio Membership|Switch Error Rate (%)|Flip Error Rate (%)|True Switch Error Rate (%)|Genotype Discordance Rate (%)|
|-|-|-|-|-|-|-|
|GRCh38|HPRC|Proband|0.408|0.186|0.029|1.114|
|GRCh38|HGSVC3|Proband|0.504|0.232|0.030|1.333|
|GRCh38|HGSVC3|Parent|0.705|0.303|0.080|1.457|
|GRCh38|HGSVC3|Non-trio|1.285|0.487|0.257|1.443|
|T2T-CHM13|HPRC|Proband|0.355|0.166|0.019|0.918|
|T2T-CHM13|HGSVC3|Proband|0.428|0.202|0.019|1.153|
|T2T-CHM13|HGSVC3|Parent|0.495|0.235|0.020|1.277|
|T2T-CHM13|HGSVC3|Non-Trio|1.110|0.474|0.130|1.253|

## Abbreviated Methods
### Variant filters used in this analysis (And string used to implement filter with bcftools)
<br>

- exclude FILTER (column in the VCF) = PASS:

      FILTER!="PASS"

- exclude variants with an alt allele of '*' after multiallelic splitting: 
        
      ALT=='*'
- exclude GT missingness rate < 5%

      F_MISSING>0.05
- exclude Hardy-Weinberg p-value < 1e−10 in any 1000G subpopulation (as calculated in the 2504 unrelated 1KGP samples):

      INFO/HWE_EUR<1e-10 && INFO/HWE_AFR<1e-10 && INFO/HWE_EAS<1e-10 && INFO/HWE_AMR<1e-10 && INFO/HWE_SAS<1e-10
- exclude sites where Mendelian Error Rate (Mendelian errors/num alleles) >= 0.05 (Note: 0.05*602 trios = 30 mendelian errors)
    
      INFO/MERR>=30
- exclude homoalellelic sites

      MAC==0
- exclude variants with a high chance of being errors as predicted by computational modeling*
      
      INFO/VQSLOD<0

Note: the SHAPEIT5 UK Biobank phasing paper excludes alternative alleles with AAscore < 0.5. This is a statistic produced by GraphTyper, which was not used to produce this dataset. The closest equivelent is the VQSLOD produced by Haplotype Caller. GraphTyper's AA score is simply the likelihood of an alternative allele truly being present in the dataset, so a cutoff of 0.5 is equivelent to 50% odds. Log odds of 1:1 is 0, so the VQSLOD log odds equivelent would be to exclude sites with a VQSLOD score of less than a cutoff of 0.

<br>

## Phasing
Phasing was performed with SHAPEIT v5.1.1, largely in accordance with the recommendations outlined in [SHAPEIT5's online tutorial](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs_200k/). For each chromosome, common variants were phased in one chunk. Rare variants were phased in chunks of 40 megabases.

Chromosome X PAR regions were phased separately from the rest of chromosome X. To phase the body of chromsome X, male samples were provided to SHAPEIT5 via the --halpoid option. The provided Ne value was reduced to 75% of the overall Ne, to account for the reduced population of haplotypes in this region of chrX. Phasing statistics were only calculated for female samples.

## Evaluation of phasing accuracy
We rely on two different methods of phasing quality evaluation:

- Using family data [as outlined in SHAPEIT5's documentation](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing).  

- Using the haplotype-phased samples present in the draft human pangenome as a ground truth.  

To perform these evaluations, we phase the data set four times:
- all 3202 samples, providing a pedigree
- all samples excluding those identified as parents in the 1KGP pedigree (2002 samples)
- all samples excluding those identified as parents of samples included in the draft pangenome (note: all shared 1KGP-pangenome samples are in the 1KGP dataset as children of trios).
- Only pangenome samples, using the first data set (with pangenome samples and their parents excluded) as a reference panel.

We then calculate different sets of performance statistics using SHAPEIT5's switch tool.  

Using family data as 'ground truth', we can:
- Evaluate switch error rate of pedigree-informed 3202 sample panel, using family data
    - Measures best-case performance of phasing when using pedigree data
- Evaluate switch error rate of 2002 sample no-parent phased panel using family data [as outlined in SHAPEIT5's documentation](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing). 
    - Most valid measure of phasing accuracy when not using pedigree data
    - Sets upper bound on error rate, as phasing is performed with ~66% of our dataset's haplotypes

Using the 39 1KGP samples present in the draft pangenome, we can use the HPRC's empirically phased haplotypes to:
- Evaluate switch error rate of the of pedigree-informed 3202 sample panel
    - Measures phasing accuracy when using a pedigree
- Evaluate switch error rate of a panel generated using samples excluding those identified as parents in the 1KGP pedigree (2002 samples)
    - Measures phasing accuracy without a pedigree
- Evaluate switch error rate when phasing a small number of samples (aka the 39 HPRC samples) using 2502 unrelated samples from the pedigree-informed 1KGP phased panel dataset as a reference (with the HPRC samples and parents excluded, of course)
    - Measures accuracy in most realistic use case
