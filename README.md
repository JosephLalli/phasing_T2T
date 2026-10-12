[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17178669.svg)](https://doi.org/10.5281/zenodo.17178669)

**[Download phased haplotype panel at this location](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/)**

**[CHM13v2.0 recombination maps used by the pipeline](resources/recombination_maps/t2t_native_scaled_maps)**

[Classic GRCh38 recombination maps used by the pipeline](resources/recombination_maps/grch38)


# Computationally phased T2T 1KGP panel

This repository contains CHM13v2-aligned 1000 Genomes Project (1KGP) variant data ([Rhie et al 2023](https://www.nature.com/articles/s41586-023-06457-y)) that have been computationally phased using SHAPEIT5 ([Hofmeister et al 2023](https://www.nature.com/articles/s41588-023-01415-w)). Phased panels (both unrelated 2504 member panels and full 3202 member panels) are available at the [T2T/HPRC aws bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/).

The previous code and analysis data are available at [Zenodo record 17178669](https://zenodo.org/record/17178669). Use the code and inputs described below to generate the current figures. Please cite our preprint, [Lalli et al. 2025](https://www.biorxiv.org/content/10.1101/2025.02.24.639687v1).

Code used to generate T2T recombination maps can be found at [1kgp_chm13_maps](https://github.com/mccoy-lab/1kgp_chm13_maps).

Read the whole-genome notebooks at [Phasing T2T notebooks](https://josephlalli.github.io/phasing_T2T/). Standalone HTML copies are included with the code archive and can be opened locally.

## Repository structure

```
├── bin/
│   └── [ GLIMPSE2 / impute5 / SHAPEIT5 static binaries ]
├── figures/ (Test-region figures)
│   ├── figure5/ (Figure 5 karyotype plots and assembled SVGs)
│   └── supplemental/
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

### Clone this repository and pull the repo image

Ensure you have Git and Docker installed on your machine.

Clone the repository and go into the phasing_T2T directory:

```bash
git clone https://github.com/JosephLalli/phasing_T2T
cd phasing_T2T
docker pull jlalli/phasing_t2t:v2.1
```

### Phase 1KGP variation and calculate summary statistics

Download the inputs listed in [resources/README.md](resources/README.md). We have created a script to do this [scripts/utility/download_resources.sh](scripts/utility/download_resources.sh). 
From the repository root, download all input files that are too large to host on github using Docker:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 \
  bash scripts/utility/download_resources.sh
```

To only download test datasets, run:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/phasing_T2T_project" \
  --workdir /phasing_T2T_project jlalli/phasing_t2t:v2.1 \
  bash scripts/utility/download_resources.sh --test
```

Then, to run the pipeline from the repository root, provide the chromosome, number of
threads, output suffix, and reference genome (in that order):

```bash
bash run_pipeline.sh chr22 12 run1 CHM13v2.0
bash run_pipeline.sh chr22 12 run1 GRCh38
```

For the chr15 and chr22 test regions, run with `chr15_test` or `chr22_test`:

```bash
bash run_pipeline.sh chr15_test 12 test_run CHM13v2.0
bash run_pipeline.sh chr15_test 12 test_run GRCh38
bash run_pipeline.sh chr22_test 12 test_run CHM13v2.0
bash run_pipeline.sh chr22_test 12 test_run GRCh38
```

`run_pipeline.sh` is a driver script that coordinates running `scripts/create_and_assess_haplotype_panels.sh` in Docker.
Put inputs at the paths listed in the resources README; results are saved in
the repository. The image and downloads require internet access.

### Produce figures from summary statistics

Processing and phasing the full panel takes a lot of compute power. To reproduce the figures without repeating those steps, we provide the intermediate files on [Zenodo](https://zenodo.org/records/23289651): 12 summary Parquet files and 176 chromosome-level files containing variant-level imputation accuracy results.

Assuming you have already cloned the Github repository, run the following commands from the repository root to download the intermediate data inputs into their appropriate folders:

```bash
zenodo_files_url="https://zenodo.org/records/23289651/files"

curl -fL "$zenodo_files_url/NG-TR68126R_concordance_files.zip?download=1" \
  -o NG-TR68126R_concordance_files.zip
unzip NG-TR68126R_concordance_files.zip

mkdir -p "intermediate_data_whole_genome"

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
  curl -fL "$zenodo_files_url/$filename?download=1" \
    -o "intermediate_data_whole_genome/$filename"
done
```

Then, you can run the driver script that runs figure generation scripts using a Docker container with all requirements installed.

```bash
bash scripts/run_make_plots.sh
```

The runner checks the input files against their supplied checksums, then runs `notebooks/notebooks_whole_genome/make_plots.ipynb` in a Docker container. Figures, tables, and the executed notebook are saved in `figures_whole_genome/`, `tables_whole_genome/`, and `notebook_runs/` inside the repository on your disk.

Figure 1 is generated as part of the [1kgp_chm13_maps](https://github.com/mccoy-lab/1kgp_chm13_maps) repository. Figure 5 is generated separately by running `scripts/figure5/Figure_5_script.R` followed by `scripts/figure5/stitch_svgs.py`.

The script uses 8 CPUs and 64 GiB of memory by default. Set `PHASING_T2T_THREADS`,
`PHASING_T2T_MEMORY`, or `PHASING_T2T_IMAGE` to change them.

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

## Inline statistics

We also report many individual statistics in our paper. These are calculated by 
`notebooks/notebooks_whole_genome/calc_figures_for_paper.ipynb`. 

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
