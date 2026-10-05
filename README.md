[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19836492.svg)](https://doi.org/10.5281/zenodo.19836492)

**[Download phased haplotype panel at this location](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/)**

**[Download CHM13v2.0 recombination maps at this location](https://doi.org/10.5281/zenodo.19601957)**


# Computationally phased T2T 1KGP panel

This repository contains CHM13v2-aligned 1000 Genomes Project (1KGP) variant data ([Rhie et al 2023](https://www.nature.com/articles/s41586-023-06457-y)) that have been computationally phased using SHAPEIT5 ([Hofmeister et al 2023](https://www.nature.com/articles/s41588-023-01415-w)). Phased panels (both unrelated 2504 member panels and full 3202 member panels) are available at the [T2T/HPRC aws bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/).

The prior Zenodo archive at [record 19836492](https://zenodo.org/record/19836492) remains a reference release. Its analysis data are not the matched reviewer inputs for the current figures. Please cite our preprint describing this work, which is available at [Lalli et al. 2025](https://www.biorxiv.org/content/10.1101/2025.02.24.639687v1).

## Repository structure

```
├── bin/
│   └── [ GLIMPSE2 / impute5 / SHAPEIT5 static binaries ]
├── figures/ (Generated publication figures)
│   ├── figure6/ (Figure 6 karyotype plots and assembled SVGs)
│   └── supplemental/
├── notebooks/
│   ├── calc_figures_for_paper.ipynb (Summary statistics)
│   ├── calc_per_variant_figures_for_paper.ipynb (Per-variant analysis)
│   └── notebooks_whole_genome/
│       ├── [ Summary and per-variant notebooks for whole-genome data ]
│       ├── make_plots.{py,ipynb} (Publication figures)
│       ├── utilities.py (Plotting and data helpers)
│       └── figure_composition.py (Vector panel composition)
├── resources/
│   ├── GIABv3.6_bedfiles/
│   │   └── [ GIAB v3.6 T2T/GRCh38 stratifications ]
│   ├── pedigrees/
│   │   └── 1kGP Pedigrees
│   ├── recombination_maps/
│   │   ├── grch38/ (HAPMAP GRCh38 recombination maps)
│   │   └── t2t_native_scaled_maps/ (T2T native recombination maps)
│   ├── sample_subsets/
│   │   └── [ Sample IDs used to, for example, remove pangenome relatives ]
│   ├── SGDP_variation/
│   │   └── [ Reference GRCh38/CHM13 Simons Genome Diversity Project variation ]
│   └── [ Additional reference files ]
├── scripts/
│   └── [Utility and analysis scripts]
├── tables/
│   └── [Generated summary tables]
├── Dockerfile (Containerized environment)
└── README.md (This file)
```

### Key Directory Functions

- **bin/**: Static executables for SHAPEIT5 phasing and validation tools
- **resources/**: Reference data including genomes, genetic maps, and validation datasets
  - **pedigrees/**: Family relationship files for trio-informed phasing
  - **sample_subsets/**: Population and family groupings for analysis
  - **recombination_maps/**: Genetic distance maps for both T2T and GRCh38 references
- **scripts/**: Analysis pipeline and utility scripts for phasing and evaluation
- **notebooks/**: Jupyter notebooks for figure generation and statistical analysis
- **figures/** & **tables/**: Generated publication outputs


## How to create this phased haplotype panel

This section provides step-by-step instructions for how to recreate the phased CHM13v2.0 haplotype panel in your working environment.

### Prerequisites

Before downloading data files, ensure you have the following tools installed:

```bash
# Required system tools
sudo apt-get update && sudo apt-get install -y \
    wget curl bcftools tabix parallel docker.io

# Python packages (if running analysis scripts)
pip install polars pandas numpy jupyter
```

### Docker Environment (Recommended)

The published Docker image supports the pipeline test workflow. The separate reported-version figure environment below is used for the revised publication figures.

**Run in Docker container:**
```bash
# Pull the container
docker pull jlalli/phasing_t2t_dep_container:v2.0

# Start an interactive session with the project mounted
docker run -it -v /path/to/your/phasing_T2T_project:/workspace/phasing_T2T_project \
    -w /workspace/phasing_T2T_project \
    jlalli/phasing_t2t_dep_container:v2.0

# Once inside the container, you can run any pipeline commands
cd scripts
./create_and_assess_haplotype_panels.sh chr22_test 12 docker_test CHM13v2.0
```

Replace `/path/to/your/phasing_T2T_project` with the absolute path to your local copy of this repository.

For the pipeline smoke test, download the external test inputs and use the
repository's configured runner:

```bash
bash scripts/utility/download_resources.sh --test
docker pull jlalli/phasing_t2t_dep_container:v2.0
cp docker.env.example docker.env
# Review the input paths and output directory in docker.env.
bash scripts/run_docker_smoke_test.sh docker.env
```

The test uses the chr15 and chr22 regions for both references. It still needs
external biological inputs; [resources/README.md](resources/README.md) lists
their canonical paths. The published image and the revised figure image are
separate environments.

### Publication figure code

The matched reviewer release provides `JosephLalli__phasing_T2T__NG-TR68126R_20261004.zip`
and `NG-TR68126R_notebook_inputs.zip` together. Extract both archives as siblings;
the input archive contains the 188 manifest-listed inputs needed by the canonical
notebook. From the directory containing the downloaded files:

```bash
unzip JosephLalli__phasing_T2T__NG-TR68126R_20261004.zip
unzip NG-TR68126R_notebook_inputs.zip
cd phasing_T2T-NG-TR68126R
bash scripts/run_publication_notebook.sh \
  --inputs ../phasing_T2T-notebook-inputs \
  --output ../phasing_T2T-notebook-reproduction-output
```

The runner verifies the release manifest and canonical notebook hash before it
executes `notebooks/notebooks_whole_genome/make_plots.ipynb`, writes outputs to
the requested directory, and records the execution receipt in `scratch/`. Building
the Docker image downloads pinned packages; after it exists, notebook execution
runs with Docker networking disabled.

The canonical publication figure workflow is the Jupytext pair
`notebooks/notebooks_whole_genome/make_plots.py` and `make_plots.ipynb`; edit
the `py:percent` file and sync it before use. The unchanged test-region
`notebooks/make_plots.ipynb` remains a historical record; the prior
whole-genome version is available in Git history:

```bash
jupytext --sync notebooks/notebooks_whole_genome/make_plots.ipynb
```

`Dockerfile.reported-versions` provides Python 3.11 with the exact reported
numpy, polars, pandas, scipy, matplotlib, seaborn, and pyarrow versions. It also
includes PyMuPDF, `rsvg-convert`, and a verified Arial installation. Build it
and run the notebook with the repository and publication-scale inputs mounted
read-only, while mounting only the figure and table destinations writable:

```bash
# This image only needs its Dockerfile and exact figure lockfile as build inputs.
figure_build_dir=$(mktemp -d)
cp Dockerfile.reported-versions requirements.figures.lock.txt "$figure_build_dir/"
docker build -f "$figure_build_dir/Dockerfile.reported-versions" \
  -t phasing-t2t-figures:reported "$figure_build_dir"
rm -r "$figure_build_dir"

# Create the writable output directories before mounting them.
mkdir -p /path/to/scratch/figures /path/to/scratch/tables /path/to/scratch/notebooks

docker run --rm --network none --read-only --tmpfs /tmp:rw,size=8g \
  --user "$(id -u):$(id -g)" \
  -e MPLCONFIGDIR=/tmp/matplotlib \
  -e XDG_CACHE_HOME=/tmp/cache \
  -e JUPYTER_DATA_DIR=/tmp/jupyter-data \
  -e JUPYTER_CONFIG_DIR=/tmp/jupyter-config \
  -e JUPYTER_RUNTIME_DIR=/tmp/jupyter-runtime \
  -e IPYTHONDIR=/tmp/ipython \
  -e PHASING_T2T_RUN_PROFILE=whole_genome \
  -v "$PWD:/work/repo:ro" \
  -v /path/to/intermediate_data_whole_genome:/work/repo/intermediate_data_whole_genome:ro \
  -v /path/to/imputation_statistics_whole_genome:/work/repo/imputation_statistics_whole_genome:ro \
  -v /path/to/scratch/figures:/work/repo/figures_whole_genome \
  -v /path/to/scratch/tables:/work/repo/tables_whole_genome \
  -v /path/to/scratch/notebooks:/scratch \
  -w /work/repo phasing-t2t-figures:reported \
  jupyter nbconvert --to notebook --execute \
    notebooks/notebooks_whole_genome/make_plots.ipynb \
    --output /scratch/make_plots.executed.ipynb \
    --ExecutePreprocessor.timeout=-1
```

Vector figure composition is implemented in
`notebooks/notebooks_whole_genome/figure_composition.py`, with compatible
exports from `utilities.py`. Figure 5 is generated by
`scripts/figure6/Figure_6_script.R` and assembled by
`scripts/figure6/stitch_svgs.py`; its text sizes and colorbar tile overlap are
defined in those sources. The historical `figure6` names are retained.

The revision notebook retains historical output filenames. Use this map
for the current main and Extended Data figures:

| Current figure | Source output under `figures_whole_genome/` |
| --- | --- |
| Figure 2 | `Figure 2.pdf` |
| Figure 3 | `Figure 4 (a,b only) with graphic.pdf` |
| Figure 4 | `Figure 5.pdf` |
| Figure 5 | `figure6/figure6.svg` from the R workflow |
| Figure 6 | `Figure 7.pdf` |
| Extended Data 1 | `supplemental/Supplemental_13.pdf` |
| Extended Data 2 | `Extended Data - chromosome detail.pdf` |
| Extended Data 3 | `supplemental/Supplemental 17.pdf` |
| Extended Data 4 | `supplemental/Supplemental 18.pdf` |
| Extended Data 5 | `supplemental/Supplemental_19_segmental_duplication_genotype_error_with_ratio.pdf` |
| Extended Data 6 | `Extended Data - LiftoverIndel.pdf` |
| Extended Data 7 | `Extended Data - ancestry imputation.pdf` |
| Extended Data 8 | `supplemental/Supplemental 23.pdf` |

Figure 1 is supplied separately. The whole-genome revision notebook has
been executed with the reported Python package versions from existing
summary and per-site inputs. This check does not rerun the upstream variant
processing or phasing pipeline. The Figure 5 R workflow was checked with
the available R environment; original R/package versions were not reported.

Panel-wide headline SER is calculated in
`notebooks/notebooks_whole_genome/calc_figures_for_paper.ipynb`: counts are pooled
within each benchmark category, and the three category SERs are averaged using
panel-membership weights of 608 probands, 1,195 parents, and 1,399 non-trio
individuals.

### 1) Clone this repository

Begin by cloning this github repository to your local hard drive and going into the phasing_T2T directory:

```bash
# Create main project directory
git clone https://github.com/JosephLalli/phasing_T2T
cd phasing_T2T
```

### 2) Obtain primary data sources

If your goal is to run the Docker smoke test or the chr15/chr22 test regions,
do not skip this section. The repository includes the scripts and region
definitions for those test cases, but the biological inputs are still external.
The fastest supported setup path is:

```bash
bash scripts/utility/download_resources.sh --test
```

That helper downloads the canonical test inputs into the repo, converts and
indexes the FASTA files (including the required `.gzi` files), and populates
the SGDP truth data needed by the smoke test.

Add `--parallel` if you want the independent download workloads to run together
across sections:

```bash
bash scripts/utility/download_resources.sh --test --parallel
```

If you are interested in replicating our work completely, the necessary data that is too large to store on github can be obtained by following the instructions below.


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

(Note that these files are only necessary for comparing how many variants are filtered by different pre-phasing filter cutoffs.)

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

Alternatively, run `scripts/utility/download_resources.sh` from the repository root to fetch the canonical runtime inputs, including reference genomes, pangenome VCFs, SGDP truth data, and the required FASTA indexes. Use `--test` to fetch only the chr15/chr22 smoke-test inputs, and add `--parallel` to launch the independent download sections together.

#### Make binaries executable 

After downloading, ensure binaries are executable:
``` bash
# Check executable permissions
chmod +x SHAPEIT5_* GLIMPSE2_* impute5_*
chmod +x scripts/*.sh
```

### 3) Run phasing scripts

As an example, to compare the performance of phasing the test region of chromosome 22, run the following commands

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

The resulting phased haplotype panels will be found in phased_panels/phased_CHM13v2.0_panel_$CHM13_suffix folder, while the equivelent subsets of the Byrska-Bishop phased GRCh38 panel can be found in phased_panels/phased_GRCh38_panel_$GRCh38_suffix.

Raw output from SHAPEIT5_switch phasing accuracy assessments can be found in SHAPEIT5_switch_output.

## To replicate analysis in paper

Information on how to generate recombination maps can be found at https://github.com/mccoy-lab/1kgp_chm13_maps. This repository contains the portion of the analysis devoted to producing and evaluating a phased T2T-CHM13 reference panel.

### Gather summary statistics

#### Option A:  Download summary statistics 

 For the current reviewer workflow, use the matched `NG-TR68126R_notebook_inputs.zip`
 together with its code archive and `scripts/run_publication_notebook.sh`. The
 older Zenodo record is a prior reference release and is not the exact input set
 for the current figure notebook.

#### Option B: Generate summary statistics from output of create_and_assess_haplotype_panels.sh

If you are interesting in replicating our work from scratch, please continue from step 3 above. These scripts gather the results of SHAPEIT5_switch and GLIMPSE2_concordance into experiment-wide dataframes.

Please note: These scripts can easily take up over 350GB of ram as written.

```bash
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

GLIMPSE2_concordance output (Both GRCh38 and CHM13v2.0) will be found in imputation_statistics/imputation_results_${CHM13_suffix}

Summary parquet files will be found in intermediate data

### Run data notebooks

The current publication notebooks are in `notebooks/notebooks_whole_genome/`:

- `calc_figures_for_paper.ipynb` calculates manuscript summary values.
- `calc_per_variant_figures_for_paper.ipynb` performs the heavier per-variant analyses.
- `make_plots.py` and its paired `.ipynb` generate the publication figures.

Use the reported-version figure workflow above to execute `make_plots.ipynb`.
The unchanged test-region `notebooks/make_plots.ipynb` remains a historical
record; the prior whole-genome version is available in Git history.
The `test` profile writes to `figures/` and `tables/`; `whole_genome` writes
to `figures_whole_genome/` and `tables_whole_genome/`.

### Run the Figure 5 R workflow

Run from the repository root in an environment with the required R packages:

```bash
PHASING_T2T_RUN_PROFILE=whole_genome Rscript scripts/figure6/Figure_6_script.R
python3 scripts/figure6/stitch_svgs.py --batch figures_whole_genome/figure6
python3 scripts/figure6/stitch_svgs.py --grid figures_whole_genome/figure6 --pdf
```

The historical `figure6.svg` and `figure6.pdf` filenames correspond to
current Figure 5. Set `PHASING_T2T_RUN_PROFILE=test` and use `figures/figure6`
for test-region inputs.



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
