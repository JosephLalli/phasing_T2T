[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7612953.svg)](https://doi.org/10.5281/zenodo.7612953)


# Computationally phased T2T 1KGP panel

# [Download panel at this location](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/)

This repository contains CHM13v2-aligned 1000 Genomes Project (1KGP) variant data ([Rhie et al 2022](https://www.biorxiv.org/content/10.1101/2022.12.01.518724v1)) that have been computationally phased using SHAPEIT5 ([Hofmeister et al 2022](https://www.biorxiv.org/content/10.1101/2022.10.19.512867v2)). Phased panels (both unrelated 2504 member panels and full 3202 member panels) are available at the [T2T/HPRC aws bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/).


A tar file containing the initial version of our phased T2T panel, along with another files containing phasing statistics has been uploaded to zenodo at [this link](https://zenodo.org/record/7612953#.Y-8VD3bMJPY). Please cite the zenodo repository if you find this panel to be useful in your work. A preprint is in the works, and once it is up we will have a proper citation for you.

## Repository structure
- resources: bed files and vcf files used to create and benchmark this panel.
  - pedigrees: 1000 Genomes Project pedigree files (and subsets of that pedigree)
  - sample_subsets: text files containing different groups of samples used for annotating or producing subsets of vcf files.
  - recombination maps
    - t2t_lifted_from_hapmap: HapMap recombination maps lifted from GRCh38 recombination map.
    - t2t_native_scaled_maps: T2T recombination maps generated from pyrho analysis of 1KGP T2T variant panel by @andrew-bortvin.
- scripts: various scripts used to evaluate performance and gather GLIMPSE2_concordance or SHAPEIT5_concordance output into single dataframes.

## Reference data sources
- 1KGP variants called against T2Tv2.0: [T2T_1KGP_calls](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202/)
- T2T reference fasta: [chm13v2.0.fa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/GCA_009914755.4/chm13v2.0.fa.gz)
- HPRC Pangenome sample variation in T2T coordinates: [hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz)
- HPRC Pangenome sample variation in GRCh38 coordinates: [hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz)
- SHAPEIT5 5.1.1 release files: [SHAPEIT5 v5.1.1 release](https://github.com/odelaneau/shapeit5/releases/tag/v5.1.1)


## Switch Error Rates

Performance was measured by comparing the phased haplotypes of 39 samples shared between the Human Pangenome Reference Consortium (HPRC)'s draft human pangenome and the 1000 Genotypes dataset.

|Method of phasing variants|Switch Error Rate|+/- 95th c.i.|Ground Truth Source|
|-|-|-|-|
|All 3202 samples, including trios, with pedigree|0.07% |0.003%|Trio concordance|
|2002 singletons or trio children, no pedigree|1.31%|0.03%|Trio concordance|
|44 HPRC samples, using the 2504 unrelated samples panel (with HPRC samples and parents removed) as a reference panel |2.29%|0.10%|HPRC concordance|



## Methods
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
      
      INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0

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
  
