#!/usr/bin/env bash
set -eo pipefail

# set up script log
## one-liner to get script location; credit to https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
basedir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
logfile=$basedir/phasing_T2T_$1_$3.log
exec 19>$logfile
export BASH_XTRACEFD=19
set -x # writes commands to logfile

chrom=$1
num_threads=$2
suffix=$3

missing_filter_cutoff=$4

## imputation options
# limit_to_syntenic_regions=$4 #'false'
limit_to_snps=$5 #'false'
filter_VQSLOD=$6
impute_5_hmm=$7
impute_only=$8

hmm_ne=$9
chrom_map=${10}


if [[ -z $chrom_map ]]
then
    chrom_map=$basedir/t2t_maps_no_filtered_regions/scaled_fixedamount/${chrom}_noMask.scaled.gmap.gz
fi

if [[ -z $hmm_ne ]]
then
    hmm_ne=135000
fi

if [[ $impute_5_hmm == 'true' ]]
then
    impute5_hmm_option="--ne $hmm_ne"
    echo $impute5_hmm_option
fi


if [[ suffix != '' ]]
then
    suffix="_$suffix"
fi


set -u

final_panel_dir=$basedir/phased_T2T_panel${suffix}
stats_dir=$basedir/phasing_stats${suffix}_custom_switch
imputation_results_dir=$basedir/imputation_results$suffix
initial_vcf_calls_folder=/mnt/ssd/lalli/phasing_T2T

ref_fasta=$basedir/chm13v2.0.fa.gz
pangenome_vcf=$basedir/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
HGSVC_vcf=$basedir/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
# HGSVC_vcf=$basedir/MC_hgsvc3-hprc_chm13_filtered_decomposed.vcf.gz

population_ids=$basedir/sample_subsets/unrelated_superpopulations.csv
chrom_chunking_coords=$basedir/regions2.txt

syntenic_site_location="$basedir/chm13v2-syntenic_to_hg38.bed"

## Identify if we are dealing with a PAR region
# GRCh38: "X:10001-2781479" "X:2781480-155701382" "X:155701383-156030895"
# T2T: "X:0-2394410" "X:2394411-153925834" "X:153925835-154259566"

PAR1_start='0'
PAR1_end='2394410'
PAR2_start='153925834'
PAR2_end='154259566'


if [[ $chrom == "chrPAR1" || $chrom == "PAR1" ]] 
then
    chrom='PAR1'
    region="chrX:$(($PAR1_start+1))-$PAR1_end"
    grch38_region="chrX:10001-2781479"
    whole_chrom=$region
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    end_chrom=$PAR1_end

elif [[ $chrom == "chrPAR2" || $chrom == "PAR2" ]]
then
    chrom='PAR2'
    region="chrX:$(($PAR2_start+1))-$PAR2_end"
    grch38_region="chrX:155701384-156030895"
    whole_chrom=$region
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    end_chrom=$PAR2_end
elif [[ $chrom == 'chrdebug' ]]
then
    chrom='debug'
    region="chr20:10000000-12000000"
    whole_chrom=$region
    grch38_region=$region
    chrom_map=$basedir/t2t_maps_no_filtered_regions/scaled_fixedamount/chr20_noMask.scaled.gmap.gz
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chr20.recalibrated.snp_indel.pass.vcf.gz
    end_chrom=12000000

elif [[ $chrom == 'chrX' ]]
then
    region="chrX:$(($PAR1_end+1))-$PAR2_start"
    whole_chrom=$region
    grch38_region="chrX:2781480-155701383"
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    end_chrom=$PAR2_start
    use_beagle='true'
    hmm_ne=$(awk "BEGIN { printf \"%.0f\", ($hmm_ne*.75) }")
else
    region=$chrom
    end_chrom=$(grep "^$chrom	" ${ref_fasta}.fai | cut -f 2)
    whole_chrom=$chrom:1-$end_chrom
    grch38_region=$chrom
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.$chrom.recalibrated.snp_indel.pass.vcf.gz
fi

chrom_working_dir=$basedir/${chrom}_working${suffix}

mkdir -p $chrom_working_dir
mkdir -p $final_panel_dir
mkdir -p $stats_dir


# move command logfile to newly created working directory
logfile=$chrom_working_dir/phasing_T2T_$1_$3.log
cp $basedir/phasing_T2T_$1_$3.log $logfile
exec 19>$logfile
export BASH_XTRACEFD=19
set -x # writes commands to logfile
rm $basedir/phasing_T2T_$1_$3.log


# Assign chrom-specific regions filename
chrom_regions=$chrom_working_dir/${chrom}_regions.txt


if [[ $chrom == *PAR* ]]
then
    echo $region > $chrom_regions
elif [[ $chrom == 'debug' ]]
then
    echo $region > $chrom_regions
else
    grep $chrom: $chrom_chunking_coords > $chrom_regions
fi

if [[ $chrom == 'chrX' ]]
then
    # Make list of males for X chrom phasing
    male_sample_list=$basedir/sample_subsets/males.txt
    cat $basedir/sample_subsets/superpopulations.samples.txt | grep " 1 " | cut -f 2 -d " " > $male_sample_list
    haploid_arg="--haploids $male_sample_list"
    impute5_haploid_arg="--haploid"
else
    haploid_arg=""
    impute5_haploid_arg=""
fi

# Define names of chrom specific vcfs that we will be phasing
vcf_to_phase=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.biallelic.bcf
fully_annotated_input_variants=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.bcf

chr_specific_reference_pangenome_variation_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.biallelic.bcf
chr_specific_reference_HGSVC_variation_biallelic=$chrom_working_dir/${chrom}_reference_HGSVC.biallelic.bcf
chr_specific_reference_pangenome_variation_trimmed_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.filtered_variants.biallelic.bcf

# Define the names of input variant files that are sample subsets
vcf_to_phase_no_parents=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unphased.noparents.bcf
vcf_to_phase_pangenome_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.pangenome.biallelic.bcf

# Define phased result file names
common_variants_phased_ped=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.common.bcf
rare_variants_phased_ped=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.rare.bcf
rare_variants_phased_ped_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.biallelic.rare.bcf
phased_panel_no_pangenome_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.no_pangenome.biallelic.bcf

common_variants_phased_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.common.pangenome_samples.biallelic.bcf
rare_variants_phased_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.biallelic.rare.pangenome_samples.bcf

vcf_phased_no_parents_common_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.biallelic.bcf
vcf_phased_no_parents_rare_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.biallelic.bcf

# output files
phased_panel_vcf_3202=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.3202.vcf.gz
phased_panel_vcf_2504=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.2504.vcf.gz
phased_panel_vcf_3202_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.vcf.gz
phased_panel_vcf_2504_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.vcf.gz

fully_annotated_input_variant_report=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv
phased_panel_vcf_2504_biallelic_variant_report=$stats_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.2504.stats.tsv

# pedigrees
pedigree=$basedir/pedigrees/1kgp.ped
bcftools_formatted_pedigree=$basedir/pedigrees/trios_only.ped

# lists of different categories of samples
no_parents=$basedir/sample_subsets/not_parents.txt
unrelated_samples=$basedir/sample_subsets/unrelated_samples.txt
pangenome_samples=$basedir/sample_subsets/pangenome_samples.txt
pangenome_and_parents=$basedir/sample_subsets/pangenome_samples_and_parents.txt
SGDP_in_1KGP=$basedir/sample_subsets/SGDP_samples_in_1KGP_numbered_format.txt

# Imputation and imputation metric gathering

GRCh38_fasta=/dev/shm/Homo_sapiens_assembly38.fasta
grch38_syntenic_site_location="$basedir/hg38.GCA_009914755.4.synNet.summary.bed.gz"

lifted_panel_folder=$basedir/liftover/lifted_panels
GRCh38_lifted_panel=$lifted_panel_folder/$(echo $(basename ${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf) | sed 's/CHM13v2.0/GRCh38.lifted_from_CHM13v2.0/' | sed 's/native_maps.//' )
GRCh38_native_panel=$basedir/phased_GRCh38_panel/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
T2T_lifted_panel=$lifted_panel_folder/$(echo $(basename $GRCh38_native_panel) | sed 's/\.GRCh38\./.CHM13v2.0.lifted_from_GRCh38./')
T2T_native_panel=${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf

SGDP_ground_truth_dir_T2T=/mnt/ssd/lalli/nf_stage/genome_refs/T2T-CHM13_v2_ncbi110/SGDP
SGDP_ground_truth_T2T=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf
SGDP_ground_truth_dir_GRCh38=$basedir/GRCh38_SGDP_full
SGDP_ground_truth_GRCh38=$chrom_working_dir/SGDP.GRCh38.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf

pangenome_ground_truth_GRCh38=$basedir/GRCh38_pangenome_variation/${chrom}_reference_pangenome.filtered_variants.biallelic.bcf
pangenome_ground_truth_T2T=$chr_specific_reference_pangenome_variation_trimmed_biallelic

T2T_lifted_panel_no_pangenome=$lifted_panel_folder/1KGP.CHM13v2.0.lifted_from_GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf
GRCh38_lifted_panel_no_pangenome=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf
T2T_native_panel_no_pangenome=$phased_panel_no_pangenome_biallelic
GRCh38_native_panel_no_pangenome=$basedir/phased_GRCh38_panel/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf


#### making symlinks to skip re-phasing. Do not use for production.

prior_working_dir=$basedir/${chrom}_working_T2T_scaled_newimpute_071524
# ln -sf $prior_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.common.{vcf.gz,bcf} $chrom_working_dir/ #$common_variants_phased_ped
# ln -sf $prior_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.native_maps.common.{vcf.gz.tbi,bcf.csi} $chrom_working_dir/ #$common_variants_phased_ped
# ln -sf $prior_working_dir/*[0-9].rare.bcf* $chrom_working_dir/
# ln -sf $prior_working_dir/*tmp*.rare.bcf* $chrom_working_dir/
# # rm $chrom_working_dir/*_tmp_pangenome.rare.bcf*
# ln -sf $prior_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.biallelic.bcf* $chrom_working_dir/ #$vcf_phased_no_parents_common_biallelic
# ln -fs $prior_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.biallelic.bcf* $chrom_working_dir/ #$vcf_phased_no_parents_rare_biallelic


#settings
default_mcmc_iteration_scheme='5b,1p,1b,1p,1b,1p,5m'
shapeit4_suggested_unlimited_resources_mcmc_iteration_scheme='10b,1p,1b,1p,1b,1p,1b,1p,10m'
mcmc_iteration_scheme=$shapeit4_suggested_unlimited_resources_mcmc_iteration_scheme  #$default_mcmc_iteration_scheme
common_pbwt_depth=8  # rare default: 2; common default: 4; high-accuracy: 8
common_pbwt_mac=5    # rare default: 2; common default: 5; shapeit4_default=2
rare_pbwt_depth=2
rare_pbwt_mac=2
pbwt_mdr=0.1         # rare default: 0.1; common default: 0.1; shapeit4_default=0.05
pbwt_modulo=0.1      # rare default: 0.1; common default: 0.1; shapeit4_sequencing_default=0.0005
window=5             # rare default: 4;   common default: 4;   1kgp paper using shapeit4: 5
rare_variant_threshold=0.001    # default: 0.001

graph_reference_missing_cutoff=0.2 #$missing_filter_cutoff #0.2

## VARIANT FILTERING, QC VARIANT SUBSETTING/FORMATTING

# Create ground truth variation vcf from pangenome
## Note: Throughout, I am comparing this vcf of genome assemblies to variant calls.
## One difference between assemblies and calls is that there is no such thing as a missing variant in an assembly.
## What is represented as a missing variant (overlapping indel) is consistently represented as a reference allele in our variant calls
## So we will convert pangenome missing variants into pangenome reference alleles.

if [ ! -s $chr_specific_reference_pangenome_variation_biallelic.csi ]
then
    echo "making reference pangenome variation"
    if [[ $chrom == 'chrX' ]] || [[ $chrom == 'PAR2' ]] || [[ $chrom == 'PAR1' ]]; then    # make missing/haploid into haploid
        bcftools view --threads 8 -r $region -s ^GRCh38 --force-samples -Ou $pangenome_vcf \
        | bcftools norm --threads 2 -Ou --atomize --atom-overlaps \. -m +snps - \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|\.,qqq,g' | sed 's,\.|,0|,g' | sed 's,|\.,|0,g' | sed 's,qqq,\.|\.,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools +setGT -Ou - -- -t a -n p \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_pangenome_variation_biallelic \
        && bcftools index $chr_specific_reference_pangenome_variation_biallelic &

    else
        bcftools view --threads 8 -r $region -s ^GRCh38  --force-samples -Ou $pangenome_vcf \
        | bcftools norm --threads 2 -Ou --atomize --atom-overlaps \. -m +snps - \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|,0|,g' | sed 's,|\.,|0,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools view --threads 1 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_pangenome_variation_biallelic \
        && bcftools index $chr_specific_reference_pangenome_variation_biallelic &
    fi
fi

if [ ! -s $chr_specific_reference_HGSVC_variation_biallelic.csi ]
then
    echo "making reference pangenome variation"
    if [[ $chrom == 'chrX' ]] || [[ $chrom == 'PAR2' ]] || [[ $chrom == 'PAR1' ]]; then    # make missing/haploid into haploid
        bcftools view --threads 8 -r $region -s ^GRCh38 --force-samples -Ou $HGSVC_vcf \
        | bcftools norm --threads 2 -Ou --atomize --atom-overlaps \. -m +snps - \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|\.,qqq,g' | sed 's,\.|,0|,g' | sed 's,|\.,|0,g' | sed 's,qqq,\.|\.,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools +setGT -Ou - -- -t a -n p \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_HGSVC_variation_biallelic \
        && bcftools index $chr_specific_reference_HGSVC_variation_biallelic &

    else
        bcftools view --threads 8 -r $region -s ^GRCh38  --force-samples -Ou $HGSVC_vcf \
        | bcftools norm --threads 2 -Ou --atomize --atom-overlaps \. -m +snps - \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|,0|,g' | sed 's,|\.,|0,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools view --threads 1 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_HGSVC_variation_biallelic \
        && bcftools index $chr_specific_reference_HGSVC_variation_biallelic &
    fi
fi

wait

# Split multiallelic sites, filter sites using criteria described above,
# convert data to bcf format, and index.
echo "creating unphased, annotated, filtered variant panel"
if [ ! -s $fully_annotated_input_variants.csi ]
then
    echo "fully_annotated_input_variants"
    bcftools norm --threads 8 -Ou -r $region -f $ref_fasta -m -any $input_vcf \
    | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +mendelian2 -Ou - --ped $bcftools_formatted_pedigree -m a -m d \
    | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,F_MISSING,HWE,MAC:1=MAC \
    | bcftools +fill-tags --threads 8 -Oz - -- -S $population_ids -t HWE \
    | tee ${fully_annotated_input_variants%%.bcf}.vcf.gz \
    | bcftools view --threads 8 -Ob - \
        > $fully_annotated_input_variants \
    && bcftools index --threads 8 -f -t ${fully_annotated_input_variants%%.bcf}.vcf.gz \
    && bcftools index --threads 8 -f $fully_annotated_input_variants
fi

if [ ! -s $vcf_to_phase.csi ]
then
    echo "vcf_to_phase"

    if [[ $chrom != 'chrX' ]]
    then
        bcftools view \
                -e "(TYPE!='snp' && (ABS(ILEN) >= 50)) || ALT=='*' || INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || INFO/MAC==0 || ( INFO/HWE_EUR<1e-10 && INFO/HWE_AFR<1e-10 && INFO/HWE_EAS<1e-10 && INFO/HWE_AMR<1e-10 && INFO/HWE_SAS<1e-10 ) || FILTER!='PASS'" \
                --threads 8 -Ou $fully_annotated_input_variants \
        | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS - > $vcf_to_phase \
        && bcftools index --threads 8 -f $vcf_to_phase
    else
        bcftools view \
                -e "(TYPE!='snp' && (ABS(ILEN) >= 50)) || ALT=='*' || INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || INFO/MAC==0 || ( INFO/HWE_EUR<1e-10 && INFO/HWE_AFR<1e-10 && INFO/HWE_EAS<1e-10 && INFO/HWE_AMR<1e-10 && INFO/HWE_SAS<1e-10 ) || FILTER!='PASS'" \
                --threads 8 -Ou $fully_annotated_input_variants \
        | bcftools +fixploidy --threads 2 -Ou - -- -f 2 \
        | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS - > $vcf_to_phase \
        && bcftools index --threads 8 -f $vcf_to_phase
    fi
fi


# While the phasing is running, create a tsv table of the to-be-phased variants and quality metrics for downstream QC
if [[ ! -s $fully_annotated_input_variant_report ]]
then
    echo "making variant report"
    if [ ! -s ${fully_annotated_input_variants%%.bcf}.vcf.gz.tbi ]
    then
        echo "fully_annotated_input_variants"
        bcftools norm --threads 8 -Ou -r $region -f $ref_fasta -m -any $input_vcf \
        | bcftools +mendelian2 -Ou - --ped $bcftools_formatted_pedigree -m a -m d \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                            -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                            -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,F_MISSING,HWE,MAC:1=MAC \
        | bcftools +fill-tags --threads 8 -Oz - -- -S $population_ids -t HWE \
        | tee ${fully_annotated_input_variants%%.bcf}.vcf.gz \
        | bcftools view --threads 8 -Ob - \
         > $fully_annotated_input_variants \
        && bcftools index --threads 8 -f -t ${fully_annotated_input_variants%%.bcf}.vcf.gz \
        && bcftools index --threads 8 -f $fully_annotated_input_variants
    fi

    docker run --user $UID -v $chrom_working_dir:$chrom_working_dir quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0 \
        gatk VariantsToTable \
            -V ${fully_annotated_input_variants%%.bcf}.vcf.gz \
            -F ID \
            -F CHROM -F POS -F ALT -F QUAL \
            -F NEGATIVE_TRAIN_SITE \
            -F VQSLOD \
            -F MERR \
            -F HWE -F HWE_EUR -F HWE_AFR -F HWE_EAS -F HWE_AMR -F HWE_SAS \
            -F FILTER \
            -F culprit \
            -F InbreedingCoeff \
            -F ExcessHet \
            -F AN -F AC -F MAF -F F_MISSING \
            -F SYNTENIC \
            --show-filtered \
            --split-multi-allelic \
            -O $fully_annotated_input_variant_report 2> /dev/null &
fi


### PHASING
## Phase full 3202 sample panel with pedigree
# (max accuracy, requires trios and therefore accuracy measurements not generalizable)
# Phase common variants
if [ ! -s $common_variants_phased_ped.csi ]
then
    ./SHAPEIT5_phase_common_static_v1.1.1 \
        --input $vcf_to_phase \
        --map $chrom_map \
        --output $common_variants_phased_ped \
        --thread $num_threads \
        --log $chrom_working_dir/$chrom.common.log \
        --filter-maf $rare_variant_threshold \
        --mcmc-iterations $mcmc_iteration_scheme \
        --pbwt-modulo $pbwt_modulo \
        --pbwt-depth $common_pbwt_depth \
        --pbwt-mac $common_pbwt_mac \
        --pbwt-mdr $pbwt_mdr \
        --region $region \
        --pbwt-window $window --hmm-window $window --hmm-ne $hmm_ne \
        $haploid_arg \
        --pedigree $pedigree && \
    bcftools index --threads 8 -f $common_variants_phased_ped
fi

### Add rare variants to previously generated common variant scaffold
i=1
for chrom_region in $(cat $chrom_regions)
do
    if [ ! -s $chrom_working_dir/$i.rare.bcf.csi ]
    then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase \
            --map $chrom_map \
            --scaffold $common_variants_phased_ped \
            --thread $num_threads \
            --pedigree $pedigree \
            --log $chrom_working_dir/$chrom.$i.rare.log \
            --pbwt-modulo $pbwt_modulo \
            --input-region $chrom_region \
            --scaffold-region $chrom_region \
            --effective-size $hmm_ne \
            $haploid_arg \
            --output $chrom_working_dir/$i.rare.bcf && \
        bcftools index -f --threads 8 $chrom_working_dir/$i.rare.bcf
    fi
    let i++
done

### Rare variant phasing had to be done in chunks due to high memory requirements.
### Stitch chunks back together, merge multiallelic sites, and index the resulting bcf.
if [ ! -s $rare_variants_phased_ped_biallelic.csi ]; then
    echo "rare_variants_phased_ped_biallelic"
    bcftools concat --threads 8 -Ou -l $chrom_working_dir/[^_].*rare.bcf \
    | bcftools norm --threads 8 -Ou -m -any --fasta $ref_fasta - \
    | bcftools view -c 1:minor --threads 8 -Ou -V other - \
    | bcftools annotate -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -Ou -x INFO/MAC --threads 2 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
    > $rare_variants_phased_ped_biallelic \
    && bcftools index --threads 8 -f $rare_variants_phased_ped_biallelic
fi

echo "making phased panel subsets"
if [ ! -s $rare_variants_phased_ped.csi ]; then
    bcftools norm --threads 8 -Ou --fasta $ref_fasta -m +any $rare_variants_phased_ped_biallelic \
    | bcftools annotate -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -Ou -x INFO/MAC --threads 2 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags -Ou --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools view --threads 8 -Ou -V other - \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_ped && \
    bcftools index --threads 8 -f $rare_variants_phased_ped
fi


if [ ! -s $phased_panel_vcf_3202.tbi ]; then
    ### Extract 2504 unrelated samples from results. This will be the official phased panel.
    echo "phased_panel_vcf_3202"
    bcftools annotate --threads 8 -Oz -x ^INFO/MAF,^INFO/MAC,^INFO/AC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT $rare_variants_phased_ped > $phased_panel_vcf_3202 \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_3202 &
fi

if [ ! -s ${phased_panel_vcf_3202_biallelic%%.vcf.gz}.bcf.csi ]; then
    echo "phased_panel_vcf_3202_biallelic"
    bcftools annotate --threads 8 -Oz -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT $rare_variants_phased_ped_biallelic > $phased_panel_vcf_3202_biallelic \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_3202_biallelic \
    && bcftools view  --threads 8 -Ob $phased_panel_vcf_3202_biallelic > ${phased_panel_vcf_3202_biallelic%%.vcf.gz}.bcf \
    && bcftools index --threads 8 ${phased_panel_vcf_3202_biallelic%%.vcf.gz}.bcf &
fi

wait

if [ ! -s ${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf.csi ]; then
    ### Note: inputation evaluation requires MAF is in the INFO field, and it's not that much of a bother/size increase
    echo "phased_panel_vcf_2504_biallelic"
    bcftools view -Ou --threads 8 -r $region -S $unrelated_samples $phased_panel_vcf_3202_biallelic \
    | bcftools view -Ou --threads 2 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ou --threads 4 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools annotate -Oz --threads 4 -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT - > $phased_panel_vcf_2504_biallelic \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_2504_biallelic \
    && bcftools view  --threads 8 -Ob $phased_panel_vcf_2504_biallelic > ${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf \
    && bcftools index --threads 8 ${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf &
fi

# while the above is going on, proceed with phasing the no-parents panel

if [ ! -s $vcf_to_phase_no_parents.csi ]; then
    ## Repeat, but with no trio parents (per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing)
    ### Remove parents from unphased vcf file
    echo "vcf_to_phase_no_parents"
    bcftools view -Ou --threads 8 -S $no_parents $vcf_to_phase \
    | bcftools view -Ob -c 1:minor --threads 8 - > $vcf_to_phase_no_parents \
    && bcftools index --threads 8 -f $vcf_to_phase_no_parents #biallelic
fi


### Perform same phasing procedure on children-and-singleton-only vcf
### Phase common variants
if [ ! -s $vcf_phased_no_parents_common_biallelic.csi ]; then
    ./SHAPEIT5_phase_common_static_v1.1.1 \
        --input $vcf_to_phase_no_parents \
        --map $chrom_map \
        --output $chrom_working_dir/noparents.bcf \
        --thread $num_threads \
        --log $chrom_working_dir/$chrom.common.noparents.log \
        --filter-maf $rare_variant_threshold \
        --mcmc-iterations $mcmc_iteration_scheme \
        --pbwt-modulo $pbwt_modulo \
        --pbwt-depth $common_pbwt_depth \
        --pbwt-mac $common_pbwt_mac \
        --pbwt-mdr $pbwt_mdr \
        --pbwt-window $window --hmm-window $window --hmm-ne $hmm_ne \
        $haploid_arg \
        --region $region && \
    bcftools view --threads 8 -Ob $chrom_working_dir/noparents.bcf > $vcf_phased_no_parents_common_biallelic \
    && bcftools index -f --threads 8 $vcf_phased_no_parents_common_biallelic #\
fi

### Phase rare variants in chunks
i=1
for chrom_region in $(cat $chrom_regions)
    do
    if [ ! -s $chrom_working_dir/${i}_tmp_noparents.rare.bcf ]; then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase_no_parents \
            --map $chrom_map \
            --scaffold $vcf_phased_no_parents_common_biallelic \
            --thread $num_threads \
            --log $chrom_working_dir/${chrom}.${i}.rare.noparents.log \
            --pbwt-modulo $pbwt_modulo \
            --input-region $chrom_region \
            --scaffold-region $chrom_region \
            --effective-size $hmm_ne \
            $haploid_arg \
            --output $chrom_working_dir/${i}_tmp_noparents.rare.bcf && \
            bcftools index -f --threads 8 $chrom_working_dir/${i}_tmp_noparents.rare.bcf
    fi
    let i++
done

wait
# once that is done, make the reports and subsets that stem from the panels generated while we were phasing a no-parents panel

### phased_panel_vcf_2504_biallelic is what is used for imputation evaluation - get stats
if [ ! -s $phased_panel_vcf_2504_biallelic_variant_report ]; then
    echo "phased_panel_vcf_2504_biallelic_variant_report"
    bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAC\t%INFO/AN\t%INFO/MAF\t%INFO/SYNTENIC\n' $phased_panel_vcf_2504_biallelic > $phased_panel_vcf_2504_biallelic_variant_report &
fi

if [ ! -s $phased_panel_vcf_2504.csi ]; then
    echo "phased_panel_vcf_2504"
    bcftools view -Ou --threads 8 -S $unrelated_samples $phased_panel_vcf_3202 \
    | bcftools view -Oz --threads 8 -c 1:minor - > $phased_panel_vcf_2504 \
    && bcftools index --threads 8 -f $phased_panel_vcf_2504 &
fi

### Concat rare variant chunks
if [ ! -s ${vcf_phased_no_parents_rare_biallelic}.csi ]; then
    bcftools concat --threads 8 -Ou -l $chrom_working_dir/*_tmp_noparents.rare.bcf \
    | bcftools view -Ou --threads 2 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ou --threads 4 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools annotate -Ob --threads 4 -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT - \
    > $vcf_phased_no_parents_rare_biallelic \
    && bcftools index --threads 8 -f  $vcf_phased_no_parents_rare_biallelic &
fi

wait

## To evaluate 1KGP T2T performance as a reference panel when phasing variants,
## phase samples present in pangenome using the phased 2504 panel as a reference
### Remove all pangenome samples and parents of pangenome samples from phased 2504 panel.
### This will be our 'reference panel'
echo "creating 'ground truth' reference panels for phasing evaluation"
if [ ! -s $phased_panel_no_pangenome_biallelic.csi ]; then
    echo "phased_panel_no_pangenome_biallelic"
    bcftools view -Ou --threads 8 -S ^$pangenome_and_parents --force-samples $phased_panel_vcf_2504_biallelic \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC  > $phased_panel_no_pangenome_biallelic \
    && bcftools index --threads 8 -f $phased_panel_no_pangenome_biallelic &
fi

if [ ! -s $chr_specific_reference_pangenome_variation_trimmed_biallelic.csi ]; then
    echo "chr_specific_reference_pangenome_variation_trimmed_biallelic"
    ### Remove pangenome reference sites that are not present in 3202 biallelic reference.
    bcftools isec -r $region --threads 8 -o $chr_specific_reference_pangenome_variation_trimmed_biallelic -Ob -n =2 -w 1 $chr_specific_reference_pangenome_variation_biallelic $phased_panel_vcf_3202_biallelic \
    && bcftools index -f --threads 8 $chr_specific_reference_pangenome_variation_trimmed_biallelic &
fi


echo "preparing references for pangenome sample phasing"
if [ ! -s $vcf_to_phase_pangenome_biallelic.csi ]; then
    echo "vcf_to_phase_pangenome_biallelic"
    ### Create unphased variant call set of samples present in the pangenomepe
    ### also perform basic prephasing filtering to identify phasable variant set
    bcftools +setGT -Ou --threads 8 $chr_specific_reference_pangenome_variation_biallelic -- -t a -n u \
    | bcftools view -c 1:minor \
        -e "ALT=='*' || F_MISSING>0.05 || INFO/MAC==0" \
        --threads 8 -Ob - > $vcf_to_phase_pangenome_biallelic \
    && bcftools index -f --threads 8 $vcf_to_phase_pangenome_biallelic &
fi

wait

if [[ ! -s $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.csi ]]; then
    bcftools view --threads 8 -S $pangenome_samples --force-samples -Ou $vcf_phased_no_parents_rare_biallelic \
    | bcftools view -Ou --threads 8 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
    > $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf \
    && bcftools index --threads 8 -f $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf &
fi

## While doing this, gather AC_AN data to later generate reference MAF tables
if [[ ! -s $chrom_working_dir/${chrom}_3202_AC_AN.tsv ]]; then
    echo "3202"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $phased_panel_vcf_3202_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_3202_AC_AN.tsv &
fi
if [[ ! -s $chrom_working_dir/${chrom}_2504_AC_AN.tsv ]]; then
    echo "2504"    
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF ${phased_panel_vcf_2504_biallelic%%.vcf.gz}.bcf \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f '%ID\t%INFO/MAC\t%INFO/AN\n' - > $chrom_working_dir/${chrom}_2504_AC_AN.tsv  &
fi
if [[ ! -s $chrom_working_dir/${chrom}_2002_AC_AN.tsv ]]; then
    echo "2002"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $vcf_phased_no_parents_rare_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_2002_AC_AN.tsv &
fi
if [[ ! -s $chrom_working_dir/${chrom}_2430_AC_AN.tsv ]]; then
    echo "2430"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $phased_panel_no_pangenome_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_2430_AC_AN.tsv &
fi

# Identify trio-private singletons
if [[ ! -s $chrom_working_dir/${chrom}_private_singletons.txt ]]; then
    mkdir -p $chrom_working_dir/tmp
    cat pedigrees/duos_and_trios.txt \
    | parallel -j $num_threads "bcftools view --force-samples -H -G -s {} -x -c 2 ${phased_panel_vcf_3202_biallelic%%.vcf.gz}.bcf | cut -f 3 > $chrom_working_dir/tmp/{#}.txt" && \
    cat $chrom_working_dir/tmp/*.txt | sort | uniq > $chrom_working_dir/${chrom}_private_singletons.txt && rm -rf $chrom_working_dir/tmp &
fi

echo "phasing pangenome samples against 2504 panel"
### Repeat phasing as above, specifying a reference panel during the common variant phasing.
if [ ! -s $common_variants_phased_pangenome_against_ref_biallelic.csi ]; then
    ./SHAPEIT5_phase_common_static_v1.1.1 \
        --input $vcf_to_phase_pangenome_biallelic \
        --reference $phased_panel_no_pangenome_biallelic \
        --map $chrom_map \
        --output $chrom_working_dir/tmp_pangenome.bcf \
        --thread $num_threads \
        --log $chrom_working_dir/$chrom.common.pangenome_with_ref_panel.log \
        --filter-maf $rare_variant_threshold \
        --mcmc-iterations $mcmc_iteration_scheme \
        --pbwt-modulo $pbwt_modulo \
        --pbwt-depth $common_pbwt_depth \
        --pbwt-mac $common_pbwt_mac \
        --pbwt-mdr $pbwt_mdr \
        --pbwt-window $window --hmm-window $window --hmm-ne $hmm_ne \
        $haploid_arg \
        --region $region \
    && bcftools view --threads 8 -Ob $chrom_working_dir/tmp_pangenome.bcf > $common_variants_phased_pangenome_against_ref_biallelic \
    && bcftools index -f --threads 8 $common_variants_phased_pangenome_against_ref_biallelic
fi


## Phase rare variants in chunks. We cannot specify a reference panel in this step.
i=1
# for chrom_region in $(cat $chrom_regions)
# do
    if [ ! -s $chrom_working_dir/${i}_tmp_pangenome.rare.bcf.csi ]; then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase_pangenome_biallelic \
            --map $chrom_map \
            --scaffold $common_variants_phased_pangenome_against_ref_biallelic \
            --thread $num_threads \
            --log $chrom_working_dir/${chrom}.${i}.rare.pangenome_with_ref_panel.log \
            --pbwt-modulo $pbwt_modulo \
            --input-region $whole_chrom \
            --scaffold-region $whole_chrom \
            --effective-size $hmm_ne \
            $haploid_arg \
            --output $chrom_working_dir/${i}_tmp_pangenome.rare.bcf && \
            bcftools index --threads 8 -f $chrom_working_dir/${i}_tmp_pangenome.rare.bcf
    fi
#     let i++
# done

### Concat rare variant 
if [ ! -s $rare_variants_phased_pangenome_against_ref_biallelic.csi ]; then
    bcftools concat --threads 8 -Ou -l $chrom_working_dir/*_tmp_pangenome.rare.bcf > $rare_variants_phased_pangenome_against_ref_biallelic && \
    bcftools index --threads 8 -f $rare_variants_phased_pangenome_against_ref_biallelic &
fi


##############################################################
# Drop men from reference sets when working with X chromosome

mkdir -p $stats_dir

# Phased panel of unrelated samples generated from 3202 panel w/ pedigree:
if [[ $chrom == 'chrX' ]]
then
    echo "removing male samples from panels for evaluation"
    if [[ ! -s ${phased_panel_vcf_3202_biallelic%.*cf.gz}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S ^$male_sample_list $phased_panel_vcf_3202_biallelic > ${phased_panel_vcf_3202_biallelic%.*cf.gz}.females_only.vcf.gz && \
        bcftools index --threads 8 -f ${phased_panel_vcf_3202_biallelic%.*cf.gz}.females_only.vcf.gz && \
        phased_panel_vcf_3202_biallelic=${phased_panel_vcf_3202_biallelic%.*cf.gz}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_3202%.*cf.gz}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S ^$male_sample_list $phased_panel_vcf_3202 > ${phased_panel_vcf_3202%.*cf.gz}.females_only.vcf.gz && \
        bcftools index --threads 8 -f ${phased_panel_vcf_3202%.*cf.gz}.females_only.vcf.gz && \
        phased_panel_vcf_3202=${phased_panel_vcf_3202%.*cf.gz}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_2504_biallelic%.*cf.gz}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S ^$male_sample_list $phased_panel_vcf_2504_biallelic > ${phased_panel_vcf_2504_biallelic%.*cf.gz}.females_only.vcf.gz && \
        bcftools index --threads 8 -f ${phased_panel_vcf_2504_biallelic%.*cf.gz}.females_only.vcf.gz && \
        phased_panel_vcf_2504_biallelic=${phased_panel_vcf_2504_biallelic%.*cf.gz}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_2504%.*cf.gz}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S ^$male_sample_list $phased_panel_vcf_2504 > ${phased_panel_vcf_2504%.*cf.gz}.females_only.vcf.gz && \
        bcftools index --threads 8 -f ${phased_panel_vcf_2504%.*cf.gz}.females_only.vcf.gz && \
        phased_panel_vcf_2504=${phased_panel_vcf_2504%.*cf.gz}.females_only.vcf.gz &
    fi

    if [[ ! -s ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list $vcf_phased_no_parents_common_biallelic > ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf && \
        bcftools index --threads 8  -f ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf && \
        vcf_phased_no_parents_common_biallelic=${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf &
    fi

    if [[ ! -s ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list ${vcf_phased_no_parents_rare_biallelic%.*cf}.bcf > ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf && \
        bcftools index --threads 8 -f ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf && \
        vcf_phased_no_parents_rare_biallelic=${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf &
    fi
    
    if [[ ! -s ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf}.females_only.bcf.csi ]]; then
        bcftools view -Ob --force-samples --threads 8 -S ^$male_sample_list ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf}.bcf > ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf}.females_only.bcf && \
        bcftools index --threads 8 -f ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf}.females_only.bcf && \
        chr_specific_reference_pangenome_variation_trimmed_biallelic=${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf}.females_only.bcf &
    fi

    wait
fi

wait

##############################################################
# Evaluate accuracy
recalc_phasing_stats='true'
if [[ "$recalc_phasing_stats" == 'true' ]]; then
    echo "evaluating phasing accuracy"

    # # Evaluate accuracy of 3202 panel via two methods:
    # # 1) by looking at within-trio phasing consistency per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing
    # #  1.1) Trio consistency when only probands + unrelated are phased (no parental genomes leak into the rest of the data)
    # ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
    #                                 --estimation $vcf_phased_no_parents_rare_biallelic \
    #                                 -P $pedigree -R $whole_chrom --singleton \
    #                                 --log $chrom_working_dir/rare_noparents_vs_trios_${chrom}.log \
    #                                 --output $stats_dir/rare_noparents_vs_trios_${chrom} 2> /dev/null &

    # #  1.2) by looking at within-trio phasing consistency - Full panel, trio + statistically phased.
    # #       Basically confirming that trio-phasing worked, best assessment of trio-phased sample accuracy
    # ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
    #                                 --estimation $phased_panel_vcf_3202_biallelic \
    #                                 -P $pedigree -R $whole_chrom --singleton \
    #                                 --log $chrom_working_dir/3202_panel_vs_trios_${chrom}.log \
    #                                 --output $stats_dir/3202_panel_vs_trios_${chrom} 2> /dev/null &

    # #  1.3) by looking at within-trio phasing consistency - 2504 panel, trio + statistically phased.
    # #       Should be the same as 1.2, but summary statistics will not include children. Parents are phased with a mix of trio-consistency and statistical phasing.
    # ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
    #                                 --estimation $phased_panel_vcf_2504_biallelic \
    #                                 -P $pedigree -R $whole_chrom --singleton \
    #                                 --log $chrom_working_dir/2504_panel_vs_trios_${chrom}.log \
    #                                 --output $stats_dir/2504_panel_vs_trios_${chrom} 2> /dev/null &

    #  2) by looking at phasing consistency with ground truth pangenome samples
    #   2.1) trio-phased 3202 panel
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HPRC_${chrom}_$4.log \
                                    --output $stats_dir/3202_panel_vs_HPRC_${chrom}_$4 2> /dev/null &


    #  2.2) evaluate phasing performance of phasing without pedigree against ground truth pangenome samples
    #        Able to compare trio-measured SER and HPRC consistency based SER
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/noparents_vs_HPRC_${chrom}_$4.log \
                                    --output $stats_dir/noparents_vs_HPRC_${chrom}_$4 2> /dev/null &

    #  3) Evaluate panel's usefulness as a reference panel
    #      3.1) Rephased pangenome samples compared to HPRC samples
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $rare_variants_phased_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/rare_pangenome_panelphased_vs_pangenome_${chrom}_$4.log \
                                    --output $stats_dir/rare_pangenome_panelphased_vs_pangenome_${chrom}_$4 2> /dev/null &
    
    #      3.2) Rephased pangenome samples compared to trio samples
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $rare_variants_phased_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/rare_pangenome_panelphased_vs_trios_${chrom}_$4.log \
                                    --output $stats_dir/rare_pangenome_panelphased_vs_trios_${chrom}_$4 2> /dev/null &

    # 4) Experiment with HGSVC samples
    #     4.1) trio-phased 3202 panel
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_HGSVC_variation_biallelic \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_${chrom}_$4.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_${chrom}_$4 2> /dev/null &

    #     4.2) only examine samples that are not part of trios
    bcftools view --threads 4 -S $basedir/sample_subsets/HGSVC_not_part_of_trio.txt -Ob -o ${chr_specific_reference_HGSVC_variation_biallelic%%.bcf}.non_trio_samples.bcf $chr_specific_reference_HGSVC_variation_biallelic && \
    bcftools index --threads 4 ${chr_specific_reference_HGSVC_variation_biallelic%%.bcf}.non_trio_samples.bcf && \
    ./SHAPEIT5_switch_static_JLL --validation ${chr_specific_reference_HGSVC_variation_biallelic%%.bcf}.non_trio_samples.bcf \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_notrios_${chrom}_$4.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_notrios_${chrom}_$4 2> /dev/null &

fi


# ################### Examine imputation ##################

# # Liftover panels if necessary
# if [[ ! -s $T2T_native_panel.csi ]] && [[ -s $phased_panel_vcf_2504_biallelic.tbi ]]; then
#         bcftools view -Ob --threads 4 $phased_panel_vcf_2504_biallelic > $T2T_native_panel && bcftools index -f --threads 4 $T2T_native_panel \
#     && ./liftover_panel.sh $T2T_native_panel $GRCh38_lifted_panel /dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta /dev/shm/chm13v2-hg38.over.chain liftover/grch38-chm13v2.sort.vcf.gz &
# elif [[ ! -s $T2T_lifted_panel.csi ]]; then
#     ./liftover_panel.sh $T2T_native_panel $GRCh38_lifted_panel /dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta /dev/shm/chm13v2-hg38.over.chain liftover/grch38-chm13v2.sort.vcf.gz &
# fi

# if [[ ! -s $GRCh38_native_panel.csi ]] && [[ -s ${GRCh38_native_panel%%.bcf}.vcf.gz.tbi ]]; then
#         bcftools view -Ob --threads 4 ${GRCh38_native_panel%%.bcf}.vcf.gz > $GRCh38_native_panel && bcftools index -f --threads 4 $GRCh38_native_panel \
#     && ./liftover_panel.sh $GRCh38_native_panel $T2T_lifted_panel /dev/shm/chm13v2.0_maskedY_rCRS.fasta /dev/shm/hg38-chm13v2.over.chain liftover/chm13v2-grch38.sort.vcf.gz &
# elif [[ ! -s $GRCh38_lifted_panel.csi ]]; then
#     ./liftover_panel.sh $GRCh38_native_panel $T2T_lifted_panel /dev/shm/chm13v2.0_maskedY_rCRS.fasta /dev/shm/hg38-chm13v2.over.chain liftover/chm13v2-grch38.sort.vcf.gz &
# fi


# wait


# if [[ $chrom == "PAR1" ]]
# then    
#     SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
#     SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chrX.recalibrated.snp_indel.pass.vcf.gz
# elif [[ $chrom == "PAR2" ]]
# then
#     SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
#     SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chrX.recalibrated.snp_indel.pass.vcf.gz
# elif [[ $chrom == "debug" ]]
# then
#     SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chr20.recalibrated.snp_indel.pass.vcf.gz
#     SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chr20.recalibrated.snp_indel.pass.vcf.gz
# else
#     SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.vcf.gz
#     SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/${chrom}.recalibrated.snp_indel.pass.vcf.gz
# fi

# if [ ! -s $SGDP_ground_truth_T2T.csi ]; then
#     echo $chrom "SGDP_ground_truth_T2T"
#     bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m -any $SGDP_variants_T2T \
#     | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
#                     -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
#                     -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
#     | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - \
#     | bcftools view -Ou  -c 1:minor -i "ALT!='*' || F_MISSING<0.05 || ABS(ILEN)<=50" - \
#     | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
#     > $SGDP_ground_truth_T2T \
#     && bcftools index --threads 8 -f $SGDP_ground_truth_T2T &
# fi
# if [ ! -s $SGDP_ground_truth_GRCh38.csi ]; then
#     echo $chrom "SGDP_ground_truth_GRCh38"
#     bcftools norm --threads 2 -Ou -r $grch38_region -f $GRCh38_fasta -m -any $SGDP_variants_GRCh38 \
#     | bcftools annotate -Ou --threads 8 -a $grch38_syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
#                     -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
#                     -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
#     | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - \
#     | bcftools view -Ou  -c 1:minor -i "ALT!='*' || F_MISSING<0.05 || ABS(ILEN)<=50" - \
#     | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
#     > $SGDP_ground_truth_GRCh38 \
#     && bcftools index --threads 8 -f $SGDP_ground_truth_GRCh38 &
# fi


# wait

# if [[ ! -s $chrom_working_dir/GRCh38_imputation_workspace/GRCh38_space_filtered/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel false SGDP
# fi
# if [[ ! -s $chrom_working_dir/T2T_imputation_workspace/T2T_space_filtered/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then    
#     ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel false SGDP
# fi
# if [[ ! -s $chrom_working_dir/GRCh38_snps_imputation_workspace/GRCh38_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel true SGDP
# fi
# if [[ ! -s $chrom_working_dir/T2T_snps_imputation_workspace/T2T_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel true SGDP
# fi

# if [[ ! -s $GRCh38_lifted_panel_no_pangenome.csi ]]; then
#     echo $chrom "GRCh38_lifted_panel_no_pangenome"
#     bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_lifted_panel > $GRCh38_lifted_panel_no_pangenome && \
#     bcftools index --threads 2 $GRCh38_lifted_panel_no_pangenome &
# fi
# if [[ ! -s $T2T_lifted_panel_no_pangenome.csi ]]; then
#     echo $chrom "T2T_lifted_panel_no_pangenome"
#     bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $T2T_lifted_panel > $T2T_lifted_panel_no_pangenome && \
#     bcftools index --threads 2 $T2T_lifted_panel_no_pangenome &
# fi
# if [[ ! -s $GRCh38_native_panel_no_pangenome.csi ]]; then
#     echo $chrom "GRCh38_native_panel_no_pangenome"
#     bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_native_panel > $GRCh38_native_panel_no_pangenome && \
#     bcftools index --threads 2 $GRCh38_native_panel_no_pangenome
# fi

# wait
# if [[ ! -s $chrom_working_dir/GRCh38_imputation_pangenome_workspace/GRCh38_space_filtered/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_pangenome_workspace      $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome false pangenome &
# fi
# if [[ ! -s $chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace/GRCh38_space_filtered_snpsOnly/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome true pangenome &
# fi
# wait
# if [[ ! -s $chrom_working_dir/T2T_imputation_pangenome_workspace/T2T_space_filtered/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then    
#     ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_pangenome_workspace            $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome false pangenome &
# fi
# if [[ ! -s $chrom_working_dir/T2T_snps_imputation_pangenome_workspace/T2T_space_filtered_snpsOnly/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
#     ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_pangenome_workspace       $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome true pangenome &
# fi

wait
# close logfile
set +x
exec 19>&-
