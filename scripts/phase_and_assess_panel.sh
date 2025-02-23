#!/usr/bin/env bash
set -eo pipefail

# set up script log
## one-liner to get script location; credit to https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
basedir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
logfile=$basedir/phasing_$4_$1_$3.log
exec 19>$logfile
export BASH_XTRACEFD=19
set -x # writes commands to logfile

chrom=$1
num_threads=$2
suffix=$3
genome=$4
atomize=$5
missing_to_ref=$6
missing_filter_cutoff=$7
trim_assemblies_to_callset=$8
filter_multiallelic_indels=$9
rephase_high_accuracy=${10}
manual_maf_for_common_rephasing=${11}
# which_pangenome=${11}


#settings
hmm_ne=135000
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

if [[ $manual_maf_for_common_rephasing == 'false' ]]
then
    minimum_common_MAC=1
elif [[ -z $manual_maf_for_common_rephasing ]]
then
    minimum_common_MAC=3
else
    minimum_common_MAC=$manual_maf_for_common_rephasing
fi

## process keywords
if [[ -z $missing_filter_cutoff ]]; then
    graph_reference_missing_cutoff=0.2
else
    graph_reference_missing_cutoff=$missing_filter_cutoff
fi

if [[ $atomize == 'true' ]]; then
    atomize_operation="bcftools norm --threads 2 -Ou --atomize --atom-overlaps . -m +snps - "
else
    atomize_operation="cat"
fi

if [[ $missing_to_ref == 'true' ]]; then
    fill_missing_to_ref(){
        sed 's,\.|\.,qqq,g' | sed 's,\.|,0|,g' | sed 's,|\.,|0,g' | sed 's,qqq,\.|\.,g'
    }

    missing_to_ref_operation=fill_missing_to_ref
else
    missing_to_ref_operation="cat"
fi


if [[ $filter_multiallelic_indels == 'true' ]]; then
    multiallelic_indel_filter="-e (N_ALT>1)&(TYPE~'indel')"
else
    multiallelic_indel_filter="-e ''"
fi


if [[ $rephase_high_accuracy == 'true' ]]; then
    non_default_shapeit_common_settings=" --mcmc-iterations $mcmc_iteration_scheme --pbwt-modulo $pbwt_modulo \
                                          --pbwt-depth $common_pbwt_depth --pbwt-mac $common_pbwt_mac --pbwt-mdr $pbwt_mdr \
                                          --pbwt-window $window --hmm-window $window --hmm-ne $hmm_ne"
    non_default_shapeit_rare_settings=" --pbwt-modulo $pbwt_modulo --effective-size $hmm_ne"
    non_default_impute5_settings="--ne $hmm_ne"
else
    non_default_shapeit_common_settings=""
    non_default_shapeit_rare_settings=""
    non_default_impute5_settings=""
fi



if [[ suffix != '' ]]
then
    suffix="_$suffix"
fi

set -u

if [[ $genome == 'GRCh38' ]]; then
    drop_reference='CHM13'
    native_maps_insert=''
    chrom_map=$basedir/GRCh38_performance_comparison/hg38_chrom_maps/${chrom}.b38.gmap.gz
    initial_vcf_calls_folder=$basedir/unphased_GRCh38_panel

    ref_fasta=$basedir/GRCh38_full_analysis_set_plus_decoy_hla.fasta
    pangenome_vcf=$basedir/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
    HGSVC_vcf=$basedir/hgsvc3-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz
    HGSVC_HPRC_vcf=$basedir/hgsvc3-hprc-2024-02-23-mc-chm13.GRCh38-vcfbub.a100k.wave.norm.vcf.gz

    chrom_chunking_coords=$basedir/GRCh38_performance_comparison/regions2.txt
    syntenic_site_location=$basedir/hg38.GCA_009914755.4.synNet.summary.bed.gz

    chrom_specific_syntenic_annotation_line_part1="CHROM,FROM,TO,SYNTENIC"
    chrom_specific_syntenic_annotation_line_part2="##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description=\"Syntenic with CHM13\">"

    PAR1_start='10000'
    PAR1_end='2781479'
    PAR2_start='155701383'
    PAR2_end='156030895'


    if [[ $chrom == *PAR* ]]; then
        phased_panel_vcf=$basedir/phased_GRCh38_panel/1KGP.GRCh38.$(echo $chrom | sed 's/chr//').recalibrated.snp_indel.pass.phased.3202.vcf.gz
        input_vcf=$initial_vcf_calls_folder/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.vcf.gz
    else
        phased_panel_vcf=$basedir/phased_GRCh38_panel/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.3202.vcf.gz
        input_vcf=$initial_vcf_calls_folder/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.recalibrated_variants.vcf.gz
    fi

elif [[ $genome == 'CHM13v2.0' ]]; then
    drop_reference='GRCh38'
    native_maps_insert='native_maps.'
    chrom_map=$basedir/t2t_maps_no_filtered_regions/scaled_fixedamount/${chrom}_noMask.scaled.gmap.gz
    initial_vcf_calls_folder=/mnt/ssd/lalli/phasing_T2T

    ref_fasta=/dev/shm/chm13v2.0_maskedY_rCRS.fasta
    pangenome_vcf=$basedir/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
    HGSVC_vcf=$basedir/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
    HGSVC_HPRC_vcf=$basedir/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz

    chrom_chunking_coords=$basedir/regions2.txt
    syntenic_site_location="$basedir/chm13v2-syntenic_to_hg38.bed"
    chrom_specific_syntenic_annotation_line_part1="CHROM,FROM,TO"
    chrom_specific_syntenic_annotation_line_part2="##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description=\"Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)\">"

    PAR1_start='0'
    PAR1_end='2394410'
    PAR2_start='153925834'
    PAR2_end='154259566'

    if [[ $chrom == *PAR* ]]; then
        input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    else
        input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.$chrom.recalibrated.snp_indel.pass.vcf.gz
    fi
else
    echo "genome must be either CHM13v2.0 or GRCh38. Sorry. No refunds."
    return 1
fi


if [[ $filter_multiallelic_indels == 'true' ]]; then
    if [[ ! -s ${pangenome_vcf%%.vcf.gz}.no_tandem_repeats.bcf ]]; then
        bcftools annotate --threads 4 -Ou -x INFO/AC,INFO/AN,INFO/AF $pangenome_vcf \
        | bcftools norm --threads 4 -m +any - \
        | bcftools view --threads 8 "-e (N_ALT>1)&(TYPE~'indel')" -Ob -W \
        -o ${pangenome_vcf%%.vcf.gz}.no_tandem_repeats.bcf && \
        pangenome_vcf=${pangenome_vcf%%.vcf.gz}.no_tandem_repeats.bcf &
    fi

    if [[ ! -s ${HGSVC_vcf%%.vcf.gz}.no_tandem_repeats.bcf ]]; then
        bcftools annotate --threads 4 -Ou -x INFO/AC,INFO/AN,INFO/AF $HGSVC_vcf \
        | bcftools norm --threads 4 -m +any - \
        | bcftools view --threads 8 "-e (N_ALT>1)&(TYPE~'indel')" -Ob -W \
        -o ${HGSVC_vcf%%.vcf.gz}.no_tandem_repeats.bcf && \
        HGSVC_vcf=${HGSVC_vcf%%.vcf.gz}.no_tandem_repeats.bcf &
    fi

    wait
fi

## Identify if we are dealing with a PAR region
# GRCh38: "X:10001-2781479" "X:2781480-155701382" "X:155701383-156030895"
# T2T: "X:0-2394410" "X:2394411-153925834" "X:153925835-154259566"


if [[ $chrom == "chrPAR1" || $chrom == "PAR1" ]] 
then
    chrom='PAR1'
    region="chrX:$(($PAR1_start+1))-$PAR1_end"
    grch38_region="chrX:10001-2781479"
    whole_chrom=$region
    end_chrom=$PAR1_end
elif [[ $chrom == "chrPAR2" || $chrom == "PAR2" ]]
then
    chrom='PAR2'
    region="chrX:$(($PAR2_start+1))-$PAR2_end"
    grch38_region="chrX:155701384-156030895"
    whole_chrom=$region
    end_chrom=$PAR2_end
elif [[ $chrom == 'chrdebug' ]]
then
    chrom='debug'
    region="chr20:10000000-12000000"
    whole_chrom=$region
    grch38_region=$region
    chrom_map=$basedir/t2t_maps_no_filtered_regions/scaled_fixedamount/chr20_noMask.scaled.gmap.gz
    end_chrom=12000000
elif [[ $chrom == 'chrX' ]]
then
    region="chrX:$(($PAR1_end+1))-$PAR2_start"
    whole_chrom=$region
    grch38_region="chrX:2781480-155701383"
    use_beagle='true'
    hmm_ne=$(awk "BEGIN { printf \"%.0f\", ($hmm_ne*.75) }")
    end_chrom=$PAR2_start
else
    region=$chrom
    end_chrom=$(grep "^$chrom	" ${ref_fasta}.fai | cut -f 2)
    whole_chrom=$chrom:1-$end_chrom
    grch38_region=$chrom
fi

chrom_working_dir=$basedir/testing_phasing_params2/${chrom}_working${suffix}

final_panel_dir=$basedir/testing_phasing_params2/phased_${genome}_panel${suffix}
stats_dir=$basedir/testing_phasing_params2/phasing_stats_${genome}${suffix}_custom_switch
imputation_results_dir=$basedir/testing_phasing_params2/imputation_results$suffix

mkdir -p $chrom_working_dir
mkdir -p $final_panel_dir
mkdir -p $stats_dir


population_ids=$basedir/sample_subsets/unrelated_superpopulations.csv


# move command logfile to newly created working directory
logfile=$chrom_working_dir/phasing_${genome}_$1_$3.log
cp $basedir/phasing_${genome}_$1_$3.log $logfile
exec 19>$logfile
export BASH_XTRACEFD=19
set -x # writes commands to logfile
rm $basedir/phasing_${genome}_$1_$3.log


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
else
    haploid_arg=""
fi

# Define names of chrom specific vcfs that we will be phasing
vcf_to_phase=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.biallelic.bcf
fully_annotated_input_variants=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.bcf

chr_specific_reference_pangenome_variation_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.biallelic.bcf
chr_specific_reference_HGSVC_variation_biallelic=$chrom_working_dir/${chrom}_reference_HGSVC.biallelic.bcf
chr_specific_reference_HGSVC_HPRC_variation_biallelic=$chrom_working_dir/${chrom}_reference_HGSVC_HPRC.biallelic.bcf
chr_specific_reference_pangenome_variation_trimmed_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.filtered_variants.biallelic.bcf

# Define the names of input variant files that are sample subsets
vcf_to_phase_pangenome_biallelic_HPRC_common=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.HPRC_pangenome_calls.common.biallelic.bcf
vcf_to_phase_pangenome_biallelic_HPRC=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.HPRC_pangenome_calls.biallelic.bcf
vcf_to_phase_pangenome_biallelic_1kgp_common=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.1kGP_pangenome_calls.common.biallelic.bcf
vcf_to_phase_pangenome_biallelic_1kgp=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.1kGP_pangenome_calls.biallelic.bcf
vcf_to_phase_pangenome_biallelic_1kgp_intersect=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.1kGP_pangenome_calls.intersect_panel.biallelic.bcf
vcf_to_phase_no_parents=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unphased.noparents.bcf


# Define phased result file names
common_variants_phased_ped=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}common.bcf
rare_variants_phased_ped=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}rare.bcf
rare_variants_phased_ped_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}biallelic.rare.bcf
phased_panel_no_pangenome_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}no_pangenome.biallelic.bcf

common_variants_phased_HPRC_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}common.HPRC_pangenome_calls.biallelic.bcf
rare_variants_phased_HPRC_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}biallelic.rare.HPRC_pangenome_calls.bcf
common_variants_phased_1kgp_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}common.1kGP_pangenome_calls.biallelic.bcf
rare_variants_phased_1kgp_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.phased.${native_maps_insert}biallelic.rare.1kGP_pangenome_calls.bcf

vcf_phased_no_parents_common_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.biallelic.bcf
vcf_phased_no_parents_rare_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.biallelic.bcf



# output files
if [[ $genome == 'GRCh38' ]]; then
    phased_panel_vcf_3202=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.3202.bcf
    phased_panel_vcf_2504=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.2504.bcf
    phased_panel_vcf_3202_biallelic=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.${native_maps_insert}biallelic.3202.bcf
    phased_panel_vcf_2504_biallelic=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.${native_maps_insert}biallelic.2504.bcf
else
    phased_panel_vcf_3202=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf
    phased_panel_vcf_2504=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.2504.bcf
    phased_panel_vcf_3202_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf
    phased_panel_vcf_2504_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf
fi
fully_annotated_input_variant_report=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv
# phased_panel_vcf_2504_biallelic_variant_report=$stats_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.2504.stats.tsv

# pedigrees
pedigree=$basedir/pedigrees/1kgp.ped
bcftools_formatted_pedigree=$basedir/pedigrees/trios_only.ped

# lists of different categories of samples
no_parents=$basedir/sample_subsets/not_parents.txt
unrelated_samples=$basedir/sample_subsets/unrelated_samples.txt
pangenome_samples=$basedir/sample_subsets/pangenome_samples.txt
pangenome_and_parents=$basedir/sample_subsets/pangenome_samples_and_parents.txt
SGDP_in_1KGP=$basedir/sample_subsets/SGDP_samples_in_1KGP_numbered_format.txt
female_samples=$basedir/sample_subsets/females.txt

# #### making symlinks to skip re-phasing. Do not use for production.

# if [[ $genome == 'CHM13v2.0' ]]; then
#     prior_working_dir=$basedir/${chrom}_working_T2T_scaled_newimpute_092424
#     prior_panel_dir=$basedir/phased_T2T_panel_T2T_scaled_newimpute_092424
#     prior_stats_dir=$basedir/phasing_stats_T2T_scaled_newimpute_092424_custom_switch
#     prior_imputation_results_dir=$basedir/imputation_results_T2T_scaled_newimpute_092424
# else
#     prior_working_dir=$basedir/${chrom}_working_GRCh38_statistics_092424
#     prior_stats_dir=$basedir/phasing_stats_${genome}_GRCh38_statistics_092424_custom_switch
#     prior_imputation_results_dir=$basedir/imputation_results_GRCh38_statistics_092424
# fi

# # ln -sf $prior_working_dir/${chrom}_private_singletons.txt $chrom_working_dir/${chrom}_private_singletons.txt

# if [[ $filter_multiallelic_indels != 'true' ]]; then
#     ln -sf $prior_working_dir/$(basename $fully_annotated_input_variants)* $chrom_working_dir/ 
#     ln -sf $prior_working_dir/$(basename $fully_annotated_input_variant_report)* $chrom_working_dir/ 
#     ln -sf $prior_working_dir/${chrom}_private_singletons.txt $chrom_working_dir/${chrom}_private_singletons.txt
    
#     if [[ $genome == 'CHM13v2.0' ]]; then

#         ln -sf $prior_working_dir/$(basename $vcf_to_phase)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $common_variants_phased_ped)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $rare_variants_phased_ped)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/*[0-9].rare.bcf* $chrom_working_dir/
#         ln -sf $prior_working_dir/*tmp*.rare.bcf* $chrom_working_dir/
#         ln -sf $prior_working_dir/$(basename $rare_variants_phased_ped_biallelic)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $phased_panel_no_pangenome_biallelic)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $fully_annotated_input_variant_report)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/${chrom}_*_AC_AN.tsv $chrom_working_dir/ 

#         ln -sf $prior_panel_dir/$(basename $phased_panel_vcf_3202)* $final_panel_dir/ 
#         ln -sf $prior_panel_dir/$(basename $phased_panel_vcf_2504)* $final_panel_dir/ 
#         ln -sf $prior_panel_dir/$(basename $phased_panel_vcf_3202_biallelic)* $final_panel_dir/ 
#         ln -sf $prior_panel_dir/$(basename $phased_panel_vcf_2504_biallelic)* $final_panel_dir/ 

#         # ln -sf $prior_stats_dir/$(basename $phased_panel_vcf_2504_biallelic_variant_report) $stats_dir/ 
    
#         ln -sf $prior_working_dir/$(basename $vcf_to_phase_no_parents)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $vcf_phased_no_parents_common_biallelic)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/$(basename $vcf_phased_no_parents_rare_biallelic)* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/noparents.bcf* $chrom_working_dir/ 
#         ln -sf $prior_working_dir/*_tmp_noparents.rare.bcf* $chrom_working_dir/

#     fi
# fi

####################

## VARIANT FILTERING, QC VARIANT SUBSETTING/FORMATTING

# Create ground truth variation vcf from pangenome
## Note: Throughout, I am comparing this vcf of genome assemblies to variant calls.
## One difference between assemblies and calls is that there is no such thing as a missing variant in an assembly.
## What is represented as a missing variant (overlapping indel) is consistently represented as a reference allele in our variant calls
## So we will convert pangenome missing variants into pangenome reference alleles.
if [ ! -s $chr_specific_reference_pangenome_variation_biallelic.csi ]
then
    echo "making reference pangenome variation"
    echo $atomize_operation
    if [[ $chrom == 'chrX' ]] || [[ $chrom == 'PAR2' ]] || [[ $chrom == 'PAR1' ]]; then    # make missing/haploid into haploid
        bcftools view --threads 8 -r $region -s ^$drop_reference -Ou $pangenome_vcf \
        | bcftools annotate -Ou -x INFO/AT - \
        | $atomize_operation \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|\.,qqq,g' | sed 's,\.|,,g' | sed 's,|\.,,g' | sed 's,qqq,\.|\.,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools +setGT -Ou - -- -t a -n p \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_pangenome_variation_biallelic \
        && bcftools index $chr_specific_reference_pangenome_variation_biallelic &
    
    else
        bcftools view --threads 8 -r $region -s ^$drop_reference -Ou $pangenome_vcf \
        | bcftools annotate -Ou -x INFO/AT - \
        | $atomize_operation \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | $missing_to_ref_operation \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_pangenome_variation_biallelic \
        && bcftools index $chr_specific_reference_pangenome_variation_biallelic &
    fi
fi

if [ ! -s $chr_specific_reference_HGSVC_variation_biallelic.csi ]
then
    echo "making HGSVC variation"
    # HGSVC generated vcfs put PARs on ChrX/ChrY instead of combining on chrX. That means only females will be appropriately heterozygous in this region. Thankfully, references are not in the female sample list, so we can drop those at the same time!
    if [[ $chrom == 'chrX' ]] || [[ $chrom == 'PAR2' ]] || [[ $chrom == 'PAR1' ]]; then    # make missing/haploid into haploid
        bcftools view --threads 8 -r $region -S $female_samples --force-samples -Ou $HGSVC_vcf 2> /dev/null \
        | bcftools annotate -Ou -x INFO/AT - \
        | $atomize_operation \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|\.,qqq,g' | sed 's,\.|,,g' | sed 's,|\.,,g' | sed 's,qqq,\.|\.,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools +setGT -Ou - -- -t a -n p \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_HGSVC_variation_biallelic \
        && bcftools index $chr_specific_reference_HGSVC_variation_biallelic &

    else
        bcftools view --threads 8 -r $region -s ^$drop_reference -Ou $HGSVC_vcf \
        | bcftools annotate -Ou -x INFO/AT - \
        | $atomize_operation \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | $missing_to_ref_operation \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools view --threads 1 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_HGSVC_variation_biallelic \
        && bcftools index $chr_specific_reference_HGSVC_variation_biallelic &
    fi
fi

if [ ! -s $chr_specific_reference_HGSVC_HPRC_variation_biallelic.csi ]; then
    #HGSVC generated vcfs put PARs on ChrX/ChrY instead of combining on chrX. That means only females will be appropriately heterozygous in this region. Thankfully, references are not in the female sample list, so we can drop those at the same time!
    if [[ $chrom == 'chrX' ]] || [[ $chrom == 'PAR2' ]] || [[ $chrom == 'PAR1' ]]; then    # make missing/haploid into haploid
        bcftools view --threads 8 -r $region -S $female_samples --force-samples -Ou $HGSVC_HPRC_vcf 2> /dev/null \
        | bcftools annotate --threads 4 -Ou -x INFO/AT - \
        | bcftools norm --threads 2 -Ou --atomize --atom-overlaps . -m +snps - \
        | bcftools norm --threads 2 -Ov -f $ref_fasta -m -any - \
        | sed 's,\.|\.,qqq,g' | sed 's,\.|,,g' | sed 's,|\.,,g' | sed 's,qqq,\.|\.,g' \
        | bcftools view -Ou -i "F_MISSING<$graph_reference_missing_cutoff" - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools +setGT -Ou - -- -t a -n p \
        | bcftools view --threads 2 -Ou -c 1:minor - \
        | bcftools sort -m 40G -T $PWD -Ob > $chr_specific_reference_HGSVC_HPRC_variation_biallelic \
        && bcftools index --threads 4 $chr_specific_reference_HGSVC_HPRC_variation_biallelic &

    else
        bcftools view --threads 8 -r $region -s ^$drop_reference -Ou $HGSVC_HPRC_vcf  \
        |  bcftools annotate --threads 4 -Ou -x INFO/AT -  \
        |  bcftools norm --threads 4 -Ou --atomize --atom-overlaps . -m +snps -  \
        |  bcftools norm --threads 4 -Ou -f $ref_fasta -m -any - \
        |  bcftools view --threads 4 -Ou -i "F_MISSING<$graph_reference_missing_cutoff" -  \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF,INFO/MISSING --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        |  bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC,MISSING:1=F_MISSING \
        |  bcftools +setGT -Ou --threads 8 - -- -t a -n p \
        |  bcftools view --threads 2 -Ou -c 1:minor -  \
        |  bcftools sort -m 40G -Ob - > $chr_specific_reference_HGSVC_HPRC_variation_biallelic \
        && bcftools index --threads 4 $chr_specific_reference_HGSVC_HPRC_variation_biallelic &
    fi
fi

# Split multiallelic sites, filter sites using criteria described above,
# convert data to bcf format, and index.
echo "annotating variants"
if [ ! -s $fully_annotated_input_variants.csi ]
then
    echo "fully_annotated_input_variants"
    bcftools view --threads 4 -Ou "$multiallelic_indel_filter" -r $region $input_vcf \
    | bcftools norm --threads 4 -Ou -f $ref_fasta -m -any - \
    | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +mendelian2 -Ou - --ped $bcftools_formatted_pedigree -m a -m d \
    | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,F_MISSING,HWE,MAC:1=MAC \
    | bcftools +fill-tags --threads 8 -Oz - -- -S $population_ids -t HWE \
    | tee ${fully_annotated_input_variants%%.bcf}.vcf.gz \
    | bcftools view --threads 8 -Ob - \
        > $fully_annotated_input_variants \
    && bcftools index --threads 8 -f -t ${fully_annotated_input_variants%%.bcf}.vcf.gz \
    && bcftools index --threads 8 -f $fully_annotated_input_variants &
fi

wait

if [ ! -s $vcf_to_phase.csi ]
then
    echo "making vcf_to_phase: $vcf_to_phase"
    if [[ $genome == 'CHM13v2.0' ]]
    then
        echo "filtering annotated variants"
        if [[ $chrom != 'chrX' ]]
        then
            bcftools view \
                    -e "(TYPE!='snp' && (ABS(ILEN) >= 50)) || ALT=='*' || INFO/VQSLOD < 0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || INFO/MAC==0 || ( INFO/HWE_EUR<1e-10 && INFO/HWE_AFR<1e-10 && INFO/HWE_EAS<1e-10 && INFO/HWE_AMR<1e-10 && INFO/HWE_SAS<1e-10 ) || FILTER!='PASS'" \
                    --threads 8 -Ou $fully_annotated_input_variants \
            | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS - > $vcf_to_phase \
            && bcftools index --threads 8 -f $vcf_to_phase
        else
            bcftools view \
                    -e "(TYPE!='snp' && (ABS(ILEN) >= 50)) || ALT=='*' || INFO/VQSLOD < 0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || INFO/MAC==0 || ( INFO/HWE_EUR<1e-10 && INFO/HWE_AFR<1e-10 && INFO/HWE_EAS<1e-10 && INFO/HWE_AMR<1e-10 && INFO/HWE_SAS<1e-10 ) || FILTER!='PASS'" \
                    --threads 8 -Ou $fully_annotated_input_variants \
            | bcftools +fixploidy --threads 2 -Ou - -- -f 2 \
            | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS - > $vcf_to_phase \
            && bcftools index --threads 8 -f $vcf_to_phase
        fi
    else  # GRCh38
        echo "extracting variants from GRCh38 panel; dropping * alleles and SVs"

        bcftools view --threads 4 -Ou -r $region "$multiallelic_indel_filter" $phased_panel_vcf \
        | bcftools +setGT -Ou --threads 8 - -- -t a -n u \
        | bcftools norm --threads 4 -Ou -f $ref_fasta -m -any - \
        | bcftools view --threads 4 -Ou -e "(TYPE!='snp' && (ABS(ILEN) >= 50)) || ALT=='*' || INFO/MAC==0" - \
        | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS \
        --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - > $vcf_to_phase \
        && bcftools index --threads 8 -f $vcf_to_phase

    fi
fi

wait

# While the phasing is running, create a tsv table of the to-be-phased variants and quality metrics for downstream QC
if [[ ! -s $fully_annotated_input_variant_report ]]
then
    echo "making variant report"
    if [ ! -s ${fully_annotated_input_variants%%.bcf}.vcf.gz.tbi ]
    then
        echo "fully_annotated_input_variants"
        bcftools view --threads 2 -Ou -r $region "$multiallelic_indel_filter" $input_vcf \
        | bcftools norm --threads 4 -Ou -f $ref_fasta -m -any - \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location --mark-sites +SYNTENIC \
                            -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
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
if [[ $genome != 'GRCh38' ]]; then
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
        | bcftools annotate -a $syntenic_site_location --mark-sites +SYNTENIC \
                            -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                            -Ou -x INFO/MAC --threads 2 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
        > $rare_variants_phased_ped_biallelic \
        && bcftools index --threads 8 -f $rare_variants_phased_ped_biallelic
    fi
else ## If GRCh38, just process phased panel
    if [ ! -s $rare_variants_phased_ped_biallelic.csi ]; then
        echo "formatting GRCh38 phased panel to mirror CHM13v2.0 panel"
        bcftools view --threads 4 -Ou -r $region "$multiallelic_indel_filter" $phased_panel_vcf \
        | bcftools norm --threads 8 -Ou -m -any --fasta $ref_fasta - \
        | bcftools +fixploidy -Ou - -- -f 2 \
        | bcftools view -c 1:minor --threads 8 -Ou -V other - \
        | bcftools annotate -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                            -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with CHM13">' \
                            -Ou -x INFO/MAC --threads 2 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags -Ou --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
        | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_ped_biallelic \
        && bcftools index --threads 8 -f $rare_variants_phased_ped_biallelic
    fi
fi


echo "making phased panel subsets"
if [ ! -s $rare_variants_phased_ped.csi ]; then
    bcftools norm --threads 8 -Ou --fasta $ref_fasta -m +any $rare_variants_phased_ped_biallelic \
    | bcftools annotate -a $syntenic_site_location --mark-sites +SYNTENIC \
                        -c "$chrom_specific_syntenic_annotation_line_part1" -H "$chrom_specific_syntenic_annotation_line_part2" \
                        -Ou -x INFO/MAC --threads 2 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags -Ou --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools view --threads 8 -Ou -V other - \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_ped && \
    bcftools index --threads 8 -f $rare_variants_phased_ped
fi


if [ ! -s $phased_panel_vcf_3202.csi ]; then
    ### Extract 2504 unrelated samples from results. This will be the official phased panel.
    echo "phased_panel_vcf_3202"
    bcftools annotate --threads 8 -Ob -x ^INFO/MAF,^INFO/MAC,^INFO/AC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT $rare_variants_phased_ped > $phased_panel_vcf_3202 \
    && bcftools index --threads 8 -f $phased_panel_vcf_3202 &
fi

if [ ! -s $phased_panel_vcf_3202_biallelic.csi ]; then
    echo "phased_panel_vcf_3202_biallelic"
    bcftools annotate --threads 8 -Oz -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT $rare_variants_phased_ped_biallelic > ${phased_panel_vcf_3202_biallelic%%.bcf}.vcf.gz \
    && bcftools index --threads 8 -f -t ${phased_panel_vcf_3202_biallelic%%.bcf}.vcf.gz \
    && bcftools view  --threads 8 -Ob ${phased_panel_vcf_3202_biallelic%%.bcf}.vcf.gz > $phased_panel_vcf_3202_biallelic \
    && bcftools index --threads 8 $phased_panel_vcf_3202_biallelic &
fi

wait

if [ ! -s $phased_panel_vcf_2504_biallelic.csi ]; then
    ### Note: inputation evaluation requires MAF is in the INFO field, and it's not that much of a bother/size increase
    echo "phased_panel_vcf_2504_biallelic"
    bcftools view -Ou --threads 8 -r $region -S $unrelated_samples $phased_panel_vcf_3202_biallelic \
    | bcftools view -Ou --threads 2 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ou --threads 4 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools annotate -Oz --threads 4 -x ^INFO/MAF,^INFO/MAC,^INFO/AN,^INFO/SYNTENIC,^FORMAT/GT - > ${phased_panel_vcf_2504_biallelic%%.bcf}.vcf.gz \
    && bcftools index --threads 8 -f -t ${phased_panel_vcf_2504_biallelic%%.bcf}.vcf.gz \
    && bcftools view  --threads 8 -Ob ${phased_panel_vcf_2504_biallelic%%.bcf}.vcf.gz > $phased_panel_vcf_2504_biallelic \
    && bcftools index --threads 8 $phased_panel_vcf_2504_biallelic &
fi

# while the above is going on, proceed with phasing the no-parents panel
if [[ $genome != 'GRCh38' ]]
then
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
fi

wait

# once that is done, make the reports and subsets that stem from the panels generated while we were phasing a no-parents panel

### phased_panel_vcf_2504_biallelic is what is used for imputation evaluation - get stats
# if [ ! -s $phased_panel_vcf_2504_biallelic_variant_report ]; then
#     echo "phased_panel_vcf_2504_biallelic_variant_report"
#     bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAC\t%INFO/AN\t%INFO/MAF\t%INFO/SYNTENIC\n' $phased_panel_vcf_2504_biallelic > $phased_panel_vcf_2504_biallelic_variant_report &
# fi

if [ ! -s $phased_panel_vcf_2504.csi ]; then
    echo "phased_panel_vcf_2504"
    bcftools view -Ou --threads 8 -S $unrelated_samples $phased_panel_vcf_3202 \
    | bcftools view -Oz --threads 8 -c 1:minor - > $phased_panel_vcf_2504 \
    && bcftools index --threads 8 -f $phased_panel_vcf_2504 &
fi

if [[ $genome != 'GRCh38' ]]
then
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
fi


## To evaluate 1KGP T2T performance as a reference panel when phasing variants,
## phase samples present in pangenome using the phased 2504 panel as a reference
### Remove all pangenome samples and parents of pangenome samples from phased 2504 panel.
### This will be our 'reference panel'
echo "creating 'ground truth' reference panels for phasing evaluation"
if [ ! -s $phased_panel_no_pangenome_biallelic.csi ]; then
    echo "phased_panel_no_pangenome_biallelic"
    bcftools view -Ou --threads 8 -S ^$pangenome_and_parents --force-samples $phased_panel_vcf_2504_biallelic 2> /dev/null \
    | bcftools view -Ou --threads 8 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC  > $phased_panel_no_pangenome_biallelic \
    && bcftools index --threads 8 -f $phased_panel_no_pangenome_biallelic &
fi

wait

if [[ $genome != 'GRCh38' ]]
then
    if [[ ! -s $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.csi ]]; then
        bcftools view --threads 8 -S $pangenome_samples --force-samples -Ou $vcf_phased_no_parents_rare_biallelic 2> /dev/null \
        | bcftools view -Ou --threads 8 -c 1:minor - \
        | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
        | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC \
        > $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf \
        && bcftools index --threads 8 -f $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf &
    fi
fi

# Identify trio-private singletons
if [[ ! -s $chrom_working_dir/${chrom}_private_singletons.txt ]]; then
    mkdir -p $chrom_working_dir/tmp
    cat pedigrees/duos_and_trios.txt \
    | parallel -j $num_threads "bcftools view --force-samples -H -G -s {} -x -c 2 $phased_panel_vcf_3202_biallelic | cut -f 3 > $chrom_working_dir/tmp/{#}.txt" && \
    cat $chrom_working_dir/tmp/*.txt | sort | uniq > $chrom_working_dir/${chrom}_private_singletons.txt && rm -rf $chrom_working_dir/tmp &
fi


if [ ! -s $chr_specific_reference_pangenome_variation_trimmed_biallelic.csi ]; then
    if [[ $trim_assemblies_to_callset == 'true' ]]; then
        echo "chr_specific_reference_pangenome_variation_trimmed_biallelic"
        ### Remove pangenome reference sites that are not present in 3202 biallelic reference.
        bcftools isec -r $region --threads 8 -o $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                -Ob -n =2 -w 1 $chr_specific_reference_pangenome_variation_biallelic $phased_panel_vcf_3202_biallelic \
        && bcftools index -f --threads 8 $chr_specific_reference_pangenome_variation_trimmed_biallelic
    else
        chr_specific_reference_pangenome_variation_trimmed_biallelic=$chr_specific_reference_pangenome_variation_biallelic
    fi
fi

wait

echo "preparing references for pangenome sample phasing"
if [ ! -s $vcf_to_phase_pangenome_biallelic_HPRC.csi ]; then
    echo "vcf_to_phase_pangenome_biallelic"
    ### Create unphased variant call set of samples present in the pangenome
    ### also perform basic prephasing filtering to identify phasable variant set
    bcftools +setGT -Ou --threads 8 $chr_specific_reference_pangenome_variation_trimmed_biallelic -- -t a -n u \
    | bcftools view -c 1:minor \
        -e "ALT=='*' || F_MISSING>0.05 || INFO/MAC==0" \
        --threads 8 -Ob -W -o $vcf_to_phase_pangenome_biallelic_HPRC - \
    && bcftools view -c $minimum_common_MAC:minor -Ob -W \
        -o $vcf_to_phase_pangenome_biallelic_HPRC_common $vcf_to_phase_pangenome_biallelic_HPRC &
fi

echo "preparing references for pangenome sample phasing"
if [ ! -s $vcf_to_phase_pangenome_biallelic_1kgp.csi ]; then
    echo "vcf_to_phase_pangenome_biallelic_1kgp"
    ### Create unphased variant call set of samples present in the pangenome
    ### also perform basic prephasing filtering to identify phasable variant set
    bcftools view --threads 2 -Ou -S $pangenome_samples --force-samples $vcf_to_phase 2> /dev/null \
    | bcftools view -c 1:minor \
        -e "ALT=='*' || F_MISSING>0.05 || INFO/MAC==0" \
        --threads 8 -Ob -W -o $vcf_to_phase_pangenome_biallelic_1kgp - &
fi

# echo "preparing references for pangenome sample phasing"
# if [ ! -s $vcf_to_phase_pangenome_biallelic_1kgp_intersect.csi ]; then
#     bcftools isec -r $region --threads 8 -o $vcf_to_phase_pangenome_biallelic_1kgp_intersect \
#             -Ob -n =2 -w 1 $vcf_to_phase_pangenome_biallelic_1kgp $phased_panel_vcf_3202_biallelic \
#     && bcftools index -f --threads 8 $vcf_to_phase_pangenome_biallelic_1kgp_intersect
# fi

wait

if (( $minimum_common_MAC > 1 )); then
    bcftools view -c $minimum_common_MAC:minor -Ob -W \
        -o $vcf_to_phase_pangenome_biallelic_1kgp_common $vcf_to_phase_pangenome_biallelic_1kgp
else
    vcf_to_phase_pangenome_biallelic_1kgp_common=$vcf_to_phase_pangenome_biallelic_1kgp
fi

wait


## While doing this, gather AC_AN data to later generate reference MAF tables
if [[ ! -s $chrom_working_dir/${chrom}_3202_AC_AN.tsv ]]; then
    echo "3202"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $phased_panel_vcf_3202_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_3202_AC_AN.tsv &
fi
if [[ ! -s $chrom_working_dir/${chrom}_2504_AC_AN.tsv ]]; then
    echo "2504"    
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $phased_panel_vcf_2504_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f '%ID\t%INFO/MAC\t%INFO/AN\n' - > $chrom_working_dir/${chrom}_2504_AC_AN.tsv  &
fi
if [[ $genome != 'GRCh38' ]]
then
    if [[ ! -s $chrom_working_dir/${chrom}_2002_AC_AN.tsv ]]; then
        echo "2002"
        bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $vcf_phased_no_parents_rare_biallelic \
        | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
        | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_2002_AC_AN.tsv &
    fi
fi
if [[ ! -s $chrom_working_dir/${chrom}_2430_AC_AN.tsv ]]; then
    echo "2430"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $phased_panel_no_pangenome_biallelic \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_2430_AC_AN.tsv &
fi
if [[ ! -s $chrom_working_dir/${chrom}_44_AC_AN.tsv ]]; then
    echo "44"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $vcf_to_phase_pangenome_biallelic_HPRC \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_44_AC_AN.tsv &
fi
if [[ ! -s $chrom_working_dir/${chrom}_39_AC_AN.tsv ]]; then
    echo "39"
    bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF $vcf_to_phase_pangenome_biallelic_1kgp \
    | bcftools +fill-tags -Ou --threads 2 - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools query -f "%ID\t%INFO/MAC\t%INFO/AN\n" - > $chrom_working_dir/${chrom}_39_AC_AN.tsv &
fi


echo "phasing pangenome samples against 2504 panel"
# ### Repeat phasing as above, specifying a reference panel during the common variant phasing.

if [ ! -s $common_variants_phased_HPRC_pangenome_against_ref_biallelic.csi ]; then
    ./SHAPEIT5_phase_common_static_v1.1.1 \
        --input $vcf_to_phase_pangenome_biallelic_HPRC_common \
        --reference $phased_panel_no_pangenome_biallelic \
        --map $chrom_map \
        --output $chrom_working_dir/tmp_HPRC_pangenome.bcf \
        --thread $num_threads \
        --log $chrom_working_dir/$chrom.common.pangenome_with_ref_panel.log \
        --filter-maf $rare_variant_threshold \
        $haploid_arg \
        --region $region \
    && bcftools view --threads 8 -Ob $chrom_working_dir/tmp_HPRC_pangenome.bcf > $common_variants_phased_HPRC_pangenome_against_ref_biallelic \
    && bcftools index -f --threads 8 $common_variants_phased_HPRC_pangenome_against_ref_biallelic
fi


## Phase rare variants in chunks. We cannot specify a reference panel in this step.
i=1
# for chrom_region in $(cat $chrom_regions)
# do
    # if [ ! -s $chrom_working_dir/${i}_tmp_pangenome.rare.bcf.csi ]; then
    if [ ! -s $rare_variants_phased_HPRC_pangenome_against_ref_biallelic.csi ]; then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase_pangenome_biallelic_HPRC \
            --map $chrom_map \
            --scaffold $common_variants_phased_HPRC_pangenome_against_ref_biallelic \
            --scaffold-region $whole_chrom \
            --input-region $whole_chrom \
            --thread $num_threads \
            $haploid_arg \
            --log $chrom_working_dir/${chrom}.${i}.rare.pangenome_with_ref_panel.log \
            --output $chrom_working_dir/${i}_tmp_pangenome_HPRC.rare.bcf && \
            bcftools index --threads 8 -f $chrom_working_dir/${i}_tmp_pangenome_HPRC.rare.bcf && \
        bcftools view --threads 8 -Ob $chrom_working_dir/${i}_tmp_pangenome_HPRC.rare.bcf > $rare_variants_phased_HPRC_pangenome_against_ref_biallelic && \
        bcftools index --threads 8 -f $rare_variants_phased_HPRC_pangenome_against_ref_biallelic

        if [[ ! -s $rare_variants_phased_HPRC_pangenome_against_ref_biallelic.csi ]]; then
            ## This can occur if there are no variants under the rare variant threshold - then shapeit5 errors with a 'no variants to phase' error
            cp $common_variants_phased_HPRC_pangenome_against_ref_biallelic $rare_variants_phased_HPRC_pangenome_against_ref_biallelic
            cp $common_variants_phased_HPRC_pangenome_against_ref_biallelic.csi $rare_variants_phased_HPRC_pangenome_against_ref_biallelic.csi
        fi

    fi
    
    wait

    echo "phasing pangenome samples against 2504 panel"
    ### Repeat phasing as above, specifying a reference panel during the common variant phasing.
    if [ ! -s $common_variants_phased_1kgp_pangenome_against_ref_biallelic.csi ]; then
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $vcf_to_phase_pangenome_biallelic_1kgp_common \
            --reference $phased_panel_no_pangenome_biallelic \
            --map $chrom_map \
            --output $chrom_working_dir/tmp_pangenome_1kgp.bcf \
            --thread $num_threads \
            --log $chrom_working_dir/$chrom.common.pangenome_with_ref_panel.log \
            --filter-maf $rare_variant_threshold \
            $haploid_arg \
            --region $region \
        && bcftools view --threads 8 -Ob $chrom_working_dir/tmp_pangenome_1kgp.bcf > $common_variants_phased_1kgp_pangenome_against_ref_biallelic \
        && bcftools index -f --threads 8 $common_variants_phased_1kgp_pangenome_against_ref_biallelic
    fi


    ## Phase rare variants in chunks. We cannot specify a reference panel in this step.
    i=1
    # for chrom_region in $(cat $chrom_regions)
    # do
                # --pbwt-modulo $pbwt_modulo \
                # --effective-size $hmm_ne \
    if [ ! -s $rare_variants_phased_1kgp_pangenome_against_ref_biallelic.csi ]; then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase_pangenome_biallelic_1kgp \
            --map $chrom_map \
            --scaffold $common_variants_phased_1kgp_pangenome_against_ref_biallelic \
            --thread $num_threads \
            --log $chrom_working_dir/${chrom}.${i}.rare.pangenome_with_ref_panel.log \
            --scaffold-region $whole_chrom \
            --input-region $whole_chrom \
            $haploid_arg \
            --output $chrom_working_dir/${i}_tmp_pangenome_1kgp.rare.bcf && \
            bcftools index --threads 8 -f $chrom_working_dir/${i}_tmp_pangenome_1kgp.rare.bcf && \
        bcftools view --threads 8 -Ob $chrom_working_dir/${i}_tmp_pangenome_1kgp.rare.bcf > $rare_variants_phased_1kgp_pangenome_against_ref_biallelic && \
        bcftools index --threads 8 -f $rare_variants_phased_1kgp_pangenome_against_ref_biallelic 

        if [[ ! -s $rare_variants_phased_1kgp_pangenome_against_ref_biallelic.csi ]]; then
            ## This can occur if there are no variants under the rare variant threshold - then shapeit5 errors with a 'no variants to phase' error
            cp $common_variants_phased_1kgp_pangenome_against_ref_biallelic $rare_variants_phased_1kgp_pangenome_against_ref_biallelic
            cp $common_variants_phased_1kgp_pangenome_against_ref_biallelic.csi $rare_variants_phased_1kgp_pangenome_against_ref_biallelic.csi
        fi
    fi
    #     let i++
    # done


### Concat rare variant 
# if [ ! -s $rare_variants_phased_pangenome_against_ref_biallelic.csi ]; then
    # bcftools concat --threads 8 -Ob -l $chrom_working_dir/*_tmp_pangenome.rare.bcf > $rare_variants_phased_pangenome_against_ref_biallelic && \
    # bcftools index --threads 8 -f $rare_variants_phased_pangenome_against_ref_biallelic &
# fi
wait

##############################################################
# Drop men from reference sets when working with X chromosome


# Phased panel of unrelated samples generated from 3202 panel w/ pedigree:
if [[ $chrom == 'chrX' ]]
then
    echo "removing male samples from panels for evaluation"
    if [[ ! -s ${phased_panel_vcf_3202_biallelic%%.bcf}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S $female_samples $phased_panel_vcf_3202_biallelic > ${phased_panel_vcf_3202_biallelic%%.bcf}.females_only.vcf.gz 2> /dev/null && \
        bcftools index --threads 8 -f ${phased_panel_vcf_3202_biallelic%%.bcf}.females_only.vcf.gz && \
        phased_panel_vcf_3202_biallelic=${phased_panel_vcf_3202_biallelic%%.bcf}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_3202%%.bcf}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S $female_samples $phased_panel_vcf_3202 > ${phased_panel_vcf_3202%%.bcf}.females_only.vcf.gz 2> /dev/null && \
        bcftools index --threads 8 -f ${phased_panel_vcf_3202%%.bcf}.females_only.vcf.gz && \
        phased_panel_vcf_3202=${phased_panel_vcf_3202%%.bcf}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_2504_biallelic%%.bcf}.females_only.vcf.gz.tbi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S $female_samples $phased_panel_vcf_2504_biallelic > ${phased_panel_vcf_2504_biallelic%%.bcf}.females_only.vcf.gz 2> /dev/null && \
        bcftools index --threads 8 -f ${phased_panel_vcf_2504_biallelic%%.bcf}.females_only.vcf.gz && \
        phased_panel_vcf_2504_biallelic=${phased_panel_vcf_2504_biallelic%%.bcf}.females_only.vcf.gz &
    fi

    if [[ ! -s ${phased_panel_vcf_2504%%.bcf}.females_only.vcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Oz -S $female_samples $phased_panel_vcf_2504 > ${phased_panel_vcf_2504%%.bcf}.females_only.vcf.gz 2> /dev/null && \
        bcftools index --threads 8 -f ${phased_panel_vcf_2504%%.bcf}.females_only.vcf.gz && \
        phased_panel_vcf_2504=${phased_panel_vcf_2504%%.bcf}.females_only.vcf.gz &
    fi

    if [[ ! -s ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.bcf}.females_only.bcf.csi ]]; then
        bcftools view -Ob --force-samples --threads 8 -S $female_samples $chr_specific_reference_pangenome_variation_trimmed_biallelic > ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.bcf}.females_only.bcf 2> /dev/null && \
        bcftools index --threads 8 -f ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.bcf}.females_only.bcf && \
        chr_specific_reference_pangenome_variation_trimmed_biallelic=${chr_specific_reference_pangenome_variation_trimmed_biallelic%.bcf}.females_only.bcf &
    fi

    if [[ $genome != 'GRCh38' ]]; then
        if [[ ! -s ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf.csi ]]; then
            bcftools view --threads 8 --force-samples -Ob -S $female_samples $vcf_phased_no_parents_common_biallelic > ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf 2> /dev/null && \
            bcftools index --threads 8  -f ${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf && \
            vcf_phased_no_parents_common_biallelic=${vcf_phased_no_parents_common_biallelic%.*cf}.females_only.bcf &
        fi

        if [[ ! -s ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf.csi ]]; then
            bcftools view --threads 8 --force-samples -Ob -S $female_samples ${vcf_phased_no_parents_rare_biallelic%.*cf}.bcf > ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf 2> /dev/null && \
            bcftools index --threads 8 -f ${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf && \
            vcf_phased_no_parents_rare_biallelic=${vcf_phased_no_parents_rare_biallelic%.*cf}.females_only.bcf &
        fi
    fi
  
    wait
fi

wait

##############################################################
# Evaluate accuracy
recalc_phasing_stats='true'
if [[ "$recalc_phasing_stats" == 'true' ]]; then
    echo "evaluating phasing accuracy"

    Evaluate accuracy of 3202 panel via two methods:
    1) by looking at within-trio phasing consistency per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing
    if [[ $genome != 'GRCh38' ]]
    then
    #  1.1) Trio consistency when only probands + unrelated are phased (no parental genomes leak into the rest of the data)
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $vcf_phased_no_parents_rare_biallelic \
                                    -P $pedigree -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/rare_noparents_vs_trios_${chrom}.log \
                                    --output $stats_dir/rare_noparents_vs_trios_${chrom} 2> /dev/null &
    fi
    #  1.2) by looking at within-trio phasing consistency - Full panel, trio + statistically phased.
    #       Basically confirming that trio-phasing worked, best assessment of trio-phased sample accuracy
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -P $pedigree -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_trios_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_trios_${chrom} 2> /dev/null &

    #  1.3) by looking at within-trio phasing consistency - 2504 panel, trio + statistically phased.
    #       Should be the same as 1.2, but summary statistics will not include children. Parents are phased with a mix of trio-consistency and statistical phasing.
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $phased_panel_vcf_2504_biallelic \
                                    -P $pedigree -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/2504_panel_vs_trios_${chrom}.log \
                                    --output $stats_dir/2504_panel_vs_trios_${chrom} 2> /dev/null &

    #  2) by looking at phasing consistency with ground truth pangenome samples
    #   2.1) trio-phased 3202 panel
    if [[ ! -s $stats_dir/3202_panel_vs_HPRC_${chrom}.calibration.switch.txt.gz ]]; then
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HPRC_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HPRC_${chrom} 2> /dev/null &
    fi

    if [[ $genome != 'GRCh38' ]]
    then
    #  2.2) evaluate phasing performance of phasing without pedigree against ground truth pangenome samples
    #        Able to compare trio-measured SER and HPRC consistency based SER
        ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                        --estimation $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf \
                                        -R $whole_chrom --singleton \
                                        --log $chrom_working_dir/noparents_vs_HPRC_${chrom}.log \
                                        --output $stats_dir/noparents_vs_HPRC_${chrom} 2> /dev/null &
    fi
    #  3) Evaluate panel's usefulness as a reference panel
    #      3.1) Rephased pangenome samples compared to HPRC samples
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $rare_variants_phased_HPRC_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/rare_HPRC_pangenome_panelphased_vs_pangenome_${chrom}.log \
                                    --output $stats_dir/rare_HPRC_pangenome_panelphased_vs_pangenome_${chrom} 2> /dev/null &
    
    #      3.2) Rephased pangenome variation compared to trio samples
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $rare_variants_phased_HPRC_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton -P $pedigree \
                                    --log $chrom_working_dir/rare_HPRC_pangenome_panelphased_vs_trios_${chrom}.log \
                                    --output $stats_dir/rare_HPRC_pangenome_panelphased_vs_trios_${chrom} 2> /dev/null &

    #      3.3) Rephased 1kgp variation from pangenome samples compared to HPRC samples
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
                                    --estimation $rare_variants_phased_1kgp_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/rare_1kgp_pangenome_panelphased_vs_pangenome_${chrom}.log \
                                    --output $stats_dir/rare_1kgp_pangenome_panelphased_vs_pangenome_${chrom} 2> /dev/null &
    
    #      3.4) Rephased 1kgp variation from pangenome samples compared to trio samples
    ./SHAPEIT5_switch_static_JLL --validation $vcf_to_phase \
                                    --estimation $rare_variants_phased_1kgp_pangenome_against_ref_biallelic \
                                    -R $whole_chrom --singleton -P $pedigree \
                                    --log $chrom_working_dir/rare_1kgp_pangenome_panelphased_vs_trios_${chrom}.log \
                                    --output $stats_dir/rare_1kgp_pangenome_panelphased_vs_trios_${chrom} 2> /dev/null &
    # 4) Experiment with HGSVC samples
    #     4.1) trio-phased 3202 panel
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_HGSVC_variation_biallelic \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_${chrom} 2> /dev/null &

    #     4.2) Examine all HGSVC/HPRC samples
    ./SHAPEIT5_switch_static_JLL --validation $chr_specific_reference_HGSVC_HPRC_variation_biallelic \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HPRC_and_HGSVC_all_samples_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HPRC_and_HGSVC_all_samples_${chrom} 2> /dev/null &

    #     4.2) Examine HGSVC/HPRC samples that were fully trio phased
    bcftools view --threads 4 -c 1:minor -r $whole_chrom --force-samples -S $basedir/sample_subsets/HGSVC_HPRC_probands.in_1KGP.txt 2> /dev/null \
                     -Ob -W -o ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_HPRC_probands.bcf \
                    $chr_specific_reference_HGSVC_HPRC_variation_biallelic && \
    ./SHAPEIT5_switch_static_JLL --validation ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_HPRC_probands.bcf \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HPRC_and_HGSVC_trio_probands_only_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HPRC_and_HGSVC_trio_probands_only_${chrom} 2> /dev/null &
    
    #     4.3) Examine subset of HGSVC samples that were fully trio phased
    bcftools view --threads 4 -c 1:minor -r $whole_chrom --force-samples -S $basedir/sample_subsets/HGSVC_probands.in_1KGP.txt 2> /dev/null \
                     -Ob -W -o ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_probands.bcf \
                    $chr_specific_reference_HGSVC_HPRC_variation_biallelic && \
    ./SHAPEIT5_switch_static_JLL --validation ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_probands.bcf \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_probands_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_probands_${chrom} 2> /dev/null &

    #     4.3) only examine HGSVC samples that were partially trio phased (parents)
    bcftools view --threads 4 -c 1:minor -r $whole_chrom --force-samples -S $basedir/sample_subsets/HGSVC_parents.in_1KGP.txt 2> /dev/null \
                     -Ob -W -o ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_parents.bcf \
                    $chr_specific_reference_HGSVC_HPRC_variation_biallelic && \
    ./SHAPEIT5_switch_static_JLL --validation ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_parents.bcf \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_parents_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_parents_${chrom} 2> /dev/null &

    #     4.4) only examine HGSVC samples that were not trio phased
    bcftools view --threads 4 -c 1:minor -r $whole_chrom --force-samples -S $basedir/sample_subsets/HGSVC_not_part_of_trio.in_1KGP.txt 2> /dev/null \
                     -Ob -W -o ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_not_part_of_trio.bcf \
                    $chr_specific_reference_HGSVC_HPRC_variation_biallelic && \
    ./SHAPEIT5_switch_static_JLL --validation ${chr_specific_reference_HGSVC_HPRC_variation_biallelic%%.bcf}.HGSVC_not_part_of_trio.bcf \
                                    --estimation $phased_panel_vcf_3202_biallelic \
                                    -R $whole_chrom --singleton \
                                    --log $chrom_working_dir/3202_panel_vs_HGSVC_not_part_of_trio_${chrom}.log \
                                    --output $stats_dir/3202_panel_vs_HGSVC_not_part_of_trio_${chrom} 2> /dev/null &


fi
wait



################### Examine imputation ##################

if [[ $genome == 'GRCh38' ]]
then
    ln -sf $chr_specific_reference_pangenome_variation_trimmed_biallelic* $basedir/GRCh38_pangenome_variation/
elif [[ $genome == 'CHM13v2.0' ]]
then
    # Imputation and imputation metric gathering

    GRCh38_fasta=/dev/shm/Homo_sapiens_assembly38.fasta
    T2T_fasta=$basedir/chm13v2.0.fa.gz

    grch38_syntenic_site_location="$basedir/hg38.GCA_009914755.4.synNet.summary.bed.gz"
    t2t_syntenic_site_location="$basedir/chm13v2-syntenic_to_hg38.bed"

    lifted_panel_folder=$basedir/liftover/lifted_panels
    GRCh38_lifted_panel=$lifted_panel_folder/$(echo $(basename $phased_panel_vcf_2504_biallelic) | sed 's/CHM13v2.0/GRCh38.lifted_from_CHM13v2.0/' | sed 's/native_maps.//' )
    GRCh38_native_panel=$basedir/phased_GRCh38_panel/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
    T2T_lifted_panel=$lifted_panel_folder/$(echo $(basename $GRCh38_native_panel) | sed 's/\.GRCh38\./.CHM13v2.0.lifted_from_GRCh38./')
    T2T_native_panel=$phased_panel_vcf_2504_biallelic

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

    # Liftover panels `if necessary
    if [[ ! -s $T2T_native_panel.csi ]] && [[ -s $phased_panel_vcf_2504_biallelic.tbi ]]; then
            bcftools view -Ob --threads 4 $phased_panel_vcf_2504_biallelic > $T2T_native_panel && bcftools index -f --threads 4 $T2T_native_panel \
        && ./liftover_panel.sh $T2T_native_panel $GRCh38_lifted_panel /dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta /dev/shm/chm13v2-hg38.over.chain liftover/grch38-chm13v2.sort.vcf.gz &
    elif [[ ! -s $T2T_lifted_panel.csi ]]; then
        ./liftover_panel.sh $T2T_native_panel $GRCh38_lifted_panel /dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta /dev/shm/chm13v2-hg38.over.chain liftover/grch38-chm13v2.sort.vcf.gz &
    fi

    if [[ ! -s $GRCh38_native_panel.csi ]] && [[ -s ${GRCh38_native_panel%%.bcf}.vcf.gz.tbi ]]; then
            bcftools view -Ob --threads 4 ${GRCh38_native_panel%%.bcf}.vcf.gz > $GRCh38_native_panel && bcftools index -f --threads 4 $GRCh38_native_panel \
        && ./liftover_panel.sh $GRCh38_native_panel $T2T_lifted_panel /dev/shm/chm13v2.0_maskedY_rCRS.fasta /dev/shm/hg38-chm13v2.over.chain liftover/chm13v2-grch38.sort.vcf.gz &
    elif [[ ! -s $GRCh38_lifted_panel.csi ]]; then
        ./liftover_panel.sh $GRCh38_native_panel $T2T_lifted_panel /dev/shm/chm13v2.0_maskedY_rCRS.fasta /dev/shm/hg38-chm13v2.over.chain liftover/chm13v2-grch38.sort.vcf.gz &
    fi


    wait


    if [[ $chrom == "PAR1" ]]
    then    
        SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
        SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chrX.recalibrated.snp_indel.pass.vcf.gz
    elif [[ $chrom == "PAR2" ]]
    then
        SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
        SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chrX.recalibrated.snp_indel.pass.vcf.gz
    elif [[ $chrom == "debug" ]]
    then
        SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.chr20.recalibrated.snp_indel.pass.vcf.gz
        SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/chr20.recalibrated.snp_indel.pass.vcf.gz
    else
        SGDP_variants_T2T=$SGDP_ground_truth_dir_T2T/SGDP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.vcf.gz
        SGDP_variants_GRCh38=$SGDP_ground_truth_dir_GRCh38/${chrom}.recalibrated.snp_indel.pass.vcf.gz
    fi

    if [ ! -s $SGDP_ground_truth_T2T.csi ]; then
        echo $chrom "SGDP_ground_truth_T2T"
        bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m -any $SGDP_variants_T2T \
        | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - 2> /dev/null \
        | bcftools view -Ou  -c 1:minor -i "ALT!='*' || F_MISSING<0.05 || ABS(ILEN)<=50" - \
        | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
        > $SGDP_ground_truth_T2T \
        && bcftools index --threads 8 -f $SGDP_ground_truth_T2T &
    fi
    if [ ! -s $SGDP_ground_truth_GRCh38.csi ]; then
        echo $chrom "SGDP_ground_truth_GRCh38"
        bcftools norm --threads 2 -Ou -r $grch38_region -f $GRCh38_fasta -m -any $SGDP_variants_GRCh38 \
        | bcftools annotate -Ou --threads 8 -a $grch38_syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                        -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - 2> /dev/null \
        | bcftools view -Ou  -c 1:minor -i "ALT!='*' || F_MISSING<0.05 || ABS(ILEN)<=50" - \
        | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
        > $SGDP_ground_truth_GRCh38 \
        && bcftools index --threads 8 -f $SGDP_ground_truth_GRCh38 &
    fi


    wait

    if [[ ! -s $chrom_working_dir/GRCh38_imputation_workspace/GRCh38_space_filtered/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel false SGDP
    fi
    if [[ ! -s $chrom_working_dir/T2T_imputation_workspace/T2T_space_filtered/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then    
        ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel false SGDP
    fi
    if [[ ! -s $chrom_working_dir/GRCh38_snps_imputation_workspace/GRCh38_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel true SGDP
    fi
    if [[ ! -s $chrom_working_dir/T2T_snps_imputation_workspace/T2T_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}/lifted_panel.common_variants.SGDP.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel true SGDP
    fi

    if [[ ! -s $GRCh38_lifted_panel_no_pangenome.csi ]]; then
        echo $chrom "GRCh38_lifted_panel_no_pangenome"
        bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_lifted_panel > $GRCh38_lifted_panel_no_pangenome && \
        bcftools index --threads 2 $GRCh38_lifted_panel_no_pangenome &
    fi
    if [[ ! -s $T2T_lifted_panel_no_pangenome.csi ]]; then
        echo $chrom "T2T_lifted_panel_no_pangenome"
        bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $T2T_lifted_panel > $T2T_lifted_panel_no_pangenome && \
        bcftools index --threads 2 $T2T_lifted_panel_no_pangenome &
    fi
    if [[ ! -s $GRCh38_native_panel_no_pangenome.csi ]]; then
        echo $chrom "GRCh38_native_panel_no_pangenome"
        bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_native_panel > $GRCh38_native_panel_no_pangenome && \
        bcftools index --threads 2 $GRCh38_native_panel_no_pangenome
    fi

    wait
    if [[ ! -s $chrom_working_dir/GRCh38_imputation_pangenome_workspace/GRCh38_space_filtered/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_pangenome_workspace      $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome false pangenome &
    fi
    if [[ ! -s $chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace/GRCh38_space_filtered_snpsOnly/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome true pangenome &
    fi
    wait
    if [[ ! -s $chrom_working_dir/T2T_imputation_pangenome_workspace/T2T_space_filtered/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then    
        ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_pangenome_workspace            $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome false pangenome &
    fi
    if [[ ! -s $chrom_working_dir/T2T_snps_imputation_pangenome_workspace/T2T_space_filtered_snpsOnly/lifted_panel.common_variants.pangenome.${chrom}/lifted_panel.common_variants.pangenome.${chrom}.rsquare.grp.txt.gz ]]; then
        ./assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_pangenome_workspace       $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome true pangenome &
    fi
fi

wait
# close logfile
set +x
exec 19>&-

wait