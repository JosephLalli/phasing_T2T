#!/usr/bin/env bash
set -e 

# command to run on my machine:
# cd into phasing_T2T dir
# num_threads=11
initial_vcf_calls_folder='/mnt/ssd/lalli/phasing_T2T'
# chrom_maps_folder=/mnt/ssd/lalli/phasing_T2T/T2T_maps_raw
# initial_vcf_calls_folder=/folder/with/chrom/vcfs
# chrom_maps_folder=/folder/with/T2T/chrom/maps
# parallel -j 8 echo "./phase_T2T.sh chr{} 11 T2T_phasing_native_maps_large_Ne" ::: $(seq 1 22)

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
hmm_ne=$6
chrom_map=$7

## imputation options
limit_to_syntenic_regions=$4 #'false'
limit_to_snps=$5 #'false'


if [[ -z $chrom_map ]]
then
    # chrom_map=$basedir/scaled_t2t_map_testing_032124/${chrom}.nomask.scaled.native.gmap.gz
    chrom_map=$basedir/t2t_maps_no_filtered_regions/${chrom}_noMask.gmap.gz
    # chrom_map=$basedir/T2T_maps_raw/${chrom}.native.gmap.gz
    # chrom_map=$basedir/T2T_maps_raw/${chrom}.scaled.gmap.tsv
    # chrom_map=$basedir/t2t_lifted_chrom_maps/${chrom}.t2t.gmap.resorted.gmap.gz
    # chrom_map=$basedir/T2T_maps_raw/${chrom}.inter.formatted.gmap.tsv

fi

if [[ -z $hmm_ne ]]
then
    hmm_ne=135000
fi


if [[ suffix != '' ]]
then
    suffix="_$suffix"
fi

##NOTE: for imputation purposes only!!!
# final_panel_dir=$basedir/phased_T2T_panel${suffix}
final_panel_dir=$basedir/phased_T2T_panel_clean_unmasked
stats_dir=$basedir/phasing_stats${suffix}
imputation_results_dir=$basedir/imputation_results$suffix

ref_fasta=$basedir/chm13v2.0.fa.gz
pangenome_vcf=$basedir/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
population_ids=$basedir/sample_subsets/unrelated_superpopulations.csv
chrom_chunking_coords=$basedir/regions.txt

syntenic_site_location="$basedir/ucsc_browser_syntentic_GCA_009914755.4_T2T-CHM13v2.0.hg38.synNet.summary.sorted.flag_fmt.bed.gz"

## Identify if we are dealing with a PAR region
# GRCh38: "X:10001-2781479" "X:2781480-155701382" "X:155701383-156030895"
# T2T: "X:0-2394410" "X:2394410-153925833" "X:153925834-154259566"

PAR1_start='0'
PAR1_end='2394410'
PAR2_start='153925834'
PAR2_end='154259566'

use_beagle='false'

if [[ $chrom == "chrPAR1" ]]
then
    chrom='PAR1'
    region="chrX:$PAR1_start-$PAR1_end"
    whole_chrom=$region
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    phased_panel_vcf=$phased_panel_folder/1KGP.CHM13v2.0.PAR1.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
    end_chrom=$PAR1_end

elif [[ $chrom == "chrPAR2" || $chrom == "PAR2" ]]
then
    chrom='PAR2'
    region="chrX:$PAR2_start-$PAR2_end"
    whole_chrom=$region
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    phased_panel_vcf=$phased_panel_folder/1KGP.CHM13v2.0.PAR2.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
    end_chrom=$PAR2_end
elif [[ $chrom == 'chrdebug' ]]
then
    chrom='debug'
    region="chr20:10000000-12000000"
    whole_chrom=$region
    chrom_map=$basedir/t2t_maps_no_filtered_regions/chr20_noMask.gmap.gz
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chr20.recalibrated.snp_indel.pass.vcf.gz
    phased_panel_vcf=$phased_panel_folder/1KGP.CHM13v2.0.chr20.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
    end_chrom=12000000

elif [[ $chrom == 'chrX' ]]
then
    region="chrX:$PAR1_end-$PAR2_start"
    whole_chrom=$region
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
    phased_panel_vcf=$phased_panel_folder/1KGP.CHM13v2.0.chrX.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
    end_chrom=$PAR2_start
    hmm_ne=$(awk "BEGIN { printf \"%.0f\", ($hmm_ne*.75) }")
    use_beagle='true'

else
    region=$chrom
    end_chrom=$(grep "^$chrom	" ${ref_fasta}.fai | cut -f 2)
    whole_chrom=$chrom:1-$end_chrom
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.$chrom.recalibrated.snp_indel.pass.vcf.gz
    phased_panel_vcf=$phased_panel_folder/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
fi

chrom_working_dir=$basedir/${chrom}_working${suffix}

mkdir -p $chrom_working_dir
mkdir -p $final_panel_dir
mkdir -p $stats_dir


### Since at this point I am focused on imputation only,
### I am setting up a new folder for each run that contains
### the 'clean unmasked' intermediate files.
# for f in $basedir/${chrom}_working_clean_unmasked/*
# do
#     rm -f $chrom_working_dir/$(basename $f)
#     # if [[ -s $chrom_working_dir/$(basename $f) ]]
#     # then
#     #     rm $chrom_working_dir/$(basename $f)
#     # fi
#     # echo $f
#     # ls $f
#     ln -s $f $chrom_working_dir/
# done


echo $chrom_working_dir

logfile=$chrom_working_dir/phasing_T2T_$1_$3.log
cp $basedir/phasing_T2T_$1_$3.log $logfile
exec 19>$logfile
export BASH_XTRACEFD=19
set -x # writes commands to logfile
rm $basedir/phasing_T2T_$1_$3.log


# Define chrom-specific reference files
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
    impute5_haploid_arg="--haploid $male_sample_list"
else
    haploid_arg=""
fi

# Define names of chrom specific vcfs that we will be phasing
vcf_to_phase=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.biallelic.bcf.gz
vcf_to_phase_multiallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.bcf.gz
fully_annotated_input_variants=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.vcf.gz

chr_specific_reference_pangenome_variation=$chrom_working_dir/${chrom}_reference_pangenome.bcf.gz
chr_specific_reference_pangenome_variation_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.biallelic.bcf.gz
chr_specific_reference_pangenome_variation_trimmed=$chrom_working_dir/${chrom}_reference_pangenome.filtered_variants.bcf.gz
chr_specific_reference_pangenome_variation_trimmed_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.filtered_variants.biallelic.bcf.gz

# Define the names of vcf files that are sample subsets
vcf_to_phase_unrelated=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unrelated.bcf.gz
vcf_to_phase_no_pangenome=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.no_pangenome.bcf.gz
vcf_to_phase_pangenome=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.pangenome.bcf.gz
vcf_to_phase_pangenome_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.pangenome.biallelic.bcf.gz

# Define phased result file names
common_variants_phased_ped=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.common.bcf
common_variants_phased_ped_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.biallelic.common.bcf
rare_variants_phased_ped=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.rare.bcf.gz
rare_variants_phased_ped_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.biallelic.rare.bcf.gz
phased_panel_no_pangenome=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.no_pangenome.bcf.gz
phased_panel_no_pangenome_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.no_pangenome.biallelic.bcf.gz

common_variants_phased_unrelated=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.common.unrelated.bcf.gz
rare_variants_phased_unrelated=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.rare.unrelated.bcf.gz

common_variants_phased_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.common.pangenome_samples.biallelic.bcf.gz
rare_variants_phased_pangenome_against_ref=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.rare.pangenome_samples.bcf.gz
rare_variants_phased_pangenome_against_ref_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.lifted_maps.biallelic.rare.pangenome_samples.bcf.gz

vcf_phased_trios_only=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.trios.bcf.gz
vcf_to_phase_no_parents=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unphased.noparents.bcf.gz
vcf_phased_no_parents_common=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.bcf.gz
vcf_phased_no_parents_common_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.biallelic.bcf.gz
vcf_phased_no_parents_rare=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.bcf.gz
vcf_phased_no_parents_rare_biallelic=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.biallelic.bcf.gz
vcf_phased_pangenome_samples=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.pangenome.bcf.gz

# output files
phased_panel_vcf_3202=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.3202.vcf.gz
phased_panel_vcf_2504=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.2504.vcf.gz
phased_panel_vcf_3202_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.biallelic.3202.vcf.gz
phased_panel_vcf_2504_biallelic=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.biallelic.2504.vcf.gz

fully_annotated_input_variant_report=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv
phased_panel_vcf_2504_biallelic_variant_report=$stats_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.2504.stats.tsv

# pedigrees
pedigree=$basedir/pedigrees/1kgp.ped
mendelian_pedigree=$basedir/pedigrees/trios_only.ped
pedigree_no_pangenome=$basedir/pedigrees/1kgp.no_pangenome.ped

# lists of different categories of samples
children_and_parents=$basedir/sample_subsets/all_1kgp_trio_samples.txt
children=$basedir/sample_subsets/children_only.txt
no_parents=$basedir/sample_subsets/not_parents.txt
unrelated_samples=$basedir/sample_subsets/unrelated_samples.txt
pangenome_samples=$basedir/sample_subsets/pangenome_samples.txt
pangenome_and_parents=$basedir/sample_subsets/pangenome_samples_and_parents.txt
parents=$basedir/sample_subsets/pangenome_parents.txt
SGDP_related_to_1KGP=$basedir/sample_subsets/SGDP_samples_in_1KGP_and_associated_trios.txt
SGDP_in_1KGP=$basedir/sample_subsets/SGDP_samples_in_1KGP_numbered_format.txt

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


## VARIANT FILTERING, QC VARIANT SUBSETTING/FORMATTING

# Create ground truth variation vcf from pangenome
if [ ! -s $chr_specific_reference_pangenome_variation_biallelic.csi ]
then
    bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m- $pangenome_vcf \
    | bcftools view -Ou --threads 8 -s ^GRCh38 - \
    | bcftools annotate -Ou --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with CHM13">' - \
                    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools sort -m 40G -Ob > $chr_specific_reference_pangenome_variation_biallelic \
    && bcftools index $chr_specific_reference_pangenome_variation_biallelic
fi

echo "creating unphased, annotated, filtered variant panel"
# Split multiallelic sites, filter sites using criteria described above,
# convert data to bcf format, and index.
if [ ! -s $vcf_to_phase.csi ]
then
    if [ ! -s $fully_annotated_input_variants.tbi ]
    then
        echo "fully_annotated_input_variants"
        bcftools norm --threads 8 -Ou -r $region -f $ref_fasta -m- $input_vcf \
        | bcftools +mendelian2 -Ou - --ped $mendelian_pedigree -m a -m d \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,F_MISSING \
        | bcftools +fill-tags --threads 8 -Ou - -- -S $population_ids -t HWE \
        | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                            -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
                            --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - > $fully_annotated_input_variants \
        && bcftools index --threads 8 -t -f $fully_annotated_input_variants
    fi

    echo "vcf_to_phase"

    if [[ $chrom != 'chrX' ]]
    then
        bcftools view \
                -e "ALT=='*' || INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || MAC==0 || INFO/HWE<1e-10 || INFO/HWE_EUR<1e-10 || INFO/HWE_AFR<1e-10 || INFO/HWE_EAS<1e-10 || INFO/HWE_AMR<1e-10 || INFO/HWE_SAS<1e-10 || FILTER!='PASS'" \
                --threads 8 -Ou $fully_annotated_input_variants \
        | bcftools annotate --threads 8 -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS - > $vcf_to_phase \
        && bcftools index --threads 8 -f $vcf_to_phase
    else
        bcftools view \
                -e "ALT=='*' || INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0 || F_MISSING>0.05 || INFO/MERR>(INFO/AN*0.05) || MAC==0 || INFO/HWE<1e-10 || INFO/HWE_EUR<1e-10 || INFO/HWE_AFR<1e-10 || INFO/HWE_EAS<1e-10 || INFO/HWE_AMR<1e-10 || INFO/HWE_SAS<1e-10 || FILTER!='PASS'" \
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
    if [ ! -s $fully_annotated_input_variants.tbi ]
    then
        echo "fully_annotated_input_variants"
        bcftools norm --threads 8 -Ou -r $region -f $ref_fasta -m- $input_vcf \
        | bcftools +mendelian2 -Ou - --ped $mendelian_pedigree -m a -m d \
        | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC,MAF,F_MISSING \
        | bcftools +fill-tags --threads 8 -Ou - -- -S $population_ids -t HWE \
        | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
                    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - > $fully_annotated_input_variants \
        && bcftools index --threads 8 -t -f $fully_annotated_input_variants
    fi

    docker run --user $UID -v $chrom_working_dir:$chrom_working_dir quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0 \
        gatk VariantsToTable \
            -V $fully_annotated_input_variants \
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
            -O $fully_annotated_input_variant_report &
fi


# PHASING
## Phase full 3202 sample panel with pedigree (max accuracy, requires trios and therefore accuracy measurements not generalizable)
### Phase common variants
if [ ! -s $common_variants_phased_ped ]
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
    if [ ! -s $chrom_working_dir/$i.rare.bcf ]
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
    | bcftools norm --threads 8 -Ou --fasta $ref_fasta - \
    | bcftools annotate -Ou --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
                    -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
    | bcftools +fill-tags -Ou --threads 8 - -- -t AN,AC,MAF \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_ped_biallelic && \
    bcftools index --threads 8 -f $rare_variants_phased_ped_biallelic
fi

echo "making phased panel subsets"
if [ ! -s $rare_variants_phased_ped.csi ]; then
    bcftools norm --threads 8 -Ou -r $region --fasta $ref_fasta -m +any $rare_variants_phased_ped_biallelic \
    | bcftools annotate -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
                        -Ou --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags -Ou --threads 8 - -- -t AN,AC,MAF \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_ped && \
    bcftools index --threads 8 -f $rare_variants_phased_ped
fi


if [ ! -s $phased_panel_vcf_3202.tbi ]; then
    ### Extract 2504 unrelated samples from results. This will be the official phased panel.
    echo "phased_panel_vcf_3202"
    bcftools view --threads 8 -Ou -V other -r $region $rare_variants_phased_ped \
    | bcftools annotate --threads 8 -Oz -x ^FORMAT/GT $rare_variants_phased_ped > $phased_panel_vcf_3202 \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_3202 &
fi

if [ ! -s $phased_panel_vcf_3202_biallelic.tbi ]; then
    echo "phased_panel_vcf_3202_biallelic"
    bcftools annotate --threads 8 -Oz -r $region -x ^FORMAT/GT $rare_variants_phased_ped_biallelic > $phased_panel_vcf_3202_biallelic \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_3202_biallelic &
fi

wait

if [ ! -s $phased_panel_vcf_2504.csi ]; then
    echo "phased_panel_vcf_2504"
    bcftools view -Oz --threads 8 -r $region -V other -S $unrelated_samples $phased_panel_vcf_3202 > $phased_panel_vcf_2504 \
    && bcftools index --threads 8 -f $phased_panel_vcf_2504
fi


if [ ! -s $phased_panel_vcf_2504_biallelic.tbi ]; then
    ### Note: inputation evaluation requires MAF is in the INFO field, and it's not that much of a bother/size increase
    echo "phased_panel_vcf_2504_biallelic"
    bcftools view -Ou --threads 8 -r $region -S $unrelated_samples $phased_panel_vcf_3202_biallelic \
    | bcftools annotate -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                        -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
                        -Ou --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags -Oz --threads 4 - -- -t AN,AC,MAF > $phased_panel_vcf_2504_biallelic \
    && bcftools index --threads 8 -f -t $phased_panel_vcf_2504_biallelic

fi

### phased_panel_vcf_2504_biallelic is what is used for imputation evaluation - get stats
if [ ! -s $phased_panel_vcf_2504_biallelic_variant_report ]; then
    echo "phased_panel_vcf_2504_biallelic_variant_report"
    if [[ $chrom == 'chrX' ]]; then
        bcftools +fixploidy -Ou --threads 4 $phased_panel_vcf_2504_biallelic -- -p bcftools_ploidy.txt -s bcftools_ploidy_sexes.txt \
        | bcftools +fill-tags -Ou --threads 4 - -- -t AN,AC,MAF  \
        | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/MAF\t%INFO/SYNTENIC\n' - > $phased_panel_vcf_2504_biallelic_variant_report &
    else
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/MAF\t%INFO/SYNTENIC\n' $phased_panel_vcf_2504_biallelic > $phased_panel_vcf_2504_biallelic_variant_report &
    fi
fi

if [ ! -s $vcf_to_phase_no_parents ]; then
## Repeat, but with no trio parents (per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing)
### Remove parents from unphased vcf file
    echo "vcf_to_phase_no_parents"

    bcftools view -Ob --threads 8 -S $no_parents $vcf_to_phase > $vcf_to_phase_no_parents \
    && bcftools index --threads 8 -f $vcf_to_phase_no_parents #biallelic
fi

### Perform same phasing procedure on children-and-singleton-only vcf
### Phase common variants
if [ ! -s $vcf_phased_no_parents_common_biallelic ]; then
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


### Concat rare variant chunks and merge multiallelic sites
if [ ! -s $vcf_phased_no_parents_rare_biallelic ]; then
    bcftools concat --threads 8 -Ou -l $chrom_working_dir/*_tmp_noparents.rare.bcf \
    | bcftools norm --threads 8 -Ob --fasta $ref_fasta - > $vcf_phased_no_parents_rare_biallelic \
    && bcftools index --threads 8 -f $vcf_phased_no_parents_rare_biallelic
fi
if [ ! -s $vcf_phased_no_parents_rare ]; then
    bcftools norm --threads 8 -Ou --fasta $ref_fasta -m +any $vcf_phased_no_parents_rare_biallelic \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $vcf_phased_no_parents_rare && \
    bcftools index --threads 8 -f $vcf_phased_no_parents_rare &
fi

if [ ! -s $vcf_phased_no_parents_common ]; then
    bcftools norm --threads 8 -Ou --fasta $ref_fasta -m +any $vcf_phased_no_parents_common_biallelic \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $vcf_phased_no_parents_common && \
    bcftools index --threads 8 -f $vcf_phased_no_parents_common &
fi


## To evaluate 1KGP T2T performance as a reference panel when phasing variants,
## phase samples present in pangenome using the phased 2504 panel as a reference
### Remove all pangenome samples and parents of pangenome samples from phased 2504 panel.
### This will be our 'reference panel'
echo "creating 'ground truth' reference panels for phasing evaluation"
if [ ! -s $phased_panel_no_pangenome_biallelic.csi ]; then
    echo "phased_panel_no_pangenome_biallelic"
    bcftools view -Ob --threads 8 -r $region -S ^$pangenome_and_parents --force-samples $phased_panel_vcf_2504_biallelic > $phased_panel_no_pangenome_biallelic \
    && bcftools index --threads 8 -f $phased_panel_no_pangenome_biallelic &
fi

if [ ! -s $phased_panel_no_pangenome.csi ]; then
    echo "phased_panel_no_pangenome"
    bcftools view -Ob --threads 8 -r $region -S ^$pangenome_and_parents --force-samples $phased_panel_vcf_2504 > $phased_panel_no_pangenome \
    && bcftools index --threads 8 -f $phased_panel_no_pangenome &
fi

if [ ! -s $chr_specific_reference_pangenome_variation_trimmed_biallelic.csi ]; then
    echo "chr_specific_reference_pangenome_variation_trimmed_biallelic"
    ### Remove sites that are not present in 3202 biallelic reference.
    bcftools isec -r $region --threads 8 -o $chr_specific_reference_pangenome_variation_trimmed_biallelic -Ob -n =2 -w 1 $chr_specific_reference_pangenome_variation_biallelic $vcf_to_phase \
    && bcftools index -f --threads 8 $chr_specific_reference_pangenome_variation_trimmed_biallelic &
fi

if [ ! -s $chr_specific_reference_pangenome_variation.csi ]; then
    echo "chr_specific_reference_pangenome_variation"
    bcftools norm --threads 8 -Ou -r $region -f $ref_fasta $pangenome_vcf \
    | bcftools view -Ou --threads 8 -s ^GRCh38 - \
    | bcftools annotate -Ou --threads 8 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    | bcftools sort -m 40G -Ob > $chr_specific_reference_pangenome_variation \
    && bcftools index -f $chr_specific_reference_pangenome_variation &
fi

if [ ! -s $vcf_to_phase_multiallelic.csi ]; then
    echo "vcf_to_phase_multiallelic"
    bcftools norm --threads 8 -r $region -Ou --fasta $ref_fasta -m +any $vcf_to_phase \
    | bcftools +fill-tags -Oz - -- -t AN,AC > $vcf_to_phase_multiallelic && \
    bcftools index --threads 8 -f $vcf_to_phase_multiallelic
fi


wait 

echo "preparing references for pangenome sample phasing"
if [ ! -s $vcf_to_phase_pangenome_biallelic.csi ]; then
    echo "vcf_to_phase_pangenome_biallelic"
    ### Create unphased variant call set of samples present in the pangenome.
    bcftools view -r $region --threads 8 -Ou -S $pangenome_samples --force-samples $vcf_to_phase \
    | bcftools view --threads 2 -Ob -i 'F_MISSING<((N_SAMPLES-1)/N_SAMPLES)' - > $vcf_to_phase_pangenome_biallelic \
    && bcftools index -f --threads 8 $vcf_to_phase_pangenome_biallelic
fi

if [ ! -s $chr_specific_reference_pangenome_variation_trimmed.csi ]; then
    echo "chr_specific_reference_pangenome_variation_trimmed"
    ### Remove sites that are not present in 3202 reference.
    bcftools isec -r $region --threads 8 -o $chr_specific_reference_pangenome_variation_trimmed -Ob -n =2 -w 1 $chr_specific_reference_pangenome_variation $vcf_to_phase_multiallelic \
    && bcftools index -f --threads 8 $chr_specific_reference_pangenome_variation_trimmed
fi

if [ ! -s $vcf_to_phase_pangenome.csi ]; then
    echo "vcf_to_phase_pangenome"  
    # SHAPEIT5 requires at least two alleles to have been called
    ### Create unphased variant call set of samples present in the pangenome.
    bcftools view --threads 8 -r $region -c 1 -Ob -S $pangenome_samples -i 'F_MISSING<((N_SAMPLES-1)/N_SAMPLES)' --force-samples $vcf_to_phase_multiallelic > $vcf_to_phase_pangenome \
    && bcftools index -f --threads 8 $vcf_to_phase_pangenome
fi

wait

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


### Phase rare variants in chunks. We cannot specify a reference panel in this step.
i=1
for chrom_region in $(cat $chrom_regions)
do
    if [ ! -s $chrom_working_dir/${i}_tmp_pangenome.rare.bcf.csi ]; then
        ./SHAPEIT5_phase_rare_static_v1.1.1 \
            --input $vcf_to_phase_pangenome_biallelic \
            --map $chrom_map \
            --scaffold $common_variants_phased_pangenome_against_ref_biallelic \
            --thread $num_threads \
            --log $chrom_working_dir/${chrom}.${i}.rare.pangenome_with_ref_panel.log \
            --pbwt-modulo $pbwt_modulo \
            --input-region $chrom_region \
            --scaffold-region $chrom_region \
            --effective-size $hmm_ne \
            $haploid_arg \
            --output $chrom_working_dir/${i}_tmp_pangenome.rare.bcf && \
            bcftools index --threads 8 -f $chrom_working_dir/${i}_tmp_pangenome.rare.bcf
    fi 
    let i++
done

### Concat rare variant 
if [ ! -s $rare_variants_phased_pangenome_against_ref_biallelic.csi ]; then
    bcftools concat --threads 8 -Ou -l $chrom_working_dir/*_tmp_pangenome.rare.bcf \
    | bcftools norm --threads 8 -Ob --fasta $ref_fasta - > $rare_variants_phased_pangenome_against_ref_biallelic && \
    bcftools index --threads 8 -f $rare_variants_phased_pangenome_against_ref_biallelic
fi

if [ ! -s $rare_variants_phased_pangenome_against_ref.csi ]; then
    bcftools norm --threads 8 -Ou --fasta $ref_fasta -m +any $rare_variants_phased_pangenome_against_ref_biallelic \
    | bcftools +setGT -Ob --threads 8 - -- -t a -n p > $rare_variants_phased_pangenome_against_ref \
    && bcftools index --threads 8 -f $rare_variants_phased_pangenome_against_ref
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

    if [[ ! -s ${vcf_phased_no_parents_common_biallelic%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list $vcf_phased_no_parents_common_biallelic > ${vcf_phased_no_parents_common_biallelic%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8  -f ${vcf_phased_no_parents_common_biallelic%.*cf.gz}.females_only.bcf.gz && \
        vcf_phased_no_parents_common_biallelic=${vcf_phased_no_parents_common_biallelic%.*cf.gz}.females_only.bcf.gz &
    fi

    if [[ ! -s ${vcf_phased_no_parents_common%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list $vcf_phased_no_parents_common > ${vcf_phased_no_parents_common%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8  -f ${vcf_phased_no_parents_common%.*cf.gz}.females_only.bcf.gz && \
        vcf_phased_no_parents_common=${vcf_phased_no_parents_common%.*cf.gz}.females_only.bcf.gz &
    fi

    if [[ ! -s ${vcf_phased_no_parents_rare_biallelic%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list ${vcf_phased_no_parents_rare_biallelic%.*cf.gz}.bcf.gz > ${vcf_phased_no_parents_rare_biallelic%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8 -f ${vcf_phased_no_parents_rare_biallelic%.*cf.gz}.females_only.bcf.gz && \
        vcf_phased_no_parents_rare_biallelic=${vcf_phased_no_parents_rare_biallelic%.*cf.gz}.females_only.bcf.gz &
    fi

    if [[ ! -s ${vcf_phased_no_parents_rare%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view --threads 8 --force-samples -Ob -S ^$male_sample_list ${vcf_phased_no_parents_rare%.*cf.gz}.bcf.gz > ${vcf_phased_no_parents_rare%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8 -f ${vcf_phased_no_parents_rare%.*cf.gz}.females_only.bcf.gz && \
        vcf_phased_no_parents_rare=${vcf_phased_no_parents_rare%.*cf.gz}.females_only.bcf.gz &
    fi
    
    if [[ ! -s ${vcf_to_phase_multiallelic%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view -Ob --force-samples --threads 8 -S ^$male_sample_list ${vcf_to_phase_multiallelic%.*cf.gz}.bcf.gz > ${vcf_to_phase_multiallelic%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8 -f ${vcf_to_phase_multiallelic%.*cf.gz}.females_only.bcf.gz && \
        vcf_to_phase_multiallelic=${vcf_to_phase_multiallelic%.*cf.gz}.females_only.bcf.gz &
    fi

    if [[ ! -s ${chr_specific_reference_pangenome_variation_trimmed%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view -Ob --force-samples --threads 8 -S ^$male_sample_list ${chr_specific_reference_pangenome_variation_trimmed%.*cf.gz}.bcf.gz > ${chr_specific_reference_pangenome_variation_trimmed%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8 -f ${chr_specific_reference_pangenome_variation_trimmed%.*cf.gz}.females_only.bcf.gz && \
        chr_specific_reference_pangenome_variation_trimmed=${chr_specific_reference_pangenome_variation_trimmed%.*cf.gz}.females_only.bcf.gz &
    fi

    if [[ ! -s ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf.gz}.females_only.bcf.gz.csi ]]; then
        bcftools view -Ob --force-samples --threads 8 -S ^$male_sample_list ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf.gz}.bcf.gz > ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf.gz}.females_only.bcf.gz && \
        bcftools index --threads 8 -f ${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf.gz}.females_only.bcf.gz && \
        chr_specific_reference_pangenome_variation_trimmed_biallelic=${chr_specific_reference_pangenome_variation_trimmed_biallelic%.*cf.gz}.females_only.bcf.gz &
    fi

    wait
fi

##############################################################
# Evaluate accuracy
echo "evaluating phasing accuracy"
# recalc_phasing_stats="false"
# if [[ "$recalc_phasing_stats" ]]; then

# # if [ ! -s "$stats_dir/rare_pangenome_panelphased_vs_pangenome_${chrom}.type.switch.txt.gz" ]; then

#     # Evaluate accuracy of 3202 panel via two methods:
#     #  1) by looking at within-trio phasing consistency per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase \
#                                     --estimation $vcf_phased_no_parents_rare_biallelic \
#                                     -P $pedigree -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/rare_noparents_vs_trios_${chrom}.log \
#                                     --output $stats_dir/rare_noparents_vs_trios_${chrom} &

#     # ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase_multiallelic \
#     #                                 --estimation $vcf_phased_no_parents_rare \
#     #                                 -P $pedigree -R $whole_chrom --singleton \
#     #                                 --log $chrom_working_dir/multiallelic_rare_noparents_vs_trios_${chrom}.log \
#     #                                 --output $stats_dir/multiallelic_rare_noparents_vs_trios_${chrom} &
#     # if [[ ! -s $chrom_working_dir/common_variants_to_phase.biallelic.bcf.gz ]]; then
#     # Separately evaluate performance of panel generated without "rare" variant phasing
#     # bcftools isec --threads 8 -r $chrom -Ob -o $chrom_working_dir/common_variants_to_phase.biallelic.bcf.gz -n =2 -w 1 $vcf_to_phase $vcf_phased_no_parents_common_biallelic \
#     # && bcftools index --threads 8 -f $chrom_working_dir/common_variants_to_phase.biallelic.bcf.gz && \
#     # ./SHAPEIT5_switch_static_v1.1.1 --validation $chrom_working_dir/common_variants_to_phase.biallelic.bcf.gz \
#     #                                 --estimation $vcf_phased_no_parents_common_biallelic \
#     #                                 -P $pedigree -R $whole_chrom --singleton \
#     #                                 --log $chrom_working_dir/common_noparents_vs_trios_${chrom}.log \
#     #                                 --output $stats_dir/common_noparents_vs_trios_${chrom} &

#     # if [[ ! -s $chrom_working_dir/common_variants_to_phase.bcf.gz ]]; then
#     # Separately evaluate performance of panel generated without "rare" variant phasing
#     # # bcftools isec --threads 8 -r $chrom -Ob -o $chrom_working_dir/common_variants_to_phase.bcf.gz -n =2 -w 1 $vcf_to_phase_multiallelic $vcf_phased_no_parents_common \
#     # # && bcftools index --threads 8 -f $chrom_working_dir/common_variants_to_phase.bcf.gz && \
#     # ./SHAPEIT5_switch_static_v1.1.1 --validation $chrom_working_dir/common_variants_to_phase.bcf.gz \
#     #                                 --estimation $vcf_phased_no_parents_common \
#     #                                 -P $pedigree -R $whole_chrom --singleton \
#     #                                 --log $chrom_working_dir/multiallelic_common_noparents_vs_trios_${chrom}.log \
#     #                                 --output $stats_dir/multiallelic_common_noparents_vs_trios_${chrom} &


#     #  1.5) by looking at within-trio phasing consistency
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase \
#                                     --estimation $phased_panel_vcf_3202_biallelic \
#                                     -P $pedigree -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/3202_panel_vs_trios_${chrom}.log \
#                                     --output $stats_dir/3202_panel_vs_trios_${chrom} &

#     #  1.6) by looking at within-trio phasing consistency -- no children evaluated
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase \
#                                     --estimation $phased_panel_vcf_2504_biallelic \
#                                     -P $pedigree -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/2504_panel_vs_trios_${chrom}.log \
#                                     --output $stats_dir/2504_panel_vs_trios_${chrom} &


#     # #  1.7) by looking at within-trio phasing consistency - multiallelics
#     # ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase_multiallelic \
#     #                                 --estimation $phased_panel_vcf_3202 \
#     #                                 -P $pedigree -R $whole_chrom --singleton \
#     #                                 --log $chrom_working_dir/multiallelic_3202_panel_vs_trios_${chrom}.log \
#     #                                 --output $stats_dir/multiallelic_3202_panel_vs_trios_${chrom} &

#     # #  1.7) by looking at within-trio phasing consistency - multiallelics
#     # ./SHAPEIT5_switch_static_v1.1.1 --validation $vcf_to_phase_multiallelic \
#     #                                 --estimation $phased_panel_vcf_2504 \
#     #                                 -P $pedigree -R $whole_chrom --singleton \
#     #                                 --log $chrom_working_dir/multiallelic_2504_panel_vs_trios_${chrom}.log \
#     #                                 --output $stats_dir/multiallelic_2504_panel_vs_trios_${chrom} &


#     #  2) evaluate phasing performance of phasing with pedigree against ground truth pangenome samples
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
#                                     --estimation $phased_panel_vcf_3202_biallelic \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/3202_panel_vs_HPRC_${chrom}.log \
#                                     --output $stats_dir/3202_panel_vs_HPRC_${chrom} &

#     #  2) evaluate phasing performance of phasing with pedigree against ground truth pangenome samples
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed \
#                                     --estimation $phased_panel_vcf_3202 \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/multiallelic_3202_panel_vs_HPRC_${chrom}.log \
#                                     --output $stats_dir/multiallelic_3202_panel_vs_HPRC_${chrom} &


#     #  3) evaluate phasing performance of phasing without pedigree against ground truth pangenome samples
#     if [[ ! -s $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.gz ]]; then
#         bcftools view --threads 8 -S $pangenome_samples --force-samples -Ob $vcf_phased_no_parents_rare_biallelic > $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.gz && \
#         bcftools index --threads 8 -f $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.gz
#     fi
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
#                                     --estimation $chrom_working_dir/phased_pangenome_noparents.biallelic.bcf.gz \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/noparents_vs_HPRC_${chrom}.log \
#                                     --output $stats_dir/noparents_vs_HPRC_${chrom} &

#     # For some reason this file in particular is giving a fit. 
#     # It does not recalculate AN/AC values when extracting samples, and that causes problems downstream. 
#     # While this is likely a problem common to all sample subsetting, it only causes problems with this file.
#     # passing --no-update to view is the solution. - AN/AC tags not needed for this.
#     if [[ ! -s $chrom_working_dir/phased_pangenome_noparents.bcf.gz ]]; then
#         bcftools view --threads 8 -S $pangenome_samples --force-samples --no-update -Ou $vcf_phased_no_parents_rare >  $chrom_working_dir/phased_pangenome_noparents.bcf.gz && \
#         bcftools index --threads 8 -f $chrom_working_dir/phased_pangenome_noparents.bcf.gz
#     fi
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed \
#                                     --estimation $chrom_working_dir/phased_pangenome_noparents.bcf.gz \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/multiallelic_noparents_vs_HPRC_${chrom}.log \
#                                     --output $stats_dir/multiallelic_noparents_vs_HPRC_${chrom} &


#     ## 4) Evaluate performance of panel phased with all samples as reference panel
#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed_biallelic \
#                                     --estimation $rare_variants_phased_pangenome_against_ref_biallelic \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/rare_pangenome_panelphased_vs_pangenome_${chrom}.log \
#                                     --output $stats_dir/rare_pangenome_panelphased_vs_pangenome_${chrom} &

#     ./SHAPEIT5_switch_static_v1.1.1 --validation $chr_specific_reference_pangenome_variation_trimmed \
#                                     --estimation $rare_variants_phased_pangenome_against_ref \
#                                     -R $whole_chrom --singleton \
#                                     --log $chrom_working_dir/multiallelic_rare_pangenome_panelphased_vs_pangenome_${chrom}.log \
#                                     --output $stats_dir/multiallelic_rare_pangenome_panelphased_vs_pangenome_${chrom} &
# fi

# #########################################################################################
# #### Begin imputation evaluation using SDGP
# echo "evaluating imputation performance - SDGP"

# ##TODO: make sure all input files have syntenic sites annotated appropriately
# ##TODO: That is SGDP references, 2504 member reference panels, and pangenome references

# basedir=
# chrom=
# suffix=
# map=
## Generic imputation settings
redo_imputation=true
limit_to_syntenic_regions=$limit_to_syntenic_regions #($4)
limit_to_snps=$limit_to_snps #($5)

chrom_working_dir=$basedir/${chrom}_working${suffix}
imputation_results_dir=$basedir/imputation_results$suffix

mkdir -p $imputation_results_dir


omni_variants=$basedir/1000G_omni2.5.hg38.t2t-chm13-v2.0.vcf.gz
omni_variants_biallelic=$chrom_working_dir/1000G_omni2.5.hg38.t2t-chm13-v2.0.${chrom}.bcf

omni_variants_biallelic_filtered=$chrom_working_dir/1000G_omni2.5.hg38.t2t-chm13-v2.0.filtered.${chrom}.bcf
omni_variants_biallelic_filtered_in_1kgp=$chrom_working_dir/1000G_omni2.5.hg38.t2t-chm13-v2.0.chr_specific.in_ref_panel.bcf

phased_panel_vcf_2504_biallelic_filtered=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.biallelic.filtered.2504.bcf
phased_panel_vcf_2504_biallelic_filtered_nogts=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.lifted_maps.biallelic.filtered.2504.no_gts.bcf

## testset-specific settings
imputation_prefix_SGDP=$imputation_results_dir/${chrom}_CHM13v2.0_1000G_omni2.5_SGDP_variants.imputed

SGDP_ground_truth_dir=/mnt/ssd/lalli/nf_stage/genome_refs/T2T-CHM13_v2_ncbi110/SGDP
SGDP_variants_biallelic=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf
SGDP_variants_biallelic_filtered=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.filtered.pass.bcf
SGDP_variants_biallelic_filtered_in_1kgp=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.filtered.in_1kgp.pass.bcf
downsampled_SGDP_variants_biallelic_filtered_in_1kgp=$chrom_working_dir/${chrom}_CHM13v2.0_1000G_omni2.5_SGDP_phased_variants.filtered.downsampled.biallelic.bcf

if [[ $chrom == "PAR1" ]]
then    
    SGDP_variants=$SGDP_ground_truth_dir/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
elif [[ $chrom == "PAR2" ]]
then
    SGDP_variants=$SGDP_ground_truth_dir/SGDP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz
elif [[ $chrom == "debug" ]]
then
    SGDP_variants=$SGDP_ground_truth_dir/SGDP.CHM13v2.0.chr20.recalibrated.snp_indel.pass.vcf.gz
else
    SGDP_variants=$SGDP_ground_truth_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.vcf.gz
fi

## Filtering settings
## at minimum, remove all * alleles, all alleles with missing fraction above 0.05, anything without a PASS or missing filter value, and any indels over 50bp long.
variant_quality_filter_string="ALT!='*' && F_MISSING<0.05 && (FILTER=='PASS' || FILTER=='.') && ((TYPE!='snp' && (ABS(ILEN) < 50)) || TYPE='snp')"

if [[ $limit_to_syntenic_regions == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && INFO/SYNTENIC=1"
fi
if [[ $limit_to_snps == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && TYPE='snp'"
fi


# echo "phased panel variants"
# if [[ ! -s $phased_panel_vcf_2504_biallelic_filtered.csi ]] || [ $redo_imputation == 'true' ]
# then
    bcftools view -Ou -i "$variant_quality_filter_string" --threads 4 -r $region $phased_panel_vcf_2504_biallelic \
    | bcftools +fill-tags --threads 4 -Ou - -- -t "AN,AC,MAF,MAC:1=MAC" \
    | bcftools view --write-index -Ob -o $phased_panel_vcf_2504_biallelic_filtered -
# fi


echo "dropping gts to increase evaluation speed"
# if [[ ! -s $phased_panel_vcf_2504_biallelic_filtered_nogts.csi ]] || [ $redo_imputation == 'true' ]
# then
    bcftools view -Ob -G $phased_panel_vcf_2504_biallelic_filtered > $phased_panel_vcf_2504_biallelic_filtered_nogts
    bcftools index $phased_panel_vcf_2504_biallelic_filtered_nogts
# fi


echo "omni variants"
## Note: I'd like the option of filtering all datasets by the same parameters, including allele frequency/allele count.
## omni vcf files do not contain GTs, and do not have valid AN/AC/MAF/MAC values.
## When trying to filter on them, I eliminate all alleles.
## There are two options: Put dummy AN/AC values to allow for filtering,
## or don't filter the omni variants.
## If all other files are filtered, then intersecting those files with the omni variants should not cause otherwise filtered variants to be reintroduced into the dataset.
## I think the easier and more robust option is to not filter this dataset.
## Below I've commented out the original code that filters the omni datset, and then placed the new code that does not filter these variants.
# if [ ! -s $omni_variants_biallelic.csi ] || [ $redo_imputation == 'true' ]; then
#     # perform qc on omni variants
#     bcftools norm --threads 8 -Ou -f $ref_fasta -m- $omni_variants \
#     | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
#                     -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with CHM13">' - \
#                     --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
#     | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC > $omni_variants_biallelic \
#     && bcftools index $omni_variants_biallelic \
#     && bcftools view --threads 4 -Ou -r $region -i "$variant_quality_filter_string" $omni_variants_biallelic > $omni_variants_biallelic_filtered \
#     && bcftools index --threads 8 -f $omni_variants_biallelic_filtered \
#     && bcftools isec --threads 4 -o $omni_variants_biallelic_filtered_in_1kgp -Ob --write-index -n=2 -w 1 \
#             -r $region $omni_variants_biallelic_filtered $phased_panel_vcf_2504_biallelic_filtered_nogts
# else
#     echo "reusing omni_variants_biallelic"
#     echo "$omni_variants_biallelic"
# fi
# if [ ! -s $omni_variants_biallelic_filtered_in_1kgp.csi ] || [ $redo_imputation == 'true' ]; then
    # perform qc on omni variants
    bcftools norm --threads 8 -Ou -r $region -f $ref_fasta -m- $omni_variants \
    | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with CHM13">' - \
                    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC > $omni_variants_biallelic \
    && bcftools index $omni_variants_biallelic \
    && bcftools isec --threads 4 -o $omni_variants_biallelic_filtered_in_1kgp -Ob --write-index -n=2 -w 1 \
            -r $region $omni_variants_biallelic $phased_panel_vcf_2504_biallelic_filtered_nogts
# else
#     echo "reusing omni_variants_biallelic_filtered_in_1kgp"
#     echo "$omni_variants_biallelic_filtered_in_1kgp"
# fi

echo "SDGP_variants"
# if [ ! -s $SGDP_variants_biallelic_filtered.csi ] || [ $redo_imputation == 'true' ]; then
#     Exclude from evaluation sites that:
#         do not pass
#         have a high missing rate
#         are SVs
    bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m- $SGDP_variants \
    | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' \
                    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - \
    | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
    > $SGDP_variants_biallelic \
    && bcftools index $SGDP_variants_biallelic \
    && bcftools view --threads 6 -c 1:minor -Ob -i "INFO/VQSLOD>0 && $variant_quality_filter_string" $SGDP_variants_biallelic \
    > $SGDP_variants_biallelic_filtered \
    && bcftools index --threads 2 -f $SGDP_variants_biallelic_filtered
# fi

echo "isec to downsample"
# if [ ! -s $SGDP_variants_biallelic_filtered_in_1kgp.csi ] || [ $redo_imputation == 'true' ]; then
    bcftools isec --threads 4 -o $SGDP_variants_biallelic_filtered_in_1kgp -Ob --write-index -n=2 -w 1 \
            -r $region $SGDP_variants_biallelic_filtered $phased_panel_vcf_2504_biallelic_filtered_nogts
# fi

# if [ ! -s $downsampled_SGDP_variants_biallelic_filtered_in_1kgp.csi ] || [ $redo_imputation == 'true' ]; then
    bcftools isec --threads 4 -o $downsampled_SGDP_variants_biallelic_filtered_in_1kgp -Ob --write-index -n=2 -w 1 \
            -r $region $SGDP_variants_biallelic_filtered_in_1kgp $omni_variants_biallelic_filtered_in_1kgp
# fi


# FYI notes: 256 samples in SGDP that aren't in 1KGP
# 2504-41=2463 samples that aren't in SGDP or related to those individuals
# Could remove one or the other in theory; removing both ensures no chance of information bleed

wait
    
echo "prepping to impute $chrom"

# echo "making syntenic"
#     bcftools view -Ou -r $region -i "INFO/SYNTENIC=1"  $phased_panel_vcf_2504_biallelic \
#     | bcftools view -Ob -c 1:minor --trim-alt-alleles - \
#     > $phased_panel_vcf_2504_biallelic_syntenic && \
#     bcftools index --threads 2 -f $phased_panel_vcf_2504_biallelic_syntenic
# echo "final downsampling"
#     bcftools isec --threads 4 -Ou -n=3 -w 1 -r $region $SGDP_variants_biallelic $omni_variants_biallelic_filtered_in_1kgp $phased_panel_vcf_2504_biallelic_syntenic \
#     | bcftools +fill-tags --threads 8 -Ou - -- -t HWE,MAF \
#     | bcftools view -c 1:minor -Ob -i "F_MISSING<0.05 || INFO/HWE>=1e-4 || INFO/MAF>0.01 || ABS(ILEN)<=50" - \
#     > $downsampled_SGDP_variants_biallelic_in_1kgp \
#     && bcftools index -f --threads 4 $downsampled_SGDP_variants_biallelic

# fi

wait

# if [ ! -s $imputation_prefix_SGDP.bcf ] &&  [ $redo_imputation != 'true' ]
# then
#     if [[ $use_beagle == 'true' ]]; then
#         if [ ! -s ${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
#         then
#             bcftools view -Oz -o ${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz $phased_panel_vcf_2504_biallelic_filtered && \
#             phased_panel_vcf_2504_biallelic_filtered=${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz && \
#             bcftools index -f -t --threads 4 $phased_panel_vcf_2504_biallelic_filtered
#         fi

#         if [ ! -s ${downsampled_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
#         then
#             bcftools view -Oz -o ${downsampled_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz $downsampled_SGDP_variants_biallelic_filtered_in_1kgp && \
#             downsampled_SGDP_variants_biallelic_filtered_in_1kgp=${downsampled_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz && \
#             bcftools index -f -t --threads 4 $downsampled_SGDP_variants_biallelic_filtered_in_1kgp
#         fi

#         java -Xmx50g -jar $basedir/beagle.22Jul22.46e.jar \
#                     gt=$downsampled_SGDP_variants_biallelic_filtered_in_1kgp \
#                     ref=$phased_panel_vcf_2504_biallelic_filtered \
#                     out=$imputation_prefix_SGDP \
#                     map=${chrom_map%.gmap.gz}.beagle.gmap.gz \
#                     chrom=$region \
#                     gp=true \
#                     nthreads=$num_threads && \
#         bcftools index --threads 4 -f -t $imputation_prefix_SGDP.vcf.gz && \
#         bcftools view -Ob $imputation_prefix_SGDP.vcf.gz > $imputation_prefix_SGDP.bcf && \
#         bcftools index --threads 4 -f $imputation_prefix_SGDP.bcf
#         phased_panel_vcf_2504_biallelic_filtered=${phased_panel_vcf_2504_biallelic_filtered%%.vcf.gz}.bcf
#         downsampled_SGDP_variants_biallelic_filtered_in_1kgp=${downsampled_SGDP_variants_biallelic_filtered_in_1kgp%%.vcf.gz}.bcf
#     else
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $downsampled_SGDP_variants_biallelic_filtered_in_1kgp \
            --reference $phased_panel_vcf_2504_biallelic_filtered \
            --map $chrom_map \
            --output $chrom_working_dir/prephased_downsampled_SGDP_variants_biallelic_filtered_in_1kgp.bcf \
            --thread $num_threads \
            --pbwt-modulo 0.02 \
            --hmm-ne $hmm_ne \
            --log $chrom_working_dir/$chrom.SGDP_prephasing.log \
            $haploid_arg \
            --region $region \
        && impute5_v1.2.0/impute5_v1.2.0_static \
                    --g $chrom_working_dir/prephased_downsampled_SGDP_variants_biallelic_filtered_in_1kgp.bcf \
                    --h $phased_panel_vcf_2504_biallelic_filtered \
                    --m $chrom_map \
                    --r $whole_chrom \
                    --buffer-region $whole_chrom \
                    --o $imputation_prefix_SGDP.bcf \
                    --l $imputation_prefix_SGDP.log \
                    --out-ap-field \
                    --contigs-fai $ref_fasta.fai \
                    --threads $num_threads \
                    $impute5_haploid_arg
#     fi
# fi


#########################################################################################
#### Evaluate imputation with GLIMPSE2
echo "evaluating imputation"

# format:
# region frequency_source truth imputed
# if [ ! -s ${imputation_results_dir}/CHM13_glimpse2_concordance_input.SGDP.txt ] || [ $(grep "$region $phased_panel_vcf_2504_biallelic_filtered $SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_SGDP.bcf" ${imputation_results_dir}/CHM13_glimpse2_concordance_input.SGDP.txt | wc -l) == 0 ]
# then 
    echo -e "$region $phased_panel_vcf_2504_biallelic_filtered $SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_SGDP.bcf" >> ${imputation_results_dir}/CHM13_glimpse2_concordance_input.SGDP.txt
    echo -e "$region $phased_panel_vcf_2504_biallelic_filtered $SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_SGDP.bcf" >  ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.SGDP.txt
# fi
r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'

GLIMPSE2_concordance \
    --gt-val \
    --bins $r2_bins \
    --threads $num_threads \
    --af-tag MAF \
    --input ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.SGDP.txt \
    --log ${chrom_working_dir}/glimpse2_concordance_r2_bins.log \
    --out-r2-per-site \
    --out-rej-sites	\
    --out-conc-sites \
    --out-disc-sites \
    --output $stats_dir/$chrom.GLIMPSE2_concordance_r2_bins.SGDP




# #########################################################################################
# #### Begin imputation evaluation using SGDP
echo "evaluating imputation performance - SGDP variants lifted from GRCh38 callset"

## testset-specific settings
imputation_prefix_lifted_SGDP=$imputation_results_dir/${chrom}_CHM13v2.0_1000G_omni2.5_lifted_SGDP_variants.imputed

# lifted_SGDP_variants_biallelic=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.bcf
lifted_SGDP_variants_biallelic_filtered=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.filtered.bcf
lifted_SGDP_variants_biallelic_filtered_in_1kgp=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.filtered.in_panel.bcf
downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.filtered.in_panel.downsampled.bcf


# Lifted variants = subset from 'real' panel that are shared between lifted over and not 
if [[ $chrom == "PAR1" ]]
then    
    lifted_SGDP_variants=$basedir/liftover/SGDP.CHM13v2.0.chrX.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.bcf
elif [[ $chrom == "PAR2" ]]
then
    lifted_SGDP_variants=$basedir/liftover/SGDP.CHM13v2.0.chrX.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.bcf
elif [[ $chrom == "debug" ]]
then
    lifted_SGDP_variants=$basedir/liftover/SGDP.CHM13v2.0.chr20.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.bcf
else
    lifted_SGDP_variants=$basedir/liftover/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.genome_neutral_variants.bcf
fi

## Filtering settings
## at minimum, remove all * alleles, all alleles with missing fraction above 0.05, anything without a PASS or missing filter value, and any indels over 50bp long.
variant_quality_filter_string="ALT!='*' && F_MISSING<0.05 && (FILTER=='PASS' || FILTER=='.') && ((TYPE!='snp' && (ABS(ILEN) < 50)) || TYPE='snp')"

if [[ $limit_to_syntenic_regions == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && INFO/SYNTENIC=1"
fi
if [[ $limit_to_snps == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && TYPE='snp'"
fi
if [[ $filter_VQSLOD == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && INFO/VQSLOD>=0"
fi

echo "lifted_SGDP_variants"
if [ ! -s $lifted_SGDP_variants_biallelic_filtered.csi ] || [ $redo_imputation == 'true' ]; then
#     Exclude from evaluation sites that:
#         do not pass
#         have a high missing rate
#         are SVs
    bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m- $lifted_SGDP_variants \
    | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)">' \
                    --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
    | bcftools view --threads 6 -c 1:minor -Ob -i "$variant_quality_filter_string" $lifted_SGDP_variants_biallelic \
    > $lifted_SGDP_variants_biallelic_filtered \
    && bcftools index --threads 2 -f $lifted_SGDP_variants_biallelic_filtered
fi

echo "isec to downsample"
if [ ! -s $lifted_SGDP_variants_biallelic_filtered_in_1kgp.csi ] || [ $redo_imputation == 'true' ]; then
    bcftools isec --threads 4 -o $lifted_SGDP_variants_biallelic_filtered_in_1kgp -Ob -n=2 -w 1 \
            -r $region $lifted_SGDP_variants_biallelic_filtered $phased_panel_vcf_2504_biallelic_filtered_nogts \
    && bcftools index --threads 2 -f $lifted_SGDP_variants_biallelic_filtered_in_1kgp
fi

if [ ! -s $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp.csi ] || [ $redo_imputation == 'true' ]; then
    bcftools isec --threads 4 -o $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp -Ob -n=2 -w 1 \
            -r $region $lifted_SGDP_variants_biallelic_filtered_in_1kgp $lifted_omni_variants_biallelic_filtered_in_1kgp \
    && bcftools index --threads 2 -f $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp
fi


# FYI notes: 256 samples in lifted_SGDP that aren't in 1KGP
# 2504-41=2463 samples that aren't in lifted_SGDP or related to those individuals
# Could remove one or the other in theory; removing both ensures no chance of information bleed

wait
    
if [ ! -s $imputation_prefix_lifted_SGDP.bcf ] || [ $redo_imputation == 'true' ]
then
    if [[ $use_beagle == 'true' ]]; then
        if [ ! -s ${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
        then
            bcftools view -Oz -o ${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz $phased_panel_vcf_2504_biallelic_filtered && \
            phased_panel_vcf_2504_biallelic_filtered=${phased_panel_vcf_2504_biallelic_filtered%%.bcf}.vcf.gz && \
            bcftools index -f -t --threads 4 $phased_panel_vcf_2504_biallelic_filtered
        fi

        if [ ! -s ${downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
        then
            bcftools view -Oz -o ${downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp && \
            downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp=${downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp%%.bcf}.vcf.gz && \
            bcftools index -f -t --threads 4 $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp
        fi

        java -Xmx50g -jar $basedir/beagle.22Jul22.46e.jar \
                    gt=$downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp \
                    ref=$phased_panel_vcf_2504_biallelic_filtered \
                    out=$imputation_prefix_lifted_SGDP \
                    map=${chrom_map%.gmap.gz}.beagle.gmap.gz \
                    chrom=$region \
                    gp=true \
                    nthreads=$num_threads && \
        bcftools index --threads 4 -f -t $imputation_prefix_lifted_SGDP.vcf.gz && \
        bcftools view -Ob $imputation_prefix_lifted_SGDP.vcf.gz > $imputation_prefix_lifted_SGDP.bcf && \
        bcftools index --threads 4 -f $imputation_prefix_lifted_SGDP.bcf
        phased_panel_vcf_2504_biallelic_filtered=${phased_panel_vcf_2504_biallelic_filtered%%.vcf.gz}.bcf
        downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp=${downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp%%.vcf.gz}.bcf
    else
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp \
            --reference $phased_panel_vcf_2504_biallelic_filtered \
            --map $chrom_map \
            --output $chrom_working_dir/prephased_downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp.bcf \
            --thread $num_threads \
            --pbwt-modulo 0.02 \
            --hmm-ne $hmm_ne \
            --log $chrom_working_dir/$chrom.lifted_SGDP_prephasing.log \
            $haploid_arg \
            --region $region \
        && impute5_v1.2.0/impute5_v1.2.0_static \
                    --g $chrom_working_dir/prephased_downsampled_lifted_SGDP_variants_biallelic_filtered_in_1kgp.bcf \
                    --h $phased_panel_vcf_2504_biallelic_filtered \
                    --m $chrom_map \
                    --r $whole_chrom \
                    --buffer-region $whole_chrom \
                    --o $imputation_prefix_lifted_SGDP.bcf \
                    --l $imputation_prefix_lifted_SGDP.log \
                    --out-ap-field \
                    --contigs-fai $ref_fasta.fai \
                    --threads $num_threads \
                    $impute5_haploid_arg
    fi
fi

#########################################################################################
#### Evaluate imputation with GLIMPSE2
echo "evaluating imputation"

# format:
# region frequency_source truth imputed
if [ ! -s ${imputation_results_dir}/CHM13_glimpse2_concordance_input.lifted_SGDP.txt ] || [ $(grep "$region $phased_panel_vcf_2504_biallelic_filtered $lifted_SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_lifted_SGDP.bcf" ${imputation_results_dir}/CHM13_glimpse2_concordance_input.lifted_SGDP.txt | wc -l) == 0 ]
then 
    echo -e "$region $phased_panel_vcf_2504_biallelic_filtered $lifted_SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_lifted_SGDP.bcf" >> ${imputation_results_dir}/CHM13_glimpse2_concordance_input.lifted_SGDP.txt
    echo -e "$region $phased_panel_vcf_2504_biallelic_filtered $lifted_SGDP_variants_biallelic_filtered_in_1kgp $imputation_prefix_lifted_SGDP.bcf" >  ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.lifted_SGDP.txt
fi
r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'

GLIMPSE2_concordance \
    --gt-val \
    --bins $r2_bins \
    --threads $num_threads \
    --af-tag MAF \
    --input ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.lifted_SGDP.txt \
    --log ${chrom_working_dir}/glimpse2_concordance_r2_bins.lifted_SGDP.log \
    --out-r2-per-site \
    --out-rej-sites	\
    --out-conc-sites \
    --out-disc-sites \
    --output $stats_dir/$chrom.GLIMPSE2_concordance_r2_bins.lifted_SGDP


#########################################################################################
#### Begin imputation evaluation using Pangenome
echo "evaluating imputation performance - pangenome"

##TODO: make sure all input files have syntenic sites annotated appropriately
##TODO: That is SGDP references, 2504 member reference panels, and pangenome references

# basedir=
# chrom=
# suffix=
# map=
## Generic imputation settings


## testset-specific settings
phased_panel_vcf_2504_no_pangenome_biallelic_filtered=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.filtered.no_pangenome.biallelic.bcf
phased_panel_vcf_2504_no_pangenome_biallelic_filtered_nogts=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.filtered.no_pangenome.biallelic.gts_only.bcf

chr_specific_reference_pangenome_variation_filtered_biallelic=$chrom_working_dir/${chrom}_reference_pangenome.biallelic.filtered.bcf
chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp=$chrom_working_dir/${chrom}_reference_pangenome.biallelic.filtered.in_1kgp.bcf
downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp=$chrom_working_dir/${chrom}_CHM13v2.0_1000G_omni2.5_pangenome_phased_variants.biallelic.filtered.downsampled.bcf

imputation_prefix_pangenome=$imputation_results_dir/${chrom}_CHM13v2.0_1000G_omni2.5_pangenome_variants.imputed


echo "phased panel variants"

# if [ ! -s $phased_panel_vcf_2504_no_pangenome_biallelic_filtered.csi ] || [ $redo_imputation == 'true' ]
# then
    bcftools view -Ou -i "$variant_quality_filter_string" --threads 4 -r $region $phased_panel_no_pangenome_biallelic \
    | bcftools +fill-tags --threads 4 -Ou - -- -t "AN,AC,MAF,MAC:1=MAC" \
    | bcftools view -Ob -o $phased_panel_vcf_2504_no_pangenome_biallelic_filtered - \
    && bcftools index -f --threads 4 $phased_panel_vcf_2504_no_pangenome_biallelic_filtered
# fi

wait

# if [ ! -s $phased_panel_vcf_2504_no_pangenome_biallelic_filtered_nogts.csi ] || [ $redo_imputation == 'true' ]
# then
    echo "dropping gts to increase evaluation speed"
        bcftools view -Ob -G $phased_panel_vcf_2504_no_pangenome_biallelic_filtered > $phased_panel_vcf_2504_no_pangenome_biallelic_filtered_nogts
        bcftools index -f $phased_panel_vcf_2504_no_pangenome_biallelic_filtered_nogts
# fi
wait

echo "pangenome_variants"
#     Exclude from evaluation sites that:
#         do not pass
#         have a high missing rate
#         are SVs
    # bcftools norm --threads 2 -Ou -r $region -f $ref_fasta -m- $chr_specific_reference_pangenome_variation \
    # | bcftools annotate -Oz --threads 8 -a $syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
    #                 -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with GRCh38">' - \
    #                 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    # | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC \
    # > $chr_specific_reference_pangenome_variation_biallelic \
    # && bcftools index $chr_specific_reference_pangenome_variation_biallelic \

# num_missing: pangenome vcf represents haploid male X chr as './G'. There are 17 male samples. 
# FMISSING calcs num samples w/ at least one missing / total samples. Note, not M_allles/total_alleles.
# The base number of expected missing variants is 17/44  38.6% + 5% = min missing of 43.6%.
# if [ ! -s $chr_specific_reference_pangenome_variation_filtered_biallelic.csi ] || [ $redo_imputation == 'true' ]
# then
    if [[ $chrom == 'chrX' ]]
    then
        chrX_pangenome_variant_quality_filter_string=$(echo $variant_quality_filter_string | sed 's/F_MISSING<0.05/F_MISSING<0.436/')
        bcftools view --threads 6 -c 1:minor -Ob -i "$chrX_pangenome_variant_quality_filter_string" $chr_specific_reference_pangenome_variation_biallelic \
        > $chr_specific_reference_pangenome_variation_filtered_biallelic \
        && bcftools index --threads 2 -f $chr_specific_reference_pangenome_variation_filtered_biallelic
    else 
        bcftools view --threads 6 -c 1:minor -Ob -i "$variant_quality_filter_string" $chr_specific_reference_pangenome_variation_biallelic \
        > $chr_specific_reference_pangenome_variation_filtered_biallelic \
        && bcftools index --threads 2 -f $chr_specific_reference_pangenome_variation_filtered_biallelic
    fi
# fi


echo "isec to downsample"
# if [ ! -s $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp.csi ] || [ $redo_imputation == 'true' ]
# then
    bcftools isec --threads 4 -o $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp -Ob -n=2 -w 1 \
            -r $region $chr_specific_reference_pangenome_variation_filtered_biallelic $phased_panel_vcf_2504_no_pangenome_biallelic_filtered_nogts \
    && bcftools index --threads 4 -f $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp
# fi

# if [ ! -s $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp.csi ] || [ $redo_imputation == 'true' ]
# then
    bcftools isec --threads 4 -o $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp -Ob -n=2 -w 1 \
            -r $region $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp $omni_variants_biallelic_filtered_in_1kgp \
    && bcftools index --threads 4 -f $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp   

# fi


# FYI notes: 256 samples in pangenome that aren't in 1KGP
# 2504-41=2463 samples that aren't in pangenome or related to those individuals
# Could remove one or the other in theory; removing both ensures no chance of information bleed

wait
    
echo "prepping to impute $chrom"

# echo "making syntenic"
#     bcftools view -Ou -r $region -i "INFO/SYNTENIC=1"  $phased_panel_vcf_2504_biallelic \
#     | bcftools view -Ob -c 1:minor --trim-alt-alleles - \
#     > $phased_panel_vcf_2504_biallelic_syntenic && \
#     bcftools index --threads 2 -f $phased_panel_vcf_2504_biallelic_syntenic
# echo "final downsampling"
#     bcftools isec --threads 4 -Ou -n=3 -w 1 -r $region $pangenome_variants_biallelic $omni_variants_biallelic_filtered_in_1kgp $phased_panel_vcf_2504_biallelic_syntenic \
#     | bcftools +fill-tags --threads 8 -Ou - -- -t HWE,MAF \
#     | bcftools view -c 1:minor -Ob -i "F_MISSING<0.05 || INFO/HWE>=1e-4 || INFO/MAF>0.01 || ABS(ILEN)<=50" - \
#     > $downsampled_pangenome_variants_biallelic_in_1kgp \
#     && bcftools index -f --threads 4 $downsampled_pangenome_variants_biallelic

# fi

wait

# if [ ! -s $imputation_prefix_pangenome.bcf.csi ] || [ $redo_imputation == 'true' ]
# then
    # if [[ $use_beagle == 'true' ]]; then
    #     if [ ! -s ${phased_panel_vcf_2504_no_pangenome_biallelic_filtered%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
    #     then
    #         bcftools view -Oz -o ${phased_panel_vcf_2504_no_pangenome_biallelic_filtered%%.bcf}.vcf.gz $phased_panel_vcf_2504_no_pangenome_biallelic_filtered && \
    #         phased_panel_vcf_2504_no_pangenome_biallelic_filtered=${phased_panel_vcf_2504_no_pangenome_biallelic_filtered%%.bcf}.vcf.gz && \
    #         bcftools index -f -t --threads 4 $phased_panel_vcf_2504_no_pangenome_biallelic_filtered
    #     fi
    #     if [ ! -s ${downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp%%.bcf}.vcf.gz.tbi ] || [ $redo_imputation == 'true' ]
    #     then
    #         bcftools view -Oz -o ${downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp%%.bcf}.vcf.gz $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp && \
    #         downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp=${downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp%%.bcf}.vcf.gz && \
    #         bcftools index -f -t --threads 4 $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp
    #     fi

    #     java -Xmx50g -jar $basedir/beagle.22Jul22.46e.jar \
    #                 gt=$downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp \
    #                 ref=$phased_panel_vcf_2504_no_pangenome_biallelic_filtered \
    #                 out=$imputation_prefix_pangenome \
    #                 map=${chrom_map%.gmap.gz}.beagle.gmap.gz \
    #                 chrom=$region \
    #                 gp=true \
    #                 nthreads=$num_threads && \
    #     bcftools index --threads 4 -f -t $imputation_prefix_pangenome.vcf.gz && \
    #     bcftools view -Ob $imputation_prefix_pangenome.vcf.gz > $imputation_prefix_pangenome.bcf && \
    #     bcftools index --threads 4 -f $imputation_prefix_pangenome.bcf
    #     phased_panel_vcf_2504_no_pangenome_biallelic_filtered=${phased_panel_vcf_2504_no_pangenome_biallelic_filtered%%.vcf.gz}.bcf
    #     downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp=${downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp%%.vcf.gz}.bcf
    # else
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $downsampled_chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp \
            --reference $phased_panel_vcf_2504_no_pangenome_biallelic_filtered \
            --map $chrom_map \
            --output $chrom_working_dir/prephased_downsampled_pangenome_variants_biallelic_filtered_in_1kgp.bcf \
            --thread $num_threads \
            --pbwt-modulo 0.02 \
            --hmm-ne $hmm_ne \
            --log $chrom_working_dir/$chrom.pangenome_prephasing.log \
            $haploid_arg \
            --region $region \
        && impute5_v1.2.0/impute5_v1.2.0_static \
                    --g $chrom_working_dir/prephased_downsampled_pangenome_variants_biallelic_filtered_in_1kgp.bcf \
                    --h $phased_panel_vcf_2504_no_pangenome_biallelic_filtered \
                    --m $chrom_map \
                    --r $whole_chrom \
                    --buffer-region $whole_chrom \
                    --o $imputation_prefix_pangenome.bcf \
                    --l $imputation_prefix_pangenome.log \
                    --out-ap-field \
                    --contigs-fai $ref_fasta.fai \
                    --threads $num_threads \
                    $impute5_haploid_arg
#     fi
# fi


#########################################################################################
#### Evaluate imputation with GLIMPSE2
echo "evaluating imputation"

# format:
# region frequency_source truth imputed
# if [ ! -s ${imputation_results_dir}/CHM13_glimpse2_concordance_input.pangenome.txt ] || [ $(grep "$region $phased_panel_vcf_2504_no_pangenome_biallelic_filtered $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp $imputation_prefix_pangenome.bcf" ${imputation_results_dir}/CHM13_glimpse2_concordance_input.pangenome.txt | wc -l) == 0 ]
# then 
    echo -e "$region $phased_panel_vcf_2504_no_pangenome_biallelic_filtered $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp $imputation_prefix_pangenome.bcf" >> ${imputation_results_dir}/CHM13_glimpse2_concordance_input.pangenome.txt
    echo -e "$region $phased_panel_vcf_2504_no_pangenome_biallelic_filtered $chr_specific_reference_pangenome_variation_filtered_biallelic_in_1kgp $imputation_prefix_pangenome.bcf" >  ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.pangenome.txt
# fi

r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'

# if [ ! -s "$stats_dir/$chrom.GLIMPSE2_concordance_r2_bins.pangenome.rsquare.grp.txt.gz" ] || [ $redo_imputation == 'true' ]
# then
    GLIMPSE2_concordance \
        --gt-val \
        --bins $r2_bins \
        --threads $num_threads \
        --af-tag MAF \
        --input ${chrom_working_dir}/per_chrom_glimpse2_concordance_input.pangenome.txt \
        --log ${chrom_working_dir}/glimpse2_concordance_r2_bins.log \
        --out-r2-per-site \
        --out-rej-sites	\
        --out-conc-sites \
        --out-disc-sites \
        --output $stats_dir/$chrom.GLIMPSE2_concordance_r2_bins.pangenome
# fi


# close logfile
set +x
exec 19>&-
