#!/usr/bin/env bash

chrom=$1
num_threads=$2
initial_vcf_calls_folder=$3

basedir=$PWD/${chrom}_working

mkdir -p $basedir

# Define reference files
input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.vcf.gz
ref_fasta=$basedir/../chm13v2.0.fa.gz
end_chrom=$(cat $ref_fasta.fai | grep $chrom | cut -f 2)
chrom_map=/mnt/data/lalli/nf_stage/genome_refs/1000G_imputation_panels/t2t_lifted_chrom_maps/${chrom}.t2t.gmap.resorted.gmap.gz
pangenome_vcf=$basedir/../hprc-v1.0-mc-chm13.vcf.gz
population_ids=$basedir/../sample_subsets/unrelated_superpopulations.csv
regions_for_rare=$basedir/../sample_subsets/regions.txt
chrom_regions=$basedir/${chrom}_regions.txt

# extract chrom regions
grep $chrom: $regions_for_rare > $chrom_regions

vcf_to_phase=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.bcf.gz
chr_specific_reference_pangenome_variation=$basedir/${chrom}_reference_pangenome.bcf.gz

# Define the names of vcf files that are sample subsets
vcf_to_phase_unrelated=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unrelated.bcf.gz
vcf_to_phase_no_pangenome=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.no_pangenome.bcf.gz
vcf_to_phase_pangenome=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.pangenome.bcf.gz

# Define phased result file names
common_variants_phased_ped=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.common.bcf
rare_variants_phased_ped=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.rare.bcf.gz
phased_panel_no_pangenome=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.no_pangenome.bcf.gz

common_variants_phased_unrelated=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.common.unrelated.bcf.gz
rare_variants_phased_unrelated=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.rare.unrelated.bcf.gz

common_variants_phased_pangenome_against_ref=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.common.pangenome_samples.bcf.gz
rare_variants_phased_pangenome_against_ref=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.phased.rare.pangenome_samples.bcf.gz

vcf_phased_trios_only=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.trios.bcf.gz
vcf_to_phase_no_parents=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.unphased.noparents.bcf.gz
vcf_phased_no_parents_common=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.common.noparents.bcf.gz
vcf_phased_no_parents_rare=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.noparents.bcf.gz
vcf_phased_pangenome_samples=$basedir/1KGP.CHM13v2.0.${chrom}.snp_indel.phasing_qual_pass.rare.pangenome.bcf.gz

#output files
phased_panel_vcf=$basedir/../phased_T2T_panel/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.bcf.gz
stats_dir=$basedir/../phasing_stats
#logs
logfile=$basedir/$chrom.common.log
rare_logfile=$basedir/$chrom.rare.log
#pedigrees
pedigree=$basedir/../pedigrees/1kgp.ped
rare_pedigree=$basedir/../pedigrees/1kgp.no_segfault_rare.ped
mendelian_pedigree=$basedir/../pedigrees/trios_only.ped
pedigree_no_pangenome=$basedir/../pedigrees/1kgp.no_pangenome.ped
#sample lists
children_and_parents=$basedir/../sample_subsets/all_1kgp_trio_samples.txt
children=$basedir/../sample_subsets/children_only.txt
no_parents=$basedir/../sample_subsets/not_parents.txt
unrelated_samples=$basedir/../sample_subsets/unrelated_samples.txt
pangenome_samples=$basedir/../sample_subsets/pangenome_samples.txt
pangenome_and_parents=$basedir/../sample_subsets/pangenome_samples_and_parents.txt
parents=$basedir/../sample_subsets/pangenome_parents.txt
#settings
default_mcmc_iteration_scheme='5b,1p,1b,1p,1b,1p,5m'
shapeit4_suggested_unlimited_resources_mcmc_iteration_scheme='10b,1p,1b,1p,1b,1p,1b,1p,10m'
mcmc_iteration_scheme=$default_mcmc_iteration_scheme
pbwt_depth=4  # default: 4; high-accuracy: 8
pbwt_mac=5    # default: 5; shapeit4_default=2
pbwt_mdr=0.1  # default: 0.1; shapeit4_default=0.05
pbwt_modulo=0.1 # default 0.1; shapeit4_sequencing_default=0.0005
window=5      # default: 4; 1kgp paper using shapeit4: 5
rare_variant_threshold=0.001    # default: 0.001

input_vcf=$basedir/../../messy_phasing_T2T/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.vcf.gz


# # Split multiallelic sites, filter sites using criteria described above, 
# # convert data to bcf format, and index.
if [ ! -s $vcf_to_phase ]
then
    bcftools norm --threads 8 -Ou -f $ref_fasta -m- $input_vcf \
    | bcftools +mendelian -Ou - --ped $mendelian_pedigree -m a -m d \
    | bcftools +fill-tags --threads 8 -Ou - -- -t AN,AC \
    | bcftools +fill-tags --threads 8 -Ou - -- -S $population_ids -t HWE \
    | bcftools view \
            -e "ALT=='*' || INFO/NEGATIVE_TRAIN_SITE || INFO/VQSLOD<0 || F_MISSING>0.05 || INFO/MERR>30 || MAC==0 || INFO/HWE_EUR<1e-10 || INFO/HWE_AFR<1e-10 || INFO/HWE_EAS<1e-10 || INFO/HWE_AMR<1e-10 || INFO/HWE_SAS<1e-10 || FILTER!='PASS'" \
            --threads 8 -Ou - \
    | bcftools annotate -Ob -x ^INFO/AC,^INFO/AN,^FORMAT/GT,^FORMAT/PS > $vcf_to_phase && \
    bcftools index --threads 8 -f $vcf_to_phase
fi


# Phase full 3202 sample panel with pedigree
if [ ! -s $common_variants_phased_ped ]
then
./SHAPEIT5_phase_common_static_v1.0.0 \
    --input $vcf_to_phase \
    --map $chrom_map \
    --output $common_variants_phased_ped \
    --thread $num_threads \
    --log $logfile \
    --filter-maf $rare_variant_threshold \
    --mcmc-iterations $mcmc_iteration_scheme \
    --pbwt-modulo=$pbwt_modulo \
    --pbwt-depth=$pbwt_depth \
    --pbwt-mac=$pbwt_mac \
    --pbwt-mdr=$pbwt_mdr \
    --region $chrom \
    --pbwt-window $window --hmm-window $window \
    --pedigree $pedigree && \
bcftools index --threads 8 -f $common_variants_phased_ped
fi

i=1
for region in $(cat $chrom_regions)
do
    if [ ! -s $basedir/$i.rare.bcf ]
    then
        ./SHAPEIT5_phase_rare_static_v1.0.0 \
            --input-plain $vcf_to_phase \
            --map $chrom_map \
            --scaffold $common_variants_phased_ped \
            --thread $num_threads \
            --pedigree $pedigree \
            --log $basedir/$chrom.$i.rare.log \
            --input-maf $rare_variant_threshold \
            --pbwt-modulo=$pbwt_modulo \
            --input-region $region \
            --scaffold-region $region \
            --output $basedir/$i.rare.bcf
    fi
    bcftools index --threads 8 $basedir/$i.rare.bcf
    let i++
done
# Convert variant positions to original values, merge multiallelics, and index
bcftools concat -Ou -l $basedir/[^_].*rare.bcf \
| bcftools norm --threads 8 -Ob --fasta $ref_fasta -m +any - > $rare_variants_phased_ped \
&& bcftools index --threads 8 -f $rare_variants_phased_ped #&& \
# rm $basedir/*.rare.bcf

# Extract phased panel from most accurate phasing set (All samples w/ pedigree)
bcftools view -Ou --threads 8 -S $unrelated_samples $rare_variants_phased_ped \
| bcftools annotate --threads 2 -Ob -x ^FORMAT/GT > $phased_panel_vcf \
&& bcftools index --threads 8 -f $phased_panel_vcf

### Phase panel with parents removed per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing
# Step 1: create panel of unphased variants with no trio parents
bcftools view -Ob --threads 8 -S $no_parents $vcf_to_phase > $vcf_to_phase_no_parents \
&& bcftools index --threads 8 -f $vcf_to_phase_no_parents

# Step 2: phase panel without parents
./SHAPEIT5_phase_common_static_v1.0.0 \
    --input $vcf_to_phase_no_parents \
    --map $chrom_map \
    --output $basedir/noparents.bcf \
    --thread $num_threads \
    --log $basedir/$chrom.common.unrelated.log \
    --filter-maf $rare_variant_threshold \
    --mcmc-iterations $mcmc_iteration_scheme \
    --pbwt-modulo=$pbwt_modulo \
    --pbwt-depth=$pbwt_depth \
    --pbwt-mac=$pbwt_mac \
    --pbwt-mdr=$pbwt_mdr \
    --pbwt-window $window --hmm-window $window \
    --region $chrom && \
bcftools view --threads 8 -Ob $basedir/noparents.bcf > $vcf_phased_no_parents_common \
&& bcftools index --threads 8 $vcf_phased_no_parents_common #\
# && rm $basedir/noparents.bcf

i=1
for region in $(cat $chrom_regions)
do
    ./SHAPEIT5_phase_rare_static_v1.0.0 \
        --input-plain $vcf_to_phase_no_parents \
        --map $chrom_map \
        --scaffold $vcf_phased_no_parents_common \
        --thread $num_threads \
        --log $basedir/$chrom.rare.noparents.log \
        --input-maf $rare_variant_threshold \
        --pbwt-modulo $pbwt_modulo \
        --input-region $region \
        --scaffold-region $region \
        --output $basedir/${i}_tmp_noparents.rare.bcf && \
        bcftools index --threads 8 $basedir/${i}_tmp_noparents.rare.bcf
    let i++
done

bcftools concat --threads 8 -Ou -l $basedir/*_tmp_noparents.rare.bcf \
| bcftools norm --threads 8 -Ob --fasta $ref_fasta -m +any - > $vcf_phased_no_parents_rare \
&& bcftools index --threads 8 -f $vcf_phased_no_parents_rare #\
# && rm $basedir/*tmpu.rare.bcf

## Phase samples in pangenome using the phased 3202 panel as a reference
# Step 1: create subset of the phased 3202 sample panel that contains no pangenome samples
# or parents of pangenome samples
bcftools view -Ob --threads 8 -S ^$pangenome_and_parents --force-samples $phased_panel_vcf > $phased_panel_no_pangenome \
&& bcftools index --threads 8 -f $phased_panel_no_pangenome

# Step 2: create ground truth phased vcf of 1KGP samples in pangenome
#   Use bcftools isec to ensure we only look at ground truth sites that were also called in this dataset.
bcftools isec --threads 8 -r $chrom -Ob -o $chr_specific_reference_pangenome_variation -n =2 -w 1 $pangenome_vcf $vcf_to_phase \
&& bcftools index --threads 8 $chr_specific_reference_pangenome_variation

# Step 3: extract pangenomic samples from unphased variants
bcftools view --threads 8 -Ob -S $pangenome_samples --force-samples $vcf_to_phase > $vcf_to_phase_pangenome \
&& bcftools index --threads 8 -f $vcf_to_phase_pangenome

# Step 4: phase the pangenome samples with a panel of unrelated samples.
# Do not include a pedigree. All parents/relatives of the pangenome samples must be removed from the panel. (We did this in step 1.)
./SHAPEIT5_phase_common_static_v1.0.0 \
    --input $vcf_to_phase_pangenome \
    --reference $phased_panel_no_pangenome \
    --map $chrom_map \
    --output $basedir/tmp_pangenome.bcf \
    --thread $num_threads \
    --log $basedir/$chrom.common.noparents.log \
    --filter-maf $rare_variant_threshold \
    --mcmc-iterations $mcmc_iteration_scheme \
    --pbwt-modulo=$pbwt_modulo \
    --pbwt-depth=$pbwt_depth \
    --pbwt-mac=$pbwt_mac \
    --pbwt-mdr=$pbwt_mdr \
    --pbwt-window $window --hmm-window $window \
    --region $chrom
bcftools view --threads 8 -Ob $basedir/tmp_pangenome.bcf > $common_variants_phased_pangenome_against_ref \
&& bcftools index --threads 8 --threads 8 $common_variants_phased_pangenome_against_ref #\
# && rm $basedir/tmp_pangenome.bcf

i=1
for region in $(cat $chrom_regions)
do
    ./SHAPEIT5_phase_rare_static_v1.0.0 \
        --input-plain $vcf_to_phase_pangenome \
        --map $chrom_map \
        --scaffold $common_variants_phased_pangenome_against_ref \
        --thread $num_threads \
        --log $basedir/${chrom}.${i}.pangenome.rare.log \
        --input-maf $rare_variant_threshold \
        --pbwt-modulo=$pbwt_modulo \
        --input-region $region \
        --scaffold-region $region \
        --output $basedir/${i}_tmp_pangenome.rare.bcf && \
        bcftools index --threads 8 $basedir/${i}_tmp_pangenome.rare.bcf
    let i++
done

bcftools concat --threads 8 -Ou -l $basedir/*_tmp_pangenome.rare.bcf \
| bcftools norm --threads 8 -Ob --fasta $ref_fasta -m +any - > $rare_variants_phased_pangenome_against_ref \
&& bcftools index --threads 8 -f $rare_variants_phased_pangenome_against_ref #\
# && rm $basedir/*tmp_pangenome.rare.bcf

###############################################################
### Evaluate accuracy
# Phased panel of unrelated samples generated from 3202 panel w/ pedigree:
# Evaluate accuracy of 3202 panel via two methods:
#  1) by looking at within-trio phasing consistency per https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_wgs/#validation-of-your-phasing
./SHAPEIT5_switch_static_v1.0.0 --validation $vcf_to_phase \
                         --estimation $vcf_phased_no_parents_rare \
                         -P $pedigree -R $chrom:1-$end_chrom --singleton \
                         --log $basedir/rare_noparents_vs_trios_${chrom}.log \
                         --output $stats_dir/rare_noparents_vs_trios_${chrom} &

## Separately evaluate performance of panel generated without "rare" variant phasing
bcftools isec --threads 8 -r $chrom -Ob -o $basedir/common_variants_to_phase.bcf.gz -n =2 -w 1 $vcf_to_phase $vcf_phased_no_parents_common \
&& bcftools index --threads 8 $basedir/common_variants_to_phase.bcf.gz && \
./SHAPEIT5_switch_static_v1.0.0 --validation $basedir/common_variants_to_phase.bcf.gz \
                         --estimation $vcf_phased_no_parents_common \
                         -P $pedigree -R $chrom:1-$end_chrom \
                         --log $basedir/common_noparents_vs_trios_${chrom}.log \
                         --output $stats_dir/common_noparents_vs_trios_${chrom} &

#  1.5) by looking at within-trio phasing consistency
./SHAPEIT5_switch_static_v1.0.0 --validation $vcf_to_phase \
                         --estimation $rare_variants_phased_ped \
                         -P $pedigree -R $chrom:1-$end_chrom --singleton \
                         --log $basedir/panel_vs_trios_${chrom}.log \
                         --output $stats_dir/panel_vs_trios_${chrom} &


#  2) evaluate phasing performance of phasing with pedigree against ground truth pangenome samples
./SHAPEIT5_switch_static_v1.0.0 --validation $chr_specific_reference_pangenome_variation \
                         --estimation $rare_variants_phased_ped \
                         -R $chrom:1-$end_chrom --singleton \
                         --log $basedir/phased_vs_HPRC_${chrom}.log \
                         --output $stats_dir/phased_vs_HPRC_${chrom} &

#  3) evaluate phasing performance of phasing without pedigree against ground truth pangenome samples
bcftools view --threads 8 -S $pangenome_samples --force-samples -Ob $vcf_phased_no_parents_rare > $basedir/phased_pangenome_noparents.bcf.gz && \
bcftools index --threads 8 $basedir/phased_pangenome_noparents.bcf.gz && \
./SHAPEIT5_switch_static_v1.0.0 --validation $chr_specific_reference_pangenome_variation \
                         --estimation $basedir/phased_pangenome_noparents.bcf.gz \
                         -R $chrom:1-$end_chrom --singleton \
                         --log $basedir/noparents_vs_HPRC_${chrom}.log \
                         --output $stats_dir/noparents_vs_HPRC_${chrom} &


## Separately evaluate performance of panel generated without "rare" variant phasing
bcftools isec --threads 8 -r $chrom -Ob -o $basedir/chr_specific_reference_pangenome_common_variants.bcf.gz -n =2 -w 1 $chr_specific_reference_pangenome_variation $common_variants_phased_pangenome_against_ref \
&& bcftools index --threads 8 $basedir/chr_specific_reference_pangenome_common_variants.bcf.gz && \
./SHAPEIT5_switch_static_v1.0.0 --validation $basedir/chr_specific_reference_pangenome_common_variants.bcf.gz \
                         --estimation $common_variants_phased_pangenome_against_ref \
                         -R $chrom:1-$end_chrom \
                         --log $basedir/common_pangenome_panelphased_vs_pangenome_${chrom}.log \
                         --output $stats_dir/common_pangenome_panelphased_vs_pangenome_${chrom} &

## Evaluate performance of panel phased with all samples as reference panel 
./SHAPEIT5_switch_static_v1.0.0 --validation $chr_specific_reference_pangenome_variation \
                         --estimation $rare_variants_phased_pangenome_against_ref \
                         -R $chrom:1-$end_chrom --singleton \
                         --log $basedir/rare_pangenome_panelphased_vs_pangenome_${chrom}.log \
                         --output $stats_dir/rare_pangenome_panelphased_vs_pangenome_${chrom} &



# Note: Shapeit5 switch output headers
# Also note: all percentages are out of 100, not 1
# mendel_solver:
#     sample.mendel.txt.gz:       writePerSample    "mendel errors per sample"
#         sample_id father_id mother_id n_errors    n_present   percent_errors(100%)
#     variant.mendel.txt.gz:      writePerVariant   "mendel errors per variant"
#         variant_id  pos minor_allele_count    n_errors    num_called_variants  (n_errors * 100 / num_called_variants)
#     variant.imbalance.txt.gz:   writeImbalance    "Duo imbalance" = more like per-parent imbalence
#         variant_id    pos num_denovo  num_one_parental_copy   num_two_parental_copies num_denovo2 num_one_paternal_copy   num_two_paternal_copy   num_one_maternal_copy   num_two_maternal_copies
#     sample.pedigree:            writePedigree
# genotype_checker:
#     sample.typing.txt.gz:       writePerSample - genotype errors (reference vs experiment)
#         sample_id   num_errors  num_checked percentage_errors
#     variant.typing.txt.gz:      writePerVariant
#         variant_id  pos   num_errors  num_checked percentage_errors
# haplotype_checker:
#     sample.switch.txt.gz: "phasing switch errors per sample"
#         sample_id   n_switch_errors  num_checked switch_error_rate
#     sample.flipsAndSwitches.txt.gz -  "phasing flip and switch errors per sample" 
#         sample_id   num_errors  n_switches  n_flips n_correct   100*(n_switches)/(n_switches + n_flips + n_correct) 100*(n_flips)/(n_switches + n_flips + n_correct) 100*(n_correct)/(n_switches + n_flips + n_correct)
#     sample.variant.switch.txt - "phasing switch errors per variant"
#         variant_id  position  n_switch_errors    n_checked   switch_error_rate
#     sample.type.switch.txt - "phasing switch errors per variant type"
#         ref_alt_status    n_switch_errors   n_checked   switch_error_rate
#     sample.frequency.switch.txt - "phasing switch errors per frequency bin"
#         minor_allele_count   n_switch_errors    n_checked   switch_error_rate
#     sample.block.switch.txt - "correct phasing blocks per sample"
#         sample    block_coordinate(?)
#     sample.calibration.switch.txt - "phasing calibration"
#         Dunno honestly