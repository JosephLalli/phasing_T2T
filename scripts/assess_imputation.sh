#!/usr/bin/env bash
set -o xtrace

num_threads=12
# reference_genome='GRCh38' # working coordinate space. T2T or GRCh38.
# limit_to_snps=false
# missing_cutoff=0.05
# suffix=''
# recovered_included='true'
# filter_VQSLOD='false'

reference_genome=$1 # working coordinate space. T2T or GRCh38.
chrom=$2
working_dir=$3
reference_dataset=$4
native_panel=$5
lifted_panel=$6
limit_to_snps=$7
ref_dataset_name=$8
# ground_truth=$6
# limit_to_snps=true #$6
# filter_VQSLOD=false #$7


missing_cutoff=0.05

# T2T sourced from *_T2T_070224_HMMfixed_scaled
## All of these should have native variant IDs

# GRCh38_native_panel=/mnt/ssd/lalli/phasing_T2T/phased_GRCh38_panel/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
# GRCh38_lifted_panel=/mnt/ssd/lalli/phasing_T2T/liftover/lifted_panels/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
# T2T_native_panel=/mnt/ssd/lalli/phasing_T2T/phased_T2T_panel_T2T_scaled_newimpute_071524/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf
# T2T_lifted_panel=/mnt/ssd/lalli/phasing_T2T/liftover/lifted_panels/1KGP.CHM13v2.0.lifted_from_GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf

# SGDP_ground_truth_GRCh38=liftover_071224/SGDP.GRCh38.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf
# SGDP_ground_truth_T2T=liftover_071224/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf


subset_110=/mnt/ssd/lalli/phasing_T2T/sample_subsets/1kgp_paper_110_sample_SGDP_subset_numeric_format.txt

variant_quality_filter_string="ALT!='*' && (FILTER=='PASS' || FILTER=='.') && ((TYPE!='snp' && (ABS(ILEN) < 50)) || TYPE='snp')"
filtered_suffix='filtered'
if [[ $limit_to_syntenic_regions == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && INFO/SYNTENIC=1"
    filtered_suffix=$filtered_suffix'_syntenic'
fi
if [[ $limit_to_snps == 'true' ]]
then
    variant_quality_filter_string="$variant_quality_filter_string && TYPE='snp'"
    filtered_suffix=$filtered_suffix'_snpsOnly'
fi
ground_truth_filter_string=$variant_quality_filter_string
if [[ $filter_VQSLOD == 'true' ]]
then
    ground_truth_filter_string="$variant_quality_filter_string && INFO/VQSLOD>=0"
    filtered_suffix=$filtered_suffix'_vqslodThreshold'
fi
ground_truth_filtered_suffix=$filtered_suffix
if [[ $sample_set == '110' ]]
then
    subsample="-S $subset_110"
    ground_truth_filtered_suffix=$filtered_suffix'_110subsample'
fi


GRCh38_PAR1='chrX:10001-2781479'
GRCh38_chrX='chrX:2781480-155701382'
GRCh38_PAR2='chrX:155701383-156030895'
T2T_PAR1='chrX:0-2394410'
T2T_chrX='chrX:2394410-153925833'
T2T_PAR2='chrX:153925834-154259566'

T2T_fasta=/dev/shm/chm13v2.0_maskedY_rCRS.fasta
GRCh38_fasta=/dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta

if [[ $chrom == 'PAR1' ]]
then
    T2T_region=$T2T_PAR1
    GRCh38_region=$GRCh38_PAR1
elif [[ $chrom == 'PAR2' ]]
then
    T2T_region=$T2T_PAR2
    GRCh38_region=$GRCh38_PAR2
elif [[ $chrom == 'chrX' ]]
then
    T2T_region=$T2T_chrX
    GRCh38_region=$GRCh38_chrX
else
    T2T_region=$chrom:0-$(grep $chrom $T2T_fasta.fai | head -n 1 | cut -f 2)
    GRCh38_region=$chrom:0-$(grep $chrom $GRCh38_fasta.fai | head -n 1 | cut -f 2)
fi


if [[ $reference_genome == 'T2T' ]]; then
    whole_chrom=$T2T_region
    ref_fasta=$T2T_fasta
    omni=1000G_omni2.5.hg38.t2t-chm13-v2.0.biallelic.vcf.gz
    chrom_map=t2t_maps_no_filtered_regions/scaled_fixedamount/${chrom}_noMask.scaled.gmap.gz
    # SGDP_ground_truth_dir=/mnt/ssd/lalli/nf_stage/genome_refs/T2T-CHM13_v2_ncbi110/SGDP
    # ref_SGDP_groundtruth=$SGDP_ground_truth_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.vcf.gz
    # ref_pangenome_ground_truth=$working_dir/../${chrom}_reference_pangenome.filtered_variants.biallelic.bcf
    # reference_dataset=$SGDP_ground_truth_T2T
    # native_panel=$T2T_native_panel
    # lifted_panel=$T2T_lifted_panel


    if [[ ! -s $lifted_panel.csi ]]; then
        ./liftover_panel.sh $native_panel $lifted_panel $ref_fasta /dev/shm/hg38-chm13v2.over.chain /mnt/data/lalli/nf_stage/genome_refs/t2t_liftover_workspace/chm13v2-grch38.sort.vcf.gz
    fi
elif [[ $reference_genome == 'GRCh38' ]]; then
    whole_chrom=$GRCh38_region
    ref_fasta=$GRCh38_fasta
    omni=1000G_omni2.5.hg38.biallelic.vcf.gz
    chrom_map=GRCh38_performance_comparison/hg38_chrom_maps/${chrom}.b38.gmap.gz
    # SGDP_ground_truth_dir=$PWD/GRCh38_SGDP_full
    # ref_SGDP_groundtruth=$SGDP_ground_truth_dir/SGDP.GRCh38.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.vcf.gz
    # ref_pangenome_ground_truth=$basedir/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
    # reference_dataset=$SGDP_ground_truth_GRCh38
    # native_panel=$GRCh38_native_panel
    # lifted_panel=$GRCh38_lifted_panel

    
    if [[ ! -s $lifted_panel.csi ]]; then
        ./liftover_panel.sh $native_panel $lifted_panel $ref_fasta /dev/shm/chm13v2-hg38.over.chain /mnt/data/lalli/nf_stage/genome_refs/t2t_liftover_workspace/grch38-chm13v2.sort.vcf.gz
    fi
else
    print "Must be either T2T or GRCh38"
fi

echo $native_panel
echo $lifted_panel
echo $reference_dataset


working_dir=$working_dir/${reference_genome}_space_${filtered_suffix}

# 0) Filter panels
    mkdir -p $working_dir
    base_groundtruth_name=${reference_dataset%%.gz}
    base_groundtruth_name=${base_groundtruth_name%%.bcf}
    base_groundtruth_name=${base_groundtruth_name%%.vcf}
    if [ ! -s $base_groundtruth_name.$ground_truth_filtered_suffix.bcf.csi ]; then
        if [[ $ref_dataset_name != 'pangenome' ]]; then
            bcftools view --threads 2 -Ou -i "$ground_truth_filter_string" $subsample $reference_dataset \
            | bcftools norm --threads 2 -Ou -f $ref_fasta - \
            | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC,INFO/NGQ0miss:1=int(COUNT(FORMAT/GQ==0)+COUNT(FORMAT/GQ=="."))' \
            | bcftools view -Ou -c 1:minor --threads 4 -e "INFO/NGQ0miss>=0.05" - \
            | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
            > $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && bcftools index -f --threads 2 $base_groundtruth_name.$ground_truth_filtered_suffix.bcf &
        else
            bcftools view --threads 2 -Ou -i "$ground_truth_filter_string" $subsample $reference_dataset \
            | bcftools norm --threads 2 -Ou -f $ref_fasta - \
            | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC' \
            | bcftools view -Ou -c 1:minor --threads 4 -e "F_MISSING>=0.05" - \
            | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
            > $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && bcftools index -f --threads 2 $base_groundtruth_name.$ground_truth_filtered_suffix.bcf &
        fi
        if [[ $reference_genome == 'GRCh38' ]]; then
            cp $base_groundtruth_name.$ground_truth_filtered_suffix.bcf* /mnt/ssd/lalli/phasing_T2T/GRCh38_pangenome_variation/
        fi
    fi
    base_native_panel_name=${native_panel%%.gz}
    base_native_panel_name=${base_native_panel_name%%.bcf}
    base_native_panel_name=${base_native_panel_name%%.vcf}
    if [ ! -s  $base_native_panel_name.$filtered_suffix.bcf.csi ]; then
        echo $chrom native panel
        bcftools view --threads 2 -Ou -i "$variant_quality_filter_string" $native_panel \
        | bcftools norm --threads 2 -Ou -f $ref_fasta - \
        | bcftools annotate -Ou --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 4 -Ob - -- -t 'AN,AC,MAF,MAC:1=MAC' \
        | bcftools view --threads 2 -Ob -c 1:minor - \
        > $base_native_panel_name.$filtered_suffix.bcf \
        && bcftools index -f --threads 2 $base_native_panel_name.$filtered_suffix.bcf &
    fi
    base_lifted_panel_name=${lifted_panel%%.gz}
    base_lifted_panel_name=${base_lifted_panel_name%%.bcf}
    base_lifted_panel_name=${base_lifted_panel_name%%.vcf}
    if [ ! -s $base_lifted_panel_name.$filtered_suffix.bcf.csi ]; then
        echo $chrom lifted panel
        bcftools view --threads 2 -Ou -i "$variant_quality_filter_string" $lifted_panel \
        | bcftools norm --threads 2 -Ou -f $ref_fasta - \
        | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC' \
        | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools view --threads 2 -Ob -c 1:minor - \
        > $base_lifted_panel_name.$filtered_suffix.bcf \
        && bcftools index -f --threads 2 $base_lifted_panel_name.$filtered_suffix.bcf &
    fi

    wait

    reference_dataset=$base_groundtruth_name.$ground_truth_filtered_suffix.bcf
    native_panel=$base_native_panel_name.$filtered_suffix.bcf
    lifted_panel=$base_lifted_panel_name.$filtered_suffix.bcf

# 1) Downsample ground truth variants
    echo "downsampling $chrom"
    ground_truth_downsampled=${reference_dataset%%.bcf}.downsampled.bcf
    ground_truth_downsampled=$working_dir/$(basename $ground_truth_downsampled)

    if [[ ! -s $ground_truth_downsampled.csi ]]; then
        if [[ $ref_dataset_name != 'pangenome' ]]; then
            bcftools isec -Ou -n=4 -w 1 $reference_dataset $omni $lifted_panel $native_panel \
            | bcftools +fill-tags -Ou - -- -t F_MISSING,MAF,HWE \
            | bcftools view --threads 2 -c 1:minor -e "INFO/NGQ0miss>=0.05 || MAF<=0.01 || HWE<=1e-10" -Ob - > $ground_truth_downsampled
            bcftools index -f --threads 2 $ground_truth_downsampled
        else
            bcftools isec -Ou -n=4 -w 1 $reference_dataset $omni $lifted_panel $native_panel \
            | bcftools +fill-tags -Ou - -- -t F_MISSING,MAF,HWE \
            | bcftools view --threads 2 -c 1:minor -e "F_MISSING>=0.05 || MAF<=0.01 || HWE<=1e-10" -Ob - > $ground_truth_downsampled
            bcftools index -f --threads 2 $ground_truth_downsampled
        fi
    fi

# 2) Impute downsampled SGDP variants
    # 2a: Native panel
    if [[ ! -s $working_dir/$reference_genome.$chrom.native.imputed.bcf.csi ]]; then
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $ground_truth_downsampled \
            --reference $native_panel \
            --map $chrom_map \
            --output ${ground_truth_downsampled%%.bcf}.native.prephased.bcf \
            --thread $num_threads \
            --pbwt-modulo 0.02 \
            --hmm-ne 1000000 \
            --log $working_dir/${ref_dataset_name}_native_prephasing.log \
            $haploid_arg \
            --region $whole_chrom \
        && impute5_v1.2.0/impute5_v1.2.0_static \
            --g ${ground_truth_downsampled%%.bcf}.native.prephased.bcf \
            --h $native_panel \
            --m $chrom_map \
            --r $whole_chrom \
            --buffer-region $whole_chrom \
            --o $working_dir/$reference_genome.$chrom.native.imputed.bcf \
            --l $working_dir/$reference_genome.$chrom.native.imputed.log \
            --out-ap-field \
            --contigs-fai $ref_fasta.fai \
            --threads $num_threads \
            $impute5_haploid_arg &
    fi

    #2b) lifted
    if [[ ! -s $working_dir/$reference_genome.$chrom.lifted.imputed.bcf.csi ]]; then
        ./SHAPEIT5_phase_common_static_v1.1.1 \
            --input $ground_truth_downsampled \
            --reference $lifted_panel \
            --map $chrom_map \
            --output ${ground_truth_downsampled%%.bcf}.lifted.prephased.bcf \
            --thread $num_threads \
            --pbwt-modulo 0.02 \
            --log $working_dir/${ref_dataset_name}_lifted_prephasing.log \
            $haploid_arg \
            --hmm-ne 1000000 \
            --region $whole_chrom \
        && impute5_v1.2.0/impute5_v1.2.0_static \
            --g ${ground_truth_downsampled%%.bcf}.lifted.prephased.bcf \
            --h $lifted_panel \
            --m $chrom_map \
            --r $whole_chrom \
            --buffer-region $whole_chrom \
            --o $working_dir/$reference_genome.$chrom.lifted.imputed.bcf \
            --l $working_dir/$reference_genome.$chrom.lifted.imputed.log \
            --out-ap-field \
            --contigs-fai $ref_fasta.fai \
            --threads $num_threads \
            $impute5_haploid_arg &
    fi

wait

# 3) Identify variants in common between imputed datasets
python3.11 get_discordant_multiallelic_sites.py $native_panel $lifted_panel $working_dir/common.$chrom.IDs.txt

# 4) Subset imputed datasets to variants in common
bcftools view --threads 4 -Ob -i "ID==@$working_dir/common.$chrom.IDs.txt" $working_dir/$reference_genome.$chrom.native.imputed.bcf > $working_dir/$reference_genome.$chrom.native.imputed.common.bcf \
&& bcftools index --threads 4 -f $working_dir/$reference_genome.$chrom.native.imputed.common.bcf &
bcftools view --threads 4 -Ob -i "ID==@$working_dir/common.$chrom.IDs.txt" $working_dir/$reference_genome.$chrom.lifted.imputed.bcf > $working_dir/$reference_genome.$chrom.lifted.imputed.common.bcf \
&& bcftools index --threads 4 -f $working_dir/$reference_genome.$chrom.lifted.imputed.common.bcf &

wait

# 5) Measure imputation accuracy of all four sets (native, lifted, native-common, lifted-common).
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir/$reference_genome.$chrom.native.imputed.bcf" >  ${working_dir}/native_panel.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir/$reference_genome.$chrom.lifted.imputed.bcf" >  ${working_dir}/lifted_panel.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir/$reference_genome.$chrom.native.imputed.common.bcf" >  ${working_dir}/native_panel.common_variants.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir/$reference_genome.$chrom.lifted.imputed.common.bcf" >  ${working_dir}/lifted_panel.common_variants.${ref_dataset_name}.$chrom.txt

r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'

for infile in ${working_dir}/native_panel.${ref_dataset_name}.$chrom.txt ${working_dir}/lifted_panel.${ref_dataset_name}.$chrom.txt ${working_dir}/native_panel.common_variants.${ref_dataset_name}.$chrom.txt ${working_dir}/lifted_panel.common_variants.${ref_dataset_name}.$chrom.txt
do
    echo $infile
    mkdir -p ${infile%%.txt}
    GLIMPSE2_concordance \
        --gt-val \
        --bins $r2_bins \
        --threads $num_threads \
        --af-tag MAF \
        --input $infile \
        --log ${infile%%.txt}.log \
        --out-r2-per-site \
        --out-rej-sites	\
        --out-conc-sites \
        --out-disc-sites \
        --output ${infile%%.txt}/$(basename ${infile%%.txt}) &
done


# cat individual chrom input files them at the end to make glimpse2 input file