#!/usr/bin/env bash
set -o xtrace
set -e

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(dirname "$script_path")"
basedir=$script_dir/..

# Track background process failures
trap 'jobs -p | xargs -r kill 2>/dev/null; exit 1' ERR EXIT


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
no_singletons=$8
ref_dataset_name=$9
num_threads=${10}

force_imputation=false

if [[ -z $num_threads ]]; then
    num_threads=12
fi


should_run() {
    local target=$1
    if [[ $force_imputation == true ]]; then
        return 0
    fi
    # If file doesn't exist or is empty, needs to run
    if [[ ! -s "$target" ]]; then
        return 0
    fi
    # If target is a VCF/BCF (or index of one), check the underlying file
    if [[ "$target" == *.vcf* || "$target" == *.bcf* ]]; then
        target="${target%.csi}"
        target="${target%.tbi}"
        if [[ ! -s "$target" ]]; then
            echo "ERROR: Expected output $target was not created or is empty." >&2
            return 1
        fi
        if [[ "$target" == *.vcf || "$target" == *.vcf.gz || "$target" == *.bcf ]]; then
            local variant_count
            variant_count=$(bcftools index -n "$target" 2>/dev/null) || return 0
            if [[ -z "$variant_count" || "$variant_count" -eq 0 ]]; then
                return 0
            fi
        fi
    fi
    # File exists and is valid
    return 1
}


did_run() {
    local target=$1
    # If file doesn't exist or is empty, report failure
    if [[ ! -s "$target" ]]; then
        # echo "ERROR: Expected output $target was not created or is empty." >&2
        return 1
    fi
    # File exists and is valid
    return 0
}


# Global array to track expanded commands for background jobs
declare -A expanded_job_cmds

log_expanded_command() {
    local cmd="$*"
    eval "$cmd" &
    local pid=$!
    expanded_job_cmds[$pid]="$cmd"
    return 0
}

wait_and_check() {
    local failed=0
    local -a failed_pids=()

    for pid in $(jobs -p); do
        if ! wait "$pid"; then
            ((failed++))
            failed_pids+=("$pid")
        fi
    done

    if (( failed > 0 )); then
        echo "ERROR: $failed background job(s) failed:" >&2
        for pid in "${failed_pids[@]}"; do
            if [[ -n "${expanded_job_cmds[$pid]:-}" ]]; then
                echo "  [$pid] ${expanded_job_cmds[$pid]}" >&2
                unset expanded_job_cmds[$pid]
            else
                # Fallback if we don't have the expanded command
                local ps_cmd
                ps_cmd=$(ps -o command= -p "$pid" 2>/dev/null || true)
                echo "  [$pid] ${ps_cmd:-<command unavailable>}" >&2
            fi
        done
        return 1
    fi
    
    # Clean up successful job entries
    for pid in "${!expanded_job_cmds[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            unset expanded_job_cmds[$pid]
        fi
    done
    return 0
}

missing_cutoff=0.05

# working dir is $basedir/example_workspace/${chrom}_working${suffix}/GRCh38_imputation_workspace

subset_110=$basedir/resources/sample_subsets/1kgp_paper_110_sample_SGDP_subset_numeric_format.txt

working_dir=$(realpath -m "$working_dir")
reference_dataset=$(realpath -m "$reference_dataset")
native_panel=$(realpath -m "$native_panel")
lifted_panel=$(realpath -m "$lifted_panel")

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
if [[ $no_singletons == 'true' ]]
then
    minMAC=1 # change - we are handling this outside of the assess_imputation script now.
             # run_and_assess_imputation.sh will generate the no singleton file and provide it here. We shouldn't touch it.
             # We will preserve naming conventions for the sake of downstream scripts.
    filtered_suffix=$filtered_suffix'_no_singletons'
else
    minMAC=1
fi
T2T_specific_filter_string=$variant_quality_filter_string
if [[ $no_singletons == 'confident' ]]
then
    minMAC=1
    filtered_suffix=$filtered_suffix'_no_coinflips'
    T2T_specific_filter_string="$T2T_specific_filter_string && (COUNT(FORMAT/PP==0.5)==0)"
else
    minMAC=1
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

T2T_fasta=$basedir/resources/chm13v2.0.fa.gz
GRCh38_fasta=$basedir/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz

GRCh38_to_t2t_chain=$basedir/resources/grch38-chm13v2.chain
t2t_to_GRCh38_chain=$basedir/resources/chm13v2-grch38.chain

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
    source_fasta=$GRCh38_fasta
    omni=$basedir/resources/1000G_omni2.5.hg38.t2t-chm13-v2.0.biallelic.vcf.gz
    chrom_map=$basedir/resources/recombination_maps/t2t_native_scaled_maps/${chrom}.t2t.scaled.gmap.gz

    # if should_run "$lifted_panel.csi"; then
        # echo "panel $lifted_panel" never created
        # exit 1
        # "$script_dir/liftover_panel.sh" \
        #     -i $native_panel \
        #     -o $lifted_panel \
        #     -t $ref_fasta \
        #     -c $GRCh38_to_t2t_chain \
        #     -r $whole_chrom \
        #     -s $source_fasta
    # fi
elif [[ $reference_genome == 'GRCh38' ]]; then
    whole_chrom=$GRCh38_region
    ref_fasta=$GRCh38_fasta
    source_fasta=$T2T_fasta
    omni=$basedir/resources/1000G_omni2.5.hg38.biallelic.vcf.gz
    chrom_map=$basedir/resources/recombination_maps/grch38/${chrom}.b38.gmap.gz

    
    # if should_run "$lifted_panel.csi"; then
    #         echo "panel $lifted_panel" never created
    #     exit 1

        # "$script_dir/liftover_panel.sh" \
            # -i $native_panel \
            # -o $lifted_panel \
            # -t $ref_fasta \
            # -c $t2t_to_GRCh38_chain \
            # -r $whole_chrom \
            # -s $source_fasta
    # fi
else
    print "Must be either T2T or GRCh38"
fi


working_dir=$working_dir/${reference_genome}_space_${filtered_suffix}

# Namespace per reference dataset (e.g., SGDP vs pangenome) to avoid collisions when sharing the same workspace
working_dir_ref=$working_dir/${ref_dataset_name}

# 0) Filter panels
    mkdir -p $working_dir
    mkdir -p $working_dir_ref
    base_groundtruth_name=${reference_dataset%%.gz}
    base_groundtruth_name=${base_groundtruth_name%%.bcf}
    base_groundtruth_name=${base_groundtruth_name%%.vcf}
    if should_run "$base_groundtruth_name.$ground_truth_filtered_suffix.bcf.csi"; then
        if [[ $ref_dataset_name != 'pangenome' ]]; then
            bcftools view --threads 2 -Ou -i "$ground_truth_filter_string" $subsample $reference_dataset \
            | bcftools norm --threads 2 -Ou -f $ref_fasta - \
            | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC,INFO/NGQ0miss:1=int(COUNT(FORMAT/GQ==0)+COUNT(FORMAT/GQ=="."))' \
            | bcftools view -Ou -c $minMAC:minor --threads 4 -e "INFO/NGQ0miss>=0.05" - \
            | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
            > $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && bcftools index -f --threads 2 $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && did_run "$base_groundtruth_name.$ground_truth_filtered_suffix.bcf.csi"
        else
            bcftools view --threads 2 -Ou -i "$ground_truth_filter_string" $subsample $reference_dataset \
            | bcftools norm --threads 2 -Ou -f $ref_fasta - \
            | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC' \
            | bcftools view -Ou -c $minMAC:minor --threads 4 -e "F_MISSING>=0.05" - \
            | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
            > $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && bcftools index -f --threads 2 $base_groundtruth_name.$ground_truth_filtered_suffix.bcf \
            && did_run "$base_groundtruth_name.$ground_truth_filtered_suffix.bcf.csi"
        fi
    fi
    
    if bcftools view -h $native_panel | grep -qE '^##FORMAT=<ID=PP\b'; then
        native_filter=$T2T_specific_filter_string
    else
        native_filter=$variant_quality_filter_string
    fi
    base_native_panel_name=${native_panel%%.gz}
    base_native_panel_name=${base_native_panel_name%%.bcf}
    base_native_panel_name=${base_native_panel_name%%.vcf}
    if should_run "$base_native_panel_name.$filtered_suffix.bcf.csi"; then
        bcftools view --threads 2 -Ou -i "$native_filter" $native_panel \
        | bcftools norm --threads 2 -Ou -f $ref_fasta - \
        | bcftools annotate -Ou --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools +fill-tags --threads 4 -Ob - -- -t 'AN,AC,MAF,MAC:1=MAC' \
        | bcftools view --threads 2 -Ob -c $minMAC:minor - \
        > $base_native_panel_name.$filtered_suffix.bcf \
        && bcftools index -f --threads 2 $base_native_panel_name.$filtered_suffix.bcf \
        && did_run "$base_native_panel_name.$filtered_suffix.bcf.csi"
    fi
    if bcftools view -h $lifted_panel | grep -qE '^##FORMAT=<ID=PP\b'; then
        lifted_filter=$T2T_specific_filter_string
    else
        lifted_filter=$variant_quality_filter_string
    fi 
    
    base_lifted_panel_name=${lifted_panel%%.gz}
    base_lifted_panel_name=${base_lifted_panel_name%%.bcf}
    base_lifted_panel_name=${base_lifted_panel_name%%.vcf}
    if should_run "$base_lifted_panel_name.$filtered_suffix.bcf.csi"; then
        bcftools view --threads 2 -Ou -i "$lifted_filter" $lifted_panel \
        | bcftools norm --threads 2 -Ou -f $ref_fasta - \
        | bcftools +fill-tags --threads 4 -Ou - -- -t 'AN,AC,MAF,MAC:1=MAC' \
        | bcftools annotate -Ob --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
        | bcftools view --threads 2 -Ob -c $minMAC:minor - \
        > $base_lifted_panel_name.$filtered_suffix.bcf \
        && bcftools index -f --threads 2 $base_lifted_panel_name.$filtered_suffix.bcf \
        && did_run "$base_lifted_panel_name.$filtered_suffix.bcf.csi"
    fi
    reference_dataset=$base_groundtruth_name.$ground_truth_filtered_suffix.bcf
    native_panel=$base_native_panel_name.$filtered_suffix.bcf
    lifted_panel=$base_lifted_panel_name.$filtered_suffix.bcf

    wait_and_check || exit 1

    did_run "$base_groundtruth_name.$ground_truth_filtered_suffix.bcf.csi" || exit 1
    did_run "$base_native_panel_name.$filtered_suffix.bcf.csi" || exit 1
    did_run "$base_lifted_panel_name.$filtered_suffix.bcf.csi" || exit 1

    wait_and_check || exit 1
wait
# 1) Downsample ground truth variants
    ground_truth_downsampled=${reference_dataset%%.bcf}.downsampled.bcf
    ground_truth_downsampled=$working_dir/$(basename $ground_truth_downsampled)
    if should_run "$ground_truth_downsampled.csi"; then
        bcftools isec -Ou -n=2 -w 1 $reference_dataset $omni \
        | bcftools +fill-tags -Ou - -- -t F_MISSING,MAF,HWE \
        | bcftools view --threads 2 -c $minMAC:minor -e "F_MISSING>=0.05 || MAF<=0.01 || HWE<=1e-10" -Ou - \
        | bcftools +setGT -Ob - -- -t a -n u > $ground_truth_downsampled \
        && bcftools index -f --threads 2 $ground_truth_downsampled \
        && did_run "$ground_truth_downsampled.csi"
    fi
    
wait_and_check || exit 1

# trap - ERR EXIT
# exit 0


# 2) Impute downsampled SGDP variants
lifted_imputed=$working_dir_ref/$reference_genome.$chrom.lifted.imputed.bcf
native_imputed=$working_dir_ref/$reference_genome.$chrom.native.imputed.bcf
    # 2a: Native panel
    if should_run "$native_imputed.csi"; then
        $basedir/bin/SHAPEIT5_phase_common_static_v1.1.1 \
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
        && $basedir/bin/impute5_v1.2.0_static \
            --g ${ground_truth_downsampled%%.bcf}.native.prephased.bcf \
            --h $native_panel \
            --m $chrom_map \
            --r $whole_chrom \
            --buffer-region $whole_chrom \
        --o $native_imputed \
        --l $working_dir_ref/$reference_genome.$chrom.native.imputed.log \
            --out-ap-field \
            --contigs-fai $ref_fasta.fai \
            --threads $num_threads \
            $impute5_haploid_arg &
    fi

    #2b) lifted
    if should_run "$lifted_imputed.csi"; then
        $basedir/bin/SHAPEIT5_phase_common_static_v1.1.1 \
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
        && $basedir/bin/impute5_v1.2.0_static \
            --g ${ground_truth_downsampled%%.bcf}.lifted.prephased.bcf \
            --h $lifted_panel \
            --m $chrom_map \
            --r $whole_chrom \
            --buffer-region $whole_chrom \
        --o $lifted_imputed \
        --l $working_dir_ref/$reference_genome.$chrom.lifted.imputed.log \
            --out-ap-field \
            --contigs-fai $ref_fasta.fai \
            --threads $num_threads \
            $impute5_haploid_arg &
    fi

wait_and_check || exit 1

did_run "$native_imputed.csi" || exit 1
did_run "$lifted_imputed.csi" || exit 1

# 2c) While the other scripts are running, create GWAS cutoff files for downstream analysis
# for info_cutoff in 0.3 0.7 0.9; do
#     for imputed_file in $native_imputed $lifted_imputed; do
#         filtered_imputed_file=${imputed_file%%.bcf}.GWAS_filtered.info_cutoff_$info_cutoff.bcf
#         if should_run "$filtered_imputed_file.csi"; then
#             bcftools +fill-tags --threads 4 -Ou $imputed_file -- -t 'AN,AC,MAF,MAC:1=MAC,HWE' \
#             | bcftools view --threads 4 -Ob -i "INFO/INFO>=$info_cutoff && HWE>=1e-6" - \
#             > $filtered_imputed_file \
#             && bcftools index --threads 4 -f $filtered_imputed_file &
#         fi
#     done
# done



# 3) Identify variants in common between imputed datasets
if should_run "$working_dir_ref/common.$chrom.IDs.txt"; then
    python3 $basedir/scripts/get_discordant_multiallelic_sites.py $native_panel $lifted_panel $working_dir_ref/common.$chrom.IDs.txt
    did_run "$working_dir_ref/common.$chrom.IDs.txt" || exit 1
fi
wait_and_check || exit 1

# 4) Subset imputed datasets to variants in common
if should_run "$working_dir_ref/$reference_genome.$chrom.native.imputed.common.bcf"; then
    bcftools view --threads 4 -Ob -i "ID==@$working_dir_ref/common.$chrom.IDs.txt" $native_imputed > $working_dir_ref/$reference_genome.$chrom.native.imputed.common.bcf \
    && bcftools index --threads 4 -f $working_dir_ref/$reference_genome.$chrom.native.imputed.common.bcf &
fi
if should_run "$working_dir_ref/$reference_genome.$chrom.lifted.imputed.common.bcf"; then
    bcftools view --threads 4 -Ob -i "ID==@$working_dir_ref/common.$chrom.IDs.txt" $lifted_imputed > $working_dir_ref/$reference_genome.$chrom.lifted.imputed.common.bcf \
    && bcftools index --threads 4 -f $working_dir_ref/$reference_genome.$chrom.lifted.imputed.common.bcf &
fi
wait_and_check || exit 1

did_run "$working_dir_ref/$reference_genome.$chrom.native.imputed.common.bcf.csi" || exit 1
did_run "$working_dir_ref/$reference_genome.$chrom.lifted.imputed.common.bcf.csi" || exit 1


# 5) Measure imputation accuracy of all four sets (native, lifted, native-common, lifted-common).
# Include filtered_suffix in output filenames to match expected file naming
filter_tag=""
if [[ $no_singletons == 'true' ]]; then
    filter_tag=".no_singletons"
elif [[ $no_singletons == 'confident' ]]; then
    filter_tag=".no_coinflips"
fi

wait_and_check

# Always regenerate just the .txt files for this dataset to avoid cross-run errors
rm -f ${working_dir}/native_panel.${ref_dataset_name}.$chrom.txt \
    ${working_dir}/lifted_panel.${ref_dataset_name}.$chrom.txt \
    ${working_dir}/native_panel.common_variants.${ref_dataset_name}.$chrom.txt \
    ${working_dir}/lifted_panel.common_variants.${ref_dataset_name}.$chrom.txt

native_imputed_base=${native_imputed%%.bcf}
lifted_imputed_base=${lifted_imputed%%.bcf}
echo -e "$whole_chrom $native_panel $reference_dataset $native_imputed" >  ${working_dir}/native_panel.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $lifted_imputed" >  ${working_dir}/lifted_panel.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $native_imputed_base.GWAS_filtered.info_cutoff_0.3.bcf" >  ${working_dir}/native_panel.GWAS_filtered.info_cutoff_0.3.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $lifted_imputed_base.GWAS_filtered.info_cutoff_0.3.bcf" >  ${working_dir}/lifted_panel.GWAS_filtered.info_cutoff_0.3.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $native_imputed_base.GWAS_filtered.info_cutoff_0.7.bcf" >  ${working_dir}/native_panel.GWAS_filtered.info_cutoff_0.7.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $lifted_imputed_base.GWAS_filtered.info_cutoff_0.7.bcf" >  ${working_dir}/lifted_panel.GWAS_filtered.info_cutoff_0.7.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $native_imputed_base.GWAS_filtered.info_cutoff_0.9.bcf" >  ${working_dir}/native_panel.GWAS_filtered.info_cutoff_0.9.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $lifted_imputed_base.GWAS_filtered.info_cutoff_0.9.bcf" >  ${working_dir}/lifted_panel.GWAS_filtered.info_cutoff_0.9.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir_ref/$reference_genome.$chrom.native.imputed.common.bcf" >  ${working_dir}/native_panel.common_variants.${ref_dataset_name}.$chrom.txt
echo -e "$whole_chrom $native_panel $reference_dataset $working_dir_ref/$reference_genome.$chrom.lifted.imputed.common.bcf" >  ${working_dir}/lifted_panel.common_variants.${ref_dataset_name}.$chrom.txt

r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'

for infile in ${working_dir}/native_panel.${ref_dataset_name}.$chrom.txt ${working_dir}/lifted_panel.${ref_dataset_name}.$chrom.txt ${working_dir}/native_panel.common_variants.${ref_dataset_name}.$chrom.txt ${working_dir}/lifted_panel.common_variants.${ref_dataset_name}.$chrom.txt
do
    if should_run "${infile%%.txt}/$(basename ${infile%%.txt}).rsquare.grp.txt.gz"; then
        echo "Assessing imputation accuracy for $infile"
    else
        continue
    fi
    echo $infile
    mkdir -p ${infile%%.txt}
    $basedir/bin/GLIMPSE2_concordance \
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

wait_and_check || exit 1

# Disable trap on successful completion
trap - ERR EXIT
