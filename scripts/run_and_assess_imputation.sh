#!/usr/bin/env bash
set -eo pipefail
set +u
# set up script log
## one-liner to get script location; credit to https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
basedir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/..

# Track background process failures
trap 'jobs -p | xargs -r kill 2>/dev/null; exit 1' ERR EXIT

usage() {
    echo "Usage: $(basename "$0") [--force] <chromosome> <threads> [suffix] [genome]" >&2
}

force_imputation=false

positional=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)
            force_imputation=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                positional+=("$1")
                shift
            done
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            exit 2
            ;;
        *)
            positional+=("$1")
            shift
            ;;
    esac
done

set -- "${positional[@]}"

if [[ $# -lt 2 ]]; then
    usage
    exit 2
fi

chrom_arg=$1
num_threads=$2
suffix_base=${3:-true_false_0.05_false_1_CHM13v2.0}
genome=${4:-CHM13v2.0}

run_label=${suffix_base:-default}

# Only set up independent logging if not called from parent script
if [[ -z "$BASH_XTRACEFD" || "$BASH_XTRACEFD" != "19" ]]; then
    logfile=$basedir/phasing_${genome}_${chrom_arg}_${run_label}.log
    exec 19>$logfile
    export BASH_XTRACEFD=19
    set -x # writes commands to logfile
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
        echo "ERROR: Expected output $target was not created or is empty." >&2
        return 1
    fi
    # If target is an index file, check the underlying VCF/BCF
    if [[ "$target" == *.vcf* || "$target" == *.bcf* ]]; then
        local target="${target%.csi}"
        target="${target%.tbi}"
        if [[ ! -s "$target" ]]; then
            echo "ERROR: Expected output $target was not created or is empty." >&2
            return 1
        fi
        if [[ "$target" == *.vcf || "$target" == *.vcf.gz || "$target" == *.bcf ]]; then
            local variant_count
            variant_count=$(bcftools index -n "$target" 2>/dev/null) || { echo "ERROR: Expected output $target exists but underlying VCF/BCF is invalid or unreadable." >&2; return 1; }
            if [[ -z "$variant_count" || "$variant_count" -eq 0 ]]; then
                echo "ERROR: Expected output $target exists but underlying VCF/BCF has zero variants." >&2
                return 1
            fi
        fi
    fi
    # File exists and is valid
    return 0
}


# Global array to track expanded commands for background jobs
declare -A expanded_job_cmds

log_expanded_command() {
    local cmd="$*"
    # echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting background job: $cmd" >&2
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

if [[ $force_imputation == true ]]; then
    echo "Force mode enabled: existing outputs will be regenerated where applicable."
fi
chrom=$chrom_arg
suffix=$suffix_base

if [[ $chrom != chr* && $chrom != PAR* ]]; then
    chrom="chr${chrom}"
fi

# exit 1
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

# default values
atomize='true'
missing_to_ref='false'
graph_reference_missing_cutoff=0.05
trim_assemblies_to_callset='false'

# Use rare variant algorithm when phasing pangenome calls
# on variants with an MAC less than this value.
# Set to 1 to skip rare variant algorithm.
minimum_common_MAC=1

# variable for applying an arbitrary bcftools filter before phasing.
# Used for experiment/testing purposes.
# for example, to filter out multiallelic indels, set to: -e '(N_ALT>1)&(TYPE~"indel")'
pre_phasing_filter="-e ''"


if [[ -n ${suffix:-} && ${suffix:0:1} != "_" ]]
then
    suffix="_$suffix"
fi

set -u

drop_reference='GRCh38'
native_maps_insert='native_maps.'
chrom_map=$basedir/resources/recombination_maps/t2t_native_scaled_maps/$(echo $chrom | cut -f 1 -d '_').t2t.scaled.gmap.gz
initial_vcf_calls_folder=$basedir/unphased_variant_calls/t2t

ref_fasta=$basedir/resources/chm13v2.0.fa.gz
pangenome_vcf=$basedir/resources/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
HGSVC_vcf=$basedir/resources/hgsvc3-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz
HGSVC_HPRC_vcf=$basedir/resources/hgsvc3-hprc-2024-02-23-mc-chm13-vcfbub.a100k.wave.norm.vcf.gz

chrom_chunking_coords=$basedir/resources/regions.txt
syntenic_site_location="$basedir/resources/chm13v2-syntenic_to_hg38.bed"
chrom_specific_syntenic_annotation_line_part1="CHROM,FROM,TO"
chrom_specific_syntenic_annotation_line_part2="##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description=\"Syntenic with GRCh38 (source: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed)\">"

PAR1_start='0'
PAR1_end='2394410'
PAR2_start='153925834'
PAR2_end='154259566'

if [[ $chrom == *PAR* ]]; then
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.mixed_ploidy.vcf.gz
else
    input_vcf=$initial_vcf_calls_folder/1KGP.CHM13v2.0.$(echo $chrom | cut -f 1 -d '_').recalibrated.snp_indel.pass.vcf.gz
fi

is_test_run=false

## Identify if we are dealing with a PAR region or a test region
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
elif [[ $chrom == 'chr22_test' ]]
then
    chrom='chr22_test'
    region="chr22:18671427-25461594"
    whole_chrom=$region
    grch38_region="chr22:18000000-25000000"
    chrom_map=$basedir/resources/recombination_maps/t2t_native_scaled_maps/chr22.t2t.scaled.gmap.gz
    end_chrom=25461594
    is_test_run=true
elif [[ $chrom == 'chr15_test' ]]
then
    chrom='chr15_test'
    region="chr15:17904139-25242443"
    whole_chrom=$region
    grch38_region="chr15:20000000-27500000"
    chrom_map=$basedir/resources/recombination_maps/t2t_native_scaled_maps/chr15.t2t.scaled.gmap.gz
    end_chrom=25242443
    is_test_run=true
elif [[ $chrom == 'chrX' ]]
then
    region="chrX:$(($PAR1_end+1))-$PAR2_start"
    whole_chrom=$region
    grch38_region="chrX:2781480-155701383"
    use_beagle='false'
    hmm_ne=$(awk "BEGIN { printf \"%.0f\", ($hmm_ne*.75) }")
    end_chrom=$PAR2_start
else
    region=$chrom
    end_chrom=$(grep "^$chrom	" ${ref_fasta}.fai | cut -f 2)
    whole_chrom=$chrom:1-$end_chrom
    grch38_region=$chrom
fi

chrom_working_dir=$basedir/working_directories/${chrom}_working${suffix}
final_panel_dir=$basedir/phased_panels/phased_${genome}_panel${suffix}
stats_dir=$basedir/SHAPEIT5_switch_output/phasing_stats_${genome}${suffix}
lifted_panel_folder=$chrom_working_dir/liftover/lifted_panels_JLL_post_fix
imputation_results_dir=$basedir/imputation_statistics/imputation_results$suffix
variant_frequency_stats_dir=$basedir/intermediate_data/variant_frequency_stats/${genome}

mkdir -p $chrom_working_dir
mkdir -p $final_panel_dir
mkdir -p $stats_dir
mkdir -p $imputation_results_dir
mkdir -p $lifted_panel_folder
mkdir -p $variant_frequency_stats_dir

population_ids=$basedir/resources/sample_subsets/unrelated_superpopulations.csv



# Assign chrom-specific regions filename
chrom_regions=$chrom_working_dir/${chrom}_regions.txt

echo $chrom


if [[ $chrom == *PAR* ]]
then
    echo $region > $chrom_regions
elif [[ $is_test_run == true ]]
then
    echo $region > $chrom_regions
    chrom=$(echo $chrom | cut -f 1 -d '_')
else
    grep $chrom: $chrom_chunking_coords > $chrom_regions
fi

if [[ $chrom == 'chrX' ]]
then
    # Make list of males for X chrom phasing
    male_sample_list=$basedir/resources/sample_subsets/males.txt
    cat $basedir/resources/sample_subsets/superpopulations.samples.txt | grep " 1 " | cut -f 2 -d " " > $male_sample_list
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
vcf_to_phase_pangenome_biallelic_1kgp=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.1kGP_pangenome_calls.biallelic.bcf
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

# pedigrees
pedigree=$basedir/resources/pedigrees/1kgp.ped
duos_and_trios=$basedir/resources/pedigrees/duos_and_trios.txt
bcftools_formatted_pedigree=$basedir/resources/pedigrees/trios_only.ped

# lists of different categories of samples
no_parents=$basedir/resources/sample_subsets/not_parents.txt
unrelated_samples=$basedir/resources/sample_subsets/unrelated_samples.txt
pangenome_samples=$basedir/resources/sample_subsets/pangenome_samples.txt
pangenome_and_parents=$basedir/resources/sample_subsets/pangenome_samples_and_parents.txt
SGDP_in_1KGP=$basedir/resources/sample_subsets/SGDP_samples_in_1KGP_numbered_format.txt
female_samples=$basedir/resources/sample_subsets/females.txt

# Imputation and imputation metric gathering

GRCh38_fasta=$basedir/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
T2T_fasta=$basedir/resources/chm13v2.0.fa.gz

GRCh38_to_t2t_chain=$basedir/resources/grch38-chm13v2.chain
t2t_to_GRCh38_chain=$basedir/resources/chm13v2-grch38.chain

# I'm sticking with filenames as present on the HPRC site for consistency
# But for our purposes, this is backwards: grch38-chm13v2 is used to lift from T2T to GRCh38

T2T_to_GRCh38_diffs=$basedir/resources/grch38-chm13v2.sort.vcf.gz
GRCh38_to_t2t_diffs=$basedir/resources/chm13v2-grch38.sort.vcf.gz

grch38_syntenic_site_location="$basedir/resources/hg38.GCA_009914755.4.synNet.summary.bed.gz"
t2t_syntenic_site_location="$basedir/resources/chm13v2-syntenic_to_hg38.bed"

grch38_chrom_working_dir=$basedir/working_directories/${chrom_arg}_working$(echo $suffix | sed 's/CHM13v2.0/GRCh38/')

GRCh38_lifted_panel=$lifted_panel_folder/$(echo $(basename $phased_panel_vcf_2504_biallelic) | sed 's/CHM13v2.0/GRCh38.lifted_from_CHM13v2.0/' | sed 's/native_maps.//' )
GRCh38_native_panel=$basedir/phased_panels/grch38/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
T2T_lifted_panel=$lifted_panel_folder/$(echo $(basename $GRCh38_native_panel) | sed 's/\.GRCh38\./.CHM13v2.0.lifted_from_GRCh38./')
T2T_native_panel=$phased_panel_vcf_2504_biallelic

SGDP_ground_truth_dir_T2T=$basedir/resources/SGDP_variation/t2t
SGDP_ground_truth_T2T=$chrom_working_dir/SGDP.CHM13v2.0.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf
SGDP_ground_truth_dir_GRCh38=$basedir/resources/SGDP_variation/grch38
SGDP_ground_truth_GRCh38=$chrom_working_dir/SGDP.GRCh38.${chrom}.recalibrated.no_1KGP_overlaps.biallelic.snp_indel.pass.bcf

pangenome_ground_truth_GRCh38=$grch38_chrom_working_dir/${chrom}_reference_pangenome.biallelic.bcf
pangenome_ground_truth_T2T=$chr_specific_reference_pangenome_variation_biallelic

T2T_lifted_panel_no_pangenome=$lifted_panel_folder/1KGP.CHM13v2.0.lifted_from_GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf
GRCh38_lifted_panel_no_pangenome=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf
T2T_native_panel_no_pangenome=$phased_panel_no_pangenome_biallelic
GRCh38_native_panel_no_pangenome=$basedir/phased_panels/grch38/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.biallelic.2504.bcf

GRCh38_lifted_panel_no_singletons=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.no_singletons.biallelic.2504.bcf
T2T_native_panel_no_singletons=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.no_singletons.biallelic.2504.bcf
GRCh38_lifted_panel_no_singletons_no_pangenome=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.no_singletons.biallelic.2504.bcf
T2T_native_panel_no_singletons_no_pangenome=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.no_singletons.biallelic.2504.bcf

GRCh38_lifted_panel_confident_singletons=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.confident_singletons.biallelic.2504.bcf
T2T_native_panel_confident_singletons=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.confident_singletons.biallelic.2504.bcf
GRCh38_lifted_panel_confident_singletons_no_pangenome=$lifted_panel_folder/1KGP.GRCh38.lifted_from_CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.confident_singletons.biallelic.2504.bcf
T2T_native_panel_confident_singletons_no_pangenome=$chrom_working_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.nopangenome.confident_singletons.biallelic.2504.bcf

if should_run $phased_panel_no_pangenome_biallelic.csi ; then
    echo "Making $phased_panel_no_pangenome_biallelic"
    bcftools view -Ou --threads 8 -S ^$pangenome_and_parents --force-samples $phased_panel_vcf_2504_biallelic 2> /dev/null \
    | bcftools view -Ou --threads 8 -c 1:minor - \
    | bcftools annotate -Ou -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF - \
    | bcftools +fill-tags -Ob --threads 8 - -- -t AN,AC,MAF,MAC:1=MAC  > $phased_panel_no_pangenome_biallelic \
    && bcftools index --threads 8 -f $phased_panel_no_pangenome_biallelic  # &
fi


# Liftover panels if necessary
T2T_native_panel=${rare_variants_phased_ped_biallelic%%.bcf}.2504.bcf
if should_run "$T2T_native_panel.csi" && [[ -s $rare_variants_phased_ped_biallelic.csi ]]; then
    echo "Creating $T2T_native_panel from $rare_variants_phased_ped_biallelic ..."
    bcftools view -Ob --threads 4 -S $unrelated_samples $rare_variants_phased_ped_biallelic > $T2T_native_panel && bcftools index -f --threads 4 $T2T_native_panel \
    && did_run $T2T_native_panel.csi
fi

if should_run "$GRCh38_lifted_panel.csi"; then
    echo "lifting $T2T_native_panel to $GRCh38_lifted_panel"
    $basedir/scripts/liftover_panel.sh -i $T2T_native_panel \
                                        -o $GRCh38_lifted_panel \
                                        -t $GRCh38_fasta \
                                        -s $T2T_fasta \
                                        -c $t2t_to_GRCh38_chain \
                                        -d $T2T_to_GRCh38_diffs \
                                        -r $whole_chrom  && did_run $GRCh38_lifted_panel &
fi

if should_run "$T2T_lifted_panel.csi"; then
    echo "lifting $GRCh38_native_panel to $T2T_lifted_panel"
    $basedir/scripts/liftover_panel.sh -i $GRCh38_native_panel \
                                        -o $T2T_lifted_panel \
                                        -t $T2T_fasta \
                                        -s $GRCh38_fasta \
                                        -c $GRCh38_to_t2t_chain \
                                        -d $GRCh38_to_t2t_diffs \
                                        -r $grch38_region && did_run $T2T_lifted_panel &
fi

wait_and_check || echo "WARNING: Some liftover jobs failed (non-fatal; lifted panels may be incomplete)" >&2

if [[ $genome == 'GRCh38' ]]; then
    phased_panel_vcf_3202_GRCh38=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.3202.bcf
    phased_panel_vcf_2504_GRCh38=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.2504.bcf
    phased_panel_vcf_3202_biallelic_GRCh38=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.${native_maps_insert}biallelic.3202.bcf
    phased_panel_vcf_2504_biallelic_GRCh38=$final_panel_dir/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.${native_maps_insert}biallelic.2504.bcf
else
    phased_panel_vcf_3202_CHM13=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf
    phased_panel_vcf_2504_CHM13=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.2504.bcf
    phased_panel_vcf_3202_biallelic_CHM13=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf
    phased_panel_vcf_2504_biallelic_CHM13=$final_panel_dir/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf
fi
fully_annotated_input_variant_report=$chrom_working_dir/1KGP.${genome}.${chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv

if should_run "$T2T_native_panel_no_singletons.csi"; then
    bcftools view -c 2:minor -Ou --threads 2 $phased_panel_vcf_3202_biallelic_CHM13 \
    | bcftools view -Ob --threads 2 -S $unrelated_samples - \
    > $T2T_native_panel_no_singletons && \
    bcftools index -f --threads 4 $T2T_native_panel_no_singletons &&
    did_run $T2T_native_panel_no_singletons \
    && echo "Created $T2T_native_panel_no_singletons, $(bcftools index -n $T2T_native_panel_no_singletons) variants"
fi

if should_run "$GRCh38_lifted_panel_no_singletons.csi"; then
    echo "lifting $T2T_native_panel_no_singletons to $GRCh38_lifted_panel_no_singletons"
    $basedir/scripts/liftover_panel.sh -i $T2T_native_panel_no_singletons \
                                        -o $GRCh38_lifted_panel_no_singletons \
                                        -t $GRCh38_fasta \
                                        -s $T2T_fasta \
                                        -c $t2t_to_GRCh38_chain \
                                        -d $T2T_to_GRCh38_diffs \
                                        -r $whole_chrom  && did_run $GRCh38_lifted_panel_no_singletons &

fi

# if should_run "$GRCh38_lifted_panel_no_singletons.csi"; then
#     echo "lifting $T2T_native_panel_no_singletons to $GRCh38_lifted_panel_no_singletons"
#     $basedir/scripts/liftover_panel_JLL.sh -i $T2T_native_panel_no_singletons \
#                                         -o $GRCh38_lifted_panel_no_singletons \
#                                         -t $GRCh38_fasta \
#                                         -s $T2T_fasta \
#                                         -c $t2t_to_GRCh38_chain \
#                                         -d $t2t_to_GRCh38_chain \
#                                         -r $whole_chrom  && did_run $GRCh38_lifted_panel_no_singletons  &
# fi

# if should_run "$GRCh38_lifted_panel_confident_singletons.csi"; then
#     bcftools view -i '(COUNT(FORMAT/PP==0.5)==0)' -Ob --threads 4 $GRCh38_lifted_panel > $GRCh38_lifted_panel_confident_singletons && bcftools index -f --threads 4 $GRCh38_lifted_panel_confident_singletons \
#     && echo "Created $GRCh38_lifted_panel_confident_singletons, $(bcftools index -n $GRCh38_lifted_panel_confident_singletons) variants" # &
# fi

# if should_run "$T2T_native_panel_confident_singletons.csi"; then
#     bcftools view -i '(COUNT(FORMAT/PP==0.5)==0)' -Ob --threads 4 $T2T_native_panel > $T2T_native_panel_confident_singletons && bcftools index -f --threads 4 $T2T_native_panel_confident_singletons \
#     && echo "Created $T2T_native_panel_confident_singletons, $(bcftools index -n $T2T_native_panel_confident_singletons) variants" # &
# fi

wait_and_check || exit 1

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


if should_run "$SGDP_ground_truth_T2T.csi"; then
    echo "Making $SGDP_ground_truth_T2T"
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
if should_run "$SGDP_ground_truth_GRCh38.csi"; then
    echo "Making $SGDP_ground_truth_GRCh38"
    bcftools norm --threads 2 -Ou -r $grch38_region -f $GRCh38_fasta -m -any $SGDP_variants_GRCh38 \
    | bcftools annotate -Ou --threads 8 -a $grch38_syntenic_site_location -c CHROM,FROM,TO,SYNTENIC --mark-sites +SYNTENIC \
                    -H '##INFO=<ID=SYNTENIC,Number=0,Type=Flag,Description="Syntenic with T2T)">' \
                    -x INFO/MAC,INFO/AN,INFO/AC,INFO/MAF --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
    | bcftools view -Ou --force-samples -S ^$SGDP_in_1KGP - 2> /dev/null \
    | bcftools view -Ou  -c 1:minor -i "ALT!='*' || F_MISSING<0.05 || ABS(ILEN)<=50" - \
    | bcftools +fill-tags --threads 8 -Ob - -- -t AN,AC,MAF,MAC:1=MAC \
    > $SGDP_ground_truth_GRCh38 \
    && bcftools index --threads 8 -f $SGDP_ground_truth_GRCh38 &
fi
wait_and_check || exit 1


### SGDP Imputation Accuracy ###

#if should_run "$chrom_working_dir/GRCh38_imputation_workspace/GRCh38_space_filtered/lifted_panel.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel false false SGDP $num_threads &
#fi
#if should_run "$chrom_working_dir/T2T_imputation_workspace/T2T_space_filtered/lifted_panel.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel false false SGDP $num_threads &
#fi

# Assess imputation accuracy when imputing with SNPs only
#if should_run "$chrom_working_dir/GRCh38_snps_imputation_workspace/GRCh38_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel true false SGDP $num_threads &
#fi
#if should_run "$chrom_working_dir/T2T_snps_imputation_workspace/T2T_space_filtered_snpsOnly/lifted_panel.common_variants.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel $T2T_lifted_panel true false SGDP $num_threads &
#fi
wait_and_check || exit 1

# Assess imputation accuracy when removing singletons
#if should_run "$chrom_working_dir/GRCh38_no_singletons_imputation_workspace/GRCh38_space_filtered_no_singletons/lifted_panel.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_singletons_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel_no_singletons false true SGDP $num_threads &
#fi
#if should_run "$chrom_working_dir/T2T_no_singletons_imputation_workspace/T2T_space_filtered_no_singletons/lifted_panel.SGDP.${chrom}.txt"; then
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_singletons_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel_no_singletons $T2T_lifted_panel false true SGDP $num_threads &
#fi
# # Assess imputation accuracy only when imputing with SNPs only, removing singletons
# if should_run "$chrom_working_dir/GRCh38_no_singletons_snps_imputation_workspace/GRCh38_space_filtered_snpsOnly_no_singletons/SGDP/GRCh38.${chrom}.native.imputed.common.bcf"; then # lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_singletons_snps_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel_no_singletons true true SGDP  $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_singletons_snps_imputation_workspace/T2T_space_filtered_snpsOnly_no_singletons/SGDP/T2T.${chrom}.native.imputed.common.bcf"; then #lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_singletons_snps_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel_no_singletons $T2T_lifted_panel true true SGDP  $num_threads &
# fi
wait_and_check || exit 1
# # Assess imputation accuracy when removing variants that contain one phasing call with a PP of 0.5 (coinflips)
# if should_run "$chrom_working_dir/GRCh38_no_coinflips_imputation_workspace/GRCh38_space_filtered_no_coinflips/lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_coinflips_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel_confident_singletons false confident SGDP $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_coinflips_imputation_workspace/T2T_space_filtered_no_coinflips/lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_coinflips_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel_confident_singletons $T2T_lifted_panel false confident SGDP  $num_threads &
# fi

# # Assess imputation accuracy only when removing coinflips and indels (snpsOnly)
# if should_run "$chrom_working_dir/GRCh38_no_coinflips_imputation_workspace/GRCh38_space_filtered_snpsOnly_no_coinflips/SGDP/GRCh38.${chrom}.native.imputed.common.bcf"; then #lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_coinflips_imputation_workspace $SGDP_ground_truth_GRCh38 $GRCh38_native_panel $GRCh38_lifted_panel_confident_singletons true confident SGDP  $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_coinflips_imputation_workspace/T2T_space_filtered_snpsOnly_no_coinflips/SGDP/T2T.${chrom}.native.imputed.common.bcf"; then #lifted_panel.SGDP.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_coinflips_imputation_workspace $SGDP_ground_truth_T2T $T2T_native_panel_confident_singletons $T2T_lifted_panel true confident SGDP  $num_threads &
# fi
# 
if should_run "$GRCh38_lifted_panel_no_pangenome.csi"; then
    echo $chrom "GRCh38_lifted_panel_no_pangenome"
    bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_lifted_panel > $GRCh38_lifted_panel_no_pangenome && \
    bcftools index -f --threads 2 $GRCh38_lifted_panel_no_pangenome &
fi
if should_run "$T2T_lifted_panel_no_pangenome.csi"; then
    echo $chrom "T2T_lifted_panel_no_pangenome"
    bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $T2T_lifted_panel > $T2T_lifted_panel_no_pangenome && \
    bcftools index -f --threads 2 $T2T_lifted_panel_no_pangenome &
fi

if should_run "$GRCh38_native_panel_no_pangenome.csi"; then
    echo $chrom "GRCh38_native_panel_no_pangenome"
    bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_native_panel > $GRCh38_native_panel_no_pangenome && \
    bcftools index -f --threads 2 $GRCh38_native_panel_no_pangenome &
fi
if should_run "$GRCh38_lifted_panel_no_singletons_no_pangenome.csi"; then
    echo $chrom "GRCh38_lifted_panel_no_singletons_no_pangenome"
    bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_lifted_panel_no_singletons > $GRCh38_lifted_panel_no_singletons_no_pangenome && \
    bcftools index -f --threads 2 $GRCh38_lifted_panel_no_singletons_no_pangenome
fi
if should_run "$T2T_native_panel_no_singletons_no_pangenome.csi"; then
    echo $chrom "T2T_native_panel_no_singletons_no_pangenome"
    bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $T2T_native_panel_no_singletons > $T2T_native_panel_no_singletons_no_pangenome && \
    bcftools index -f --threads 2 $T2T_native_panel_no_singletons_no_pangenome
fi

# if should_run "$GRCh38_lifted_panel_confident_singletons_no_pangenome.csi"; then
#     echo $chrom "GRCh38_lifted_panel_confident_singletons_no_pangenome"
#     bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $GRCh38_lifted_panel_confident_singletons > $GRCh38_lifted_panel_confident_singletons_no_pangenome && \
#     bcftools index --threads 2 $GRCh38_lifted_panel_confident_singletons_no_pangenome &
# fi
# if should_run "$T2T_native_panel_confident_singletons_no_pangenome.csi"; then
#     echo $chrom "T2T_native_panel_confident_singletons_no_pangenome"
#     bcftools view -Ob --threads 6 --force-samples -S ^$pangenome_and_parents $T2T_native_panel_confident_singletons > $T2T_native_panel_confident_singletons_no_pangenome && \
#     bcftools index --threads 2 $T2T_native_panel_confident_singletons_no_pangenome &
# fi

wait_and_check || exit 1
wait
# ### Pangenome Imputation Accuracy ###
# # Assess imputation accuracy
#if should_run "$chrom_working_dir/GRCh38_imputation_pangenome_workspace/GRCh38_space_filtered/pangenome/GRCh38.${chrom}.native.imputed.common.bcf"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_imputation_pangenome_workspace      $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome false false pangenome  $num_threads &
#fi
#if should_run "$chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace/GRCh38_space_filtered_snpsOnly/pangenome/GRCh38.${chrom}.native.imputed.common.bcf"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_snps_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_pangenome true false pangenome  $num_threads &
#fi
## Assess imputation accuracy when imputing with SNPs only
#if should_run "$chrom_working_dir/T2T_imputation_pangenome_workspace/T2T_space_filtered/pangenome/T2T.${chrom}.native.imputed.common.bcf"; then    
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_imputation_pangenome_workspace            $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome false false pangenome  $num_threads &
#fi
#if should_run "$chrom_working_dir/T2T_snps_imputation_pangenome_workspace/T2T_space_filtered_snpsOnly/pangenome/T2T.${chrom}.native.imputed.common.bcf"; then
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_snps_imputation_pangenome_workspace       $pangenome_ground_truth_T2T    $T2T_native_panel_no_pangenome    $T2T_lifted_panel_no_pangenome true false pangenome  $num_threads &
#fi
## wait_and_check || exit 1

# Assess imputation accuracy when removing singletons
#if should_run "$chrom_working_dir/GRCh38_no_singletons_imputation_pangenome_workspace/GRCh38_space_filtered_no_singletons/pangenome/GRCh38.${chrom}.native.imputed.common.bcf"; then
    $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_singletons_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_singletons_no_pangenome false true pangenome  $num_threads &
#fi
#if should_run "$chrom_working_dir/T2T_no_singletons_imputation_pangenome_workspace/T2T_space_filtered_no_singletons/pangenome/T2T.${chrom}.native.imputed.common.bcf"; then
    $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_singletons_imputation_pangenome_workspace $pangenome_ground_truth_T2T $T2T_native_panel_no_singletons_no_pangenome $T2T_lifted_panel_no_pangenome false true pangenome $num_threads &
#fi
# # Assess imputation accuracy only when removing singletons and indels (snpsOnly)
# if should_run "$chrom_working_dir/GRCh38_no_singletons_snps_imputation_pangenome_workspace/GRCh38_space_filtered_snpsOnly_no_singletons/pangenome/GRCh38.${chrom}.native.imputed.common.bcf"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_singletons_snps_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_no_singletons_no_pangenome true true pangenome  $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_singletons_snps_imputation_pangenome_workspace/T2T_space_filtered_snpsOnly_no_singletons/pangenome/T2T.${chrom}.native.imputed.common.bcf"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_singletons_snps_imputation_pangenome_workspace $pangenome_ground_truth_T2T $T2T_native_panel_no_singletons_no_pangenome $T2T_lifted_panel_no_pangenome true true pangenome  $num_threads &
# fi

wait_and_check || exit 1
# # Assess imputation accuracy when removing variants that contain one phasing call with a PP of 0.5 (coinflips)
# if should_run "$chrom_working_dir/GRCh38_no_coinflips_imputation_pangenome_workspace/GRCh38_space_filtered_no_coinflips/lifted_panel.pangenome.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_coinflips_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_confident_singletons_no_pangenome false confident pangenome  $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_coinflips_imputation_pangenome_workspace/T2T_space_filtered_no_coinflips/lifted_panel.pangenome.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_coinflips_imputation_pangenome_workspace $pangenome_ground_truth_T2T $T2T_native_panel_confident_singletons_no_pangenome $T2T_lifted_panel_no_pangenome false confident pangenome  $num_threads &
# fi
# # Assess imputation accuracy only when removing singletons and indels (snpsOnly)
# if should_run "$chrom_working_dir/GRCh38_no_coinflips_imputation_pangenome_workspace/GRCh38_space_filtered_snpsOnly_no_coinflips/lifted_panel.pangenome.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh GRCh38 $chrom $chrom_working_dir/GRCh38_no_coinflips_imputation_pangenome_workspace $pangenome_ground_truth_GRCh38 $GRCh38_native_panel_no_pangenome $GRCh38_lifted_panel_confident_singletons_no_pangenome true confident pangenome  $num_threads &
# fi
# if should_run "$chrom_working_dir/T2T_no_coinflips_imputation_pangenome_workspace/T2T_space_filtered_snpsOnly_no_coinflips/lifted_panel.pangenome.${chrom}.txt"; then
#     $basedir/scripts/assess_imputation.sh T2T $chrom $chrom_working_dir/T2T_no_coinflips_imputation_pangenome_workspace $pangenome_ground_truth_T2T $T2T_native_panel_confident_singletons_no_pangenome $T2T_lifted_panel_no_pangenome true confident pangenome  $num_threads &
# fi


wait_and_check || exit 1

# Disable trap on successful completion
trap - ERR EXIT

# close logfile
set +x
unset BASH_XTRACEFD
exec 19>&- || true

echo "Phasing and imputation assessment for $chrom on $genome panel completed."
