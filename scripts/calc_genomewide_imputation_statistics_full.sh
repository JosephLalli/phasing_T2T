#!/usr/bin/env bash
set -euo pipefail

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(dirname "$script_path")"
basedir="$(dirname "$script_dir")"
source $script_dir/parameters.sh

# Track background process failures
trap 'jobs -p | xargs -r kill 2>/dev/null; exit 1' ERR EXIT

test_run=false
force_imputation=false


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


# Global array to track expanded commands for background jobs
declare -A expanded_job_cmds

log_expanded_command() {
    local cmd="$*"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting background job: $cmd" >&2
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

GRCh38_run_suffix="${1:-true_false_0.05_false_1_GRCh38}"
CHM13v2_run_suffix="${2:-true_false_0.05_false_1_CHM13v2.0}"
num_threads="${3:-16}"
test_run="${4:-false}"

subset_folder=$basedir/resources/sample_subsets


# Use new imputation_statistics layout
outfolder=$basedir/imputation_statistics/imputation_results_${CHM13v2_run_suffix}
logfolder=$outfolder/logs
mkdir -p $outfolder
mkdir -p $logfolder
mkdir -p $outfolder/ancestry_specific
mkdir -p $outfolder/ancestry_specific/logs


for genome in GRCh38 CHM13v2.0; do
    if [[ ! -s "$outfolder/${genome}_syntenic-nonsyntenic_overall.tsv" || \
          ! -s "$outfolder/${genome}_syntenic-nonsyntenic_vartype.tsv" ]]; then
        if [[ $genome == 'GRCh38' ]]; then
            run_suffix=$GRCh38_run_suffix
        else
            run_suffix=$CHM13v2_run_suffix
        fi
        echo "Identifying binned syntenic/nonsyntenic variants for $genome"
        echo """
            Did not find $outfolder/${genome}_syntenic-nonsyntenic_overall.tsv and/or $outfolder/${genome}_syntenic-nonsyntenic_vartype.tsv
            Need to run:
            $script_dir/create_syn_nonsyn_bins.sh $genome $run_suffix $outfolder $test_run
        """
        exit 1
    fi
done

wait_and_check || exit 1

for dataset in SGDP pangenome
do 
    for genomic_variants in T2T GRCh38 T2T_snps T2T_no_singletons  GRCh38 GRCh38_snps GRCh38_no_singletons #T2T_no_singletons_snps GRCh38_no_singletons_snps T2T_no_coinflips GRCh38_no_coinflips
    do 
        for run in native_panel "native_panel.common_variants" lifted_panel "lifted_panel.common_variants"
        do
            if [[ $genomic_variants != *GRCh38* ]]; then
                genome='T2T'
                group_genome='CHM13v2.0'
                suffix=$CHM13v2_run_suffix
            else
                genome=GRCh38
                group_genome='GRCh38'
                suffix=$CHM13v2_run_suffix
            fi

            echo "$outfolder/$genomic_variants.$dataset.$run.txt"
            rm -f $outfolder/$genomic_variants.$dataset.$run.txt 

            cat $basedir/working_directories/*_working_${suffix}/${genomic_variants}_imputation*_workspace/*/$run.$dataset.*.txt | sort > $outfolder/$genomic_variants.$dataset.$run.txt
            n_contigs=$(wc -l < $outfolder/$genomic_variants.$dataset.$run.txt)
            if [[ "$n_contigs" == "0" ]]; then
                echo "WARNING: No files found for $outfolder/$genomic_variants.$dataset.$run.txt, skipping." >&2
                continue
            fi
            if [[ "$n_contigs" != "25" ]]; then
                echo "NOTE: Found $n_contigs contigs (expected 25) for $outfolder/$genomic_variants.$dataset.$run.txt — proceeding anyway." >&2
            fi


            run_prefix=$outfolder/$genomic_variants.$dataset.$run
            log_prefix=$outfolder/logs/$genomic_variants.$dataset.$run
            if should_run "$run_prefix.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz"; then
                $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --bins $r2_bins \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $run_prefix.txt \
                    --log $log_prefix.glimpse2_concordance_r2_bins.log \
                    --out-r2-per-site \
                    --out-rej-sites	\
                    --out-disc-sites \
                    --output $run_prefix.glimpse2_concordance_r2_bins &
            fi
            if should_run "$run_prefix.glimpse2_concordance_syntenic_maf_overall_bins.rsquare.grp.txt.gz"; then
                $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --groups $outfolder/${group_genome}_syntenic-nonsyntenic_overall.tsv \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $run_prefix.txt \
                    --log $log_prefix.glimpse2_concordance_syntenic_maf_overall_bins.log \
                    --output $run_prefix.glimpse2_concordance_syntenic_maf_overall_bins &
            fi
            if should_run "$run_prefix.glimpse2_concordance_syntenic_maf_vartype_bins.rsquare.grp.txt.gz"; then
                $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --groups $outfolder/${group_genome}_syntenic-nonsyntenic_vartype.tsv \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $run_prefix.txt \
                    --log $log_prefix.glimpse2_concordance_syntenic_maf_vartype_bins.log \
                    --output $run_prefix.glimpse2_concordance_syntenic_maf_vartype_bins &
            fi
        done
    done
    wait_and_check || exit 1
done

### Per-ancestry concordance runs take far less time than the whole sample concordance runs.
### So we'll run the whole sample commands 4 at a time
### But run the per-ancestry concordance runs 16 at a time.

for dataset in SGDP
do 
    for genomic_variants in T2T GRCh38 #T2T_snps T2T_no_singletons T2T_no_singletons_snps GRCh38 GRCh38_snps GRCh38_no_singletons GRCh38_no_singletons_snps #T2T_no_coinflips GRCh38_no_coinflips
    do
        for run in native_panel "native_panel.common_variants" lifted_panel "lifted_panel.common_variants"
        do
            input_location="$outfolder/$genomic_variants.$dataset.$run.txt"
            run_prefix="$outfolder/ancestry_specific/$genomic_variants.$dataset.$run"
            log_prefix="$outfolder/ancestry_specific/logs/$genomic_variants.$dataset.$run"

            if [[ $genome != *GRCh38* ]]; then
                genome="CHM13v2.0"
                suffix=$CHM13v2_run_suffix
            else
                genome=GRCh38
                suffix=$CHM13v2_run_suffix # keep both genomes in same folder for comparison's sake
            fi

            if [[ $dataset == *SGDP* ]]
            then
                for ancestry_samples in $basedir/resources/sample_subsets/CHM13_SGDP*_samples.txt
                do
                    ancestry=$(echo $(basename $ancestry_samples) | cut -f 3 -d '_')
                    if should_run "$run_prefix.$ancestry.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz"; then
                        $basedir/bin/GLIMPSE2_concordance \
                            --samples $ancestry_samples \
                            --gt-val \
                            --bins $r2_bins \
                            --threads $num_threads \
                            --af-tag MAF \
                            --input $input_location \
                            --log $log_prefix.$ancestry.glimpse2_concordance_r2_bins.log \
                            --output $run_prefix.$ancestry.glimpse2_concordance_r2_bins &
                    fi
                done
            fi
        done
    wait_and_check || exit 1
    done
done

# Disable trap on successful completion
trap - ERR EXIT
