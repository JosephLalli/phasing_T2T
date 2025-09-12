#!/usr/bin/env bash
set -euo pipefail

script_path="$(realpath "${BASH_SOURCE[0]}")"
script_dir="$(dirname "$script_path")"
basedir="$(dirname "$script_dir")"
source $script_dir/parameters.sh

test_run=false


GRCh38_run_suffix=$1
CHM13v2_run_suffix=$2
num_threads=$3
test_run=$4

subset_folder=$basedir/resources/sample_subsets


# Use new imputation_statistics layout
outfolder=$basedir/imputation_statistics/imputation_results_${CHM13v2_run_suffix}
logfolder=$outfolder/logs
mkdir -p $outfolder
mkdir -p $logfolder


for alt_genome in GRCh38 CHM13v2.0; do
    if [[ ! -s $outfolder/${alt_genome}_syntenic-nonsyntenic_overall.tsv ]]; then
        if [[ $alt_genome == 'GRCh38' ]]; then
            run_suffix=$GRCh38_run_suffix
        else
            run_suffix=$CHM13v2_run_suffix
        fi
        echo "Identifying binned syntenic/nonsyntenic variants for $alt_genome"
        $script_dir/create_syn_nonsyn_bins.sh $alt_genome $run_suffix $outfolder $test_run
    fi
done

wait

for dataset in SGDP pangenome
do 
    for genomic_variants in T2T T2T_snps GRCh38 GRCh38_snps 
    do 
        for run in native_panel "native_panel.common_variants" lifted_panel "lifted_panel.common_variants"
        do 
            if [[ $genomic_variants != *GRCh38* ]]; then
                genome="CHM13v2.0"
                alt_genome='T2T'
                suffix=$CHM13v2_run_suffix
            else
                genome=GRCh38
                alt_genome=GRCh38
                suffix=$CHM13v2_run_suffix
            fi
            

            echo "$outfolder/$genomic_variants.$dataset.$run.txt"
            if [[ ! -s $outfolder/$genomic_variants.$dataset.$run.txt ]]; then
                cat $basedir/working_directories/*_working_${suffix}/${genomic_variants}_imputation*_workspace/*/$run.$dataset.*.txt | sort \
                    > $outfolder/$genomic_variants.$dataset.$run.txt
            fi

            syntenic_nonsyntenic_variant_grouping=$outfolder/${genome}_syntenic-nonsyntenic_overall.tsv
            syntenic_nonsyntenic_vartype_grouping=$outfolder/${genome}_syntenic-nonsyntenic_vartype.tsv

            if [[ ! -s $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz ]]; then
                $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --bins $r2_bins \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $outfolder/$genomic_variants.$dataset.$run.txt \
                    --log $logfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_r2_bins.log \
                    --out-r2-per-site \
                    --out-rej-sites	\
                    --out-disc-sites \
                    --output $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_r2_bins &
            fi
            if [[ ! -s $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_overall_bins.rsquare.grp.txt.gz ]]; then
               $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --groups $syntenic_nonsyntenic_variant_grouping \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $outfolder/$genomic_variants.$dataset.$run.txt \
                    --log $logfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_overall_bins.log \
                    --output $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_overall_bins &
            fi
            if [[ ! -s $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_vartype_bins.rsquare.grp.txt.gz ]]; then
                $basedir/bin/GLIMPSE2_concordance \
                    --gt-val \
                    --groups $syntenic_nonsyntenic_vartype_grouping \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $outfolder/$genomic_variants.$dataset.$run.txt \
                    --log $logfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_vartype_bins.log \
                    --output $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_vartype_bins &
            fi
        done
        wait
    done
done

### Per-ancestry concordance runs take far less time than the whole sample concordance runs.
### So we'll run the whole sample commands 4 at a time
### But run the per-ancestry concordance runs 16 at a time.

for dataset in SGDP pangenome
do 
    for genomic_variants in T2T T2T_snps GRCh38 GRCh38_snps 
    do 
        for run in native_panel "native_panel.common_variants" lifted_panel "lifted_panel.common_variants"
        do 
            if [[ $genomic_variants != *GRCh38* ]]; then
                genome="CHM13v2.0"
                alt_genome='T2T'
                suffix=$CHM13v2_run_suffix
            else
                genome=GRCh38
                alt_genome=GRCh38
                suffix=$CHM13v2_run_suffix
            fi

            if [[ $dataset == *SGDP* ]]
            then
                for ancestry_samples in $basedir/resources/sample_subsets/CHM13_SGDP*_samples.txt
                do
                    ancestry=$(echo $(basename $ancestry_samples) | cut -f 3 -d '_')
                    mkdir -p $outfolder/ancestry_specific
                    if [[ ! -s $outfolder/ancestry_specific/${genomic_variants}.$dataset.$run.$ancestry.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz ]]; then
                        $basedir/bin/GLIMPSE2_concordance \
                            --samples $ancestry_samples \
                            --gt-val \
                            --bins $r2_bins \
                            --threads $num_threads \
                            --af-tag MAF \
                            --input $outfolder/$genomic_variants.$dataset.$run.txt \
                            --log $logfolder/${genomic_variants}.$dataset.$run.$ancestry.glimpse2_concordance_r2_bins.log \
                            --output $outfolder/ancestry_specific/${genomic_variants}.$dataset.$run.$ancestry.glimpse2_concordance_r2_bins &
                    fi
                done
            fi
        done
    done
    wait
done