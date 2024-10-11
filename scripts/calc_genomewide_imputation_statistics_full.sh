#!/usr/bin/env bash
set -euo pipefail
source parameters.sh

GRCh38_suffix=$1
T2T_suffix=$2
run_suffix=$3
num_threads=$4

basedir=$PWD
subset_folder=$basedir/sample_subsets


outfolder=$PWD/genomewide_imputation_evaluation_${run_suffix}
logfolder=$PWD/genomewide_imputation_evaluation_${run_suffix}/working
mkdir -p $outfolder
mkdir -p $logfolder



for dataset in SGDP pangenome
do 
    for genomic_variants in T2T T2T_snps GRCh38 GRCh38_snps 
    do 
        for run in native_panel "native_panel.common_variants" lifted_panel "lifted_panel.common_variants"
        do 
            if [[ $genomic_variants != *GRCh38* ]]; then
                genome="CHM13"
                alt_genome='T2T'
                suffix=$T2T_suffix
            else
                genome=GRCh38
                alt_genome=GRCh38
                suffix=$T2T_suffix
            fi

            if [[ ! -s $outfolder/$genomic_variants.$dataset.$run.txt ]]; then
                cat *_working*${suffix}/${genomic_variants}_imputation*_workspace/*/$run.$dataset.*.txt | sort \
                    > $outfolder/$genomic_variants.$dataset.$run.txt
            fi

            syntenic_nonsyntenic_variant_grouping=$outfolder/${alt_genome}_syntenic-nonsyntenic_overall.tsv
            syntenic_nonsyntenic_vartype_grouping=$outfolder/${alt_genome}_syntenic-nonsyntenic_vartype.tsv
            
            if [[ ! -s $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz ]]; then
                GLIMPSE2_concordance \
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
               GLIMPSE2_concordance \
                    --gt-val \
                    --groups $syntenic_nonsyntenic_variant_grouping \
                    --threads $num_threads \
                    --af-tag MAF \
                    --input $outfolder/$genomic_variants.$dataset.$run.txt \
                    --log $logfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_overall_bins.log \
                    --output $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_overall_bins &
            fi
            if [[ ! -s $outfolder/${genomic_variants}.$dataset.$run.glimpse2_concordance_syntenic_maf_vartype_bins.rsquare.grp.txt.gz ]]; then
                GLIMPSE2_concordance \
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
                genome="CHM13"
                alt_genome='T2T'
                suffix=$T2T_suffix
            else
                genome=GRCh38
                alt_genome=GRCh38
                suffix=$T2T_suffix
            fi

            if [[ $dataset == *SGDP* ]]
            then
                for ancestry_samples in sample_subsets/CHM13_SGDP*_samples.txt
                do
                    ancestry=$(echo $(basename $ancestry_samples) | cut -f 3 -d '_')
                    mkdir -p $outfolder/ancestry_specific
                    if [[ ! -s $outfolder/ancestry_specific/${genomic_variants}.$dataset.$run.$ancestry.glimpse2_concordance_r2_bins.rsquare.grp.txt.gz ]]; then
                        GLIMPSE2_concordance \
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