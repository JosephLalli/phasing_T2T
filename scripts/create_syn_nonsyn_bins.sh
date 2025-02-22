#!/usr/bin/env bash
set -euo pipefail

genome=$1
genome_suffix=$2
outdir_suffix=$3

basedir=$PWD

outfolder=$basedir/genomewide_imputation_evaluation_${outdir_suffix}
mkdir -p $outfolder

r2_bins="0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5"

for chrom in $(seq 1 14) X $(seq 15 22) PAR1 PAR2
do
    if [[ $chrom != *PAR* ]]; then
        chrom=chr$chrom
    fi
    
    chrom_working_dir=${chrom}_working_${genome_suffix}

    if [[ $genome == 'GRCh38' ]]; then
        phased_panel_location=$basedir/phased_GRCh38_panel
        phased_panel_vcf_2504_biallelic=$phased_panel_location/1KGP.GRCh38.${chrom}.recalibrated.snp_indel.pass.phased.biallelic.2504.bcf
    else
        phased_panel_location=$basedir/phased_T2T_panel_${genome_suffix}
        phased_panel_vcf_2504_biallelic=$phased_panel_location/1KGP.CHM13v2.0.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf
    fi
    ## Loop through each r2 bin number
    ## If bin number is 0, then continue. At start, bin is defined as start=0, end=0.00021.
    ## For each bin, use bcftools view to stream variants from the reference panel that have MAFs within the bin.
    ## Use bcftool query to record the chrom/pos/ref/alt of these variatns, and their label:
    ## SYNTENIC (or NONSYNTENIC) and bin (labeled $start_bin-$end_bin)
    ## Do this in parallel. My machine can handle running this for all bins for one chromosome at a time.
    for end_bin in $r2_bins; do
        if [[ $end_bin == 0 ]]; then
            start_bin=$end_bin
            continue
        fi
        bcftools view --threads 6 -i "INFO/MAF>$start_bin && INFO/MAF<=$end_bin" $phased_panel_vcf_2504_biallelic \
            | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYNTENIC\_$start_bin-$end_bin\n" - \
            | sed -e 's/\t1_/\tSYNTENIC_/' | sed -e 's/\t\._/\tNONSYNTENIC_/' > $chrom_working_dir/syntenic-nonsyntenic_overall.$end_bin.tsv &

        start_bin=$end_bin
    done

    wait

    # Repeat, this time including the type (SNP/INDEL) in the variant bin label.
    ## (Is this wasteful of CPU resources? Hell yes! But it's also the method the requires
    ##  the least custom coding on my part, which seems wise.)
    for end_bin in $r2_bins; do
        if [[ $end_bin == 0 ]]; then
            start_bin=$end_bin
            continue
        fi
        bcftools view --threads 6 -i "INFO/MAF>$start_bin && INFO/MAF<=$end_bin" $phased_panel_vcf_2504_biallelic \
            | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYNTENIC\_%TYPE\_$start_bin-$end_bin\n" - \
            | sed -e 's/\t1_/\tSYNTENIC_/' | sed -e 's/\t\._/\tNONSYNTENIC_/' > $chrom_working_dir/syntenic-nonsyntenic_vartype.$end_bin.tsv &

        start_bin=$end_bin
    done

    wait
done

# Once done, concat and sort all per-bin, per-contig labels into one whole-genome label tsv
cat *working*_${genome_suffix}/syntenic-nonsyntenic_vartype.*.tsv | sort -k1,1 -k2,2n > genomewide_imputation_evaluation_${outdir_suffix}/${genome}_syntenic-nonsyntenic_vartype.tsv
cat *working*_${genome_suffix}/syntenic-nonsyntenic_overall.*.tsv | sort -k1,1 -k2,2n > genomewide_imputation_evaluation_${outdir_suffix}/${genome}_syntenic-nonsyntenic_overall.tsv
wait

# Remove intermediate files
rm -f *working*_$genome_suffix/syntenic-nonsyntenic_overall.*.tsv
rm -f *working*_$genome_suffix/syntenic-nonsyntenic_vartype.*.tsv