#!/usr/bin/env bash
set -eo pipefail

n_singleton_jobs=$1
suffix=$2

basedir=$PWD

for i in $(seq 1 14) X $(seq 15 22); do
if [[ ! -s $basedir/chr${i}_working_${suffix}/chr${i}_private_singletons.txt ]]; then
    mkdir -p $basedir/chr${i}_working_${suffix}/tmp
    cat $basedir/pedigrees/duos_and_trios.txt     | parallel -j $n_singleton_jobs "bcftools view --force-samples -H -G -s {} -x -c 2 $basedir/phased_T2T_panel_${suffix}/1KGP.CHM13v2.0.chr${i}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf | cut -f 3 > $basedir/chr${i}_working_${suffix}/tmp/{#}.txt" &&     cat $basedir/chr${i}_working_${suffix}/tmp/*.txt | sort | uniq > $basedir/chr${i}_working_${suffix}/chr${i}_private_singletons.txt && rm -rf $basedir/chr${i}_working_${suffix}/tmp
fi

for i in PAR1 PAR2; do
if [[ ! -s $basedir/${i}_working_${suffix}/${i}_private_singletons.txt ]]; then
    mkdir -p $basedir/${i}_working_${suffix}/tmp
    cat $basedir/pedigrees/duos_and_trios.txt     | parallel -j  $n_singleton_jobs "bcftools view --force-samples -H -G -s {} -x -c 2 $basedir/phased_T2T_panel_${suffix}/1KGP.CHM13v2.0.${i}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf | cut -f 3 > $basedir/${i}_working_${suffix}/tmp/{#}.txt" &&     cat $basedir/${i}_working_${suffix}/tmp/*.txt | sort | uniq > $basedir/${i}_working_${suffix}/${i}_private_singletons.txt && rm -rf $basedir/${i}_working_${suffix}/tmp
fi
