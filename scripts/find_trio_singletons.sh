#!/usr/bin/env bash
set -eo pipefail

# Default values
genome=""
pedigree=""
infile=""
outdir="intermediate_data/variant_frequency_stats"
n_singleton_jobs=2

while getopts "c:g:p:i:o:j:" opt; do
    case $opt in
        c) chrom="$OPTARG" ;;
        g) genome="$OPTARG" ;;
        p) pedigree="$OPTARG" ;;
        i) infile="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        j) n_singleton_jobs="$OPTARG" ;;
        *) echo "Usage: $0 -g genome -p pedigree -i infile -o outdir [-j n_jobs]" >&2; exit 1 ;;
    esac
done

# Check required
if [ -z "$genome" ] || [ -z "$pedigree" ] || [ -z "$infile" ] || [ -z "$outdir" ]; then
    echo "Missing required options: -g, -p, -i, -o" >&2
    exit 1
fi

if [[ $genome == 'CHM13v2.0' ]]; then
    short_genome='t2t'
else
    short_genome='grch38'
fi

basedir=$PWD/..

mkdir -p $outdir/${genome}

if [[ ! -s $outdir/${genome}/${chrom}_private_singletons.txt ]]; then
    mkdir -p $outdir/../singleton_tmp/${genome}
    # /${short_genome}/1KGP.${genome}.${chrom}.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf
    cat $pedigree \
        | parallel -j $n_singleton_jobs "bcftools view --force-samples -H -G -s {} -x -c 2 \
        $infile | \
        cut -f 3 > $outdir/../singleton_tmp/${genome}/{#}.txt" \
    && \
    cat $outdir/../singleton_tmp/${genome}/*.txt | sort | uniq \
        > $outdir/${genome}/${chrom}_private_singletons.txt \
    && \
    rm -rf $outdir/../singleton_tmp/${genome}
fi

