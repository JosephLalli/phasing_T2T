#!/usr/bin/env bash
set -o xtrace
set -x

# Simple CLI: support flags and backward-compatible positional args
usage() {
        cat <<'USAGE'
Usage: $(basename "$0") -i <input.vcf/BCF> -o <output.bcf> -t <target_fasta> -c <chain> -s <src_fasta> [-r <region>] [-h]

Options:
    -i INPUT       input VCF/BCF to liftover (required)
    -o OUTPUT      output BCF path (required)
    -t TARGET      target fasta (required)
    -c CHAIN       chain file for liftover (required)
    -s SRC_FASTA   source fasta used by liftover (required)
    -r REGION      optional region to limit processing (e.g. chr1:100000-200000)
    -h             show this help and exit
USAGE
}

region=""
while getopts ":hi:o:t:c:s:r:" opt; do
    case ${opt} in
        h )
            usage
            exit 0
            ;;
        i ) in=${OPTARG} ;;
        o ) out=${OPTARG} ;;
        t ) target_fasta=${OPTARG} ;;
        c ) chain=${OPTARG} ;;
        s ) in_fasta=${OPTARG} ;;
        r ) region=${OPTARG} ;;
        \: ) echo "Missing argument for -${OPTARG}" >&2; usage; exit 2 ;;
        \? ) echo "Unknown option: -${OPTARG}" >&2; usage; exit 2 ;;
    esac
done

shift $((OPTIND -1))

# Backwards-compatible positional fallback
if [ -z "${in:-}" ] && [ "$#" -ge 1 ]; then
    in=$1
fi
if [ -z "${out:-}" ] && [ "$#" -ge 2 ]; then
    out=$2
fi
if [ -z "${target_fasta:-}" ] && [ "$#" -ge 3 ]; then
    target_fasta=$3
fi
if [ -z "${chain:-}" ] && [ "$#" -ge 4 ]; then
    chain=$4
fi
if [ -z "${in_fasta:-}" ] && [ "$#" -ge 6 ]; then
    in_fasta=$6
fi

# Validate required args (will exit with message if missing)
: ${in:?"Missing input (--input or -i). See -h for help."}
: ${out:?"Missing output (--output or -o). See -h for help."}
: ${target_fasta:?"Missing target fasta (--target or -t). See -h for help."}
: ${chain:?"Missing chain file (--chain or -c). See -h for help."}
: ${in_fasta:?"Missing source fasta (--src-fasta or -s). See -h for help."}


## this code strips out all forms of vcf/bcf/vcf.gz/bcf.gz suffix
## strip a gz suffix if present, then strip everything after the last period.
root_name=$(dirname $out)/$(basename ${out%.gz})
root_name=${root_name%.*}
# tmpfile=$(dirname $out)/temp.$(basename $root_name).bcf
# tmpfolder=$(dirname $out)/tmp
# mkdir -p $tmpfolder

# python3 liftover_indels.py $in $vcf_of_differences $tmpfile $chain $target_fasta \
bcftools +liftover -Ou $in --regions $region -- --chain $chain --src-fasta-ref $in_fasta --fasta-ref $target_fasta \
                    --write-src --write-fail --fix-tags \
                    --reject $root_name.unlifted.bcf -Ob \
| bcftools norm -Ou --threads 4 -f $target_fasta -m -any - \
| bcftools annotate -Ou --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
| bcftools +fill-tags -Ou --threads 2 - -- -t "AN,AC,MAF,MAC:1=MAC" \
| bcftools sort -Ou -m 8G -T $tmpfolder - \
| bcftools view -Ob --threads 8 -l1 - > $out \
&& bcftools index --threads 2 -f $out # \
# && rm $tmpfile

wait

echo "Lifted $(basename $in) to $(basename $out)"
