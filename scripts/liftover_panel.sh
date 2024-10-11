#!/usr/bin/env bash
set -o xtrace
set -x

in=$1
out=$2
target_fasta=$3
chain=$4
vcf_of_differences=$5


## this code strips out all forms of vcf/bcf/vcf.gz/bcf.gz suffix
## strip a gz suffix if present, then strip everything after the last period.
root_name=$(dirname $out)/$(basename ${out%.gz})
root_name=${root_name%.*}
tmpfile=$(dirname $out)/temp.$(basename $root_name).bcf

echo "liftover_indels.py $in $vcf_of_differences $tmpfile $chain $target_fasta"

python3 liftover_indels.py $in $vcf_of_differences $tmpfile $chain $target_fasta \
&& bcftools norm -Ou --threads 4 -f $target_fasta -m -any $tmpfile \
| bcftools annotate -Ou --threads 4 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
| bcftools +fill-tags -Ou --threads 2 - -- -t "AN,AC,MAF,MAC:1=MAC" \
| bcftools sort -Ou -m 8G -T /tmp/ - \
| bcftools view -Ob --threads 8 -l1 - > $out \
&& bcftools index -f --threads 2 -f $out \
&& rm $tmpfile

wait

echo "Lifted $(basename $in) to $(basename $out)"
