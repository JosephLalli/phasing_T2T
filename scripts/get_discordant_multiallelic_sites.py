import sys
import pandas as pd
from pysam import VariantFile

native=sys.argv[1]
lifted=sys.argv[2]
out=sys.argv[3]

print(native)
print(lifted)
print(out)
if out is None:
    out='native_lifted_same_calls.IDs.txt'

native = VariantFile(native)
lifted = VariantFile(lifted)

def iter_vars(vf):
    for var in vf:
        assert len(var.alts) == 1
        yield (var.chrom, var.pos, var.ref, var.alts[0], var.id, var.info['MAF'])

print("loading variants")
native = pd.DataFrame(iter_vars(native), columns=['chrom','pos','ref','alt','ID','MAF']).sort_values(['chrom','pos','ref','alt'])
lifted = pd.DataFrame(iter_vars(lifted), columns=['chrom','pos','ref','alt','ID','MAF']).sort_values(['chrom','pos','ref','alt'])

native.ref += ','
lifted.ref += ','
native.alt += ','
lifted.alt += ','

native['alt'] = native.ref + native.alt
lifted['alt'] = lifted.ref + lifted.alt

print ('finding multiallelic variants')
native = native[['chrom','pos','ID','MAF']].merge(native[['chrom','pos','alt']]
                                                    .groupby(['chrom','pos'])
                                                    .sum()
                                                    .alt.str[:-1].str.split(',')
                                                    .apply(lambda x: tuple(sorted(list(set(x)))))
                                                    .reset_index(),
                                                on=['chrom','pos'], how='outer')

lifted = lifted[['chrom','pos','ID','MAF']].merge(lifted[['chrom','pos','alt']]
                                                    .groupby(['chrom','pos'])
                                                    .sum()
                                                    .alt.str[:-1].str.split(',')
                                                    .apply(lambda x: tuple(sorted(list(set(x)))))
                                                    .reset_index(),
                                                on=['chrom','pos'], how='outer')

print('identifying variants identical between the two files')
vars = native.merge(lifted, on=['chrom','pos','ID'], how='outer', suffixes=('_native','_lifted'))

samesies = vars.loc[vars.alt_native == vars.alt_lifted]
samesies[['ID','MAF_native','MAF_lifted','alt_native','alt_lifted']].to_csv(out, index=False, header=False, sep='\t')
print(f'variant IDs exported to {out}')

bad_out='.'.join(out.split('.')[:-1])+'mismatched_num_alleles.txt'
vars.loc[vars.alt_native != vars.alt_lifted, ['ID','MAF_native','MAF_lifted','alt_native','alt_lifted']].dropna().to_csv(bad_out, index=False, header=True, sep='\t')
