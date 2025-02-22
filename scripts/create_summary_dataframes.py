import polars as pl
import pandas as pd
import numpy as np
import os
import glob
from copy import copy
from datetime import datetime

StartTime = datetime.now() 
def print_runtime():
    time = datetime.now() - StartTime
    s = time.seconds % 60
    m = (time.seconds - s)//60
    return f"{m}m {s}s"

run_version='092424'

summary_statistics_folder = f'global_phasing_stats_summary_{run_version}'
if not os.path.exists(summary_statistics_folder):
    os.mkdir(summary_statistics_folder)


contigs = ['chr' + str(x) for x in list(range(1,23))] + ['PAR1', 'chrX','PAR2']
t2t_run_suffix = f'T2T_scaled_newimpute_{run_version}'
grch38_run_suffix = f'GRCh38_statistics_{run_version}'
t2t_phasing_stats_folder = f'phasing_stats_{t2t_run_suffix}'
grch_phasing_stats_folder = f'phasing_stats_GRCh38_{grch38_run_suffix}'

ground_truths = ['trios','HPRC_samples','HGSVC3']
methods = ['phased_with_parents_and_pedigree', 'phased_without_parents_or_pedigree', 'phased_with_parents_pedigree_children_excluded_during_measurement', 'HPRC_variation_phased_with_reference_panel']

common_sample_columns = ['sample_id','chrom', 'genome']
common_variant_columns = ['variant_id','position','chrom', 'genome']
expt_columns = ['method_of_phasing','ground_truth_data_source']

sample_metadata = (pl.read_csv('sample_subsets/superpopulations.samples.txt', separator=' ', infer_schema_length=1000,
                                                                             new_columns=['trio_id', 'sample_id','father','mother','sex','population','superpopulation'])
                     .with_columns(sex = pl.col('sex').cast(pl.String()).replace('1','Male').replace('2','Female')))
pangenome_only_samples = ['HG01123','HG02109','HG02486','HG02559','NA21309']

catChrom = pl.Enum(['chr' + str(x) for x in list(range(1,23))] + ['chrX','chrY','PAR1', 'PAR2'])
catMethod = pl.Enum(methods)
catGroundTruth = pl.Enum(ground_truths)
catGenomes = pl.Enum(['GRCh38','CHM13v2.0'])
catSamples = pl.Enum(sorted(list(set(pl.Series(sample_metadata.select('sample_id')).to_list() + pangenome_only_samples))))
catSyn = pl.Enum(['Syntenic','Nonsyntenic'])
catVarType = pl.Enum(['SNP','Indel','SNPs + Indels'])

sample_metadata = sample_metadata.with_columns(sample_id = pl.col('sample_id').cast(catSamples))

all_variant_annotations_dtypes = pl.Schema( {'ID':pl.String(),
                                        'CHROM': pl.String(),
                                        'POS': pl.UInt32(),
                                        'ALT':pl.String(),
                                        'REF':pl.String(),
                                        'NEGATIVE_TRAIN_SITE': pl.Boolean(),
                                        'QUAL': pl.Float32(),
                                        'VQSLOD': pl.Float32(),
                                        'MERR': pl.UInt16(),
                                        'HWE': pl.Float32(),
                                        'HWE_EUR': pl.Float32(),
                                        'HWE_AFR': pl.Float32(),
                                        'HWE_EAS': pl.Float32(),
                                        'HWE_AMR': pl.Float32(),
                                        'HWE_SAS': pl.Float32(),
                                        'FILTER': pl.String(),    # good for categorical if all known in advance
                                        'culprit': pl.String(),   # good for categorical if all known in advance
                                        'InbreedingCoeff': pl.Float32(),
                                        'ExcessHet': pl.Float32(),
                                        'MAF': pl.Float32(),
                                        'AN': pl.UInt16(),
                                        'AC': pl.UInt16(),
                                        'F_MISSING': pl.Float32(),
                                        'SYNTENIC': pl.Boolean()})

headers = {'flips_and_switches'           :['sample_id','n_pure_switch_errors','n_flips','num_correct_phased_hets','pure_switch_error_rate','flip_error_rate','accurately_phased_rate'],
           'switch_errors'                :['sample_id', 'n_switch_errors',  'n_checked', 'switch_error_rate'],
           'frequency'                    :['minor_allele_count', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'type'                         :['is_SNP', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'variant_genotype_concordance' :['variant_id', 'position', 'n_gt_errors', 'n_gt_checked', 'gt_error_rate'],
           'block_length'                 :['sample_id','switch_error_site'],
           'variant_switch'               :['variant_id', 'position', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'sample_genotype_concordance'  :['sample_id','n_gt_errors','n_gt_checked','gt_error_rate']}

schemas = {'flips_and_switches'           :{'sample_id':pl.String(), 'n_pure_switch_errors':pl.UInt32(),'n_flips':pl.UInt32(),'num_correct_phased_hets':pl.UInt32(),'pure_switch_error_rate':pl.Float32(),'flip_error_rate':pl.Float32(),'accurately_phased_rate':pl.Float32()},
           'switch_errors'                :{'sample_id':pl.String(), 'n_switch_errors':pl.UInt32(),  'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'frequency'                    :{'minor_allele_count':pl.UInt16(), 'n_switch_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'type'                         :{'is_SNP':pl.Boolean(), 'n_gt_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'gt_error_rate':pl.UInt32()},
           'variant_genotype_concordance' :{'variant_id':pl.String(), 'position':pl.UInt32(), 'n_gt_errors':pl.UInt32(), 'n_gt_checked':pl.UInt32(), 'gt_error_rate':pl.Float32()},
           'block_length'                 :{'sample_id':pl.String(),'switch_error_site':pl.UInt32()},
           'variant_switch'               :{'variant_id':pl.String(), 'position':pl.UInt32(), 'n_switch_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'sample_genotype_concordance'  :{'sample_id':pl.String(), 'n_gt_errors':pl.UInt32(),  'n_gt_checked':pl.UInt32(), 'gt_error_rate':pl.Float32()}}

analysis_suffixes = {'trios':       {'phased_without_parents_or_pedigree'        :'rare_noparents_vs_trios',
                                     'phased_with_parents_and_pedigree'          :'3202_panel_vs_trios',
                                     'phased_with_parents_pedigree_children_excluded_during_measurement':'2504_panel_vs_trios'}, #For GRCh38, this is just trio validation of the reference phased set.
                     'HPRC_samples':{'phased_without_parents_or_pedigree'        :'noparents_vs_HPRC',
                                     'phased_with_parents_and_pedigree'          :'3202_panel_vs_HPRC',
                                     'HPRC_variation_phased_with_reference_panel':'rare_pangenome_panelphased_vs_pangenome'}}

stat_prefixes =     {'flips_and_switches':'flipsAndSwitches',
                     'switch_errors':'sample.switch',
                     'frequency':'frequency',
                     'type':'type',
                     'variant_switch':'variant.switch',
                     'sample_mendel':'sample.mendel',
                     'variant_mendel':'variant.mendel',
                     'calibration': 'calibration',
                     'sample_genotype_concordance': 'sample.typing',
                     'variant_genotype_concordance':'variant.typing',
                     'block_length':'block.switch'}

r2_bin_names = ('singleton',
                '0.00021-0.00042',
                '0.00042-0.00064',
                '0.00064-0.001',
                '0.001-0.0016',
                '0.0016-0.0022',
                '0.0022-0.003',
                '0.003-0.004',
                '0.004-0.0054',
                '0.0054-0.0072',
                '0.0072-0.0094',
                '0.0094-0.0126',
                '0.0126-0.0172',
                '0.0172-0.0244',
                '0.0244-0.037',
                '0.037-0.06',
                '0.06-0.102',
                '0.102-0.166',
                '0.166-0.255',
                '0.255-0.372',
                '0.372-0.5')
r2_bins='0 0.00021 0.00042 0.00064 0.001 0.0016 0.0022 0.003 0.004 0.0054 0.0072 0.0094 0.0126 0.0172 0.0244 0.0369 0.0601 0.1018 0.1661 0.2556 0.3724 0.5'
r2_bins = tuple([float(x) for x in r2_bins.split(' ')])

out_file=list()

def group_phasing_statistics(df, columns_to_groupby):
    df_cols = df.columns
    assert all([c in df_cols for c in columns_to_groupby])

    if 'N50' in df_cols:
        means = df.group_by(columns_to_groupby).mean()
        N50_L50 = means.select(columns_to_groupby+['N50','L50'])
        df = df.group_by(columns_to_groupby).sum()
        df = df.drop(['N50','L50']).join(N50_L50, on=columns_to_groupby)
    else:
        df = df.group_by(columns_to_groupby).sum()

    expressions = dict()
    if 'n_gt_errors' in df_cols:
        expressions.update({'gt_error_rate':(pl.col('n_gt_errors')/pl.col('n_gt_checked')) * 100})
    if 'n_mendel_errors' in df_cols:
        expressions.update({'mendel_error_rate':(pl.col('n_mendel_errors')/pl.col('n_variant_calls')) * 100})
    if 'n_pure_switch_errors' in df_cols:
        # See definitions below.
        expressions.update({'flip_error_rate':(pl.col('n_flips')/(pl.col('n_checked')-1))*100,
                            'pure_switch_error_rate':(pl.col('n_pure_switch_errors')/pl.col('n_checked')) * 100,
                            'accurately_phased_rate':(pl.col('num_correct_phased_hets')/pl.col('n_checked'))*100})
    if 'n_switch_errors' in df_cols:
        expressions.update({'switch_error_rate':(pl.col('n_switch_errors')/pl.col('n_checked'))*100})

    return df.with_columns(**expressions)

print(f'Gathering alt frequencies for each phased panel, dropping trio-private singletons ({print_runtime()})')
for genome, run_suffix in [('GRCh38', grch38_run_suffix),#]:
                           ('CHM13v2.0', t2t_run_suffix)]:
    for contig in contigs:
        chrom_working_dir = f"{contig}_working_{run_suffix}"
        with open (f'{chrom_working_dir}/{contig}_private_singletons.txt', 'r') as f:
            private_variants = {x.strip() for x in f}
        an_acs = ((['phased_with_parents_and_pedigree','phased_with_parents_pedigree_children_excluded_during_measurement'], f"{chrom_working_dir}/{contig}_3202_AC_AN.tsv"),
                  (['HPRC_variation_phased_with_reference_panel'], f"{chrom_working_dir}/{contig}_2430_AC_AN.tsv"),
                  (['phased_without_parents_or_pedigree'], f"{chrom_working_dir}/{contig}_2002_AC_AN.tsv"))
        for methods, tsv in an_acs:
            # if os.path.exists(tsv):
            try:
                anac = pl.scan_csv(tsv, separator='\t', has_header=False, new_columns=['variant_id','contextual_MAC','contextual_AN'],
                                   schema_overrides={'variant_id':pl.String(),
                                                     'contextual_MAC':pl.UInt16(),
                                                     'contextual_AN':pl.UInt16()})
            except:
                if ('2002' not in tsv) and (genome!='GRCh38'):
                    print(f'{tsv} missing')
                continue
            for method in methods:
                if method in ('HPRC_variation_phased_with_reference_panel', 'phased_without_parents_or_pedigree'):
                    out_file.append(anac.filter(~pl.col('variant_id').is_in(private_variants))
                                        .with_columns(method_of_phasing=pl.lit(method).cast(catMethod),
                                                      genome=pl.lit(genome).cast(catGenomes)))
                else:
                    out_file.append(anac.with_columns(method_of_phasing=pl.lit(method).cast(catMethod),
                                                      genome=pl.lit(genome).cast(catGenomes)))

print (f'Collecting each reference panel\'s minor allele frequency statistics ({print_runtime()})')
contextual_variants = pl.concat(out_file)#.collect()
# contextual_variants.write_parquet(f'{summary_statistics_folder}/contextual_variants.parquet')


# Gather per variant statistics calculated by bcftools query
variant_data = list()
for chrom in contigs:
    for genome, phasing_stats_folder, workdir_suffix in [('GRCh38',grch_phasing_stats_folder, grch38_run_suffix),
                                                         ('CHM13v2.0', t2t_phasing_stats_folder, t2t_run_suffix)]:
        full_panel_variant_data = f'{chrom}_working_{workdir_suffix}/1KGP.{genome}.{chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv'
        if os.path.exists(full_panel_variant_data):
            contig_variant_data = pl.scan_csv(full_panel_variant_data, separator='\t', null_values=['.','NA'],
                                                schema_overrides = all_variant_annotations_dtypes
                                                ).with_columns(genome=pl.lit(genome).cast(catGenomes),
                                                            CHROM=pl.col('CHROM').cast(catChrom),
                                                            SYNTENIC=pl.col('SYNTENIC').fill_null(False).cast(pl.Boolean()),
                                                            NEGATIVE_TRAIN_SITE=pl.col('NEGATIVE_TRAIN_SITE').fill_null(False).cast(pl.Boolean()),
                                                ).rename({'AN':'AN_original_panel',
                                                          'AC':'AC_original_panel',
                                                          'MAF':'MAF_original_panel'})
        else:
            continue
        variant_data.append(contig_variant_data)

print (f"Collecting variant statistics for all variants in study ({print_runtime()})")
variant_data = pl.concat(variant_data).unique(subset=['ID','SYNTENIC','genome'], keep='first')#.collect()
# variant_data = variant_data.filter( # Due to what appears to be a quirk of how multiallelic SVs were called and represented in the original GRCh38 dataset,
#                                     # on rare occasions (<30 in the whole callset) splitting, left-aligning, and normalizing indels results in
#                                     # alleles at different sites (see chr21:45918284 and chr21:45918304 as an example) 
#                                     # are left-aligned to the same alt allele and same site. That results in duplicate variant IDs.
#                                     # I do not know how to merge these duplicates into one variant using bcftools, 
#                                     # since the allele's INFO data are so different.
#                                     # Thankfully these sites are all SVs, which are beyond the scope of this paper. 
#                                     # Therefore I will just drop any duplicate variant IDs after collecting all data.write_parquet(f'{summary_statistics_folder}/all_variant_data.parquet')
#                                     ~variant_data.select(['ID','SYNTENIC','genome']).is_duplicated())
# print (f"Writing variant statistics to disk ({print_runtime()})")
# variant_data.write_parquet(f'{summary_statistics_folder}/all_variant_data.parquet')


print (f'loading SHAPEIT5_switch statistics! ({print_runtime()})')
total_phasing_stats_datasets = {stat: list() for stat in schemas.keys()}
del total_phasing_stats_datasets['variant_genotype_concordance'] #merged with variant_switch and gathered as one per-variant df
for contig in contigs:
    per_contig_phasing_stats_datasets = {stat: list() for stat in schemas.keys()}
    for statistic, schema in schemas.items():
        for genome, phasing_stats_folder in [('GRCh38',grch_phasing_stats_folder),('CHM13v2.0', t2t_phasing_stats_folder)]:
            for verification_sample in analysis_suffixes.keys():
                for phasing_method, prefix in analysis_suffixes[verification_sample].items():
                    if genome == 'GRCh38' and 'noparents' in prefix:
                        continue
                    try:
                        file = glob.glob(phasing_stats_folder + f'/{prefix}_{contig}.*{stat_prefixes[statistic]}*')[0]
                    except IndexError:
                        print (phasing_stats_folder + f'/{prefix}_{contig}.*{stat_prefixes[statistic]}*')
                        raise
                    ### note: there are some hacks involved here to allow for setting column names and dtypes
                    ### while scanning .gz files that use non-standard separator values
                    contig_file = pl.scan_csv(file, separator=' ', 
                                                null_values=['.','NA'], has_header=False
                                             ).cast({f'column_{i+1}':x for i,x in enumerate(schema.values())}
                                             ).rename({f'column_{i+1}':x for i,x in enumerate(schema.keys())}
                                             ).with_columns(
                                                method_of_phasing=pl.lit(phasing_method).cast(catMethod),
                                                ground_truth_data_source=pl.lit(verification_sample).cast(catGroundTruth),
                                                chrom=pl.lit(contig).cast(catChrom),
                                                genome=pl.lit(genome).cast(catGenomes))
                    if 'sample_id' in schema.keys():
                        contig_file = contig_file.with_columns(sample_id=pl.col('sample_id').cast(catSamples))
                    if 'variant_id' in schema.keys():
                        if "PAR" in contig:
                            position_chr = 'chrX'
                        else:
                            position_chr = contig
                        contig_file = contig_file.with_columns(variant_id = 'chr' + pl.col('variant_id').str.slice(3).str.to_uppercase(),
                                                               chrom=pl.lit(position_chr).cast(catChrom))
                    if statistic == 'variant_switch':
                        contig_file = contig_file.with_columns(ref_len=pl.col("variant_id").str.split(by="_").list.get(2).str.len_chars().cast(pl.UInt16,strict=False).fill_null(255),
                                                               alt_len=pl.col("variant_id").str.split(by="_").list.get(3).str.len_chars().cast(pl.UInt16,strict=False).fill_null(255)
                                                ).with_columns(type = pl.when(pl.col("ref_len") == 1,
                                                                              pl.col("alt_len") == 1)
                                                                        .then(pl.lit('SNP'))
                                                                        .otherwise(pl.lit('Indel')).cast(catVarType)
                                                 ).filter(pl.col('ref_len') < 50, pl.col('alt_len') < 50
                                                 ).drop(['switch_error_rate','ref_len','alt_len'])
                    elif statistic == 'sample_genotype_concordance':
                        contig_file = contig_file.with_columns(n_gt_errors = pl.when(pl.col('n_gt_errors')==0)
                                                                               .then(pl.lit(float("nan")))
                                                                               .otherwise(pl.col('n_gt_errors')))
                    per_contig_phasing_stats_datasets[statistic].append(contig_file)
        per_contig_phasing_stats_datasets[statistic] = pl.concat(per_contig_phasing_stats_datasets[statistic])
        if statistic not in ('variant_switch','variant_genotype_concordance'):
            total_phasing_stats_datasets[statistic].append(per_contig_phasing_stats_datasets[statistic])
    per_variant_data = (per_contig_phasing_stats_datasets['variant_switch'].join(per_contig_phasing_stats_datasets['variant_genotype_concordance'], 
                                                                                                                    on=['variant_id','position','chrom','ground_truth_data_source','genome','method_of_phasing'], 
                                                                                                                    how='left', validate='1:1').drop(['gt_error_rate']))
    total_phasing_stats_datasets['variant_switch'].append(per_variant_data)

for stat, datasets in total_phasing_stats_datasets.items():
    print (f'Collecting statistics on {stat}... ({print_runtime()})')
    dataset = pl.concat(datasets)
    if stat != 'variant_switch':
        total_phasing_stats_datasets[stat] = dataset.collect()
        total_phasing_stats_datasets[stat].write_parquet(f'{summary_statistics_folder}/{stat}.parquet')
    else:
        total_phasing_stats_datasets[stat] = (dataset.join(contextual_variants, on=['variant_id', 'method_of_phasing', 'genome'], how='left', validate='m:1')
                                                     .join(variant_data.select(['ID','genome','SYNTENIC'])
                                                                       .rename({'ID':'variant_id','SYNTENIC':'syntenic'}), 
                                                           on=['variant_id','genome'], how='left', validate='m:1'))
        total_phasing_stats_datasets[stat] = total_phasing_stats_datasets[stat].collect(streaming=True)
        total_phasing_stats_datasets[stat].write_parquet(f'{summary_statistics_folder}/variants.parquet')

print (f'Writing a second, abbreviated per-variant switch statistics to disk ({print_runtime()})')
total_phasing_stats_datasets['variant_switch'].filter(pl.col('contextual_AN') > 0).write_parquet(f'{summary_statistics_folder}/variants.nozeros.parquet')


print (f'Binning variants ({print_runtime()})')
MAF_performance_variant_df = total_phasing_stats_datasets['variant_switch'].select(['variant_id','n_switch_errors','n_checked','n_gt_errors', 'n_gt_checked','method_of_phasing','ground_truth_data_source','genome', 'contextual_MAC', 'contextual_AN', 'syntenic', 'type']).rename({'contextual_MAC':'MAC', 'contextual_AN':'AN'})


# num with MAC==0: 
# │ HPRC_samples             ┆ HPRC_variation_phased_with_reference_panel ┆ GRCh38    ┆ 31882   │
# │ trios                    ┆ phased_without_parents_or_pedigree         ┆ CHM13v2.0 ┆ 3734053 │ <- case number one
# num with MAC is null:
# │ HPRC_samples             ┆ HPRC_variation_phased_with_reference_panel ┆ GRCh38    ┆ 828045   │
# │ trios                    ┆ phased_without_parents_or_pedigree         ┆ CHM13v2.0 ┆ 10894327 │
# │ HPRC_samples             ┆ HPRC_variation_phased_with_reference_panel ┆ CHM13v2.0 ┆ 711405   │
# │ HPRC_samples             ┆ phased_without_parents_or_pedigree         ┆ CHM13v2.0 ┆ 712131   │

# cases where contextual_MAC == 0:
# 1) Site is present in verification dataset (full 3202 panel when trio phasing) but not in analysis dataset (eg, using 3202 panel to evaluation no-reference phasing of 2003 no-parents panel)
    # -drop these. They have an undefined switch error rate an a 100% GT accuracy rate - always accurately "called" as reference. Uninformative.
# 2) Site is only in the 5 samples that are in the HPRC variation but not in the 1KGP
    # - not sure what to do. The 'contextual AC/AN' is going to be the number of alleles/haplotypes available during phasing.
    # I think we drop them. The point of using MAF bins is to allow comparison between wildly different situations -
    #     eg, rephasing 44 samples, natively phasing 3202 samples, natively phasing 2002 samples. 
    #     if an individual vairant does not have a reference MAF, then it cannot be compared to other situations.
    #     we are talking about a small number of inherently rare variants. Let it go.
# cases where contextual_mac==null:
# 1) Site is in private trio, dropped to avoid falsely high singleton rate
    # - drop these, explain in methods
    # - yes, even in HPRC_samples. many of these seem to be private trios in the referenece panel, but some are common in the HPRC dataset.
    #      - this is due to the panel being miscalled.

MAF_performance_variant_df = (MAF_performance_variant_df.filter(~pl.col('MAC').is_null())
                                                        .filter(~(pl.col('MAC') == 0), ~(pl.col('AN')==0))
                                                        .with_columns(MAF = ((pl.col('MAC') / pl.col('AN'))).cast(pl.Float32())))
MAF_performance_variant_df = MAF_performance_variant_df.with_columns(rounded_MAF = pl.col('MAF').cut(r2_bins[1:-1], labels=r2_bin_names))
MAF_performance_variant_df = MAF_performance_variant_df.with_columns(MAF = pl.col('MAF') * 100)

MAF_performance_variant_df.write_parquet(f'{summary_statistics_folder}/MAF_performance_variants.parquet')

MAF_performance_variant_df = MAF_performance_variant_df.drop(['variant_id','MAF'])

print(f'Calculating MAF binned statistics ({print_runtime()})')
fig_ab_data = (MAF_performance_variant_df.group_by(['type','rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source'])
                                         .sum()
                                         .with_columns(syntenic=pl.lit('All'),
                                                       switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                       gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                       MAF=(pl.col('MAC')/pl.col('AN')) * 100))
fig_ab_data = pl.concat([fig_ab_data, MAF_performance_variant_df.group_by(['rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source'])
                                         .sum()
                                         .with_columns(syntenic=pl.lit('All'),
                                                       switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                       gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                       MAF=((pl.col('MAC')/pl.col('AN')) * 100),
                                                       type=pl.lit('SNPs + Indels').cast(catVarType))
                                         .select(fig_ab_data.columns)])
fig_cd_data = (MAF_performance_variant_df.group_by(['type','rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source', 'syntenic'])
                                         .sum()
                                         .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                       gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                       MAF=(pl.col('MAC')/pl.col('AN')) * 100))
fig_cd_data = pl.concat([fig_cd_data, MAF_performance_variant_df.group_by(['rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source', 'syntenic'])
                                         .sum()
                                         .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                       gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                       MAF=((pl.col('MAC')/pl.col('AN')) * 100),
                                                       type=pl.lit('SNPs + Indels').cast(catVarType))
                                         .select(fig_cd_data.columns)])

fig_ab_data.write_parquet(f'{summary_statistics_folder}/fig_ab_data.parquet')
fig_cd_data.write_parquet(f'{summary_statistics_folder}/fig_cd_data.parquet')

### Annotate variants with cytoband region
chm13_fai = 'chm13v2.0.fa.gz.fai'
grch38_fai = 'GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai'
chm13_cytobands = 'chm13v2.0_cytobands_allchrs.bed'
grch38_cytobands = 'grch38_cytobands_allchrs.bed'
fais={'CHM13v2.0':chm13_fai,
      'GRCh38':grch38_fai}

def create_actual_ends(df, fais=fais):
    ends = {k:read_fai(fai) for k, fai in fais.items()}
    out = list()
    # note: ends are 1-based, inclusive. so end - 1 (0 base) + 1 (exclusive) = unchanged value.
    for (genome, contig), d in df.groupby(['genome','chrom']):
        end = ends[genome][contig]
        d.loc[d.end==d.end.iloc[-1], 'end'] = end
        out.append(d)
    return pd.concat(out).reset_index(drop=True)

def create_end_intervals(fais=fais, window_size=int(1e6)):
    ends = {k:read_fai(fai) for k, fai in fais.items()}
    for genome, contigs in ends.items():
        for contig, end in contigs.items():
            ends[genome][contig] = pd.Interval(end-end%window_size, end)
    return ends

def read_fai(fai):
    with open(fai, 'r') as f:
        fai = [x.strip().split('\t') for x in f.readlines()]
    ends = {x[0]:int(x[1]) for x in fai}
    return ends

def add_start_and_end_columns_to_intervals(intervals, fais):
    intervals = intervals.reset_index()
    intervals['start'] = intervals.intervals.apply(lambda x: x.left).astype(int)
    intervals['end'] = intervals.intervals.apply(lambda x: x.right).astype(int)
    df = create_actual_ends(intervals, fais)
    return df.reset_index()

def get_fai_in_bed_format(fai):
    fai = pd.read_csv(fai, sep='\t', header=None, names=['chrom','end', 'cum','linewidtha','linewidthb'])
    fai['start'] = 0
    return fai[['chrom','start','end']]

def get_cytobands_in_bed_format(cytobands):
    return pd.read_csv(cytobands, sep='\t', header=None, names=['chrom','start','end','name','color']).drop(columns='color').dropna()

def make_per_window_interval_df(variants, fais, window_size):
    intervals = pd.cut(variants.set_index(['genome','chrom']).position, range(0, variants.position.max()+window_size, window_size))
    intervals.name='intervals'
    intervals = add_start_and_end_columns_to_intervals(intervals, fais)
    return intervals

def make_per_section_interval_df(variants, intervals):
    for col in ['chrom','start','end','name','genome']:
        assert col in intervals.columns
    intervals = pl.from_pandas(intervals).with_columns(chrom=pl.col('chrom').cast(variants.schema['chrom']),
                                                       genome=pl.col('genome').cast(variants.schema['genome']),
                                                       start=pl.col('start').cast(pl.UInt32()),
                                                       end=pl.col('end').cast(pl.UInt32())
                                         ).dropna()
    if 'name' in intervals.columns:
        intervals = intervals.with_columns(name=pl.col('name').cast(pl.Categorical()))


    per_variant_interval_dfs = list()
    for (genome, chrom), df in variants[['genome','chrom','position']].groupby(['genome','chrom']):
        print ((genome, chrom))
        contig_intervals = intervals.loc[(intervals.genome==genome) & (intervals.chrom==chrom)]
        contig_interval_names = contig_intervals.name.to_list()
        contig_intervals_list = [0] + contig_intervals.end.to_list()
        per_variant_interval = pd.cut(df.position, contig_intervals_list).drop(columns='position')
        per_variant_interval_names = pd.cut(df.position, contig_intervals_list, labels=contig_interval_names).drop(columns='position')
        df = df.assign(intervals=per_variant_interval.astype('object'), 
                        interval_names=per_variant_interval_names).drop(columns='position')
        per_variant_interval_dfs.append(df.merge(contig_intervals[['start','end','name']].rename(columns={'name':'interval_names'}), on='interval_names').reset_index(drop=True))
    intervals = pd.concat(per_variant_interval_dfs).reset_index(drop=True)
    return intervals

def variants_to_compressed_ideogram_data(variants: pl.Dataframe, fais=fais, window_size=int(1e6), intervals=None: pd.DataFrame):
    assert ('position' in variants.columns)
    print ('making intervals')
    variants = variants.sort(['genome','chrom','position','method_of_phasing','ground_truth_data_source'])
    if intervals is None:
        intervals = make_per_window_interval_df(variants, fais, window_size)
    else:
        intervals1 = make_per_section_interval_df(variants, intervals)
    print ('summing across intervals')
    compressed_ideogram = (variants.drop('position')
                                    .group_by(['method_of_phasing','ground_truth_data_source',
                                                'genome', 'chrom', intervals1.intervals, intervals1.interval_names])
                                    .sum())

    print ('dropping interval duplicates')
    intervals2 = intervals.drop_duplicates().reset_index(drop=True)

    print ('merging interval subs with interval metadata')
    compressed_ideogram = compressed_ideogram.merge(intervals[['genome','chrom','interval_names','start','end']], on=['genome','chrom','interval_names'], how='left')
    print('prepping interval sums for export')
    compressed_ideogram = (compressed_ideogram.drop(columns='intervals')
                                                .set_index(['genome','method_of_phasing','ground_truth_data_source','chrom']))
    return compressed_ideogram


cytoband_coords = get_cytobands_in_bed_format(chm13_cytobands)
cytoband_coords['genome'] = 'CHM13v2.0'
cytoband_coords2 = get_cytobands_in_bed_format(grch38_cytobands)
cytoband_coords2['genome'] = 'GRCh38'
cytoband_coords = pd.concat((cytoband_coords, cytoband_coords2))
all_cytobands = variants_to_compressed_ideogram_data(var_nodup, intervals=cytoband_coords)
all_cytobands.to_parquet('all_cytobands.parquet')


### gather chrom and sample specific data
del MAF_performance_variant_df
del total_phasing_stats_datasets['variant_switch']

print (f'Calculating per contig and per sample phasing statistics ({print_runtime()})')
# for samples that weren't analyzed, gt errors are reported as nan and gt_checked is reported as total gts. I want the reverse- 0 errors, 0 checked.
total_phasing_stats_datasets['sample_genotype_concordance'] = total_phasing_stats_datasets['sample_genotype_concordance'].with_columns(n_gt_errors = pl.when(pl.col('n_gt_errors')==0)
                                                                              .then(pl.lit(float("nan")))
                                                                              .otherwise(pl.col('n_gt_errors'))
                                                                             ).with_columns(n_gt_checked = pl.when(pl.col('n_gt_errors').is_nan())
                                                                                                             .then(pl.lit(0))
                                                                                                             .otherwise(pl.col('n_gt_checked')),
                                                                                            n_gt_errors = pl.col('n_gt_errors').fill_nan(0)
                                                                             ).with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked')
                                                                             )
all_chroms = total_phasing_stats_datasets['switch_errors'].join(total_phasing_stats_datasets['flips_and_switches'], 
                                                                on=common_sample_columns+expt_columns, how='full')
all_chroms = all_chroms.drop([x for x in all_chroms.columns if '_right' in x]).join(total_phasing_stats_datasets['sample_genotype_concordance'],
                                                                                    on=common_sample_columns+expt_columns, how='full',
                                                                             ).join(sample_metadata.select(['sample_id','population','superpopulation','sex']),
                                                                                    on='sample_id', how='left')
all_chroms = all_chroms.drop([x for x in all_chroms.columns if '_right' in x])
all_samples = group_phasing_statistics(all_chroms, ['sample_id','genome','population','superpopulation','sex'] + expt_columns).drop('chrom')

HPRC_samples = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='HPRC_samples')).select('sample_id').unique().to_numpy().flatten()
trio_samples = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='trios')).select('sample_id').unique().to_numpy().flatten()

all_chroms = all_chroms.filter(((pl.col('ground_truth_data_source') == 'trios') & (pl.col('sample_id').is_in(trio_samples))) | ((pl.col('ground_truth_data_source') == 'HPRC_samples') & (pl.col('sample_id').is_in(HPRC_samples))))
all_samples = all_samples.filter(((pl.col('ground_truth_data_source') == 'trios') & (pl.col('sample_id').is_in(trio_samples))) | ((pl.col('ground_truth_data_source') == 'HPRC_samples') & (pl.col('sample_id').is_in(HPRC_samples))))
males = all_samples.filter(pl.col('sex')=='Male').select('sample_id').to_numpy().flatten()

###############
###Calc N50/L50
print ("calculating per-contig and per-sample N50/L50 - can take a minute")
print (f"organizing switch error block data ({print_runtime()})")
chrom_ends=dict()
chrom_length=dict()
chrom_starts=dict()
chm13_fai = '/dev/shm/chm13v2.0_maskedY_rCRS.fasta.fai'
grch38_fai = '/dev/shm/GRCh38_full_analysis_set_plus_decoy_hla.fasta.fai'

def get_chrom_ends(faidx):
    chrom_ends = dict()
    with open(faidx) as f:
        for line in f:
            chrom, end, *_ = line.split()
            chrom_ends[chrom] = int(end)
    return chrom_ends


chrom_ends['CHM13v2.0'] = get_chrom_ends(chm13_fai)
chrom_ends['GRCh38'] = get_chrom_ends(grch38_fai)
chrom_ends['GRCh38']['PAR1'] = 2781479 + 1    #convert bedfile 1-inclusive to python 0-exclusive
chrom_ends['GRCh38']['chrX'] = 155701382 + 1
chrom_ends['GRCh38']['PAR2'] = 156030895 + 1
chrom_ends['CHM13v2.0']['PAR1'] = 2394410 + 1
chrom_ends['CHM13v2.0']['chrX'] = 153925833 + 1
chrom_ends['CHM13v2.0']['PAR2'] = 154259566 + 1
chrom_length['CHM13v2.0'] = get_chrom_ends(chm13_fai)
chrom_length['GRCh38'] = get_chrom_ends(grch38_fai)
chrom_length['GRCh38']['PAR1'] = chrom_ends['GRCh38']['PAR1']
chrom_length['GRCh38']['chrX'] = chrom_ends['GRCh38']['chrX']-chrom_length['GRCh38']['PAR1']
chrom_length['GRCh38']['PAR2'] = chrom_ends['GRCh38']['PAR2']-chrom_length['GRCh38']['chrX']
chrom_length['CHM13v2.0']['PAR1'] = chrom_ends['CHM13v2.0']['PAR1'] = 2394410
chrom_length['CHM13v2.0']['chrX'] = chrom_ends['CHM13v2.0']['chrX']-chrom_length['CHM13v2.0']['PAR1']
chrom_length['CHM13v2.0']['PAR2'] = chrom_ends['CHM13v2.0']['PAR2']-chrom_length['CHM13v2.0']['chrX']
chrom_starts = {'GRCh38':dict(), 'CHM13v2.0':dict()}
chrom_starts['GRCh38']['PAR1'] = 0
chrom_starts['GRCh38']['chrX'] = chrom_ends['GRCh38']['PAR1']
chrom_starts['GRCh38']['PAR2'] = chrom_ends['GRCh38']['PAR1']
chrom_starts['CHM13v2.0']['PAR1'] = 0
chrom_starts['CHM13v2.0']['chrX'] = chrom_ends['CHM13v2.0']['PAR1']
chrom_starts['CHM13v2.0']['PAR2'] = chrom_ends['CHM13v2.0']['PAR1']

# Assemble dataframe of per-chrom ends and lengths in pandas - I'm still learning polars, and this is complicated. (Though I'm probably doing this in a very inefficient way. At least it works.)
chrom_ends = pd.concat((pd.DataFrame.from_dict(chrom_ends['GRCh38'], orient='index', columns=['end_chrom'])
                                    .reset_index().rename(columns={'index':'chrom'})
                                    .assign(genome='GRCh38'),
                        pd.DataFrame.from_dict(chrom_ends['CHM13v2.0'], orient='index', columns=['end_chrom'])
                                    .reset_index().rename(columns={'index':'chrom'})
                                    .assign(genome='CHM13v2.0')))
chrom_lengths = pd.concat((pd.DataFrame.from_dict(chrom_length['GRCh38'], orient='index', columns=['chrom_len'])
                                       .reset_index().rename(columns={'index':'chrom'})
                                       .assign(genome='GRCh38'),
                           pd.DataFrame.from_dict(chrom_length['CHM13v2.0'], orient='index', columns=['chrom_len'])
                                       .reset_index().rename(columns={'index':'chrom'})
                                       .assign(genome='CHM13v2.0')))
chrom_ends = chrom_ends.loc[chrom_ends.chrom.isin(contigs)]
chrom_ends = chrom_ends.merge(chrom_lengths, on=['genome','chrom'], how='left')

# Convert chrom ends/lengths df to polars and set data types
chrom_ends = pl.from_pandas(chrom_ends).with_columns(chrom=pl.col('chrom').cast(catChrom),
                                                     genome=pl.col('genome').cast(catGenomes),
                                                     end_chrom=pl.col('end_chrom').cast(pl.UInt32()),
                                                     chrom_len=pl.col('chrom_len').cast(pl.UInt32()))
block = total_phasing_stats_datasets['block_length'].unique(maintain_order=True)
# block length is probably a misnomer - it actually records switch error locations. 
# We need to get the per-block length. Note that the first and last switch error site 
# are not actually switch errors, just the first/last het site. (It's a quirk/error in shapeit5.)
# So we'll set first and last site from [loc of first het, loc of last het] to [0, len of contig].
    # First, identify sites that are the minimum site for that sample and contig
id_columns=['sample_id','method_of_phasing','ground_truth_data_source','chrom','genome']
block = (block.join(block.group_by(id_columns).min().rename({'switch_error_site':'min'}), on=id_columns)
                .join(block.group_by(id_columns).max().rename({'switch_error_site':'max'}), on=id_columns)
                .join(chrom_ends, on=['chrom','genome']))
block = (block.with_columns(switch_error_site=pl.when(pl.col('switch_error_site')==pl.col('min'), ~pl.col('chrom').is_in(('PAR1','PAR2','chrX')))
                                                .then(0)
                                                .when(pl.col('genome')=='GRCh38', pl.col('switch_error_site')==pl.col('min'), pl.col('chrom').is_in(('PAR1','PAR2','chrX')))
                                                .then(pl.col('chrom').replace_strict(chrom_starts['GRCh38'], default=pl.col('switch_error_site')))
                                                .when(pl.col('genome')=='CHM13v2.0', pl.col('switch_error_site')==pl.col('min'), pl.col('chrom').is_in(('PAR1','PAR2','chrX')))
                                                .then(pl.col('chrom').replace_strict(chrom_starts['CHM13v2.0'], default=pl.col('switch_error_site')))
                                                .when(pl.col('switch_error_site')==pl.col('max'))
                                                .then(pl.col('end_chrom'))
                                                .otherwise(pl.col('switch_error_site')))
                .drop(['min','max','end_chrom']))

# Don't include male chrXs in N50 calculations - they're phased by default, so it would inflate the N50 of males
block = block.filter(~((pl.col('chrom')=='chrX') & (pl.col('sample_id').is_in(males))))
                
sorting_cols = ['genome','method_of_phasing','ground_truth_data_source','sample_id','chrom','chrom_len','switch_error_site']
sample_id_columns = ['genome','method_of_phasing','ground_truth_data_source','sample_id']
print (f"calculating per contig N50 ({print_runtime()})")
per_contig_N50 = (block.sort(sorting_cols)                           # Sort data to make sure blocks are in sequential order
                       .group_by(id_columns+['chrom_len'])           # Grouby unique contigs, add chrom_len to prevent any aggregation so we can use it downstream.
                       .agg(blocksize=pl.col('switch_error_site')    # Generate sorted list of blocksizes per contig by:
                                        .diff()                      # Calculating size of each block
                                        .fill_null(0)                # n_blocks = n_switch_sites + 1, the first site is always null when using diff. Fill null with 0 to ensure it is ignored when calculating N50. Dropping null doesn't always work the way I expect it to, filling null with 0 does.
                                        .sort(descending=True),      # Sort from largest to smallest block. 
                            L50=pl.col('switch_error_site')              # L50: How many blocks (largest to smallest) needed to get to 50%+1 of contig length
                                  .diff()                                # First, repeat code to get blocksize. Polars should see that the code is repeating and avoid duplicate calculations.
                                  .drop_nulls()                          #
                                  .sort(descending=True)                 # With per-contig list of block lengths, sorted largest to smallest...
                                  .cum_sum()                             # take cumulative sum to convert to list of cumulative block sizes
                                  .search_sorted(pl.col('chrom_len')     # L50 = the number of contigs necessary to get to contig length/2
                                                   .first()/2))          # This code finds the index you would insert (contig length/2) to keep the list of cumsum lengths sorted. To do this, you would find the number of blocks necessary to get to (contig_length/2), and insert (contig_len/2) afterwards, at the location of the next block. So the insertion index is the number of blocks necessary to get to the first block whose cumulative sum is over contig_length/50 - i.e., the L50.
                       .with_columns(N50=pl.col('blocksize')             # The N50 is just the size of the block needed to get to 50%+1 of the contig_length
                                           .list.get(pl.col('L50'))))     # As we've just shown, given a list of blocks sorted from biggest to smallest, that block is located at the L50's index.
print (f"Calculating per sample N50 ({print_runtime()})")
per_sample_N50 = (per_contig_N50.group_by(sample_id_columns)             # Group data by sample.
                                .agg(blocksize=pl.col('blocksize')       # blocksize becomes list of lists, 
                                                 .flatten()              # flatten to get one list of contig blocks in sample
                                                 .sort(descending=True), # Then sort largest to smallest
                                     L50=pl.col('blocksize')             # Do the same thing as above to get the sample L50
                                           .flatten()
                                           .sort(descending=True)
                                           .cum_sum()
                                           .search_sorted(pl.col('chrom_len') # This time sum chrom lengths to get genome length
                                                            .unique()
                                                            .sum()/2))
                                .with_columns(N50=pl.col('blocksize')
                                                    .list.get(pl.col('L50')))
                                .drop('blocksize'))

per_contig_N50 = (per_contig_N50.drop(['blocksize','chrom_len'])
                               .filter(((pl.col('ground_truth_data_source') == 'trios') & 
                                        (pl.col('sample_id').is_in(trio_samples))) | 
                                       ((pl.col('ground_truth_data_source') == 'HPRC_samples') & 
                                        (pl.col('sample_id').is_in(HPRC_samples))))) # Filter out samples where we got a blank entry in the SHAPEIT5 report
per_sample_N50 = (per_sample_N50.filter(((pl.col('ground_truth_data_source') == 'trios') & 
                                        (pl.col('sample_id').is_in(trio_samples))) | 
                                       ((pl.col('ground_truth_data_source') == 'HPRC_samples') & 
                                        (pl.col('sample_id').is_in(HPRC_samples)))))

all_chroms = all_chroms.join(per_contig_N50, on=['genome','method_of_phasing','ground_truth_data_source','sample_id','chrom'], how='left')
all_samples = all_samples.join(per_sample_N50, on=['genome','method_of_phasing','ground_truth_data_source','sample_id'], how='left')

all_ancestries = group_phasing_statistics(all_chroms, ['genome','population','superpopulation'] + expt_columns).drop(['sex','sample_id'])
all_methods = group_phasing_statistics(all_samples, ['genome'] + expt_columns).drop(['sample_id','sex','population','superpopulation'])



print (f'saving summary stats to disk ({print_runtime()})')
all_chroms.write_parquet(f'{summary_statistics_folder}/chroms.parquet')
all_samples.write_parquet(f'{summary_statistics_folder}/samples.parquet')
all_ancestries.write_parquet(f'{summary_statistics_folder}/ancestries.parquet')
all_methods.write_csv(f'{summary_statistics_folder}/methods.tsv')

print (f'Loading imputation stats ({print_runtime()})')

###Finally, load imputation stats

per_sample_imputation_performance = list()
per_variant_category_imputation_performance = list()
per_genome_imputation_performance = list()
per_ancestry_imputation_performance = list()


ancestries = ['all','WestEurasia','SouthAsia','EastAsia','Africa','America','Oceania','CentralAsiaSiberia']
glimpse2_headers={'syntenic_maf_overall_bins': {'error.grp':  ('index','group_name','num_variants','mean_AF','num_AA_mismatches','num_AA','num_Aa_mismatches','num_Aa','num_aa_mismatches','num_aa','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches'),
                                                'error.spl':  ('var_category_id','id','sample_name','num_AA','num_Aa','num_aa','num_calls_filtered_out','num_AA_matches','num_Aa_matches','num_aa_matches','num_AA_mismatches','num_Aa_mismatches','num_aa_mismatches','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches','non_reference_discordance_percent','best_gt_rsquared','imputed_ds_rsquared'),
                                                'rsquare.grp':('group_name','num_variants','mean_AF','best_gt_rsquared','imputed_ds_rsquared')},
                  'syntenic_maf_vartype_bins': {'error.grp':  ('index','group_name','num_variants','mean_AF','num_AA_mismatches','num_AA','num_Aa_mismatches','num_Aa','num_aa_mismatches','num_aa','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches'),
                                                'error.spl':  ('var_category_id','id','sample_name','num_AA','num_Aa','num_aa','num_calls_filtered_out','num_AA_matches','num_Aa_matches','num_aa_matches','num_AA_mismatches','num_Aa_mismatches','num_aa_mismatches','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches','non_reference_discordance_percent','best_gt_rsquared','imputed_ds_rsquared'),
                                                'rsquare.grp':('group_name','num_variants','mean_AF','best_gt_rsquared','imputed_ds_rsquared')},
                  'r2_bins':                   {'error.grp':  ('var_category_id','group_name','num_variants','mean_AF','num_AA','num_Aa','num_aa','num_calls_filtered_out','num_AA_matches','num_Aa_matches','num_aa_matches','num_AA_mismatches','num_Aa_mismatches','num_aa_mismatches','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches','best_gt_rsquared','imputed_ds_rsquared'),
                                                'error.spl':  ('var_category_id','id','sample_name','num_AA','num_Aa','num_aa','num_calls_filtered_out','num_AA_matches','num_Aa_matches','num_aa_matches','num_AA_mismatches','num_Aa_mismatches','num_aa_mismatches','percent_AA_mismatches','percent_Aa_mismatches','percent_aa_mismatches','non_reference_discordance_percent','best_gt_rsquared','imputed_ds_rsquared'),
                                                'rsquare.grp':('group_name','num_variants','mean_AF','best_gt_rsquared','imputed_ds_rsquared')}}

for genome in ('T2T','GRCh38'):
    for variant_set in ('', '_snps'):
        for ground_truth in ('SGDP','pangenome'):
            for panel in ('native_panel','native_panel.common_variants','lifted_panel','lifted_panel.common_variants'):
                for ancestry in ancestries:
                    if ancestry == 'all':
                        base_report_name = f'genomewide_imputation_evaluation_{run_version}/{genome}{variant_set}.{ground_truth}.{panel}.glimpse2_concordance'
                    else:
                        base_report_name = f'genomewide_imputation_evaluation_{run_version}/ancestry_specific/{genome}{variant_set}.{ground_truth}.{panel}.{ancestry}.glimpse2_concordance'
                    for bin_grouping in glimpse2_headers.keys():
                        if (ancestry != 'all') and (bin_grouping != 'r2_bins'):
                            continue
                        if (ancestry != 'all') and (ground_truth == 'pangenome'):
                            # Pangenome contains only 44 highly diverse samples, not enough samples per ancestry to discern anything meaningful
                            continue 
                        bin_grouping_list = list()                    
                        for i, report_name in enumerate(glimpse2_headers[bin_grouping].keys()):
                            try:
                                report = pd.read_csv(base_report_name + '_' + bin_grouping + '.' + report_name + '.txt.gz', sep=' ',
                                                        header=None, comment='#', names=glimpse2_headers[bin_grouping][report_name]
                                                        ).assign(genome=genome, 
                                                                 ancestry=ancestry,
                                                                 dataset=ground_truth+variant_set,
                                                                 panel=panel)
                            except:
                                print(base_report_name + '_' + bin_grouping + '.' + report_name + '.txt.gz')
                                continue
                            if 'var_category_id' in report.columns:
                                report['var_category_id'] = report.var_category_id.map({'GCsV':'All variants','GCsVAF':'All variants',
                                                                                        'GCsI':'Indels','GCsIAF':'Indels',
                                                                                        'GCsS':'SNPs','GCsSAF':'SNPs'})
                            if 'sample_name' in report.columns:
                                per_sample_imputation_performance.append(report)
                                continue
                             
                            if bin_grouping == 'r2_bins':
                                report['variant_bin'] = report.group_name.apply(lambda x: r2_bin_names[x])
                                report['Synteny'] = 'All'
                            else:
                                report['variant_bin'] = report.group_name.str.split('_').str[-1]
                            if "syntenic" in bin_grouping:
                                report['Synteny'] = report.group_name.str.split('_').str[0].map({'SYNTENIC':'Syntenic','NONSYNTENIC':'Nonsyntenic'})#.astype(catSyn)
                            if 'vartype' in bin_grouping:
                                report['var_category_id'] = report.group_name.str.split('_').str[1].map({'INDEL':'Indels','SNP':'SNPs'})#.astype(catVarType)
                            elif 'var_category_id' not in report.columns:
                                report['var_category_id'] = 'All variants'
                            
                            bin_grouping_list.append(report.drop(columns='group_name'))

                    per_var_cat = bin_grouping_list[0].merge(bin_grouping_list[1], on=['ancestry','genome','variant_bin','Synteny','var_category_id','num_variants', 'dataset', 'panel'], how='outer', suffixes=('','_y'))
                    if 'index' in per_var_cat.columns:
                        per_var_cat = per_var_cat.drop(columns=['index'])
                    per_var_cat = per_var_cat.drop(columns=[c for c in per_var_cat.columns if c[-2:] == '_y'])
                    per_variant_category_imputation_performance.append(per_var_cat)

per_variant_category_imputation_performance = pd.concat(per_variant_category_imputation_performance)
per_sample_imputation_performance = pd.concat(per_sample_imputation_performance)

per_variant_category_imputation_performance['non_reference_discordance_percent'] = ((per_variant_category_imputation_performance.num_Aa_mismatches + per_variant_category_imputation_performance.num_aa_mismatches) / 
                                                                                    (per_variant_category_imputation_performance['num_Aa']+per_variant_category_imputation_performance['num_aa']))*100

per_variant_category_imputation_performance.to_parquet(f'{summary_statistics_folder}/per_variant_category_imputation_performance.parquet')
per_sample_imputation_performance.to_parquet(f'{summary_statistics_folder}/per_sample_imputation_performance.parquet')

print ('Done!')
print (f"Total runtime: {datetime.now()}")