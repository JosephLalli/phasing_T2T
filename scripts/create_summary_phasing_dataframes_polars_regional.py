#!/usr/bin/env python3

import polars as pl
import polars.selectors as cs
import pandas as pd
import os
import glob
from datetime import datetime

# When running on test datasets chr15 and chr22, uses these resources on our server:
# user(s)   sys(s)    elapsed    cpu(%)  maxRSS(GB)   io_in(MB)    io_out(MB)   exit
# 2088.63   281.44    0:34.32    6905%   29.072       0.00         815.77       0

# Get current time to measure runtime
StartTime = datetime.now() 
def print_runtime():
    time = datetime.now() - StartTime
    s = time.seconds % 60
    m = (time.seconds - s)//60
    return f"{m}m {s}s"

# Add CLI to accept run suffixes for CHM13 and GRCh38
import argparse
parser = argparse.ArgumentParser(description='Create summary phasing dataframes')
parser.add_argument('--CHM13_run_suffix', default='CHM13v2.0', help='suffix used for CHM13 run (default: CHM13v2.0)')
parser.add_argument('--GRCh38_run_suffix', default='GRCh38', help='suffix used for GRCh38 run (default: GRCh38)')
parser.add_argument('--test', action='store_true', help='process a test with fewer contigs (chr15 and chr22 only)')
args = parser.parse_args()
test_run = args.test
t2t_suffix = args.CHM13_run_suffix
grch38_suffix = args.GRCh38_run_suffix

# Ensure we are in the project root directory
if os.getcwd().split('/')[-1] == 'phasing_T2T_project':
    pass
elif os.getcwd().split('/')[-1] == 'scripts':
    os.chdir('..')

# -------- Define folder and file locations --------
# File locations
summary_statistics_folder = 'intermediate_data'
t2t_phasing_stats_folder = f'SHAPEIT5_switch_output/phasing_stats_CHM13v2.0_{t2t_suffix}'
grch_phasing_stats_folder = f'SHAPEIT5_switch_output/phasing_stats_GRCh38_{grch38_suffix}'
imputation_statistics_folder = f'imputation_statistics/imputation_results_{t2t_suffix}'

regional_bedfiles = {'CHM13v2.0': {'in_STRs':   'resources/GIABv3.6_bedfiles/CHM13_AllTandemRepeats.bed.gz',
                                   'in_segdups':'resources/GIABv3.6_bedfiles/CHM13_segdups.bed.gz',
                                   'in_platinum_STRs':'resources/CHM13_PlatinumTandemRepeats.bed.gz'},
                        'GRCh38': {'in_STRs':   'resources/GIABv3.6_bedfiles/GRCh38_AllTandemRepeats.bed.gz',
                                   'in_segdups':'resources/GIABv3.6_bedfiles/GRCh38_segdups.bed.gz',
                                   'in_platinum_STRs':'resources/GRCh38_PlatinumTandemRepeats.bed.gz'}}


# -------- Load individuals metadata --------

sample_metadata = (pl.read_csv('resources/sample_subsets/superpopulations.samples.txt', separator=' ', infer_schema_length=1000,
                                                                             new_columns=['trio_id', 'sample_id','father','mother','sex','population','superpopulation'])
                     .with_columns(sex = pl.col('sex').cast(pl.String()).replace('1','Male').replace('2','Female')))
non_trio_HGSVC_samples = pd.read_csv('resources/sample_subsets/HGSVC_not_part_of_trio.txt', header=None, names=['sample_id'])
HGSVC_HPRC_samples = pd.read_csv('resources/sample_subsets/HGSVC_HPRC_samples.txt', header=None, names=['sample_id']).sample_id.to_list()
HGSVC_samples = pd.read_csv('resources/sample_subsets/HGSVC_samples_in_1KGP.txt', header=None, names=['sample_id']).sample_id.to_list()
trio_phased_samples=pd.read_csv('resources/sample_subsets/children.txt', header=None, names=['sample_id']).sample_id.to_list()


# -------- Define script-wide constants --------

pangenome_only_samples = ['HG01123','HG02109','HG02486','HG02559','NA21309']

ground_truths = ['trios','HPRC_samples', 'HGSVC_samples', 'HPRC_HGSVC_all_samples','HPRC_HGSVC_probands', 
                 'HGSVC_probands', 'HGSVC_parents','HGSVC_samples_nontrios_only']
methods = ['phased_with_parents_and_pedigree', 'phased_without_parents_or_pedigree','1kgp_variation_phased_with_reference_panel']
common_sample_columns = ['sample_id','chrom', 'genome']
common_variant_columns = ['variant_id','position','chrom', 'genome']
expt_columns = ['method_of_phasing','ground_truth_data_source']


# -------- Polars enum datatypes and schemas --------
# Whole genome contig definition used in paper
contigs = ['chr' + str(x) for x in list(range(1,23))] + ['PAR1', 'chrX', 'PAR2']
catChrom = pl.Enum(contigs + ['chrY'])

if test_run:
    # Test contig set
    contigs = ['chr15','chr22']

catMethod = pl.Enum(methods)
catGroundTruth = pl.Enum(ground_truths)
catGenomes = pl.Enum(['GRCh38','CHM13v2.0'])
catSamples = pl.Enum(sorted(list(set(pl.Series(sample_metadata.select('sample_id')).to_list() + pangenome_only_samples))))
catSyn = pl.Enum(['All','Syntenic','Nonsyntenic']+
                 [region for region in regional_bedfiles['GRCh38'].keys()]+
                 ['not_'+region for region in regional_bedfiles['GRCh38'].keys()]+
                 ['multiallelic', 'biallelic'])
catVarType = pl.Enum(['SNP','Indel','SNPs + Indels'])

# catSamples is defined, convert sample_metadata to the enum
sample_metadata = sample_metadata.with_columns(sample_id = pl.col('sample_id').cast(catSamples))

all_variant_annotations_dtypes = pl.Schema( {'ID':pl.String(),
                                             'CHROM': pl.String(),
                                             'POS': pl.UInt32(),
                                             'ALT': pl.String(),
                                             'REF': pl.String(),
                                             'NEGATIVE_TRAIN_SITE': pl.String(),
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
                                             'SYNTENIC': pl.String()})

headers = {'flips_and_switches'           :['sample_id','n_total_switch_errors','n_flips', 'n_consecutive_flips', 'n_true_switch_errors', 'num_correct_phased_hets', 'n_total_hets', 'flip_and_switches_SER','flip_error_rate','true_switch_error_rate','accurately_phased_rate'],
           'switch_errors'                :['sample_id', 'n_switch_errors',  'n_checked', 'switch_error_rate'],
           'frequency'                    :['minor_allele_count', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'type'                         :['is_SNP', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'variant_genotype_concordance' :['variant_id', 'position', 'n_gt_errors', 'n_gt_checked', 'gt_error_rate'],
           'block_length'                 :['sample_id','switch_error_site'],
           'variant_switch'               :['variant_id', 'position', 'n_switch_errors', 'n_checked', 'switch_error_rate'],
           'sample_genotype_concordance'  :['sample_id','n_gt_errors','n_gt_checked','gt_error_rate']}

schemas = {'flips_and_switches'           :{'sample_id':pl.String(), 'n_total_switch_errors':pl.UInt32(),'n_flips':pl.UInt32(),'n_consecutive_flips':pl.UInt32(),'n_true_switch_errors':pl.UInt32(),'num_correct_phased_hets':pl.UInt32(),'n_total_hets':pl.UInt32(),'flip_and_switches_SER':pl.Float32(),'flip_error_rate':pl.Float32(),'true_switch_error_rate':pl.Float32(),'accurately_phased_rate':pl.Float32()},
           'switch_errors'                :{'sample_id':pl.String(), 'n_switch_errors':pl.UInt32(),  'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'frequency'                    :{'minor_allele_count':pl.UInt16(), 'n_switch_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'type'                         :{'is_SNP':pl.Boolean(), 'n_gt_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'gt_error_rate':pl.UInt32()},
           'variant_genotype_concordance' :{'variant_id':pl.String(), 'position':pl.UInt32(), 'n_gt_errors':pl.UInt32(), 'n_gt_checked':pl.UInt32(), 'gt_error_rate':pl.Float32()},
           'block_length'                 :{'sample_id':pl.String(),'switch_error_site':pl.Int32()},
           'variant_switch'               :{'variant_id':pl.String(), 'position':pl.UInt32(), 'n_switch_errors':pl.UInt32(), 'n_checked':pl.UInt32(), 'switch_error_rate':pl.Float32()},
           'sample_genotype_concordance'  :{'sample_id':pl.String(), 'n_gt_errors':pl.UInt32(),  'n_gt_checked':pl.UInt32(), 'gt_error_rate':pl.Float32()}}

analysis_suffixes = {'trios':        {'phased_without_parents_or_pedigree'        :'rare_noparents_vs_trios'},
                     'HPRC_samples': {'phased_without_parents_or_pedigree'        :'noparents_vs_HPRC',
                                      'phased_with_parents_and_pedigree'          :'3202_panel_vs_HPRC',
                                      '1kgp_variation_phased_with_reference_panel':'rare_1kgp_pangenome_panelphased_vs_pangenome'},
                     'HGSVC_samples':{'phased_with_parents_and_pedigree'          :'3202_panel_vs_HGSVC'},
                     'HPRC_HGSVC_all_samples':{'phased_with_parents_and_pedigree'     :'3202_panel_vs_HPRC_and_HGSVC_all_samples'},
                     'HPRC_HGSVC_probands':{'phased_with_parents_and_pedigree'     :'3202_panel_vs_HPRC_and_HGSVC_trio_probands_only'},
                     'HGSVC_probands':{'phased_with_parents_and_pedigree'     :'3202_panel_vs_HGSVC_probands'},
                     'HGSVC_parents':{'phased_with_parents_and_pedigree'     :'3202_panel_vs_HGSVC_parents'},
                     'HGSVC_samples_nontrios_only':{'phased_with_parents_and_pedigree' :'3202_panel_vs_HGSVC_not_part_of_trio'}}

# Shapeit5 puts out a set of statistic files with $prefix.$suffix.tsv format
# where suffix is one of the following:
shapeit5_stat_suffixes = {'flips_and_switches':'flipsAndSwitches',
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

# I've tried to define in advance bins that are roughly logorimically spaced 
# with roughly equal numbers of variants in each bin.
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


# -------- Define convenience functions --------

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
    if 'n_true_switch_errors' in df_cols:
        # See definitions below.
        expressions.update({'flip_error_rate':(pl.col('n_flips')/(pl.col('n_total_hets')))*100,
                            'true_switch_error_rate':(pl.col('n_true_switch_errors')/pl.col('n_total_hets')) * 100,
                            'flip_and_switches_SER':(pl.col('n_total_switch_errors')/pl.col('n_total_hets')) * 100,
                            'accurately_phased_rate':(pl.col('num_correct_phased_hets')/pl.col('n_total_hets'))*100})
    if 'n_switch_errors' in df_cols:
        expressions.update({'switch_error_rate':(pl.col('n_switch_errors')/pl.col('n_checked'))*100})

    return df.with_columns(**expressions)

def add_bedfile_column(variant_df: pl.DataFrame, bedfile: str, name: str = 'in_STR', sorted: bool = False) -> pl.DataFrame:
    bed_df = pl.read_csv(bedfile, separator='\t', columns=['chrom','start','end'], 
                         schema_overrides={'chrom': catChrom, 'start': pl.UInt32, 'end': pl.UInt32}
                         ).unique().sort(by=["chrom", "start"])
    if not sorted:
        variant_df = variant_df.sort(by=["chrom", "position"])

    variant_df = variant_df.join_asof(
            bed_df,
            left_on="position",
            right_on="start",
            by="chrom",
            strategy="backward",   # nearest earlier or equal BED start
            suffix="_bed",        # avoid col-name clashes
        ).with_columns(
            **{name:(pl.col("start").is_not_null() & 
                     (pl.col("position") < pl.col("end"))).fill_null(False)}
                    # Check if variant position falls within the BED region: start <= position < end
        ).drop("start", "end")
    return variant_df


def add_bedfiles(variant_df: pl.DataFrame, bedfiles: dict[str, str]) -> pl.DataFrame:
    variant_df = variant_df.with_columns(position = pl.col("position").cast(pl.UInt32))
    for name, bedfile in bedfiles.items():
        variant_df = add_bedfile_column(variant_df, bedfile, name=name, sorted=False)
    return variant_df


def sum_an_acs(files):
    anac = list()
    for tsv in files:
        anac.append(pl.scan_csv(tsv, separator='\t', has_header=False, new_columns=['variant_id','contextual_MAC','contextual_AN'],
                                        schema_overrides={'variant_id':pl.String(),
                                                            'contextual_MAC':pl.UInt16(),
                                                            'contextual_AN':pl.UInt16()}))
    anac = pl.concat(anac).group_by('variant_id').sum().with_columns(contextual_MAC=pl.col('contextual_MAC').cast(pl.UInt16()),
                                                                     contextual_AN=pl.col('contextual_AN').cast(pl.UInt16()))
    return anac


# ------- start of script -------
variant_frequency_data=list()

print(f'Gathering alt frequencies for each phased panel, dropping trio-private singletons ({print_runtime()})')
for genome, run_suffix in [('GRCh38', grch38_suffix),
                           ('CHM13v2.0', t2t_suffix)]:
    for contig in contigs:
        bash_script_output_dir = f"{summary_statistics_folder}/variant_frequency_stats/{genome}"
        with open (f'{bash_script_output_dir}/{contig}_private_singletons.txt', 'r') as f:
            private_variants = {x.strip() for x in f}

        an_acs = ((['phased_with_parents_and_pedigree'], f"{bash_script_output_dir}/{contig}_3202_AC_AN.tsv"),
                  (['1kgp_variation_phased_with_reference_panel'], f"{bash_script_output_dir}/{contig}_2430_AC_AN.tsv"),
                  (['phased_without_parents_or_pedigree'], f"{bash_script_output_dir}/{contig}_2002_AC_AN.tsv"))
        for methods, tsv in an_acs:
            if type(tsv) == list:
                anac = sum_an_acs(tsv)
            else:
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
                if method in ('phased_without_parents_or_pedigree'):
                    variant_frequency_data.append(anac.filter(~pl.col('variant_id').is_in(private_variants))
                                        .with_columns(is_private_variant=pl.col('variant_id').is_in(private_variants),
                                                      method_of_phasing=pl.lit(method).cast(catMethod),
                                                      genome=pl.lit(genome).cast(catGenomes),
                                                      contextual_AN=pl.col('contextual_AN').max()),)
                else:
                    variant_frequency_data.append(anac.with_columns(is_private_variant=pl.col('variant_id').is_in(private_variants),
                                                      method_of_phasing=pl.lit(method).cast(catMethod),
                                                      genome=pl.lit(genome).cast(catGenomes),
                                                      contextual_AN=pl.col('contextual_AN').max()))
contextual_variants = pl.concat(variant_frequency_data)


# Gather per variant statistics calculated by bcftools query
print (f'Collecting each reference panel\'s minor allele frequency statistics ({print_runtime()})')
variant_data_lazy = list()
for chrom in contigs:
    for genome, run_suffix in [('GRCh38', grch38_suffix),
                               ('CHM13v2.0', t2t_suffix)]:
        contig_variant_data = f"{summary_statistics_folder}/variant_frequency_stats/{genome}/1KGP.{genome}.{chrom}.snp_indel.phasing_qual_pass.fully_annotated.tsv"
        if os.path.exists(contig_variant_data):
            # collect each file as a lazyframe
            contig_variant_data = pl.scan_csv(contig_variant_data, separator='\t', null_values=['.','NA'],
                                                schema_overrides = all_variant_annotations_dtypes
                                                ).drop('REF', strict=False
                                                # Cast genome and chrom to enum categories
                                                # and since fully_annotated.tsv sometimes uses true, sometimes use 1 to indicate true,
                                                # use boolean logic to convert to consistent boolean values
                                                ).with_columns(genome=pl.lit(genome).cast(catGenomes),
                                                               CHROM=pl.col('CHROM').cast(catChrom),
                                                               SYNTENIC=((pl.col('SYNTENIC') == "true") | (pl.col('SYNTENIC') == "1")).cast(pl.Boolean()).fill_null(False),
                                                               NEGATIVE_TRAIN_SITE=pl.when(pl.col('NEGATIVE_TRAIN_SITE').eq("true"))
                                                                                                .then(True)
                                                                                                .when(pl.col('NEGATIVE_TRAIN_SITE').eq("1"))
                                                                                                .then(True)
                                                                                                .otherwise(False)
                                                                                                .cast(pl.Boolean()),
                                                ).rename({'AN':'AN_original_panel',
                                                          'AC':'AC_original_panel',
                                                          'MAF':'MAF_original_panel',
                                                          'SYNTENIC':'Syntenic'})

            # Once annotated variant is scanned in, annotate variants which are multiallelic cause problems.
            # I'd like to track these separately.
            # I'm specifically interested in sites with many different potential indel alts (eg STRs).
            # VCF representation of these STRs is complex.
            # To cover all my bases, I'll classify a variant as multiallelic if there is more than one variant at the same position,
            # Or a variant is a deletion that overlaps one of those variants
            # This requires first identifying variants which start at the same location,
            # then identifying variants whose start-end range overlaps one of these variants.
            # First: identify positions with multiple variants
            multiallelic_pos = (contig_variant_data.group_by("POS")
                                                   .len()
                                                   .filter(pl.col("len") > 1)
                                                   .select(pl.col("POS").alias("multi_pos"))
                                                   .unique()
                                                   .sort(["multi_pos"])) # identify positions w/ more than one possible alt allele
            # second, extract variant length from ID field
            parts = pl.col("ID").str.split_exact("_", 3)
            variant_len_expr = (parts.struct.field("field_2").str.len_bytes()
                               + parts.struct.field("field_3").str.len_bytes()
                               ).cast(pl.Int64).fill_null(0)
            # and identify start/end position of all variants
            contig_variant_data = (
                contig_variant_data
                .with_columns(
                    start=(pl.col("POS") - variant_len_expr).cast(pl.Int64),
                    end=(pl.col("POS") + variant_len_expr).cast(pl.Int64),
                )
                .sort("end")
                # use join_asof to identify variants whose start-end range overlaps a multiallelic position
                .join_asof(
                    multiallelic_pos,
                    left_on="end",
                    right_on="multi_pos",
                    strategy="backward",
                )
                .with_columns(
                    multiallelic = (
                        (pl.col("multi_pos").is_not_null()
                        & (pl.col("multi_pos") >= pl.col("start"))
                    ).cast(pl.Boolean()).fill_null(False)  # fill nulls with False
                ))
                .drop("multi_pos", "start", "end")
                )
        # If there wasn't any data for this config, skip it
        else:
            continue
        variant_data_lazy.append(contig_variant_data)

print (f"Collecting variant statistics for all variants in study ({print_runtime()})")
variant_data_lazy = pl.concat(variant_data_lazy).unique(subset=['ID','Syntenic','genome'], keep='first')
print (f"Writing variant statistics to disk ({print_runtime()})")
# raise Exception
# since variant data can be larger than memory, i'd like to process these files and then sink them all into a parquet file
variant_data_lazy.sink_parquet(f'{summary_statistics_folder}/bcftools_query_variant_data.parquet')

# Now we have a variant data df on disk in parquet format with the most efficient data types
# But this just contains variant annotations and allele frequencies.
# We still need to collect each variant's phasing statistics from SHAPEIT5_switch output files,
# First, for each kind of shapeit5_switch output file, 
# let's collect the whole study's worth of data into a single dataframe per output file type
print (f'loading SHAPEIT5_switch statistics! ({print_runtime()})')
total_phasing_stats_datasets = {stat: list() for stat in schemas.keys()}
for contig in contigs:
    per_contig_phasing_stats_datasets = {stat: list() for stat in schemas.keys()}
    for statistic, schema in schemas.items():
        for genome, phasing_stats_folder in [('GRCh38',grch_phasing_stats_folder),('CHM13v2.0', t2t_phasing_stats_folder)]:
            for verification_sample in analysis_suffixes.keys():
                for phasing_method, prefix in analysis_suffixes[verification_sample].items():
                    if genome == 'GRCh38' and 'noparents' in prefix:
                        continue
                    file = f'{phasing_stats_folder}/{prefix}_{contig}.{shapeit5_stat_suffixes[statistic]}*'
                    try:
                        file = glob.glob(f'{phasing_stats_folder}/{prefix}_{contig}.{shapeit5_stat_suffixes[statistic]}*')
                        if len(file) == 1:  ##, f"multiple files found for {phasing_stats_folder}/{prefix}_{contig}_0.2.*{shapeit5_stat_suffixes[statistic]}*:\n{file}"
                            file = file[0]
                        elif len(file) == 0:
                            file = glob.glob(f'{phasing_stats_folder}/{prefix}_{contig}_CHM13v2.0.{shapeit5_stat_suffixes[statistic]}*')
                            if len(file) == 1:
                                file = file[0]
                            else:
                                raise IndexError
                        else:
                            raise IndexError
                    except IndexError:
                        print (file)
                        raise FileNotFoundError(f"no file found for {file}")
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
        
        # I want to merge these two per-variant datasets together before appending to the study-wide dataframes
        if statistic not in ('variant_switch','variant_genotype_concordance'):
            total_phasing_stats_datasets[statistic].append(per_contig_phasing_stats_datasets[statistic])
    
    # I do that merging here
    merged_variant_switch_and_genotype_data = (per_contig_phasing_stats_datasets['variant_switch'].join(per_contig_phasing_stats_datasets['variant_genotype_concordance'], 
                                                                                                        on=['variant_id','position','chrom','ground_truth_data_source','genome','method_of_phasing'], 
                                                                                                        how='left', validate='1:1').drop(['gt_error_rate']))
    total_phasing_stats_datasets['variant_switch'].append(merged_variant_switch_and_genotype_data)

# Write results to disk.
# For per variant results, add per-variant metadata generated with bcftools query above
bcftools_query_variant_data = pl.scan_parquet(f'{summary_statistics_folder}/bcftools_query_variant_data.parquet')


# we never saved anything to this key. Instead we merged all the data with variant_switch 
# and gathered as one per-variant df
del total_phasing_stats_datasets['variant_genotype_concordance'] 

# For the other datasets, save to disk
for stat, datasets in total_phasing_stats_datasets.items():
    print (f'Collecting statistics on {stat}... ({print_runtime()})')
    dataset = pl.concat(datasets)
    if stat != 'variant_switch':
        total_phasing_stats_datasets[stat] = dataset.collect()
        total_phasing_stats_datasets[stat].write_parquet(f'{summary_statistics_folder}/{stat}.parquet')
    else:
        # add the variant frequency of the reference panel used to create the measured phased haplotypes that shapeit5 reported on
        total_phasing_stats_datasets[stat] = (dataset.join(contextual_variants, on=['variant_id', 'method_of_phasing', 'genome'], how='left', validate='m:1') 
                                                     # and whether each variant is syntenic or multiallelic
                                                     .join(bcftools_query_variant_data.select(['ID','genome','Syntenic','multiallelic'])
                                                                       .rename({'ID':'variant_id'}), 
                                                           on=['variant_id','genome'], how='left', validate='m:1'))
        total_phasing_stats_datasets[stat] = total_phasing_stats_datasets[stat].collect(engine='streaming')
        total_phasing_stats_datasets[stat].write_parquet(f'{summary_statistics_folder}/variants.parquet')

print (f'Annotating variants by region...')
total_phasing_stats_datasets['variant_switch'] = pl.concat([add_bedfiles(
                                                            total_phasing_stats_datasets['variant_switch'].filter(pl.col('genome')==genome), 
                                                            regions) 
                                                            for genome, regions in regional_bedfiles.items()])

print (f'Writing a second, abbreviated per-variant switch statistics to disk ({print_runtime()})')
total_phasing_stats_datasets['variant_switch'].filter(pl.col('contextual_AN') > 0).write_parquet(f'{summary_statistics_folder}/variants.nozeros.parquet')

print (f'Binning variants ({print_runtime()})')
regional_cols = list(set([r for r in regional_bedfiles['GRCh38'].keys()] + list(regional_bedfiles['CHM13v2.0'].keys())))
MAF_performance_variant_df = total_phasing_stats_datasets['variant_switch'].select(['variant_id','n_switch_errors','n_checked','n_gt_errors', 
                                                                                    'n_gt_checked','method_of_phasing','ground_truth_data_source',
                                                                                    'genome', 'contextual_MAC', 'contextual_AN', 'Syntenic', 
                                                                                    'type','multiallelic'] + regional_cols
                                                                                    ).rename({'contextual_MAC':'MAC', 'contextual_AN':'AN'})

MAF_performance_variant_df = (MAF_performance_variant_df.filter(~pl.col('MAC').is_null())
                                                        .filter(~(pl.col('MAC') == 0), ~(pl.col('AN')==0))
                                                        .with_columns(MAF = ((pl.col('MAC') / pl.col('AN'))).cast(pl.Float32())))
MAF_performance_variant_df = MAF_performance_variant_df.with_columns(rounded_MAF = pl.col('MAF').cut(r2_bins[1:-1], labels=r2_bin_names))
MAF_performance_variant_df = MAF_performance_variant_df.with_columns(MAF = pl.col('MAF') * 100)

MAF_performance_variant_df.write_parquet(f'{summary_statistics_folder}/MAF_performance_variants.parquet')


print(f'Calculating MAF binned statistics ({print_runtime()})')

def groupby_region(df, colname):
    """
    Group by region for MAF binning analysis
    """
    if colname == 'Syntenic':
        non = 'Nonsyntenic'
    elif colname == 'multiallelic':
        non = 'biallelic'
    else:
        non = 'not_' + colname

    df = (df.rename({colname:'region'})
            .with_columns(region=pl.when(pl.col('region'))
                                .then(pl.lit(colname)
                                .cast(catSyn))
                                .otherwise(pl.lit(non).cast(catSyn))))
    binned_df = (df.group_by(['type','rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source', 'region'])
                                            .sum()
                                            .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                        gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                        MAF=(pl.col('MAC')/pl.col('AN')) * 100))
    binned_df = pl.concat([binned_df, df.group_by(['rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source', 'region'])
                                            .sum()
                                            .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                        gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                        MAF=((pl.col('MAC')/pl.col('AN')) * 100),
                                                        type=pl.lit('SNPs + Indels').cast(catVarType))
                                            .select(binned_df.columns)])
    return binned_df

regions = [r for r in regional_bedfiles['GRCh38'].keys()]+['Syntenic','multiallelic']

MAF_bins = (MAF_performance_variant_df.drop(['variant_id','MAF']).group_by(['type','rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source'])
                                        .sum()
                                        .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                    gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                    MAF=(pl.col('MAC')/pl.col('AN')) * 100))
MAF_bins = pl.concat([MAF_bins, MAF_performance_variant_df.drop(['variant_id','MAF']).group_by(['rounded_MAF', 'genome','method_of_phasing', 'ground_truth_data_source'])
                                        .sum()
                                        .with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked')*100,
                                                    gt_error_rate=pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                    MAF=((pl.col('MAC')/pl.col('AN')) * 100),
                                                    type=pl.lit('SNPs + Indels').cast(catVarType))
                                        .select(MAF_bins.columns)]).with_columns(region=pl.lit('All').cast(catSyn)).drop(regions, strict=False)

MAF_bins = pl.concat([MAF_bins] +
                     [groupby_region(MAF_performance_variant_df.drop(['variant_id','MAF']), region).drop(regions, strict=False).select(MAF_bins.columns) for region in regions])

MAF_bins.write_parquet(f'{summary_statistics_folder}/binned_maf_data.parquet')

MAF_bins = pl.read_parquet(f'{summary_statistics_folder}/binned_maf_data.parquet')
MAF_performance_variant_df= pl.read_parquet(f'{summary_statistics_folder}/MAF_performance_variants.parquet')
### Annotate variants with cytoband region
chm13_fai = 'resources/chm13v2.0.fa.gz.fai'
grch38_fai = 'resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai'
chm13_cytobands = 'resources/chm13v2.0_cytobands_allchrs.w_header.bed'
grch38_cytobands = 'resources/grch38_cytobands_allchrs.w_header.bed'
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

def create_end_interval_dict(fais=fais, window_size=int(1e6)):
    ends = {k:read_fai(fai) for k, fai in fais.items()}
    for genome, contigs in ends.items():
        for contig, end in contigs.items():
            ends[genome][contig] = pd.Interval(end-end%window_size, end)
    return ends

def create_end_intervals(fais=fais, window_size=int(1e6)):
    ends = {k:read_fai(fai) for k, fai in fais.items()}
    interval_df = list()
    for genome, contigs in ends.items():
        for contig, end in contigs.items():
            for i, window in enumerate(list(range(0, end, window_size))):
                if end - window < window_size:
                    window_end=end
                else:
                    window_end=window+window_size
                interval_df.append({'genome':genome,'chrom':contig, 'start':window, 'end':window_end, 'name':f'{contig}:{window}-{window_end}'})
    ends = pd.DataFrame(interval_df)
            # ends[genome][contig] = list(range(0, end, window_size)) + [end]
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

def add_intervals_to_variants(variants, intervals):
    for col in ['chrom','start','end','name','genome']:
        assert col in intervals.columns
    intervals = pl.from_pandas(intervals).filter(pl.col('chrom').is_in(variants.schema['chrom'].categories.to_list())
                                        ).with_columns(chrom=pl.col('chrom').cast(variants.schema['chrom']),
                                                        genome=pl.col('genome').cast(variants.schema['genome']),
                                                        start=pl.col('start').cast(variants.schema['position']),
                                                        end=pl.col('end').cast(variants.schema['position'])
                                         ).drop_nulls()
    interval_names = {g: dict() for g in intervals['genome'].unique()}
    interval_breaks = {g: dict() for g in intervals['genome'].unique()}
    ends = intervals.group_by(['chrom','genome']).agg(pl.col('start'), pl.col('name'))

    for d in ends.to_dicts():
        interval_breaks[d['genome']][d['chrom']] = d['start'][1:]
        interval_names[d['genome']][d['chrom']] = d['name']
    
    def get_interval_names(df, interval_names, interval_breaks):
        genome = df['genome'][0]
        chrom = df['chrom'][0]
        names = interval_names[genome][chrom]
        breaks = interval_breaks[genome][chrom]
        return {'breaks':breaks, 'labels':names}

    with pl.StringCache():
        if 'name' in intervals.columns:
            intervals = intervals.with_columns(name=pl.col('name').cast(pl.Categorical()))
        
        variants = (variants.group_by('genome','chrom')
                            .map_groups(lambda df: df.with_columns(interval=pl.col('position').cut(**get_interval_names(df, 
                                                                                                                        interval_names, 
                                                                                                                        interval_breaks)))
                            ).join(intervals.rename({'name':'interval'}), on=['genome','chrom','interval'], how='left')
                    )
    return variants

def variants_to_compressed_ideogram_data(variants, fais=fais, window_size=int(1e6), intervals=None):
    assert ('position' in variants.columns)
    print ('making intervals')
    # variants = variants.sort(['genome','chrom','position','method_of_phasing','ground_truth_data_source'])
    if intervals is None:
        intervals = create_end_intervals(fais, window_size)
        intervals = intervals.loc[intervals.chrom.isin(contigs)]
    # intervals['name'] = intervals['name'].astyp
    #     intervals = make_per_window_interval_df(variants, fais, window_size)
    # else:
    #     intervals = make_per_section_interval_df(variants, intervals)
    print (f'annotating variants with interval info ({print_runtime()})')
    variants = add_intervals_to_variants(variants, intervals)
    print (f'summing across intervals ({print_runtime()})')
    compressed_ideogram = (variants.drop('position')
                                    .group_by(['method_of_phasing','ground_truth_data_source',
                                                'genome', 'chrom', 'interval','start','end'])
                                    .sum())
    
    # print (f'dropping interval duplicates')
    # l = len(intervals)
    # intervals2 = intervals.drop_duplicates().reset_index(drop=True)
    # print (f'dropped interval duplicates (there are {l-len(intervals2)} duplicates)')
    
    # print ('merging interval subs with interval metadata')
    # compressed_ideogram = compressed_ideogram.merge(intervals[['genome','chrom','interval_names','start','end']], on=['genome','chrom','interval_names'], how='left')
    # print('prepping interval sums for export')
    # compressed_ideogram = (compressed_ideogram.drop(columns='intervals')
    #                                             .set_index(['genome','method_of_phasing','ground_truth_data_source','chrom']))
    return compressed_ideogram


def sum_by_interval(df, intervals, keys, return_lazy=True, slop=int(0)):
    result = list()
    df = df.lazy()
    for chrom, start, end, name in intervals:
        start = max(0, start-slop)
        end = end+slop
        result.append(df.filter(pl.col('chrom')==chrom, pl.col('position')>=start, pl.col('position')<end)
                        .drop(['chrom','position'])
                        .sum()
                        .with_columns(pl.lit(name).alias('Syndrome'), 
                                      pl.lit(start).cast(schemas['variant_switch']['position']).alias('start'),
                                      pl.lit(end).cast(schemas['variant_switch']['position']).alias('end'),
                                      **keys))
    if return_lazy:
        return result
    else:
        return [pl.concat(result).collect()]

def add_zero_start_position(inDF, id_cols=['chrom','ground_truth_data_source','genome','method_of_phasing']):
    out = list()
    fill_backwards=[cs.numeric().fill_null(0), (~cs.numeric()).backward_fill()]
    for _, df in inDF.group_by(id_cols):
        df = df.sort('position')
        if df.head(1).select('position').item() > 0:
            df = df.shift()
        out.append(df.select(fill_backwards))
    return pl.concat(out).sort(id_cols+['position'])

def rolling_stats(df, window_size=10000, smoothing_factor=10, offset=None):
    if offset is None:
        offset = -window_size//2
    fill_forwards=[cs.numeric().fill_null(0), (~cs.numeric()).forward_fill()]

    df = df.with_columns(position=pl.col('position').cast(pl.Int32())
          ).group_by_dynamic(group_by=['chrom','ground_truth_data_source','genome','method_of_phasing'],
                                every=f'{window_size//smoothing_factor}i',
                                index_column=pl.col('position'), 
                                period=f'{window_size//smoothing_factor}i', 
                                offset=f'{offset}i'
          ).agg(n_switch_errors=pl.sum('n_switch_errors'),
                n_checked=pl.sum('n_checked'),
                n_gt_errors=pl.sum('n_gt_errors'),
                n_gt_checked=pl.sum('n_gt_checked'),
                num_variants=pl.len())
    df = (add_zero_start_position(df)
            .upsample(every=f'{window_size//smoothing_factor}i', 
                        time_column='position',
                        group_by=['chrom','ground_truth_data_source','genome','method_of_phasing'])
              .select(fill_forwards)
              .rolling(group_by=['chrom','ground_truth_data_source','genome','method_of_phasing'],
                       index_column=pl.col('position'), 
                       period=f'{window_size}i', 
                       offset=f'{offset}i')
                  .agg(n_switch_errors=pl.sum('n_switch_errors'),
                       n_checked=pl.sum('n_checked'),
                       n_gt_errors=pl.sum('n_gt_errors'),
                       n_gt_checked=pl.sum('n_gt_checked'),
                       variants_per_1k=pl.mean('num_variants')/(window_size//1000),
                       )
               .rename({'position':'start'})
               .with_columns(end=pl.col('start')+window_size))
    return df

def make_r_compatible(df: pl.DataFrame | pl.LazyFrame):
    # Cast all categorical/enum columns to plain strings to avoid Arrow dictionary arrays
    return df.with_columns(cs.by_dtype(pl.Categorical, pl.Enum).cast(pl.Utf8()))

print (f'Gathering per-cnv stats ({print_runtime()})')

decipher_cnvs = pd.read_csv('resources/decipher_syndromes.txt', sep='\t')
decipher_cnvs.loc[decipher_cnvs.end_grch38 < decipher_cnvs.start_grch38, ['start_grch38','end_grch38']] = decipher_cnvs.loc[decipher_cnvs.end_grch38 < decipher_cnvs.start_grch38, ['end_grch38','start_grch38']].values
decipher_cnvs.loc[decipher_cnvs.end_chm13 < decipher_cnvs.start_chm13, ['start_chm13','end_chm13']] = decipher_cnvs.loc[decipher_cnvs.end_chm13 < decipher_cnvs.start_chm13, ['end_chm13','start_chm13']].values

chm13_intervals = list(decipher_cnvs[['chrom','start_chm13','end_chm13', 'Syndrome']].itertuples(index=False, name=None))
grch38_intervals = list(decipher_cnvs[['chrom','start_grch38','end_grch38', 'Syndrome']].itertuples(index=False, name=None))
interval_dict = {'CHM13v2.0':chm13_intervals, 'GRCh38':grch38_intervals}

var_info = pl.col('variant_id').str.split('_')
MAF_performance_variant_df = MAF_performance_variant_df.with_columns(chrom = var_info.list[0].cast(catChrom),
                                                                     position = var_info.list[1].cast(all_variant_annotations_dtypes['POS']))

ideogram_variants = MAF_performance_variant_df.drop(['variant_id','Syntenic','type','rounded_MAF'] + regional_cols)

per_genome_cnv_regions = list()
per_genome_cnv_regions_with_slop = list()
slop=int(1e6)
for (genome, method_of_phasing, ground_truth_data_source), df in ideogram_variants.group_by(['genome','method_of_phasing','ground_truth_data_source']):
    keys = {'genome':pl.lit(genome),
            'method_of_phasing':pl.lit(method_of_phasing),
            'ground_truth_data_source':pl.lit(ground_truth_data_source)}
    per_genome_cnv_regions.extend(sum_by_interval(df, interval_dict[genome], keys=keys, return_lazy=True))
    per_genome_cnv_regions_with_slop.extend(sum_by_interval(df, interval_dict[genome], keys=keys, return_lazy=True, slop=slop))
per_genome_cnv_regions = pl.concat(per_genome_cnv_regions
                                  ).collect(
                                  ).with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked') * 100)
per_genome_cnv_regions_with_slop = pl.concat(per_genome_cnv_regions_with_slop
                                  ).collect(
                                  ).with_columns(switch_error_rate=pl.col('n_switch_errors')/pl.col('n_checked') * 100)
per_genome_cnv_regions.write_parquet(f'{summary_statistics_folder}/per_genome_cnv_regions.parquet')
per_genome_cnv_regions_with_slop.write_parquet(f'{summary_statistics_folder}/per_genome_cnv_regions_with_slop.parquet')


cytoband_coords = get_cytobands_in_bed_format(chm13_cytobands)
cytoband_coords['genome'] = 'CHM13v2.0'
cytoband_coords2 = get_cytobands_in_bed_format(grch38_cytobands)
cytoband_coords2['genome'] = 'GRCh38'
cytoband_coords = pd.concat((cytoband_coords, cytoband_coords2))

print (f'Calculating per cytoband stats ({print_runtime()})')
all_cytobands = variants_to_compressed_ideogram_data(ideogram_variants, intervals=cytoband_coords)
all_cytobands = make_r_compatible(all_cytobands)
all_cytobands.write_parquet(f'{summary_statistics_folder}/per_cytoband_variant_data.parquet')

print (f'Calculating per 100,000 stats ({print_runtime()})')
compressed_ideogram_100k = variants_to_compressed_ideogram_data(ideogram_variants, window_size=100000)
compressed_ideogram_100k = make_r_compatible(compressed_ideogram_100k)
compressed_ideogram_100k.write_parquet(f'{summary_statistics_folder}/compressed_ideogram_100000_window.parquet')

print (f'Calculating per 1,000,000 stats ({print_runtime()})')
compressed_ideogram_1m = variants_to_compressed_ideogram_data(ideogram_variants, window_size=1000000)
compressed_ideogram_1m = make_r_compatible(compressed_ideogram_1m)
compressed_ideogram_1m.write_parquet(f'{summary_statistics_folder}/compressed_ideogram_1000000_window.parquet')

print (f'Calculating per 250,000 stats ({print_runtime()})')
rolling_stats_250k = rolling_stats(ideogram_variants, window_size=250000, smoothing_factor=25)
rolling_stats_250k = make_r_compatible(rolling_stats_250k)
rolling_stats_250k.write_parquet(f'{summary_statistics_folder}/rolling_stats_250k_window.parquet')

print (f'Calculating per 500,000 stats ({print_runtime()})')
rolling_stats_500k = rolling_stats(ideogram_variants, window_size=500000, smoothing_factor=50)
rolling_stats_500k = make_r_compatible(rolling_stats_500k)
rolling_stats_500k.write_parquet(f'{summary_statistics_folder}/rolling_stats_500k_window.parquet')

print (f'Calculating per 1,000,000 stats ({print_runtime()})')
rolling_stats_1m = rolling_stats(ideogram_variants, window_size=1000000, smoothing_factor=100)
rolling_stats_1m = make_r_compatible(rolling_stats_1m)
rolling_stats_1m.write_parquet(f'{summary_statistics_folder}/rolling_stats_1m_window.parquet')



### gather chrom and sample specific data
del MAF_performance_variant_df
del total_phasing_stats_datasets['variant_switch']
del ideogram_variants

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
all_chroms = all_chroms.drop([x for x in all_chroms.columns if '_right' in x]
                            ).with_columns(gt_error_rate=(pl.col('n_gt_errors')/pl.col('n_gt_checked')) * 100,
                                           flip_error_rate=(pl.col('n_flips')/(pl.col('n_total_hets')))*100,
                                           switch_error_rate=(pl.col('n_switch_errors')/pl.col('n_checked'))*100,
                                           flip_and_switches_SER=(pl.col('n_total_switch_errors')/pl.col('n_total_hets')) * 100,
                                           accurately_phased_rate=(pl.col('num_correct_phased_hets')/pl.col('n_total_hets'))*100,
                                           true_switch_error_rate=(pl.col('n_true_switch_errors')/pl.col('n_total_hets')) * 100,
                                           trio_phased=(pl.col('sample_id').is_in(trio_phased_samples)))

all_samples = group_phasing_statistics(all_chroms, ['sample_id','genome','population','superpopulation','sex','trio_phased'] + expt_columns).drop('chrom')


HGSVC_samples_nontrios_only = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='HGSVC_samples_nontrios_only')).select('sample_id').unique().to_numpy().flatten()
HGSVC_samples = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='HGSVC_samples')).select('sample_id').unique().to_numpy().flatten()
HPRC_samples = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='HPRC_samples')).select('sample_id').unique().to_numpy().flatten()
trio_samples = all_samples.filter((~pl.col('switch_error_rate').is_nan()), (pl.col('ground_truth_data_source')=='trios')).select('sample_id').unique().to_numpy().flatten()

all_chroms = all_chroms#.filter(((pl.col('ground_truth_data_source') == 'trios') & (pl.col('sample_id').is_in(trio_samples))) | ((pl.col('ground_truth_data_source') == 'HPRC_samples') & (pl.col('sample_id').is_in(HPRC_samples))) | ((pl.col('ground_truth_data_source') == 'HGSVC_samples') & (pl.col('sample_id').is_in(HGSVC_samples))) | ((pl.col('ground_truth_data_source') == 'HGSVC_samples_nontrios_only') & (pl.col('sample_id').is_in(HGSVC_samples_nontrios_only))))
all_samples = all_samples#.filter(((pl.col('ground_truth_data_source') == 'trios') & (pl.col('sample_id').is_in(trio_samples))) | ((pl.col('ground_truth_data_source') == 'HPRC_samples') & (pl.col('sample_id').is_in(HPRC_samples))) | ((pl.col('ground_truth_data_source') == 'HGSVC_samples') & (pl.col('sample_id').is_in(HGSVC_samples))) | ((pl.col('ground_truth_data_source') == 'HGSVC_samples_nontrios_only') & (pl.col('sample_id').is_in(HGSVC_samples_nontrios_only))))
males = all_samples.filter(pl.col('sex')=='Male').select('sample_id').to_numpy().flatten()



###############
###Calc N50/L50
print ("calculating per-contig and per-sample N50/L50 - can take a minute")
print (f"organizing switch error block data ({print_runtime()})")
chrom_ends=dict()
chrom_length=dict()
chrom_starts=dict()

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
chrom_length['CHM13v2.0']['PAR1'] = chrom_ends['CHM13v2.0']['PAR1'] - 1
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
                                                     end_chrom=pl.col('end_chrom').cast(pl.Int32()),
                                                     chrom_len=pl.col('chrom_len').cast(pl.Int32()))
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
                                        .cast(pl.Int64)              # Ensure signed integer type for list operations
                                        .sort(descending=True),      # Sort from largest to smallest block. 
                            L50=pl.col('switch_error_site')              # L50: How many blocks (largest to smallest) needed to get to 50%+1 of contig length
                                  .diff()                                # First, repeat code to get blocksize. Polars should see that the code is repeating and avoid duplicate calculations.
                                  .drop_nulls()                          #
                                  .sort(descending=True)                 # With per-contig list of block lengths, sorted largest to smallest...
                                  .cast(pl.Int64)
                                  .cum_sum()                             # take cumulative sum to convert to list of cumulative block sizes
                                  .cast(pl.Int64)                        # Cast cumsum result before search_sorted
                                  .search_sorted(pl.col('chrom_len')     # L50 = the number of contigs necessary to get to contig length/2
                                                   .first()/2))          # This code finds the index you would insert (contig length/2) to keep the list of cumsum lengths sorted. To do this, you would find the number of blocks necessary to get to (contig_length/2), and insert (contig_len/2) afterwards, at the location of the next block. So the insertion index is the number of blocks necessary to get to the first block whose cumulative sum is over contig_length/50 - i.e., the L50.
                       .with_columns(N50=pl.col('blocksize')             # The N50 is just the size of the block needed to get to 50%+1 of the contig_length
                                           .list.get(pl.col('L50').list.first())))     # As we've just shown, given a list of blocks sorted from biggest to smallest, that block is located at the L50's index.
print (f"Calculating per sample N50 ({print_runtime()})")
per_sample_N50 = (per_contig_N50.group_by(sample_id_columns)             # Group data by sample.
                                .agg(blocksize=pl.col('blocksize')       # blocksize becomes list of lists, 
                                                 .flatten()              # flatten to get one list of contig blocks in sample
                                                 .sort(descending=True), # Then sort largest to smallest
                                     L50=pl.col('blocksize')             # Do the same thing as above to get the sample L50
                                           .flatten()
                                           .sort(descending=True)
                                           .cast(pl.Int64)
                                           .cum_sum()
                                           .cast(pl.Int64)                   # Cast cumsum result before search_sorted
                                           .search_sorted(pl.col('chrom_len') # This time sum chrom lengths to get genome length
                                                            .unique()
                                                            .sum()/2))
                                .with_columns(N50=pl.col('blocksize')
                                                    .list.get(pl.col('L50').list.first()))
                                .drop('blocksize'))

per_contig_N50 = (per_contig_N50.drop(['blocksize','chrom_len']))

all_chroms = all_chroms.join(per_contig_N50, on=['genome','method_of_phasing','ground_truth_data_source','sample_id','chrom'], how='left')
all_samples = all_samples.join(per_sample_N50, on=['genome','method_of_phasing','ground_truth_data_source','sample_id'], how='left')

all_ancestries = group_phasing_statistics(all_chroms, ['genome','population','superpopulation'] + expt_columns).drop(['sex','sample_id'])
all_methods = group_phasing_statistics(all_samples, ['genome'] + expt_columns).drop(['sample_id','sex','population','superpopulation'])



print (f'saving summary stats to disk ({print_runtime()})')
all_chroms.write_parquet(f'{summary_statistics_folder}/chroms.parquet')
all_samples.write_parquet(f'{summary_statistics_folder}/samples.parquet')
all_ancestries.write_parquet(f'{summary_statistics_folder}/ancestries.parquet')
all_methods.write_parquet(f'{summary_statistics_folder}/methods.parquet')

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
                        base_report_name = f'{imputation_statistics_folder}/{genome}{variant_set}.{ground_truth}.{panel}.glimpse2_concordance'
                    else:
                        base_report_name = f'{imputation_statistics_folder}/ancestry_specific/{genome}{variant_set}.{ground_truth}.{panel}.{ancestry}.glimpse2_concordance'
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
                            if ("syntenic" in bin_grouping) or ('Syntenic' in bin_grouping):
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

print (f'Calculating variant filtering statistics {print_runtime()}')

grch38_ids = pl.scan_parquet(f"{summary_statistics_folder}/MAF_performance_variants.parquet").select(['variant_id','genome']).filter(pl.col('genome')=='GRCh38').select('variant_id').unique().collect().to_series().implode()
chm13_ids = pl.scan_parquet(f"{summary_statistics_folder}/MAF_performance_variants.parquet").select(['variant_id','genome']).filter(pl.col('genome')=='CHM13v2.0').select('variant_id').unique().collect().to_series().implode()
SV_cutoff = 50
summary = pl.scan_parquet(f"{summary_statistics_folder}/bcftools_query_variant_data.parquet"
   ).with_columns(MAC=pl.min_horizontal(pl.col('AC_original_panel'), pl.col('AN_original_panel')-pl.col('AC_original_panel')),
   ).with_columns(VQSLOD_filter=(pl.col('VQSLOD')<0).cast(pl.Boolean),
                  MERR_filter=(pl.col('MERR')>pl.col('AN_original_panel')*0.05).cast(pl.Boolean),
                  HWE_pop_filter=((pl.col('HWE_EUR')<1e-10) & (pl.col('HWE_AFR')<1e-10) & (pl.col('HWE_EAS')<1e-10) & (pl.col('HWE_AMR')<1e-10) & (pl.col('HWE_SAS')<1e-10)).cast(pl.Boolean),
                  MAC_filter=(pl.col('MAC')==0).cast(pl.Boolean),
                  AC_filter=(pl.col('AC_original_panel')<=1),
                  f_missing_filter=(pl.col('F_MISSING')>0.05).cast(pl.Boolean),
                  var_len_filter=((pl.col('ID').str.split('_').list[2].str.len_chars().cast(pl.Int64)
                                   - 
                                   pl.col('ID').str.split('_').list[3].str.len_chars().cast(pl.Int64)
                                  ).abs() >= SV_cutoff-1).cast(pl.Boolean),
                                  # minus 1 because bcftools ILEN does not count the anchor base
                  alt_star_filter=(pl.col('ALT') == '*').cast(pl.Boolean),
                  singleton=(pl.col('MAC')==1).cast(pl.Boolean),
                  pass_filter = (pl.col('FILTER') != 'PASS').cast(pl.Boolean)
   ).select(['ID','genome','Syntenic','singleton', 'VQSLOD_filter','MERR_filter','HWE_pop_filter',
               'MAC_filter','AC_filter','f_missing_filter','var_len_filter','alt_star_filter','pass_filter']
   ).with_columns(GRCh38_filtered= ~pl.col('ID').is_in(grch38_ids),
                   CHM13_filtered = ~pl.col('ID').is_in(chm13_ids),
                   GRCh38_criteria_fail=(pl.col('AC_filter')|
                                            pl.col('f_missing_filter')|pl.col('pass_filter')|pl.col('HWE_pop_filter')|pl.col('MERR_filter')|pl.col('var_len_filter')|pl.col('alt_star_filter')).cast(pl.Boolean),
                   CHM13_criteria_fail =(pl.col('MAC_filter')|pl.col('VQSLOD_filter')|#pl.col('neg_train_site_filter')|
                                            pl.col('f_missing_filter')|pl.col('pass_filter')|pl.col('HWE_pop_filter')|pl.col('MERR_filter')|pl.col('var_len_filter')|pl.col('alt_star_filter')).cast(pl.Boolean)
   ).group_by(['genome','Syntenic','singleton', 'VQSLOD_filter','MERR_filter','HWE_pop_filter','MAC_filter','AC_filter',
                     'f_missing_filter','var_len_filter','alt_star_filter','pass_filter','CHM13_filtered','GRCh38_filtered','GRCh38_criteria_fail','CHM13_criteria_fail']
   ).len().collect()

summary.write_parquet(f"{summary_statistics_folder}/filter_summary_stats.parquet")
print ('Done!')
print (f"Total runtime: {print_runtime()}")
