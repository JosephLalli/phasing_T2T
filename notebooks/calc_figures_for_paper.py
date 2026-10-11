# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import polars as plnotebooks
from pandas.api.types import CategoricalDtype
import seaborn as sns
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import importlib

plt.rcParams["figure.dpi"] = 120


summary_statistics_folder = 'intermediate_data'
os.chdir('..') # Change to the parent directory

# %%
binned_maf_data = pd.read_parquet(f"{summary_statistics_folder}/binned_maf_data.parquet")

# compressed_ideogram = pd.read_parquet(f'{summary_statistics_folder}/compressed_ideogram_100000_window.parquet')
all_chroms = pd.read_parquet(f'{summary_statistics_folder}/chroms.parquet')
all_samples = pd.read_parquet(f'{summary_statistics_folder}/samples.parquet')
all_ancestries = pd.read_parquet(f'{summary_statistics_folder}/ancestries.parquet')
all_methods = pd.read_parquet(f'{summary_statistics_folder}/methods.parquet')

try:
    per_variant_category_imputation_performance = pd.read_parquet(f'{summary_statistics_folder}/per_variant_category_imputation_performance.parquet')
except FileNotFoundError:
    per_variant_category_imputation_performance = None
    print("WARNING: per_variant_category_imputation_performance.parquet not found (imputation data unavailable)")

try:
    per_sample_imputation_performance = pd.read_parquet(f'{summary_statistics_folder}/per_sample_imputation_performance.parquet')
except FileNotFoundError:
    per_sample_imputation_performance = None
    print("WARNING: per_sample_imputation_performance.parquet not found (imputation data unavailable)")

rolling_window_ideogram_250k = pd.read_parquet(f'{summary_statistics_folder}/rolling_stats_250k_window.parquet')

# %%
binned_maf_data

# %%
#Post hoc fixes - one off events that are tricky to deal with when compiling stats

# there is one variant (I have not investigated, but likely a denovo SNP) that is a singleton variant in the HPRC dataset. 
# That variant passed all filters, but makes for a very odd figure (correctly phased: 0% error rate, incorrect: 100% error rate)
# We will trim it here
binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')&
                    (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')&
                    (binned_maf_data.rounded_MAF == 'singleton'), ['n_switch_errors','n_checked','switch_error_rate']] = np.nan

# Isolate syntenic vs nonsyntenic regions
regional_binned_maf_data = binned_maf_data.copy()
binned_maf_data = (binned_maf_data.loc[binned_maf_data.region.isin(['syntenic','nonsyntenic','All'])]
                                  .replace('syntenic','Syntenic')
                                  .replace('nonsyntenic','Nonsyntenic')
                                  .rename(columns={'region':'Syntenic'}))

# Nonsyntenic GRCh38 regions are simply GRCh38 errors, and we are not interested in the performance in these regions
binned_maf_data = binned_maf_data.loc[~((binned_maf_data.genome=='GRCh38')&(binned_maf_data.Syntenic=='Nonsyntenic'))]

# NA20355 is a duo member, and as such was not part of duohmm trio correction for the GRCh38 2022 effort.
all_chroms.loc[(all_chroms.genome=='GRCh38')&(all_chroms['sample_id']=='NA20355'), 'trio_phased'] = False
all_samples.loc[(all_samples.genome=='GRCh38')&(all_samples['sample_id']=='NA20355'), 'trio_phased'] = False
all_chroms = all_chroms.loc[~((all_chroms.genome=='GRCh38')&(all_chroms['sample_id']=='NA20355')&(all_chroms.ground_truth_data_source.isin(['HPRC_HGSVC_probands','HGSVC_probands'])))]
all_samples = all_samples.loc[~((all_samples.genome=='GRCh38')&(all_samples['sample_id']=='NA20355')&
                                (all_samples.ground_truth_data_source.isin(['HPRC_HGSVC_probands','HGSVC_probands'])))]

# very high variance stats when n_checked is less than 10. Replace switch error rate at these sites with np.nan
binned_maf_data.loc[(binned_maf_data.n_checked>10), 'switch_error_rate'] = binned_maf_data.loc[(binned_maf_data.n_checked>10), 'switch_error_rate'].replace(0, np.nan)

# bins with a switch error rate of 0 make for insane figures when plotted on a log scale. Replace with np.nan, and maybe note
binned_maf_data['switch_error_rate'] = binned_maf_data['switch_error_rate'].replace(0, np.nan)

# A reviewer suggests combining chrX, PAR1, and PAR2 statistics. I will do that, and set male chrX switch error rates to np.nan.
## Calculate percent discordance in imputation performance
if per_variant_category_imputation_performance is not None and 'perc_discordance' not in per_variant_category_imputation_performance.columns:
    per_variant_category_imputation_performance = per_variant_category_imputation_performance.loc[per_variant_category_imputation_performance.num_AA > 0]
    if max(per_variant_category_imputation_performance.mean_AF) < 1:
        per_variant_category_imputation_performance.mean_AF *= 100
    per_variant_category_imputation_performance['num_alt_mismatches'] = per_variant_category_imputation_performance.num_Aa_mismatches + per_variant_category_imputation_performance.num_aa_mismatches
    per_variant_category_imputation_performance['num_alt_variants'] = per_variant_category_imputation_performance.num_Aa + per_variant_category_imputation_performance.num_aa
    per_variant_category_imputation_performance['perc_discordance'] = per_variant_category_imputation_performance['num_alt_mismatches']/per_variant_category_imputation_performance['num_alt_variants']*100
    per_variant_category_imputation_performance['perc_concordance'] = 100-per_variant_category_imputation_performance['perc_discordance']

# imputed df needs to have CHM13v2.0 as genome name instead of T2T
if per_variant_category_imputation_performance is not None:
    per_variant_category_imputation_performance = per_variant_category_imputation_performance.replace('T2T','CHM13v2.0')
if per_sample_imputation_performance is not None:
    per_sample_imputation_performance = per_sample_imputation_performance.replace('T2T','CHM13v2.0')

# A supplemental calculates SER by % of chromosome that is new. we calc that here.
if 'percent_new' not in all_chroms.columns:
    syntenic_bed = pd.read_csv('resources/hg38.GCA_009914755.4.synNet.summary.bed.gz', sep='\t', header=None)
    syntenic_bed.columns = ['chrom','start','end','status']
    starts = syntenic_bed.groupby('chrom').first().reset_index()
    syntenic_bed['length'] = syntenic_bed.end-syntenic_bed.start
    syntenic_bed['dif'] = syntenic_bed.end.diff()
    syntenic_bed['nonsyn_len'] = syntenic_bed.dif - syntenic_bed.length
    syntenic_bed.loc[syntenic_bed.nonsyn_len < 0, 'nonsyn_len'] = syntenic_bed.loc[syntenic_bed.nonsyn_len < 0, 'start']
    syntenic_bed.loc[syntenic_bed.nonsyn_len.isna(), 'nonsyn_len'] = syntenic_bed.loc[syntenic_bed.nonsyn_len.isna(), 'start']
    rough_length = syntenic_bed.groupby('chrom').end.max().reset_index()
    chr_nonsyn_lengths = syntenic_bed.groupby('chrom').sum(numeric_only=True).reset_index()[['chrom','nonsyn_len']]
    chr_nonsyn_lengths = chr_nonsyn_lengths.merge(rough_length, on='chrom')
    chr_nonsyn_lengths['percent_new'] = chr_nonsyn_lengths.nonsyn_len/chr_nonsyn_lengths.end * 100
    all_chroms = all_chroms.merge(chr_nonsyn_lengths[['chrom', 'percent_new']], on=['chrom'],how='left')

all_chroms_chrX = all_chroms.loc[all_chroms.chrom.isin(('chrX','PAR1','PAR2'))].assign(chrom='chrX')
all_chroms_chrX = all_chroms_chrX.groupby(['chrom','genome','sample_id','method_of_phasing','ground_truth_data_source','population', 'superpopulation', 'sex', 'trio_phased'],observed=True).sum(numeric_only=True).reset_index()
all_chroms_chrX.loc[all_chroms_chrX.sex=='Male', ['n_switch_errors','n_total_switch_errors','n_flips','n_consecutive_flips','n_true_switch_errors','num_correct_phased_hets',
                                                  'n_total_hets','flip_and_switches_SER', 'flip_error_rate', 'true_switch_error_rate','accurately_phased_rate']] = np.nan

all_chroms_chrX['switch_error_rate'] = all_chroms_chrX.n_switch_errors/all_chroms_chrX.n_checked * 100
all_chroms_chrX['gt_error_rate'] = all_chroms_chrX.n_gt_errors/all_chroms_chrX.n_gt_checked * 100
all_chroms_chrX['flip_and_switches_SER'] = all_chroms_chrX.n_total_switch_errors/all_chroms_chrX.n_total_hets * 100
all_chroms_chrX['flip_error_rate'] = all_chroms_chrX.n_flips/all_chroms_chrX.n_total_hets * 100
all_chroms_chrX['true_switch_error_rate'] = all_chroms_chrX.n_true_switch_errors/all_chroms_chrX.n_total_hets * 100

all_chroms_chrX = pd.concat((all_chroms_chrX, all_chroms.loc[~all_chroms.chrom.isin(('chrX','PAR1','PAR2'))])).drop(columns=['L50','N50','percent_new'])


binned_maf_data = binned_maf_data.sort_values('type', key=lambda col: col.apply(lambda x: {'Indel':0, 'SNP':1, 'SNPs + Indels':2}.get(x,x)))
binned_maf_data = binned_maf_data.sort_values('genome', key=lambda col: col.apply(lambda x: {'GRCh38':0, 'CHM13v2.0':1}.get(x,x)))

if binned_maf_data['switch_error_rate'].max() < 1:
    binned_maf_data['switch_error_rate'] *= 100
if binned_maf_data['gt_error_rate'].max() < 1:
    binned_maf_data['gt_error_rate'] *= 100
if binned_maf_data['MAF'].max() < 1:
    binned_maf_data['MAF'] *= 100

# a handful of variants are miscategorized as singletons. There are no true GRCh38 singletons. We can fix this.
binned_maf_data = binned_maf_data.loc[~((binned_maf_data.rounded_MAF=='singleton')&(binned_maf_data.genome=='GRCh38')&(binned_maf_data.method_of_phasing=='phased_with_parents_and_pedigree')&(binned_maf_data.ground_truth_data_source!='trios'))]



# %%
oneKGP_meta = pd.read_csv('resources/1000_genomes_meta.tsv', sep='\t')
overlaps = oneKGP_meta.loc[oneKGP_meta['Sample name'].isin(all_samples.sample_id.unique()) & oneKGP_meta['Data collections'].str.contains('Simons Genome'), 'Sample name']
overlaps.to_csv('resources/sample_subsets/SGDP_samples_in_1KGP.txt', index=False, header=False)

# %%
pedigree = pd.read_csv('resources/pedigrees/1kGP.3202_samples.pedigree_info.txt', sep=' ')

# %% [markdown]
# ### Variant Filtering Stats 
# Stats in paragraph after "Full length, telomere to telomere maps of human genetic recombination rates are largely consistent with deCODE recombination rates"

# %%
sum_stats_var_filters = pd.read_parquet(f"{summary_statistics_folder}/filter_summary_stats.parquet")

grch38_filter_columns=['MERR_filter', 'HWE_pop_filter','AC_filter',
       'f_missing_filter', 'var_len_filter', 'alt_star_filter', 'pass_filter']
chm13_filter_columns=[#'neg_train_site_filter',
       'VQSLOD_filter', 'MERR_filter', 'HWE_pop_filter', 'MAC_filter',
       'f_missing_filter', 'var_len_filter', 'alt_star_filter', 'pass_filter']

sum_stats_var_filters['gatk_filters'] = sum_stats_var_filters.VQSLOD_filter #| sum_stats_var_filters.neg_train_site_filter
sum_stats_var_filters['GRCh38_criteria_pass_if_singletons_passed'] = ~sum_stats_var_filters[[x for x in grch38_filter_columns if x != 'AC_filter'] + ['MAC_filter']].any(axis=1)

num_singletons_gained = sum_stats_var_filters.loc[ (sum_stats_var_filters.genome=='CHM13v2.0') &
                                                   ~sum_stats_var_filters.CHM13_criteria_fail &
                                                    sum_stats_var_filters.singleton &
                                                    sum_stats_var_filters.GRCh38_criteria_fail &
                                                    sum_stats_var_filters.GRCh38_criteria_pass_if_singletons_passed].len.sum()
num_variants_lost_to_gatk = sum_stats_var_filters.loc[ (sum_stats_var_filters.genome=='CHM13v2.0') &
                                                        sum_stats_var_filters.CHM13_criteria_fail &
                                                       ~sum_stats_var_filters.GRCh38_criteria_fail &
                                                        sum_stats_var_filters.gatk_filters &
                                                        sum_stats_var_filters.GRCh38_criteria_pass_if_singletons_passed].len.sum()

print(f"CHM13 variants lost to gatk-based filters:\n{num_variants_lost_to_gatk}")
print(f"CHM13 variants gained due to ability to include singletons:\n{num_singletons_gained}")

# %% [markdown]
# The resulting panel contains 6404 haplotypes and XXX unique variants, including XXX phased variants in novel regions of the CHM13v2.0 genome.

# %%
number_unique_CHM13_variants = sum_stats_var_filters.loc[~sum_stats_var_filters.CHM13_criteria_fail & (sum_stats_var_filters.genome=='CHM13v2.0')].len.sum()
number_unique_syntenic_CHM13_variants = sum_stats_var_filters.loc[~sum_stats_var_filters.Syntenic & ~sum_stats_var_filters.CHM13_criteria_fail & (sum_stats_var_filters.genome=='CHM13v2.0')].len.sum()
print(f"Unique variants: {number_unique_CHM13_variants}")
print(f"Unique novel variants: {number_unique_syntenic_CHM13_variants}")

# %% [markdown]
# We observed that samples that underwent Mendelian pre-phasing (Figure 2a, c) had lower SERs than samples that underwent statistical phasing (Figure 2b, d; XXX% vs XXX%).

# %%
trio_vs_nontrio = all_samples.loc[all_samples.ground_truth_data_source.isin(['HPRC_samples', 'HGSVC_samples_nontrios_only']) & 
                                    (all_samples.method_of_phasing == 'phased_with_parents_and_pedigree')
                                    # (all_samples.genome=='GRCh38')
                                    ].groupby(['genome', 'ground_truth_data_source'], observed=True).sum(numeric_only=True)

trio_vs_nontrio['switch_error_rate'] = trio_vs_nontrio.n_switch_errors/trio_vs_nontrio.n_checked * 100
trio_vs_nontrio['gt_error_rate'] = trio_vs_nontrio.n_gt_errors/trio_vs_nontrio.n_gt_checked * 100
trio_vs_nontrio['true_switch_error_rate'] = trio_vs_nontrio.n_true_switch_errors/trio_vs_nontrio.n_total_hets * 100
trio_vs_nontrio[['switch_error_rate', 'gt_error_rate', 'true_switch_error_rate']]

# %% [markdown]
# ### Number of reference samples that are probands/parents/unrelateds

# %%
pedigree = pedigree.replace('0', np.nan)
num_that_are_children_and_parents = pedigree.loc[pedigree.sampleID.isin(pedigree.fatherID) | pedigree.sampleID.isin(pedigree.motherID) & ~(pedigree.fatherID.isna() | pedigree.motherID.isna())]

# %%
children_that_are_parents = pedigree.loc[(pedigree.sampleID.isin(pedigree.fatherID) | pedigree.sampleID.isin(pedigree.motherID)) & (~pedigree.fatherID.isna() | ~pedigree.motherID.isna())].sampleID
pedigree['trio_status'] = np.nan
complex_pedigrees = pedigree.loc[pedigree.sampleID.isin(children_that_are_parents) | pedigree.motherID.isin(children_that_are_parents) | pedigree.fatherID.isin(children_that_are_parents)]
simple_trios = pedigree.loc[~pedigree.index.isin(complex_pedigrees.index) & (~pedigree.fatherID.isna() & ~pedigree.motherID.isna())]
simple_duos = pedigree.loc[~pedigree.index.isin(complex_pedigrees.index) & (~pedigree.fatherID.isna() ^ ~pedigree.motherID.isna())]
relatives = pd.concat((simple_duos, simple_trios, complex_pedigrees))
relatives = relatives.fatherID.dropna().tolist() + relatives.motherID.dropna().tolist() + relatives.sampleID.tolist()
relatives = pedigree.loc[pedigree.sampleID.isin(relatives)]
no_relatives = pedigree.loc[~pedigree.sampleID.isin(relatives.sampleID)]


# %%
proband_or_not_trio = pedigree.loc[~(pedigree.sampleID.isin(pedigree.fatherID) | pedigree.sampleID.isin(pedigree.motherID))]


# %%
def get_parents(df):
    return df.fatherID.dropna().to_list() + df.motherID.dropna().to_list()

complex_people = [x for x in get_parents(complex_pedigrees) if x in complex_pedigrees.sampleID.tolist()]

complex_trio_parents = [x for x in get_parents(complex_pedigrees.loc[~(complex_pedigrees.fatherID.isna() | complex_pedigrees.motherID.isna())]) if x not in complex_pedigrees.sampleID.tolist()]
complex_duo_parents = [x for x in get_parents(complex_pedigrees.loc[(complex_pedigrees.fatherID.isna() ^ complex_pedigrees.motherID.isna())]) if x not in complex_pedigrees.sampleID.tolist()]
complex_trio_children = [x for x in complex_pedigrees.loc[~(complex_pedigrees.fatherID.isna() | complex_pedigrees.motherID.isna())].sampleID.to_list() if x not in get_parents(complex_pedigrees)]
complex_duo_children = [x for x in complex_pedigrees.loc[(complex_pedigrees.fatherID.isna() ^ complex_pedigrees.motherID.isna())].sampleID.to_list() if x not in get_parents(complex_pedigrees)]
simple_pedigrees = pedigree.loc[~pedigree.index.isin(complex_pedigrees.index)]
simple_trios = simple_pedigrees.loc[~(simple_pedigrees.fatherID.isna() | simple_pedigrees.motherID.isna())]
trio_parents = get_parents(simple_trios) + complex_trio_parents
trio_children = simple_trios.sampleID.to_list() + complex_trio_children
simple_duos = simple_pedigrees.loc[(simple_pedigrees.fatherID.isna() ^ simple_pedigrees.motherID.isna())]
duo_parents = get_parents(simple_duos) + complex_duo_parents
duo_children = simple_duos.sampleID.to_list() + complex_duo_children

trio_dads = simple_trios.fatherID.dropna().to_list() + [x for x in complex_trio_parents if x in complex_pedigrees.fatherID.tolist()]
trio_moms = simple_trios.motherID.dropna().to_list() + [x for x in complex_trio_parents if x in complex_pedigrees.motherID.tolist()]

len(set(trio_moms)), len(set(trio_dads)), len(set(trio_parents)), len(set(duo_parents)), len(set(trio_children)), len(set(duo_children)), len(complex_people)

## Describes the parent who is in two duo/trios
pedigree.loc[(pedigree.sampleID=='NA19238')|(pedigree.fatherID=='NA19238')|(pedigree.motherID=='NA19238')]

# %% [markdown]
# ### Export per-sample statistics for supplemental table 1

# %%
samples_for_export = all_samples.loc[(all_samples.method_of_phasing=='phased_with_parents_and_pedigree')&(all_samples.ground_truth_data_source.isin(('HPRC_samples','HGSVC_samples')))]
sample_export_columns = ['sample_id', 'panel_genome', 'ground_truth_data_source', 'family_status','population', 'superpopulation', 'sex', 'n_switch_errors', 'n_checked', 'switch_error_rate',
       'n_total_switch_errors', 'n_flips', 'n_consecutive_flips',
       'n_true_switch_errors', 'num_correct_phased_hets', 'n_total_hets',
       'flip_and_switches_SER', 'flip_error_rate', 'true_switch_error_rate',
       'accurately_phased_rate', 'n_gt_errors', 'n_gt_checked',
       'gt_error_rate', 'N50']

source_accuracy_columns = ['panel_genome', 'ground_truth_data_source', 'family_status', 'n_switch_errors', 'n_checked', 'switch_error_rate',
       'n_total_switch_errors', 'n_flips', 'flip_error_rate',
       'n_true_switch_errors', 'n_total_hets','true_switch_error_rate',
       'n_gt_errors', 'n_gt_checked',
       'gt_error_rate', 'N50']

samples_for_export = samples_for_export.sort_values(['genome','ground_truth_data_source', 'sample_id']).rename(columns={'genome':'panel_genome'})
samples_for_export['family_status'] = 'non_trio'
samples_for_export.loc[samples_for_export.trio_phased,'family_status'] = 'trio_proband'
samples_for_export.loc[samples_for_export.sample_id.isin(pedigree.fatherID.to_list()+pedigree.motherID.to_list()),'family_status'] = 'trio_parent'
source_accuracy = samples_for_export[sample_export_columns].groupby(['panel_genome','ground_truth_data_source','family_status'], observed=True).sum(numeric_only=True).reset_index()
source_accuracy['switch_error_rate'] = source_accuracy.n_switch_errors/source_accuracy['n_checked'] * 100
source_accuracy['gt_error_rate'] = source_accuracy['n_gt_errors']/source_accuracy['n_gt_checked'] * 100
source_accuracy['true_switch_error_rate'] = source_accuracy.n_true_switch_errors/source_accuracy.n_total_hets * 100
source_accuracy['flip_error_rate'] = source_accuracy.n_flips/source_accuracy.n_total_hets * 100
os.makedirs('figures/supplemental/tables', exist_ok=True)
source_accuracy[source_accuracy_columns].to_csv('figures/supplemental/tables/supplemental_ground_truth_source_accuracy_data_unformatted.tsv', sep='\t', index=False)
samples_for_export[sample_export_columns].to_csv('figures/supplemental/tables/supplemental_per_sample_data_unformatted.tsv', sep='\t', index=False)
source_accuracy.T


# %%
def group_phasing_statistics(df, columns_to_groupby):
    df_cols = df.columns
    assert all([c in df_cols for c in columns_to_groupby])

    df = df.groupby(columns_to_groupby, observed=True).sum(numeric_only=True).reset_index()

    if 'n_gt_errors' in df_cols:
        df['gt_error_rate'] = df.n_gt_errors/df.n_gt_checked * 100
    if 'n_true_switch_errors' in df_cols:
        # See definitions below.
        df['flip_error_rate'] = df.n_flips/df.n_checked*100
        df['true_switch_error_rate'] = df.n_true_switch_errors/df.n_checked*100
    if 'n_switch_errors' in df_cols:
        df['switch_error_rate'] = df.n_switch_errors/df.n_checked*100

    return df

chrX = all_chroms.loc[(all_chroms.sex=='Female')&(all_chroms.chrom=='chrX') & (all_chroms.method_of_phasing=='phased_with_parents_and_pedigree') & all_chroms.ground_truth_data_source.isin(['HPRC_samples','HGSVC_samples_nontrios_only',])]

group_phasing_statistics(chrX, ['genome','ground_truth_data_source'])[['genome','ground_truth_data_source','flip_error_rate', 'true_switch_error_rate']]

# %%
### Data generated in make_figures.ipynb

error_rate_corr = pd.read_csv('figures/correlation_error_rates_and_gt_error_rate.csv')
error_rate_corr

# %%
HPRC_samples_phased_with_panel = binned_maf_data.loc[(binned_maf_data.ground_truth_data_source=='HPRC_samples')&(binned_maf_data.method_of_phasing == '1kgp_variation_phased_with_reference_panel')]
HPRC_samples_phased_with_panel = group_phasing_statistics(HPRC_samples_phased_with_panel, ['type','genome','Syntenic']).drop(columns=['MAC','AN','MAF'])
error_rate_corr.to_csv('figures/panel_phased_HPRC_samples_summary_data.csv', index=False)

# %% [markdown]
# ## Get number of probands/parents in 1KGP

# %%
parents = len(trio_parents) + len(duo_parents)
probands = len(trio_children) + len(duo_children)
statistical = 3202-parents-probands
parents, probands, statistical

# %%
number_of_datasets_per_sample = all_samples.loc[(all_samples.genome=='GRCh38') &
                                (all_samples['method_of_phasing']=='phased_with_parents_and_pedigree') &
                                (all_samples.ground_truth_data_source.isin(('HPRC_samples','HGSVC_samples')))
                                ].groupby(['sample_id'], observed=True).size().sort_values()
print ("These samples are in both the HPRC and HGSVC datasets. We will use the HPRC genomes.")
print (number_of_datasets_per_sample.loc[number_of_datasets_per_sample>1].index.tolist())

# %%
all_methods.loc[(all_methods.method_of_phasing == 'phased_with_parents_and_pedigree')& (all_methods.ground_truth_data_source=='HPRC_samples')]

# %%
calc_ancestries = all_ancestries.loc[(all_ancestries.method_of_phasing == 'phased_with_parents_and_pedigree') & 
                                     (all_ancestries.ground_truth_data_source.isin(('HPRC_samples','HGSVC_samples')))]

calc_ancestries = calc_ancestries.groupby(['superpopulation','genome'], observed=True).sum(numeric_only=True)
calc_ancestries['switch_error_rate'] = calc_ancestries.n_switch_errors/calc_ancestries.n_checked*100
calc_ancestries['flip_error_rate'] = calc_ancestries.n_flips/calc_ancestries.n_checked*100
calc_ancestries['true_switch_error_rate'] = calc_ancestries.n_true_switch_errors/calc_ancestries.n_checked*100
calc_ancestries = calc_ancestries.reset_index()
calc_ancestries.loc[calc_ancestries.genome=='GRCh38', 'switch_error_ratio'] = calc_ancestries.loc[calc_ancestries.genome=='GRCh38', 'switch_error_rate'] / calc_ancestries.loc[calc_ancestries.genome=='GRCh38', 'switch_error_rate'].min()
calc_ancestries.loc[calc_ancestries.genome=='CHM13v2.0', 'switch_error_ratio'] = calc_ancestries.loc[calc_ancestries.genome=='CHM13v2.0', 'switch_error_rate'] / calc_ancestries.loc[calc_ancestries.genome=='CHM13v2.0', 'switch_error_rate'].min()
ax = sns.barplot(data=calc_ancestries.reset_index(), x='superpopulation', hue='genome', y='switch_error_ratio')
calc_ancestries.sort_values(['genome','superpopulation'])[['genome','superpopulation', 'switch_error_ratio']]

# %%
chrX = all_chroms.loc[
    (all_chroms.sex == 'Female') &
    (all_chroms.chrom == 'chrX') &
    (all_chroms.method_of_phasing == 'phased_with_parents_and_pedigree') &
    all_chroms.ground_truth_data_source.isin(['HPRC_samples'])
]

group_phasing_statistics(chrX, ['genome'])[['genome', 'flip_error_rate', 'true_switch_error_rate']]

# %%
if per_variant_category_imputation_performance is not None:
    # ── Reproduce the make_plots.ipynb wrangling to check for bugs ──

    # Step 1: Reproduce the per-variant wrangling (cells 7-8 in make_plots)
    pv = per_variant_category_imputation_performance.copy()

    separate_panel_and_subset = pv.assign(
        info_cutoff = pv.panel.str.extract(r'info_cutoff_(\d+\.\d+)', expand=False).fillna('0'),
        panel       = pv.panel.str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True)
    ).copy()

    separate_panel_and_subset = separate_panel_and_subset.assign(**{'Variant Subset':'All panel variants', 'Reference Panel':'GRCh38'})
    separate_panel_and_subset.loc[separate_panel_and_subset.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
    separate_panel_and_subset['panel'] = (separate_panel_and_subset['panel']
        .str.replace('.common_variants','')
        .str.replace('native_panel','Native')
        .str.replace('lifted_panel','Lifted'))
    separate_panel_and_subset.loc[
        ((separate_panel_and_subset.genome=='GRCh38') & (separate_panel_and_subset.panel=='Lifted')) |
        ((separate_panel_and_subset.genome=='CHM13v2.0') & (separate_panel_and_subset.panel=='Native')),
        'Reference Panel'] = 'CHM13v2.0'

    separate_panel_and_subset['panel_filter'] = (separate_panel_and_subset['dataset']
        .str.replace('SGDP', '')
        .str.replace('pangenome', '')
        .str.split('_', n=1).str[1].fillna('No Filter'))
    separate_panel_and_subset['dataset'] = (separate_panel_and_subset['dataset']
        .str.split('_', n=1).str[0].replace('pangenome','Pangenome'))
    separate_panel_and_subset = separate_panel_and_subset.rename(columns={'var_category_id':'Variant Type'})
    if separate_panel_and_subset.mean_AF.max() < 1:
        separate_panel_and_subset.mean_AF *= 100

    print("=== Wrangled per-variant imputation data ===")
    print(f"Unique 'Reference Panel': {separate_panel_and_subset['Reference Panel'].unique()}")
    print(f"Unique 'panel': {separate_panel_and_subset['panel'].unique()}")
    print(f"Unique 'panel_filter': {separate_panel_and_subset['panel_filter'].unique()}")
    print(f"Unique 'dataset': {separate_panel_and_subset['dataset'].unique()}")
    print(f"Unique 'Variant Subset': {separate_panel_and_subset['Variant Subset'].unique()}")
    print(f"Unique 'Variant Type': {separate_panel_and_subset['Variant Type'].unique()}")

    # Step 2: Reproduce per-sample wrangling
    ps = per_sample_imputation_performance.copy()
    # NOTE: make_plots.ipynb does NOT filter out syntenic files for per_sample the same way 
    # But calc_figures_for_paper.ipynb already replaced T2T->CHM13v2.0. make_plots.ipynb does NOT do this for per_sample.
    # Let's check:
    ps_filtered = ps.loc[~ps.file.str.contains('yntenic')]
    print(f"\n=== Per-sample data ===")
    print(f"Before syntenic filter: {ps.shape[0]}, After: {ps_filtered.shape[0]}")
    print(f"Unique genome values in per_sample: {ps_filtered['genome'].unique()}")

    separate_sample_panel_and_subset = ps_filtered.assign(
        info_cutoff = ps_filtered.panel.str.extract(r'info_cutoff_(\d+\.\d+)', expand=False).fillna('0'),
        panel       = ps_filtered.panel.str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True)
    ).copy()
    separate_sample_panel_and_subset = separate_sample_panel_and_subset.assign(**{'Variant Subset':'All panel variants', 'Reference Panel':'GRCh38'})
    separate_sample_panel_and_subset.loc[separate_sample_panel_and_subset.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
    separate_sample_panel_and_subset['panel'] = (separate_sample_panel_and_subset['panel']
        .str.replace('.common_variants','')
        .str.replace('native_panel','Native')
        .str.replace('lifted_panel','Lifted'))

    # BUG CHECK: make_plots.ipynb uses genome=='T2T' here, but paper_inline_statistics already replaced T2T->CHM13v2.0
    # This would cause the filter to MISS CHM13 samples
    separate_sample_panel_and_subset.loc[
        ((separate_sample_panel_and_subset.genome=='GRCh38') & (separate_sample_panel_and_subset.panel=='Lifted')) |
        ((separate_sample_panel_and_subset.genome=='T2T') & (separate_sample_panel_and_subset.panel=='Native')),
        'Reference Panel'] = 'CHM13v2.0'

    n_chm13_refpanel_T2T = (separate_sample_panel_and_subset['Reference Panel'] == 'CHM13v2.0').sum()
    print(f"\nWith genome=='T2T' filter: {n_chm13_refpanel_T2T} rows get Reference Panel = 'CHM13v2.0'")

    # Try with CHM13v2.0 instead
    separate_sample_panel_and_subset.loc[
        ((separate_sample_panel_and_subset.genome=='GRCh38') & (separate_sample_panel_and_subset.panel=='Lifted')) |
        ((separate_sample_panel_and_subset.genome=='CHM13v2.0') & (separate_sample_panel_and_subset.panel=='Native')),
        'Reference Panel'] = 'CHM13v2.0'

    n_chm13_refpanel_fixed = (separate_sample_panel_and_subset['Reference Panel'] == 'CHM13v2.0').sum()
    print(f"With genome=='CHM13v2.0' filter: {n_chm13_refpanel_fixed} rows get Reference Panel = 'CHM13v2.0'")

    # Check: what unique Reference Panel values exist?
    print(f"\nFinal unique 'Reference Panel': {separate_sample_panel_and_subset['Reference Panel'].unique()}")

    # BUG CHECK: The singleton analysis in make_plots uses Reference Panel == 'Lalli et al. CHM13v2.0 panel'
    # Does that value ever exist?
    n_lalli_match = (separate_sample_panel_and_subset['Reference Panel'] == 'Lalli et al. CHM13v2.0 panel').sum()
    print(f"\n⚠️  Rows matching 'Lalli et al. CHM13v2.0 panel': {n_lalli_match}")
    print("   (The singleton analysis in make_plots.ipynb filters for this value!)")
else:
    print('Skipping: per_variant_category_imputation_performance not available (no imputation data)')

# %%
if per_variant_category_imputation_performance is not None:
    # ── Reproduce the per-variant singleton Δr² analysis (correctly) ──
    # This mirrors make_plots.ipynb cells 160-165 but we already have the wrangled data.

    # Filter: SGDP, Native panel, CHM13v2.0, All panel variants
    analysis_data = separate_panel_and_subset[
        (separate_panel_and_subset['dataset'] == 'SGDP') &
        (separate_panel_and_subset['panel'] == 'Native') &
        (separate_panel_and_subset['Reference Panel'] == 'CHM13v2.0') &
        (separate_panel_and_subset['Variant Subset'] == 'All panel variants') &
        (separate_panel_and_subset['panel_filter'].isin(['No Filter', 'no_singletons']))
    ]

    print(f"Analysis data shape: {analysis_data.shape}")
    print(f"Panel filters: {analysis_data['panel_filter'].unique()}")
    print(f"Variant types: {analysis_data['Variant Type'].unique()}")

    # ── SNPs Δr² ──
    snp_data = analysis_data[analysis_data['Variant Type'] == 'SNPs'].copy()
    snp_pivot = snp_data.pivot_table(
        index=['variant_bin'],
        columns='panel_filter',
        values='imputed_ds_rsquared',
        aggfunc='mean'
    )

    if 'No Filter' in snp_pivot.columns and 'no_singletons' in snp_pivot.columns:
        snp_pivot['delta_r2'] = snp_pivot['No Filter'] - snp_pivot['no_singletons']
        print("\n=== SNPs: Δr² (with singletons - without singletons) ===")
        print(f"Mean Δr²:    {snp_pivot['delta_r2'].mean():.6f}")
        print(f"Median Δr²:  {snp_pivot['delta_r2'].median():.6f}")
        print(f"Min Δr²:     {snp_pivot['delta_r2'].min():.6f}")
        print(f"Max Δr²:     {snp_pivot['delta_r2'].max():.6f}")
        print(f"Bins analyzed: {len(snp_pivot)}")
    else:
        print("WARNING: Missing panel_filter columns for SNPs")
        print(f"Available: {snp_pivot.columns.tolist()}")

    # ── Indels Δr² ──
    indel_data = analysis_data[analysis_data['Variant Type'] == 'Indels'].copy()
    indel_pivot = indel_data.pivot_table(
        index=['mean_AF'],
        columns='panel_filter',
        values='imputed_ds_rsquared',
        aggfunc='mean'
    )

    if 'No Filter' in indel_pivot.columns and 'no_singletons' in indel_pivot.columns:
        indel_pivot['delta_r2'] = indel_pivot['No Filter'] - indel_pivot['no_singletons']
        print("\n=== Indels: Δr² (with singletons - without singletons) ===")
        print(f"Mean Δr²:    {indel_pivot['delta_r2'].mean():.6f}")
        print(f"Median Δr²:  {indel_pivot['delta_r2'].median():.6f}")
        print(f"Min Δr²:     {indel_pivot['delta_r2'].min():.6f}")
        print(f"Max Δr²:     {indel_pivot['delta_r2'].max():.6f}")
        print(f"Bins analyzed: {len(indel_pivot)}")
    else:
        print("WARNING: Missing panel_filter columns for Indels")
        print(f"Available: {indel_pivot.columns.tolist()}")

    # ── Check make_plots.ipynb hardcoded claims ──
    print("\n=== Checking hardcoded claims in make_plots.ipynb summary ===")
    print(f"Claim: 'Small positive effect on average (+0.003)' for SNPs")
    print(f"  Actual: {snp_pivot['delta_r2'].mean():+.6f}  {'✅ MATCH' if abs(snp_pivot['delta_r2'].mean() - 0.003) < 0.001 else '❌ MISMATCH'}")
    print(f"Claim: 'Near-zero effect on average (+0.0005)' for Indels")
    print(f"  Actual: {indel_pivot['delta_r2'].mean():+.6f}  {'✅ MATCH' if abs(indel_pivot['delta_r2'].mean() - 0.0005) < 0.001 else '❌ MISMATCH'}")
    print(f"Claim: 'Very few bins show non-trivial harm (2% SNPs, 3% indels)' [Δr² < -0.01]")
    snp_harm_pct = (snp_pivot['delta_r2'] < -0.01).sum() / len(snp_pivot) * 100
    indel_harm_pct = (indel_pivot['delta_r2'] < -0.01).sum() / len(indel_pivot) * 100
    print(f"  SNPs:   {snp_harm_pct:.1f}% bins harmed  {'✅' if abs(snp_harm_pct - 2) < 2 else '❌'}")
    print(f"  Indels: {indel_harm_pct:.1f}% bins harmed  {'✅' if abs(indel_harm_pct - 3) < 2 else '❌'}")
else:
    print('Skipping: per_variant_category_imputation_performance not available (no imputation data)')

# %%
# ── Per-sample Δr² analysis (with CORRECTED filter values) ──
# make_plots.ipynb bug: filters for Reference Panel == 'Lalli et al. CHM13v2.0 panel'
# which matches 0 rows. Correct value is 'CHM13v2.0'.

# Fix the Syntenic column assignment (as make_plots does in cell 167)
separate_sample_panel_and_subset['Syntenic'] = 'unseen'
separate_sample_panel_and_subset.loc[separate_sample_panel_and_subset.file.str.contains('syntenic'), 'Syntenic'] = 'Syntenic'
separate_sample_panel_and_subset.loc[~separate_sample_panel_and_subset.file.str.contains('syntenic'), 'Syntenic'] = 'All'

# Add panel_filter + dataset columns
separate_sample_panel_and_subset['panel_filter'] = (separate_sample_panel_and_subset['dataset']
    .str.replace('SGDP', '')
    .str.replace('pangenome', '')
    .str.split('_', n=1).str[1].fillna('No Filter'))
separate_sample_panel_and_subset['dataset'] = (separate_sample_panel_and_subset['dataset']
    .str.split('_', n=1).str[0].replace('pangenome','Pangenome'))
separate_sample_panel_and_subset = separate_sample_panel_and_subset.rename(columns={'var_category_id':'Variant Type'})

# Apply CORRECT filter (matching what make_plots SHOULD do)
sample_data = separate_sample_panel_and_subset[
    (separate_sample_panel_and_subset['dataset'] == 'SGDP') &
    (separate_sample_panel_and_subset['panel_filter'].isin(['No Filter', 'no_singletons'])) &
    (separate_sample_panel_and_subset['ancestry'] == 'all') &
    (separate_sample_panel_and_subset['Variant Subset'] == 'All panel variants') &
    (separate_sample_panel_and_subset['Syntenic'] == 'All') &
    (separate_sample_panel_and_subset['panel'] == 'Native') &
    (separate_sample_panel_and_subset['Reference Panel'] == 'CHM13v2.0')  # ← CORRECTED
].copy()

print(f"Per-sample data shape (corrected): {sample_data.shape}")
print(f"Unique samples: {sample_data['sample_name'].nunique()}")
print(f"Panel filters: {sample_data['panel_filter'].unique()}")
print(f"Variant types: {sample_data['Variant Type'].unique()}")

if len(sample_data) == 0:
    print("⚠️  No data matched even with corrected filter!")
else:
    # Calculate per-sample Δr²
    sample_pivot = sample_data.pivot_table(
        index=['sample_name', 'Variant Type'],
        columns='panel_filter',
        values='imputed_ds_rsquared',
        aggfunc='mean'
    )
    
    if 'No Filter' in sample_pivot.columns and 'no_singletons' in sample_pivot.columns:
        sample_pivot['delta_r2'] = sample_pivot['No Filter'] - sample_pivot['no_singletons']
        
        print(f"\n=== Per-sample Δr² (all variant types combined) ===")
        print(f"Median: {sample_pivot['delta_r2'].median():+.6f}")
        print(f"Mean:   {sample_pivot['delta_r2'].mean():+.6f}")
        print(f"Std:    {sample_pivot['delta_r2'].std():.6f}")
        print(f"Positive: {(sample_pivot['delta_r2'] > 0).sum()}/{len(sample_pivot)}")
        
        for vtype in sample_data['Variant Type'].unique():
            vtype_data = sample_pivot.xs(vtype, level='Variant Type')
            print(f"\n  {vtype}:")
            print(f"    Median Δr²: {vtype_data['delta_r2'].median():+.6f}")
            print(f"    Mean Δr²:   {vtype_data['delta_r2'].mean():+.6f}")
            print(f"    Positive:   {(vtype_data['delta_r2'] > 0).sum()}/{len(vtype_data)}")
        
        # Check hardcoded claims
        print("\n=== Checking hardcoded per-sample claims ===")
        snp_samples = sample_pivot.xs('SNPs', level='Variant Type')
        print(f"Claim: 'Consistently negative (0/N positive)' for SNPs per-sample")
        print(f"  Actual: {(snp_samples['delta_r2'] > 0).sum()}/{len(snp_samples)} positive")
        n_near_zero = (sample_pivot['delta_r2'].abs() < 0.001).sum()
        print(f"Claim: '|Δr²| < 0.001 for 98% of samples'")
        print(f"  Actual: {n_near_zero/len(sample_pivot)*100:.1f}%")
        n_harmed = (sample_pivot['delta_r2'] < -0.01).sum()
        print(f"Claim: 'NO samples show non-trivial harm (Δr² < -0.01)'")
        print(f"  Actual: {n_harmed} samples harmed  {'✅' if n_harmed == 0 else '❌'}")
    else:
        print(f"Missing columns in pivot. Available: {sample_pivot.columns.tolist()}")

# Now show what make_plots.ipynb would get with the BUGGY filter
print("\n" + "="*60)
print("BUG IMPACT: What make_plots.ipynb actually filters:")
print("="*60)
buggy_data = separate_sample_panel_and_subset[
    (separate_sample_panel_and_subset['dataset'] == 'SGDP') &
    (separate_sample_panel_and_subset['panel_filter'].isin(['No Filter', 'no_singletons'])) &
    (separate_sample_panel_and_subset['ancestry'] == 'all') &
    (separate_sample_panel_and_subset['Variant Subset'] == 'All panel variants') &
    (separate_sample_panel_and_subset['Syntenic'] == 'All') &
    (separate_sample_panel_and_subset['panel'] == 'Native') &
    (separate_sample_panel_and_subset['Reference Panel'] == 'Lalli et al. CHM13v2.0 panel')  # ← BUGGY
]
print(f"Rows matching buggy filter: {len(buggy_data)}")
print(f"→ Per-sample singleton impact analysis produces EMPTY results in make_plots.ipynb!")

# %%
if per_variant_category_imputation_performance is not None:
    # ═══════════════════════════════════════════════════════════════════
    # FULL END-TO-END VERIFICATION of the make_plots.ipynb fixes
    # Reproduces cells 6→7→8→9 (wrangling) then 160→163→165→167→169→173
    # using the FIXED code, from raw parquet data.
    # ═══════════════════════════════════════════════════════════════════
    import pandas as pd, numpy as np

    # ── Reload raw data (clean slate) ──
    _pv = pd.read_parquet(f"{summary_statistics_folder}/per_variant_category_imputation_performance.parquet")
    _ps = pd.read_parquet(f"{summary_statistics_folder}/per_sample_imputation_performance.parquet")

    # ── Cell 6 equivalent: compute perc_discordance, replace T2T ──
    _pv = _pv.loc[_pv.num_AA > 0]
    if _pv.mean_AF.max() < 1:
        _pv.mean_AF *= 100
    _pv['num_alt_mismatches'] = _pv.num_Aa_mismatches + _pv.num_aa_mismatches
    _pv['num_alt_variants'] = _pv.num_Aa + _pv.num_aa
    _pv['perc_discordance'] = _pv['num_alt_mismatches'] / _pv['num_alt_variants'] * 100
    _pv['perc_concordance'] = 100 - _pv['perc_discordance']
    _pv = _pv.replace('T2T', 'CHM13v2.0')
    _ps = _ps.replace('T2T', 'CHM13v2.0')

    # ── Cell 7: per-variant wrangling ──
    _sep_pv = _pv.assign(
        info_cutoff=_pv.panel.str.extract(r'info_cutoff_(\d+\.\d+)', expand=False).fillna('0'),
        panel=_pv.panel.str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True)
    ).copy()
    _sep_pv = _sep_pv.assign(**{'Variant Subset': 'All panel variants', 'Reference Panel': 'GRCh38'})
    _sep_pv.loc[_sep_pv.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
    _sep_pv['panel'] = _sep_pv['panel'].str.replace('.common_variants', '').str.replace('native_panel', 'Native').str.replace('lifted_panel', 'Lifted')
    _sep_pv.loc[
        ((_sep_pv.genome == 'GRCh38') & (_sep_pv.panel == 'Lifted')) |
        ((_sep_pv.genome == 'CHM13v2.0') & (_sep_pv.panel == 'Native')),
        'Reference Panel'] = 'CHM13v2.0'
    _sep_pv['panel_filter'] = _sep_pv['dataset'].str.replace('SGDP', '').str.replace('pangenome', '').str.split('_', n=1).str[1].fillna('No Filter')
    _sep_pv['dataset'] = _sep_pv['dataset'].str.split('_', n=1).str[0].replace('pangenome', 'Pangenome')
    _sep_pv = _sep_pv.rename(columns={'var_category_id': 'Variant Type'})
    if _sep_pv.mean_AF.max() < 1:
        _sep_pv.mean_AF *= 100

    # ── Cell 8 (FIXED): per-sample wrangling ──
    _ps = _ps.loc[~_ps.file.str.contains('yntenic')]
    _sep_ps = _ps.assign(
        info_cutoff=_ps.panel.str.extract(r'info_cutoff_(\d+\.\d+)', expand=False).fillna('0'),
        panel=_ps.panel.str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True)
    ).copy()
    _sep_ps = _sep_ps.assign(**{'Variant Subset': 'All panel variants', 'Reference Panel': 'GRCh38'})
    _sep_ps.loc[_sep_ps.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
    _sep_ps['panel'] = _sep_ps['panel'].str.replace('.common_variants', '').str.replace('native_panel', 'Native').str.replace('lifted_panel', 'Lifted')
    # ★ FIX 1: genome=='CHM13v2.0' instead of genome=='T2T'
    _sep_ps.loc[
        ((_sep_ps.genome == 'GRCh38') & (_sep_ps.panel == 'Lifted')) |
        ((_sep_ps.genome == 'CHM13v2.0') & (_sep_ps.panel == 'Native')),
        'Reference Panel'] = 'CHM13v2.0'
    _sep_ps['panel_filter'] = _sep_ps['dataset'].str.replace('SGDP', '').str.replace('pangenome', '').str.split('_', n=1).str[1].fillna('No Filter')
    _sep_ps['dataset'] = _sep_ps['dataset'].str.split('_', n=1).str[0].replace('pangenome', 'Pangenome')
    _sep_ps['Synteny'] = _sep_ps['file'].str.contains('syntenic').map({True: 'Syntenic', False: 'Nonsyntenic'})
    _sep_ps = _sep_ps.rename(columns={'var_category_id': 'Variant Type'})

    # ── Cell 160: per-variant analysis data filter ──
    _analysis = _sep_pv[
        (_sep_pv['dataset'] == 'SGDP') &
        (_sep_pv['panel'] == 'Native') &
        (_sep_pv['Reference Panel'] == 'CHM13v2.0') &
        (_sep_pv['ancestry'] == 'all') &
        (_sep_pv['Variant Subset'] == 'All panel variants') &
        (_sep_pv['panel_filter'].isin(['No Filter', 'no_singletons']))
    ].copy().drop_duplicates()

    # ── Cell 163: SNP Δr² ──
    _snp = _analysis[_analysis['Variant Type'] == 'SNPs']
    _snp_piv = _snp.pivot_table(index='variant_bin', columns='panel_filter', values='imputed_ds_rsquared', aggfunc='mean')
    _snp_piv['delta_r2'] = _snp_piv['No Filter'] - _snp_piv['no_singletons']
    snp_delta_mean = _snp_piv['delta_r2'].mean()

    # ── Cell 165: Indel Δr² ──
    _indel = _analysis[_analysis['Variant Type'] == 'Indels']
    _indel_piv = _indel.pivot_table(index='mean_AF', columns='panel_filter', values='imputed_ds_rsquared', aggfunc='mean')
    _indel_piv['delta_r2'] = _indel_piv['No Filter'] - _indel_piv['no_singletons']
    indel_delta_mean = _indel_piv['delta_r2'].mean()

    # ── Cell 167 (FIXED): per-sample filter ──
    _sep_ps['Syntenic'] = _sep_ps.Synteny.astype('str')
    _sep_ps.loc[~_sep_ps.file.str.contains('syntenic'), 'Syntenic'] = 'All'

    # ★ FIX 2: Reference Panel == 'CHM13v2.0' instead of 'Lalli et al. CHM13v2.0 panel'
    _sdata = _sep_ps[
        (_sep_ps['dataset'] == 'SGDP') &
        (_sep_ps.panel_filter.isin(['No Filter', 'no_singletons'])) &
        (_sep_ps.ancestry == 'all') &
        (_sep_ps['Variant Subset'] == 'All panel variants') &
        (_sep_ps['Syntenic'] == 'All') &
        (_sep_ps['panel'] == 'Native') &
        (_sep_ps['Reference Panel'] == 'CHM13v2.0')
    ].copy()

    # ── Cell 169: per-sample Δr² pivot ──
    _sample_piv = _sdata.pivot_table(
        index=['sample_name', 'Variant Type'],
        columns='panel_filter',
        values='imputed_ds_rsquared',
        aggfunc='mean'
    )
    _sample_piv['delta_r2'] = _sample_piv['No Filter'] - _sample_piv['no_singletons']

    # ── Cell 173: summary (FIXED — all computed, no hardcoded values) ──
    _snp_s = _sample_piv.xs('SNPs', level='Variant Type')
    _indel_s = _sample_piv.xs('Indels', level='Variant Type')
    _all_s = _sample_piv.xs('All variants', level='Variant Type')
    _n_snp_harmed = (_snp_piv['delta_r2'] < -0.01).sum()
    _n_indel_harmed = (_indel_piv['delta_r2'] < -0.01).sum()
    _n_sample_harmed = (_sample_piv['delta_r2'] < -0.01).sum()
    _pct_near_zero = (_sample_piv['delta_r2'].abs() < 0.001).sum() / len(_sample_piv) * 100

    print("=" * 70)
    print("  VERIFICATION: Fixed make_plots.ipynb singleton analysis")
    print("=" * 70)
    print()
    print("── Per-variant Δr² ──")
    print(f"  SNP mean:   {snp_delta_mean:+.6f}   (bins: {len(_snp_piv)})")
    print(f"  Indel mean: {indel_delta_mean:+.6f}   (bins: {len(_indel_piv)})")
    print()
    print("── Per-sample filter (BUG FIX CHECK) ──")
    print(f"  Rows matched by FIXED filter:  {len(_sdata)}  ← should be ~1536")
    print(f"  Unique samples:                {_sdata['sample_name'].nunique()}  ← should be 256")
    print(f"  Variant types:                 {sorted(_sdata['Variant Type'].unique())}")
    print()
    print("── Per-sample Δr² ──")
    print(f"  SNPs:   median={_snp_s['delta_r2'].median():+.6f}  positive={(_snp_s['delta_r2']>0).sum()}/{len(_snp_s)}")
    print(f"  Indels: median={_indel_s['delta_r2'].median():+.6f}  positive={(_indel_s['delta_r2']>0).sum()}/{len(_indel_s)}")
    print(f"  All:    median={_all_s['delta_r2'].median():+.6f}  positive={(_all_s['delta_r2']>0).sum()}/{len(_all_s)}")
    print()
    print("── Harm signal ──")
    print(f"  SNP bins harmed (Δr²<-0.01):    {_n_snp_harmed}/{len(_snp_piv)} ({_n_snp_harmed/len(_snp_piv)*100:.1f}%)")
    print(f"  Indel bins harmed:               {_n_indel_harmed}/{len(_indel_piv)} ({_n_indel_harmed/len(_indel_piv)*100:.1f}%)")
    print(f"  Samples harmed (Δr²<-0.01):     {_n_sample_harmed}/{len(_sample_piv)} ({_n_sample_harmed/len(_sample_piv)*100:.1f}%)")
    print(f"  |Δr²| < 0.001:                  {_pct_near_zero:.1f}% of sample×variant_type combos")
    print()

    # ── Sanity checks ──
    errors = []
    if len(_sdata) == 0:
        errors.append("FAIL: Per-sample filter matched 0 rows!")
    if _sdata['sample_name'].nunique() != 256:
        errors.append(f"WARN: Expected 256 samples, got {_sdata['sample_name'].nunique()}")
    if len(_sdata) < 1500:
        errors.append(f"WARN: Expected ~1536 rows, got {len(_sdata)}")
    if _n_sample_harmed > 0:
        errors.append(f"WARN: {_n_sample_harmed} samples show Δr² < -0.01")

    if errors:
        print("⚠️  ISSUES:")
        for e in errors:
            print(f"   {e}")
    else:
        print("✅ ALL CHECKS PASSED")
        print("   • Per-sample filter matches 1536 rows (256 samples × 3 variant types × 2 filters)")
        print("   • No samples show non-trivial harm")
        print("   • Summary cell will now print computed values instead of hardcoded ones")
    print()
    print("=" * 70)
else:
    print('Skipping: per_variant_category_imputation_performance not available (no imputation data)')
