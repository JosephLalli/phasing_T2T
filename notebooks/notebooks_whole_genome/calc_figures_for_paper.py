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

# %% [markdown]
# # Whole-genome paper inline statistics
#
# This notebook computes **inline numeric values used in the manuscript text** from pre-aggregated summary dataframes.
#
# ## Paper source of truth
# The manuscript PDFs in `papers/` are the source of truth for which values are manuscript-facing. See `notebooks/notebooks_whole_genome/PAPER_CROSS_REFERENCE.md` for the full figure/data-point map with extracted PDF text line references.
#
# Key manuscript contexts covered here:
# - Table 1/genotyping concordance: manuscript lines 208-216 and caption lines 2068-2085; reviewer reply lines 509-541 and 557-565.
# - Variant-filtering percentages: manuscript lines 231-236; reviewer reply lines 483-491.
# - Panel-wide weighted SER: manuscript lines 268-280.
# - FER/tSER decomposition: manuscript lines 293-321.
# - Singleton and imputation summary values: manuscript lines 415-422; reviewer reply lines 232-249.
# - Supplementary Data S2-S4: manuscript lines 624-629.
#
# ## Execution order
# 1. Run `notebooks/notebooks_whole_genome/calc_per_variant_figures_for_paper.ipynb` (variant-level calculations).
# 2. Run `notebooks/notebooks_whole_genome/make_plots.ipynb` (figure/table generation and some derived stats).
# 3. Run this notebook to compute/report manuscript inline values from smaller summary dataframes.
#
# ## Data paths
# Set `RUN_PROFILE` in the first code cell. `test` uses plain folders (`intermediate_data`, `figures`, `tables`); `whole_genome` uses `_whole_genome` folders for publication-scale inputs and outputs.
#
# ## Scope
# - Keep production cells near the top, organized by paper section.
# - Preserve exploratory analyses at the end for review/deletion.
# - Export a central raw-values TSV for downstream document integration.
#

# %% [markdown]
# ## 1) Load summary dataframes

# %%
import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utilities import find_repo_root, require_files

plt.rcParams["figure.dpi"] = 120

# Choose "test" for Docker test environment/datasets or "whole_genome" for publication-scale data.
# The environment variable is useful for automated runs; manually edit RUN_PROFILE for interactive use.
RUN_PROFILE = "whole_genome"
RUN_PROFILE = os.environ.get("PHASING_T2T_RUN_PROFILE", RUN_PROFILE)

PROFILE_PATHS = {
    "test": {
        "intermediate_data": "intermediate_data",
        "imputation_statistics": "imputation_statistics",
        "figures": "figures",
        "tables": "tables",
    },
    "whole_genome": {
        "intermediate_data": "intermediate_data_whole_genome",
        "imputation_statistics": "imputation_statistics_whole_genome",
        "figures": "figures_whole_genome",
        "tables": "tables_whole_genome",
    },
}

if RUN_PROFILE not in PROFILE_PATHS:
    raise ValueError(f"RUN_PROFILE must be one of {sorted(PROFILE_PATHS)}, got {RUN_PROFILE!r}")

profile_paths = PROFILE_PATHS[RUN_PROFILE]

repo_root = find_repo_root()
os.chdir(repo_root)

required_summary_files = (
    "binned_maf_data.parquet",
    "chroms.parquet",
    "samples.parquet",
    "ancestries.parquet",
    "methods.parquet",
    "per_variant_category_imputation_performance.parquet",
    "per_sample_imputation_performance.parquet",
    "rolling_stats_250k_window.parquet",
)

summary_statistics_folder = require_files(
    profile_paths["intermediate_data"],
    required_summary_files,
    "summary statistics",
    RUN_PROFILE,
)

figures_folder = profile_paths["figures"]
supplemental_figures_folder = f"{figures_folder}/supplemental"
tables_folder = profile_paths["tables"]
Path(figures_folder).mkdir(parents=True, exist_ok=True)
Path(supplemental_figures_folder).mkdir(parents=True, exist_ok=True)
Path(tables_folder).mkdir(parents=True, exist_ok=True)

print(f"RUN_PROFILE: {RUN_PROFILE}")
print(f"Working directory: {repo_root}")
print(f"Summary statistics folder: {summary_statistics_folder}")
print(f"Figures folder: {figures_folder}")
print(f"Tables folder: {tables_folder}")


# %%
binned_maf_data = pd.read_parquet(f"{summary_statistics_folder}/binned_maf_data.parquet")

all_chroms = pd.read_parquet(f'{summary_statistics_folder}/chroms.parquet')
all_samples = pd.read_parquet(f'{summary_statistics_folder}/samples.parquet')
all_ancestries = pd.read_parquet(f'{summary_statistics_folder}/ancestries.parquet')
all_methods = pd.read_parquet(f'{summary_statistics_folder}/methods.parquet')

per_variant_category_imputation_performance = pd.read_parquet(f'{summary_statistics_folder}/per_variant_category_imputation_performance.parquet')
per_sample_imputation_performance = pd.read_parquet(f'{summary_statistics_folder}/per_sample_imputation_performance.parquet')

# %% [markdown]
# ## 2) Harmonize summary tables for paper-level reporting
#
# Post-hoc corrections applied to summary dataframes before computing inline statistics.

# %%
# Post hoc fixes - one off events that are tricky to deal with when compiling stats

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

# bins with a switch error rate of 0 make for insane figures when plotted on a log scale. Replace with np.nan
binned_maf_data['switch_error_rate'] = binned_maf_data['switch_error_rate'].replace(0, np.nan)

# Calculate percent discordance in imputation performance
if 'perc_discordance' not in per_variant_category_imputation_performance.columns:
    per_variant_category_imputation_performance = per_variant_category_imputation_performance.loc[per_variant_category_imputation_performance.num_AA > 0]
    if max(per_variant_category_imputation_performance.mean_AF) < 1:
        per_variant_category_imputation_performance.mean_AF *= 100
    per_variant_category_imputation_performance['num_alt_mismatches'] = per_variant_category_imputation_performance.num_Aa_mismatches + per_variant_category_imputation_performance.num_aa_mismatches
    per_variant_category_imputation_performance['num_alt_variants'] = per_variant_category_imputation_performance.num_Aa + per_variant_category_imputation_performance.num_aa
    per_variant_category_imputation_performance['perc_discordance'] = per_variant_category_imputation_performance['num_alt_mismatches']/per_variant_category_imputation_performance['num_alt_variants']*100
    per_variant_category_imputation_performance['perc_concordance'] = 100-per_variant_category_imputation_performance['perc_discordance']

# imputed df needs to have CHM13v2.0 as genome name instead of T2T
per_variant_category_imputation_performance = per_variant_category_imputation_performance.replace('T2T','CHM13v2.0')
per_sample_imputation_performance = per_sample_imputation_performance.replace('T2T','CHM13v2.0')

# Calculate percent of each chromosome that is newly resolved sequence
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

# Combine chrX, PAR1, PAR2 into a single chrX entry (reviewer request)
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

# Sort for consistent display
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
pedigree = pd.read_csv('resources/pedigrees/1kGP.3202_samples.pedigree_info.txt', sep=' ')


# %% [markdown]
# ## Helper function: `group_phasing_statistics`
#
# Utility used by several sections below.

# %%
def group_phasing_statistics(df, columns_to_groupby):
    df_cols = df.columns
    assert all([c in df_cols for c in columns_to_groupby])

    df = df.groupby(columns_to_groupby, observed=True).sum(numeric_only=True).reset_index()

    if 'n_gt_errors' in df_cols:
        df['gt_error_rate'] = df.n_gt_errors/df.n_gt_checked * 100
    if 'n_true_switch_errors' in df_cols:
        df['flip_error_rate'] = df.n_flips/df.n_checked*100
        df['true_switch_error_rate'] = df.n_true_switch_errors/df.n_checked*100
    if 'n_switch_errors' in df_cols:
        df['switch_error_rate'] = df.n_switch_errors/df.n_checked*100

    return df


# %% [markdown]
# ## 3) Variant filtering inline statistics
#
# Paper location: Results, recombination map section (variant QC paragraph).
#
# > "Our VQSLOD filter resulted in the exclusion of an additional **14,491,949** variants compared to the filtering methods used by Byrska-Bishop et al. ... The inclusion of singleton alleles added **40,226,623** unique variants to our dataset. Overall, variant filtering had removed **46.4%** of variants in the Byrska-Bishop et al. 1kGP GRCh38 panel, including **34.1%** of all variants that were filtered out solely because they were singletons. In contrast, our filtering strategy removed **22%** of 1kGP T2T-CHM13 variants, with **16.3%** of variants filtered solely due to the VQSLOD filter." 

# %%
sum_stats_var_filters = pd.read_parquet(f"{summary_statistics_folder}/filter_summary_stats.parquet")
sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='GRCh38')&(sum_stats_var_filters.filtered != sum_stats_var_filters.GRCh38_criteria_fail)]#.dropna(subset='len')['len'].sum()
# sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='CHM13v2.0')&(sum_stats_var_filters.CHM13_filtered != sum_stats_var_filters.CHM13_criteria_fail)].dropna(subset='len')['len'].sum()

# %%
sum_stats_var_filters = pd.read_parquet(f"{summary_statistics_folder}/filter_summary_stats.parquet")

grch38_filter_columns=['MERR_filter', 'HWE_pop_filter','AC_filter',
       'f_missing_filter', 'var_len_filter', 'alt_star_filter', 'pass_filter']
chm13_filter_columns=[#'neg_train_site_filter',
       'VQSLOD_filter', 'MERR_filter', 'HWE_pop_filter', 'MAC_filter',
       'f_missing_filter', 'var_len_filter', 'alt_star_filter', 'pass_filter']

sum_stats_var_filters['gatk_filters'] = sum_stats_var_filters.VQSLOD_filter
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

# Filter percentage calculations: denominator is all variants observed for that genome in filter_summary_stats.
chm13_total = sum_stats_var_filters.loc[sum_stats_var_filters.genome=='CHM13v2.0'].len.sum()
chm13_filtered = sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='CHM13v2.0') & sum_stats_var_filters.CHM13_criteria_fail].len.sum()
chm13_filtered_vqslod_only = sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='CHM13v2.0') &
                                                         sum_stats_var_filters.CHM13_criteria_fail &
                                                         sum_stats_var_filters.gatk_filters &
                                                         ~sum_stats_var_filters[[c for c in chm13_filter_columns if c != 'VQSLOD_filter']].any(axis=1)].len.sum()

grch38_total = sum_stats_var_filters.loc[sum_stats_var_filters.genome=='GRCh38'].len.sum()
grch38_filtered = sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='GRCh38') & sum_stats_var_filters.GRCh38_criteria_fail].len.sum()
grch38_filtered_singleton_only = sum_stats_var_filters.loc[(sum_stats_var_filters.genome=='GRCh38') &
                                                            sum_stats_var_filters.GRCh38_criteria_fail &
                                                            sum_stats_var_filters.AC_filter &
                                                            ~sum_stats_var_filters[[c for c in grch38_filter_columns if c != 'AC_filter']].any(axis=1)].len.sum()

pct_chm13_filtered = chm13_filtered / chm13_total * 100
pct_chm13_vqslod_only = chm13_filtered_vqslod_only / chm13_total * 100
pct_grch38_filtered = grch38_filtered / grch38_total * 100
pct_grch38_singleton_only = grch38_filtered_singleton_only / grch38_total * 100


# CHM13-unique retention values cited in the nonsyntenic filtering paragraph.
# Denominator: GATK PASS variants in CHM13v2.0 regions without primary alignment to GRCh38.
chm13_unique_gatk_passing_mask = (
    (sum_stats_var_filters.genome == 'CHM13v2.0') &
    ~sum_stats_var_filters.Syntenic &
    ~sum_stats_var_filters.pass_filter
)
chm13_unique_gatk_passing = sum_stats_var_filters.loc[chm13_unique_gatk_passing_mask].len.sum()
chm13_unique_vqslod_only = sum_stats_var_filters.loc[
    chm13_unique_gatk_passing_mask &
    sum_stats_var_filters.VQSLOD_filter &
    ~sum_stats_var_filters[[c for c in chm13_filter_columns if c != 'VQSLOD_filter']].any(axis=1)
].len.sum()
chm13_unique_retained = sum_stats_var_filters.loc[
    chm13_unique_gatk_passing_mask &
    ~sum_stats_var_filters.CHM13_criteria_fail
].len.sum()
pct_chm13_unique_vqslod_only_of_gatk_passing = chm13_unique_vqslod_only / chm13_unique_gatk_passing * 100
pct_chm13_unique_retained_of_gatk_passing = chm13_unique_retained / chm13_unique_gatk_passing * 100

print(f"CHM13 variants lost to GATK-based filters: {num_variants_lost_to_gatk:,}")
print(f"  Selected manuscript value: {num_variants_lost_to_gatk:,}")
print(f"CHM13 singleton variants gained: {num_singletons_gained:,}")
print(f"  Selected manuscript value: {num_singletons_gained:,}")
print(f"GRCh38 % variants filtered: {pct_grch38_filtered:.1f}%")
print(f"  Selected manuscript value: {pct_grch38_filtered:.1f}%")
print(f"GRCh38 % filtered solely as singletons: {pct_grch38_singleton_only:.1f}%")
print(f"  Selected manuscript value: {pct_grch38_singleton_only:.1f}%")
print(f"CHM13 % variants filtered: {pct_chm13_filtered:.1f}%")
print(f"  Selected manuscript value: {pct_chm13_filtered:.1f}%")
print(f"CHM13 % filtered solely by VQSLOD: {pct_chm13_vqslod_only:.1f}%")
print(f"  Selected manuscript value: {pct_chm13_vqslod_only:.1f}%")
print(f"CHM13-unique GATK PASS variants removed solely by VQSLOD: {pct_chm13_unique_vqslod_only_of_gatk_passing:.1f}% ({chm13_unique_vqslod_only:,}/{chm13_unique_gatk_passing:,})")
print(f"  Selected manuscript value: {pct_chm13_unique_vqslod_only_of_gatk_passing:.1f}%")
print(f"CHM13-unique GATK PASS variants retained after all CHM13 filters: {pct_chm13_unique_retained_of_gatk_passing:.1f}% ({chm13_unique_retained:,}/{chm13_unique_gatk_passing:,})")
print(f"  Selected manuscript value: {pct_chm13_unique_retained_of_gatk_passing:.1f}%")

# %% [markdown]
# ## 4) Panel composition inline statistics
#
# Paper location: Results, haplotype panel section.
#
# Selected after reconciliation:
#
# > "The resulting panel contains **6,404 haplotypes** and **101,152,794** unique variants, including **305,914** high-confidence phased variants in previously unresolved regions of the CHM13v2.0 genome."

# %%
number_unique_CHM13_variants = sum_stats_var_filters.loc[~sum_stats_var_filters.CHM13_criteria_fail & (sum_stats_var_filters.genome=='CHM13v2.0')].len.sum()
number_unique_syntenic_CHM13_variants = sum_stats_var_filters.loc[~sum_stats_var_filters.Syntenic & ~sum_stats_var_filters.CHM13_criteria_fail & (sum_stats_var_filters.genome=='CHM13v2.0')].len.sum()

print(f"Unique CHM13 variants in panel: {number_unique_CHM13_variants:,}")
print(f"  Selected manuscript value: {number_unique_CHM13_variants:,}")
print(f"Unique novel-region variants: {number_unique_syntenic_CHM13_variants:,}")
print(f"  Selected manuscript value: {number_unique_syntenic_CHM13_variants:,}")

# %% [markdown]
# ## 5) Family-structure phasing accuracy inline statistics
#
# Paper location: Results, phasing accuracy section.
#
# > "We observed that samples that underwent Mendelian pre-phasing (Figure 2a, c) had lower SERs than samples that underwent statistical phasing (Figure 2b, d)." 

# %%
trio_vs_nontrio = all_samples.loc[all_samples.ground_truth_data_source.isin(['HPRC_samples', 'HGSVC_samples_nontrios_only']) &
                                    (all_samples.method_of_phasing == 'phased_with_parents_and_pedigree')
                                    ].groupby(['genome', 'ground_truth_data_source'], observed=True).sum(numeric_only=True)

trio_vs_nontrio['switch_error_rate'] = trio_vs_nontrio.n_switch_errors/trio_vs_nontrio.n_checked * 100
trio_vs_nontrio['gt_error_rate'] = trio_vs_nontrio.n_gt_errors/trio_vs_nontrio.n_gt_checked * 100
trio_vs_nontrio['true_switch_error_rate'] = trio_vs_nontrio.n_true_switch_errors/trio_vs_nontrio.n_total_hets * 100
trio_vs_nontrio[['switch_error_rate', 'gt_error_rate', 'true_switch_error_rate']]

# %% [markdown]
# ## 6) 1kGP pedigree composition inline statistics
#
# Paper location: Results, phasing accuracy section.
#
# > "we expect probands (**19%** of 1kGP samples) ... trio parents (**37%** of 1kGP) ... samples that are not part of a trio (**44%**)." 

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

raw_pedigree_parents = len(trio_parents) + len(duo_parents)
raw_pedigree_probands = len(trio_children) + len(duo_children)
raw_pedigree_statistical = 3202 - raw_pedigree_parents - raw_pedigree_probands

# Target categories used for Figure 2/panel-wide weighting. These intentionally use the
# manuscript category definition rather than treating overlapping proband/parent roles as exclusive.
target_probands = 608
target_parents = 1195
target_nontrios = 1399

target_panel_counts = {
    'probands': target_probands,
    'parents': target_parents,
    'nontrios': target_nontrios,
}

# Backwards-compatible names used by the central export.
probands = target_probands
parents = target_parents
statistical = target_nontrios

print("Raw pedigree-derived mutually exclusive counts before reconciliation:")
print(f"  Parents: {raw_pedigree_parents} ({raw_pedigree_parents/3202*100:.0f}% of 1kGP)")
print(f"  Probands: {raw_pedigree_probands} ({raw_pedigree_probands/3202*100:.0f}% of 1kGP)")
print(f"  Statistical: {raw_pedigree_statistical} ({raw_pedigree_statistical/3202*100:.0f}% of 1kGP)")
print("\nTarget panel categories used for Figure 2 and panel-wide weighting:")
print(f"  Probands: {target_probands} ({target_probands/3202*100:.0f}% of 1kGP)")
print(f"  Parents: {target_parents} ({target_parents/3202*100:.0f}% of 1kGP)")
print(f"  Non-trio/statistically phased: {target_nontrios} ({target_nontrios/3202*100:.0f}% of 1kGP)")
print(f"Trio parents to withhold for out-of-panel test: {len(set(trio_parents + duo_parents))}")
print(f"  Paper expects: 1,803")

# %% [markdown]
# ## 7) Panel-wide weighted SER
#
# Paper location: Results, phasing accuracy section.
#
# > "we estimated a panel-wide SER of **0.90%** for the GRCh38 panel and **0.74%** for the T2T-CHM13 panel...chromosome X...T2T-CHM13: **0.42%** SER...GRCh38: **0.91%** SER." 

# %%
# Panel-wide SER: weighted category estimate using target 1kGP panel counts.
# Evaluated categories: probands (HPRC+HGSVC probands), parents, and HGSVC non-trio samples.
panel_wide_category_definitions = pd.DataFrame([
    {
        'panel_category': 'proband',
        'ground_truth_data_source': 'HPRC_HGSVC_probands',
        'target_count': target_probands,
    },
    {
        'panel_category': 'parent',
        'ground_truth_data_source': 'HGSVC_parents',
        'target_count': target_parents,
    },
    {
        'panel_category': 'nontrio',
        'ground_truth_data_source': 'HGSVC_samples_nontrios_only',
        'target_count': target_nontrios,
    },
])

panel_wide_category_rows = []
for _, category in panel_wide_category_definitions.iterrows():
    category_samples = all_samples.loc[
        (all_samples.method_of_phasing == 'phased_with_parents_and_pedigree') &
        (all_samples.ground_truth_data_source == category.ground_truth_data_source)
    ]
    category_stats = group_phasing_statistics(category_samples, ['genome'])
    category_stats['panel_category'] = category.panel_category
    category_stats['ground_truth_data_source'] = category.ground_truth_data_source
    category_stats['target_count'] = int(category.target_count)
    category_stats['evaluated_sample_count'] = category_samples.groupby('genome', observed=True).sample_id.nunique().reindex(category_stats.genome).to_numpy()
    panel_wide_category_rows.append(category_stats)

panel_wide_category_ser = pd.concat(panel_wide_category_rows, ignore_index=True)
panel_wide_ser = (
    panel_wide_category_ser
    .assign(weighted_switch_error_rate=lambda df: df.switch_error_rate * df.target_count)
    .groupby('genome', observed=True)
    .agg(target_count=('target_count', 'sum'),
         weighted_switch_error_rate=('weighted_switch_error_rate', 'sum'))
    .reset_index()
)
panel_wide_ser['switch_error_rate'] = panel_wide_ser.weighted_switch_error_rate / panel_wide_ser.target_count

print("Panel-wide SER category inputs:")
display(panel_wide_category_ser[['genome', 'panel_category', 'ground_truth_data_source', 'target_count', 'evaluated_sample_count', 'switch_error_rate']])
print("\nWeighted panel-wide SER:")
for _, row in panel_wide_ser.iterrows():
    print(f"{row.genome} panel-wide SER: {row.switch_error_rate:.2f}%")
print(f"  Selected manuscript values: GRCh38 0.90%, CHM13 0.74%")

# ChrX SER (females only, by ground truth source)
chrX_ser = all_chroms_chrX.loc[(all_chroms_chrX.sex=='Female') &
                                (all_chroms_chrX.chrom=='chrX') &
                                (all_chroms_chrX.method_of_phasing=='phased_with_parents_and_pedigree') &
                                all_chroms_chrX.ground_truth_data_source.isin(['HPRC_samples','HGSVC_samples_nontrios_only'])]
chrX_ser_grouped = group_phasing_statistics(chrX_ser, ['genome', 'ground_truth_data_source'])
print()
print("ChrX SER by ground truth source:")
for _, row in chrX_ser_grouped.sort_values(['ground_truth_data_source','genome']).iterrows():
    print(f"  {row.ground_truth_data_source:40s} {row.genome}: SER={row.switch_error_rate:.2f}%")
print(f"  Paper expects: GRCh38 0.91%, CHM13 0.42%")

# %% [markdown]
# ## 8) Switch error decomposition (FER, tSER)
#
# Paper location: Results, switch error decomposition section.
#
# > "SER r²=**0.31**, FER r²=**0.35**, tSER r²=**0.006**"
# >
# > "HPRC probands from the 1kGP GRCh38 panel: FER **0.186%**, tSER **0.029%**"
# > "HPRC probands from the T2T-CHM13 panel: FER **0.166%**, tSER **0.018%**"
# >
# > "tSERs were **0.13%** for statistically phased ... **0.018%** for Mendelian pre-phased in T2T-CHM13"
# > "compared to **0.26%** ... **0.029%** in GRCh38"
# >
# > "ten-fold reduction in tSERs on chromosome X: GRCh38 **0.26%** vs T2T-CHM13 **0.016%**" 

# %%
# Error rate correlations (generated in make_plots.ipynb)
error_rate_corr = pd.read_csv(f'{figures_folder}/correlation_error_rates_and_gt_error_rate.csv')
print("Error rate correlations with genotyping error rate:")
print(error_rate_corr.to_string(index=False))
print(f"  Paper expects: SER r²=0.31, FER r²=0.35, tSER r²=0.006")

# %%
# HPRC proband FER and tSER
hprc_proband_stats = all_samples.loc[(all_samples.ground_truth_data_source == 'HPRC_samples') &
                                      (all_samples.method_of_phasing == 'phased_with_parents_and_pedigree')]
hprc_proband_grouped = group_phasing_statistics(hprc_proband_stats, ['genome'])
print("HPRC proband FER and tSER:")
for _, row in hprc_proband_grouped.iterrows():
    print(f"  {row.genome}: FER={row.flip_error_rate:.3f}%, tSER={row.true_switch_error_rate:.3f}%")
print(f"  Paper expects: GRCh38 FER=0.186%, tSER=0.029%; CHM13 FER=0.166%, tSER=0.018%")

# %%
# Statistical vs Mendelian tSER
stat_vs_mend = all_samples.loc[(all_samples.method_of_phasing == 'phased_with_parents_and_pedigree') &
                                all_samples.ground_truth_data_source.isin(['HPRC_samples','HGSVC_samples_nontrios_only'])]
stat_vs_mend_grouped = group_phasing_statistics(stat_vs_mend, ['genome','ground_truth_data_source'])
print("Statistical vs Mendelian tSER:")
for _, row in stat_vs_mend_grouped.iterrows():
    label = 'Mendelian' if row.ground_truth_data_source == 'HPRC_samples' else 'Statistical'
    print(f"  {row.genome} ({label}): tSER={row.true_switch_error_rate:.3f}%")
print(f"  Paper expects: CHM13 stat=0.13%, mend=0.018%; GRCh38 stat=0.26%, mend=0.029%")

# %%
# ChrX FER and tSER (HPRC Mendelian-phased probands)
chrX = all_chroms_chrX.loc[(all_chroms_chrX.sex=='Female') &
                            (all_chroms_chrX.chrom=='chrX') &
                            (all_chroms_chrX.method_of_phasing=='phased_with_parents_and_pedigree') &
                            all_chroms_chrX.ground_truth_data_source.isin(['HPRC_samples', 'HGSVC_samples_nontrios_only'])]
chrX_grouped = group_phasing_statistics(chrX, ['genome', 'ground_truth_data_source'])
print("ChrX FER and tSER (HPRC Mendelian-phased probands):")
for _, row in chrX_grouped.loc[chrX_grouped.ground_truth_data_source=='HPRC_samples'].iterrows():
    print(f"  {row.genome}: FER={row.flip_error_rate:.3f}%, tSER={row.true_switch_error_rate:.3f}%")
print("ChrX FER and tSER (HGSVC Statistically-phased probands):")
for _, row in chrX_grouped.loc[chrX_grouped.ground_truth_data_source=='HGSVC_samples_nontrios_only'].iterrows():
    print(f"  {row.genome}: FER={row.flip_error_rate:.3f}%, tSER={row.true_switch_error_rate:.3f}%")
print(f"  Paper expects: GRCh38 tSER=0.26%, CHM13 tSER=0.016% (ten-fold reduction)")

# %% [markdown]
# ## 8b) Summary and derived statistics
#
# Paper location: Summary, Introduction, and Discussion.
#
# > "Among HPRC samples with published long-read-based assemblies, alignment to T2T-CHM13 resulted in **17%** fewer assembly-discordant genotypes and **13%** fewer phasing switch errors"
# >
# > "true phasing error rates are **50%** lower in the T2T-CHM13 panel overall, and **ten-fold** lower on chromosome X"

# %%
pv_stats

# %%
# Load HPRC-only per-variant concordance values for cross-notebook summary sentence
per_variant_tsv_path = f'{tables_folder}/paper_per_variant_inline_statistics_reference.tsv'
if os.path.exists(per_variant_tsv_path):
    pv_stats = pd.read_csv(per_variant_tsv_path, sep='\t')
    chm13_hprc_disc = pv_stats.loc[pv_stats.Stat_Code == 'C_All_All', 'Value_CHM13'].values[0]
    grch38_hprc_disc = pv_stats.loc[pv_stats.Stat_Code == 'C_All_All', 'Value_GRCh38'].values[0]
    pct_fewer_discordant = (1 - chm13_hprc_disc / grch38_hprc_disc) * 100
else:
    chm13_hprc_disc = np.nan
    grch38_hprc_disc = np.nan
    pct_fewer_discordant = np.nan
    print(f"WARNING: {per_variant_tsv_path} not found. Run calc_per_variant_figures_for_paper.ipynb first.")

# Panel-wide SER values used in the phasing accuracy paragraph
grch38_pw_ser = panel_wide_ser.loc[panel_wide_ser.genome == 'GRCh38', 'switch_error_rate'].values[0]
chm13_pw_ser = panel_wide_ser.loc[panel_wide_ser.genome == 'CHM13v2.0', 'switch_error_rate'].values[0]

# Percent fewer switch errors among HPRC samples with published long-read-based assemblies
grch38_hprc_ser = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'GRCh38', 'switch_error_rate'].values[0]
chm13_hprc_ser = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'CHM13v2.0', 'switch_error_rate'].values[0]
pct_fewer_switch_errors = (1 - chm13_hprc_ser / grch38_hprc_ser) * 100

# tSER improvement (from HPRC proband tSER computed above)
grch38_proband_tser = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'GRCh38', 'true_switch_error_rate'].values[0]
chm13_proband_tser = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'CHM13v2.0', 'true_switch_error_rate'].values[0]
pct_fewer_tser = (1 - chm13_proband_tser / grch38_proband_tser) * 100

# ChrX tSER fold-change
grch38_chrx_tser = chrX_grouped.loc[chrX_grouped.genome == 'GRCh38', 'true_switch_error_rate'].values[0]
chm13_chrx_tser = chrX_grouped.loc[chrX_grouped.genome == 'CHM13v2.0', 'true_switch_error_rate'].values[0]
chrx_tser_fold_change = grch38_chrx_tser / chm13_chrx_tser

print(f"Percent fewer assembly-discordant genotypes (HPRC samples; CHM13 vs GRCh38): {pct_fewer_discordant:.0f}%")
print(f"Percent fewer switch errors (HPRC samples; SER): {pct_fewer_switch_errors:.0f}%")
print(f"Percent lower tSER (HPRC probands): {pct_fewer_tser:.0f}%")
print(f"ChrX tSER fold-change: {chrx_tser_fold_change:.0f}-fold")


# %%
# Load HGSVC-statistically-phased-only per-variant concordance values for cross-notebook summary sentence
hgsvc_stat_grouped = group_phasing_statistics(all_samples.loc[(all_samples.method_of_phasing=='phased_with_parents_and_pedigree')
                & (all_samples.ground_truth_data_source=='HGSVC_samples_nontrios_only')], ['genome'])
# Percent fewer switch errors among HPRC samples with published long-read-based assemblies
grch38_hgsvc_stat_ser = hgsvc_stat_grouped.loc[hgsvc_stat_grouped.genome == 'GRCh38', 'switch_error_rate'].values[0]
chm13_hgsvc_stat_ser = hgsvc_stat_grouped.loc[hgsvc_stat_grouped.genome == 'CHM13v2.0', 'switch_error_rate'].values[0]
pct_fewer_switch_errors = (1 - chm13_hgsvc_stat_ser / grch38_hgsvc_stat_ser) * 100

# tSER improvement (from HPRC proband tSER computed above)
grch38_hgsvc_stat_tser = hgsvc_stat_grouped.loc[hgsvc_stat_grouped.genome == 'GRCh38', 'true_switch_error_rate'].values[0]
chm13_hgsvc_stat_tser = hgsvc_stat_grouped.loc[hgsvc_stat_grouped.genome == 'CHM13v2.0', 'true_switch_error_rate'].values[0]
pct_fewer_tser = (1 - chm13_hgsvc_stat_tser / grch38_hgsvc_stat_tser) * 100

# ChrX tSER fold-change
grch38_chrx_tser = chrX_grouped.loc[chrX_grouped.genome == 'GRCh38', 'true_switch_error_rate'].values[0]
chm13_chrx_tser = chrX_grouped.loc[chrX_grouped.genome == 'CHM13v2.0', 'true_switch_error_rate'].values[0]
chrx_tser_fold_change = grch38_chrx_tser / chm13_chrx_tser

print(f"Percent fewer assembly-discordant genotypes (HPRC samples; CHM13 vs GRCh38): {pct_fewer_discordant:.0f}%")
print(f"Percent fewer switch errors (HPRC samples; SER): {pct_fewer_switch_errors:.0f}%")
print(f"Percent lower tSER (HPRC probands): {grch38_hgsvc_stat_tser:.3f}% → {chm13_hgsvc_stat_tser:.3f}% ({pct_fewer_tser:.0f}%)")
print(f"ChrX tSER fold-change: {chrx_tser_fold_change:.0f}-fold")


# %%
all_samples.ground_truth_data_source.unique()

# %%
tser=all_samples.loc[(all_samples.method_of_phasing == 'phased_with_parents_and_pedigree')&(all_samples.ground_truth_data_source.isin(['HPRC_samples','HGSVC_samples_nontrios_only']))].pivot(index=['sample_id', 'ground_truth_data_source'], columns='genome', values='true_switch_error_rate').reset_index()
tser['percent_reduction_tser'] = ((1 - tser['CHM13v2.0'] / tser['GRCh38']) * 100).abs()
tser.sort_values('CHM13v2.0', ascending=True)#[['GRCh38', 'CHM13v2.0', 'percent_reduction_tser']]
mean_HPRC_GRCh38=tser.loc[tser.ground_truth_data_source=='HPRC_samples'].GRCh38.mean()
tser.loc[(tser['CHM13v2.0'] >= mean_HPRC_GRCh38)&(tser.ground_truth_data_source=='HGSVC_samples_nontrios_only')][['sample_id', 'ground_truth_data_source', 'GRCh38', 'CHM13v2.0', 'percent_reduction_tser']].sort_values('CHM13v2.0',ascending=True)
tser.loc[tser.ground_truth_data_source=='HPRC_samples'].GRCh38 > 0.049756451557635105
compare_stats = all_samples.loc[all_samples.ground_truth_data_source.isin(('HPRC_samples','HGSVC_probands','HGSVC_parents','HGSVC_samples_nontrios_only')) &
                                (all_samples.method_of_phasing=='phased_with_parents_and_pedigree')
                               ].rename(columns={'true_switch_error_rate':'True SER', 'flip_error_rate':'Flip rate', 'flip_and_switches_SER':'Total SER', 'genome':'Genome'}
                               ).melt(value_vars=['True SER','Flip rate','Total SER'],
                                       id_vars=['sample_id', 'Genome', 'population', 'superpopulation', 'sex','trio_phased', 'method_of_phasing', 'ground_truth_data_source'],
                                      var_name='Measure of Phasing Error',
                                      value_name='Error Rate').replace('CHM13v2.0', 'T2T-CHM13')
compare_stats

# %% [markdown]
# ## 8c) Plug-and-play manuscript text
#
# The sentences below use the exact manuscript sentence structure with computed values inserted.
# Copy-paste these into the paper and use the surrounding context to ctrl-F to the correct section.

# %%
# --- Variant filtering (manuscript lines 231-236) ---
print("=== VARIANT FILTERING ===")
print(f"Our VQSLOD filter resulted in the exclusion of an additional {num_variants_lost_to_gatk:,} variants "
      f"compared to the filtering methods used by Byrska-Bishop et al. to produce the GRCh38 1kGP panel.")
print()
print(f"The inclusion of singleton alleles added {num_singletons_gained:,} unique variants to our dataset.")
print()
print(f"Overall, variant filtering had removed {pct_grch38_filtered:.1f}% of variants in the "
      f"Byrska-Bishop et al. 1kGP GRCh38 panel, including {pct_grch38_singleton_only:.1f}% of all variants "
      f"that were filtered out solely because they were singletons.")
print()
print(f"In contrast, our filtering strategy removed {pct_chm13_filtered:.0f}% of 1kGP T2T-CHM13 variants, "
      f"with {pct_chm13_vqslod_only:.1f}% of variants filtered solely due to the VQSLOD filter "
      f"(Supplementary Figure 1).")
print()
print(f"Our stricter VQSLOD filter removes {pct_chm13_unique_vqslod_only_of_gatk_passing:.1f}% of GATK-passing "
      f"variants in T2T-CHM13-unique regions. When all filters are considered, only "
      f"{pct_chm13_unique_retained_of_gatk_passing:.1f}% of GATK-passing variants in these regions are retained "
      f"(Supplementary Figure 1).")
print()

# --- Panel composition (manuscript lines 237-244) ---
print("=== PANEL COMPOSITION ===")
print(f"The resulting panel contains 6,404 haplotypes and {number_unique_CHM13_variants:,} unique variants, "
      f"including {number_unique_syntenic_CHM13_variants:,} high-confidence phased variants in previously "
      f"unresolved regions of the T2T-CHM13 genome.")
print()

# --- Pedigree composition (manuscript lines 268-280) ---
print("=== PEDIGREE COMPOSITION ===")
print(f"As a result, we expect probands ({target_probands/3202*100:.0f}% of 1kGP samples) and trio parents "
      f"({target_parents/3202*100:.0f}% of 1kGP) to be more accurately phased than samples that are not "
      f"part of a trio ({target_nontrios/3202*100:.0f}%).")
print()

# --- Panel-wide SER (manuscript lines 268-280) ---
print("=== PANEL-WIDE SER ===")
print(f"Weighting the observed proband/parent/non-trio SER by panel membership "
      f"({target_probands:,} probands, {target_parents:,} parents, and {target_nontrios:,} non-trio "
      f"individuals) allowed us to estimate a {target_probands + target_parents + target_nontrios:,}-member "
      f"panel-wide SER of {grch38_pw_ser:.2f}% for the GRCh38 panel and {chm13_pw_ser:.2f}% for the "
      f"T2T-CHM13 panel.")
print()

# --- ChrX SER (manuscript lines 282-290) ---
print("=== ChrX SER ===")
chrx_ser_chm13 = chrX_ser_grouped.loc[chrX_ser_grouped.genome == 'CHM13v2.0', 'switch_error_rate'].values[0]
chrx_ser_grch38 = chrX_ser_grouped.loc[chrX_ser_grouped.genome == 'GRCh38', 'switch_error_rate'].values[0]
print(f"In statistically-phased HGSVC non-trio samples, chromosome X SER decreased from "
      f"{chrx_ser_grch38:.2f}% in the GRCh38 panel to {chrx_ser_chm13:.2f}% in the T2T-CHM13 panel.")
print()

# --- Switch error decomposition (manuscript lines 293-321) ---
print("=== SWITCH ERROR DECOMPOSITION ===")

# Error rate correlations
corr_rows_mend = error_rate_corr.loc[error_rate_corr['Mendelian-based error correction'].astype(str).eq('True')]
ser_r2 = corr_rows_mend.loc[corr_rows_mend['Log error rate'] == 'Log(Switch Error Rate)'].iloc[0]
fer_r2 = corr_rows_mend.loc[corr_rows_mend['Log error rate'] == 'Log(Flip Error Rate)'].iloc[0]
tser_r2 = corr_rows_mend.loc[corr_rows_mend['Log error rate'] == 'Log(True Switch Error Rate)'].iloc[0]
print(f"In our dataset, the majority of switch errors were associated with flips. FERs were correlated with "
      f"genotyping error rates, while tSERs were uncorrelated with genotyping error rates "
      f"(SER r2 = {ser_r2['CHM13v2.0']:.2f}, FER r2 = {fer_r2['CHM13v2.0']:.2f}, "
      f"tSER = {tser_r2['CHM13v2.0']:.3f}, Pearson correlation, Mendelian pre-phased samples), "
      f"supporting the interpretation that flips stem from genotyping errors (Extended Data Figure 3).")
print()

# HPRC proband FER and tSER
grch38_fer = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'GRCh38', 'flip_error_rate'].values[0]
chm13_fer = hprc_proband_grouped.loc[hprc_proband_grouped.genome == 'CHM13v2.0', 'flip_error_rate'].values[0]
print(f"The measured flip error rate for HPRC probands from the 1kGP GRCh38 panel was {grch38_fer:.3f}%, "
      f"while the tSER was {grch38_proband_tser:.3f}%. HPRC probands from the T2T-CHM13 panel had a flip "
      f"error rate of {chm13_fer:.3f}% and a tSER of {chm13_proband_tser:.3f}% (Figure 3b).")
print()

# Statistical vs Mendelian tSER
stat_mend = stat_vs_mend_grouped.reset_index() if isinstance(stat_vs_mend_grouped.index, pd.MultiIndex) else stat_vs_mend_grouped
chm13_stat_tser = stat_mend.loc[(stat_mend.genome == 'CHM13v2.0') & (stat_mend.ground_truth_data_source == 'HGSVC_samples_nontrios_only'), 'true_switch_error_rate'].values[0]
chm13_mend_tser = stat_mend.loc[(stat_mend.genome == 'CHM13v2.0') & (stat_mend.ground_truth_data_source == 'HPRC_samples'), 'true_switch_error_rate'].values[0]
grch38_stat_tser = stat_mend.loc[(stat_mend.genome == 'GRCh38') & (stat_mend.ground_truth_data_source == 'HGSVC_samples_nontrios_only'), 'true_switch_error_rate'].values[0]
grch38_mend_tser = stat_mend.loc[(stat_mend.genome == 'GRCh38') & (stat_mend.ground_truth_data_source == 'HPRC_samples'), 'true_switch_error_rate'].values[0]
print(f"Specifically, tSERs were {chm13_stat_tser:.2f}% for statistically-phased samples and "
      f"{chm13_mend_tser:.3f}% for Mendelian pre-phased samples in the T2T-CHM13 panel, compared to "
      f"{grch38_stat_tser:.2f}% and {grch38_mend_tser:.3f}% in the GRCh38 panel.")
print()

# ChrX tSER
print(f"When stratified by chromosome, we observed a {chrx_tser_fold_change:.0f}-fold reduction in tSERs "
      f"on the chromosome X (GRCh38: {grch38_chrx_tser:.2f}% tSER; T2T-CHM13: {chm13_chrx_tser:.3f}% tSER).")
print()

# --- Summary / Introduction / Discussion (manuscript Summary and lines 119-125) ---
print("=== SUMMARY / INTRODUCTION / DISCUSSION ===")
print(f"Among HPRC samples with published long-read-based assemblies, alignment to T2T-CHM13 resulted in "
      f"{pct_fewer_discordant:.0f}% fewer assembly-discordant genotypes and "
      f"{pct_fewer_switch_errors:.0f}% fewer phasing switch errors")
print()
print(f"Among HPRC samples with published long-read-based assemblies, true phasing error rates are "
      f"{pct_fewer_tser:.0f}% lower in the T2T-CHM13 panel genome-wide and "
      f"{chrx_tser_fold_change:.0f}-fold lower on chromosome X")


# %% [markdown]
# ## 9) Export per-sample statistics for supplemental table 1

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
os.makedirs(f'{supplemental_figures_folder}/tables', exist_ok=True)
source_accuracy[source_accuracy_columns].to_csv(f'{supplemental_figures_folder}/tables/supplemental_ground_truth_source_accuracy_data_unformatted.tsv', sep='\t', index=False)
samples_for_export[sample_export_columns].to_csv(f'{supplemental_figures_folder}/tables/supplemental_per_sample_data_unformatted.tsv', sep='\t', index=False)
source_accuracy.T

# %% [markdown]
# ## 10) HPRC panel-phased summary export

# %%
HPRC_samples_phased_with_panel = binned_maf_data.loc[(binned_maf_data.ground_truth_data_source=='HPRC_samples')&(binned_maf_data.method_of_phasing == '1kgp_variation_phased_with_reference_panel')]
HPRC_samples_phased_with_panel = group_phasing_statistics(HPRC_samples_phased_with_panel, ['type','genome','Syntenic']).drop(columns=['MAC','AN','MAF'])
error_rate_corr.to_csv(f'{figures_folder}/panel_phased_HPRC_samples_summary_data.csv', index=False)

# %% [markdown]
# ## 11) Central export of inline statistics (raw numeric values)
#
# This table is intended as the long-term source for manuscript value injection.
#
# Paper refs covered in this export include Table 1/genotyping concordance (manuscript lines 208-216), variant filtering (manuscript lines 231-236), panel-wide SER (manuscript lines 268-280), FER/tSER (manuscript lines 293-321), and singleton/imputation summary values (manuscript lines 415-422).
#

# %%
paper_inline_rows = []

def add_stat(section, statistic_description, value_chm13=np.nan, value_grch38=np.nan,
             percent_change=np.nan, paper_text='', paper_location='',
             source_notebook='calc_figures_for_paper', notes=''):
    paper_inline_rows.append({
        'Section': section,
        'Statistic_Description': statistic_description,
        'Value_CHM13': value_chm13,
        'Value_GRCh38': value_grch38,
        'Percent_Change': percent_change,
        'Paper_Text': paper_text,
        'Paper_Location': paper_location,
        'Source_Notebook': source_notebook,
        'Notes': notes,
    })

# -- Variant filtering and panel composition --
variant_filter_note = 'Selected denominator: filtered variants / total variants for the same genome in filter_summary_stats.'
add_stat('Variant filtering', 'CHM13 variants lost to GATK-based filters',
         value_chm13=float(num_variants_lost_to_gatk),
         paper_text=f'Our VQSLOD filter resulted in the exclusion of an additional {num_variants_lost_to_gatk:,} variants compared to the filtering methods used by Byrska-Bishop et al. to produce the GRCh38 1kGP panel.',
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)
add_stat('Variant filtering', 'CHM13 singleton variants gained',
         value_chm13=float(num_singletons_gained),
         paper_text=f'The inclusion of singleton alleles added {num_singletons_gained:,} unique variants to our dataset.',
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)
add_stat('Variant filtering', 'GRCh38 % variants filtered',
         value_grch38=pct_grch38_filtered,
         paper_text=f'Overall, variant filtering had removed {pct_grch38_filtered:.1f}% of variants in the Byrska-Bishop et al. 1kGP GRCh38 panel, including {pct_grch38_singleton_only:.1f}% of all variants that were filtered out solely because they were singletons.',
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)
add_stat('Variant filtering', 'GRCh38 % filtered solely as singletons',
         value_grch38=pct_grch38_singleton_only,
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)
add_stat('Variant filtering', 'CHM13 % variants filtered',
         value_chm13=pct_chm13_filtered,
         paper_text=f'In contrast, our filtering strategy removed {pct_chm13_filtered:.0f}% of 1kGP T2T-CHM13 variants, with {pct_chm13_vqslod_only:.1f}% of variants filtered solely due to the VQSLOD filter (Supplementary Figure 1).',
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)
add_stat('Variant filtering', 'CHM13 % filtered solely by VQSLOD',
         value_chm13=pct_chm13_vqslod_only,
         paper_location='Results: variant QC paragraph',
         notes=variant_filter_note)

chm13_unique_filter_note = (
    'Denominator: GATK PASS (pass_filter=False) variants in CHM13v2.0 nonsyntenic regions. '
    'VQSLOD-only uses the CHM13 filter-criteria list; AC_filter is not a CHM13 criterion.'
)
add_stat('CHM13-unique filtering', 'CHM13-unique GATK PASS variants (denominator)',
         value_chm13=float(chm13_unique_gatk_passing),
         paper_location='Results: CHM13-unique filtering paragraph',
         notes=chm13_unique_filter_note)
add_stat('CHM13-unique filtering', 'CHM13-unique % GATK PASS variants removed solely by VQSLOD',
         value_chm13=pct_chm13_unique_vqslod_only_of_gatk_passing,
         paper_text=f'Our stricter VQSLOD filter removes {pct_chm13_unique_vqslod_only_of_gatk_passing:.1f}% of GATK-passing variants in T2T-CHM13-unique regions. When all filters are considered, only {pct_chm13_unique_retained_of_gatk_passing:.1f}% of GATK-passing variants in these regions are retained (Supplementary Figure 1).',
         paper_location='Results: CHM13-unique filtering paragraph',
         notes=f'{chm13_unique_vqslod_only:,}/{chm13_unique_gatk_passing:,}; {chm13_unique_filter_note}')
add_stat('CHM13-unique filtering', 'CHM13-unique % GATK PASS variants retained after all CHM13 filters',
         value_chm13=pct_chm13_unique_retained_of_gatk_passing,
         paper_location='Results: CHM13-unique filtering paragraph',
         notes=f'{chm13_unique_retained:,}/{chm13_unique_gatk_passing:,}; {chm13_unique_filter_note}')

add_stat('Panel composition', 'Unique CHM13 variants in panel',
         value_chm13=float(number_unique_CHM13_variants),
         paper_text=f'The resulting panel contains 6,404 haplotypes and {number_unique_CHM13_variants:,} unique variants, including {number_unique_syntenic_CHM13_variants:,} high-confidence phased variants in previously unresolved regions of the T2T-CHM13 genome.',
         paper_location='Results: panel composition paragraph',
         notes='Selected value after reconciliation; count is CHM13 criteria-pass variants in filter_summary_stats.')
add_stat('Panel composition', 'Unique CHM13 novel-region variants in panel',
         value_chm13=float(number_unique_syntenic_CHM13_variants),
         paper_location='Results: panel composition paragraph',
         notes='Selected value after reconciliation; nonsyntenic CHM13 criteria-pass variants.')

# -- Pedigree composition --
pedigree_note = 'Target 1kGP category count used for Figure 2/panel-wide weighting; sample roles can overlap in the raw pedigree.'
add_stat('Pedigree composition', 'Number of trio/duo probands',
         value_chm13=float(target_probands),
         paper_text=f'As a result, we expect probands ({target_probands/3202*100:.0f}% of 1kGP samples) and trio parents ({target_parents/3202*100:.0f}% of 1kGP) to be more accurately phased than samples that are not part of a trio ({target_nontrios/3202*100:.0f}%).',
         paper_location='Results: phasing accuracy paragraph',
         notes=pedigree_note)
add_stat('Pedigree composition', 'Number of trio/duo parents',
         value_chm13=float(target_parents),
         paper_location='Results: phasing accuracy paragraph',
         notes=pedigree_note)
add_stat('Pedigree composition', 'Number of statistically phased samples',
         value_chm13=float(target_nontrios),
         paper_location='Results: phasing accuracy paragraph',
         notes=pedigree_note)

# -- Panel-wide SER --
for _, row in panel_wide_category_ser.iterrows():
    add_stat('Panel-wide SER components',
             f"Panel-wide SER component ({row.genome}; {row.panel_category})",
             value_chm13=float(row.switch_error_rate) if 'CHM13' in str(row.genome) else np.nan,
             value_grch38=float(row.switch_error_rate) if 'GRCh38' in str(row.genome) else np.nan,
             paper_location='Results: phasing accuracy paragraph / Figure 2 estimated panel-wide row',
             notes=(f"ground_truth_data_source={row.ground_truth_data_source}; "
                    f"target_count={int(row.target_count)}; evaluated_sample_count={int(row.evaluated_sample_count)}"))

for _, row in panel_wide_ser.iterrows():
    add_stat('Panel-wide SER', f'Weighted panel-wide SER ({row.genome})',
             value_chm13=float(row.switch_error_rate) if 'CHM13' in str(row.genome) else np.nan,
             value_grch38=float(row.switch_error_rate) if 'GRCh38' in str(row.genome) else np.nan,
             paper_text=f'Weighting the observed proband/parent/non-trio SER by panel membership ({target_probands:,} probands, {target_parents:,} parents, and {target_nontrios:,} non-trio individuals) allowed us to estimate a {target_probands + target_parents + target_nontrios:,}-member panel-wide SER of {grch38_pw_ser:.2f}% for the GRCh38 panel and {chm13_pw_ser:.2f}% for the T2T-CHM13 panel.' if 'GRCh38' in str(row.genome) else '',
             paper_location='Results: phasing accuracy paragraph',
             notes='Weighted average of proband, parent, and non-trio category SER using target counts.')

# -- ChrX SER --
for _, row in chrX_ser_grouped.iterrows():
    g = row.genome
    add_stat('ChrX accuracy', f'ChrX SER ({g})',
             value_chm13=float(row.switch_error_rate) if 'CHM13' in str(g) else np.nan,
             value_grch38=float(row.switch_error_rate) if 'GRCh38' in str(g) else np.nan,
             paper_text=f'In statistically-phased HGSVC non-trio samples, chromosome X SER decreased from {chrX_ser_grouped.loc[chrX_ser_grouped.genome == "GRCh38", "switch_error_rate"].values[0]:.2f}% in the GRCh38 panel to {chrX_ser_grouped.loc[chrX_ser_grouped.genome == "CHM13v2.0", "switch_error_rate"].values[0]:.2f}% in the T2T-CHM13 panel.' if 'GRCh38' in str(g) else '',
             paper_location='Results: phasing accuracy paragraph')

# -- Switch error decomposition --
for _, row in hprc_proband_grouped.iterrows():
    g = row.genome
    add_stat('Switch error decomposition', f'HPRC proband SER ({g})',
             value_chm13=float(row.switch_error_rate) if 'CHM13' in str(g) else np.nan,
             value_grch38=float(row.switch_error_rate) if 'GRCh38' in str(g) else np.nan,
             paper_location='Summary; Introduction; Discussion',
             notes='HPRC-only switch error rate used for the headline phasing switch-error reduction.')
    add_stat('Switch error decomposition', f'HPRC proband FER ({g})',
             value_chm13=float(row.flip_error_rate) if 'CHM13' in str(g) else np.nan,
             value_grch38=float(row.flip_error_rate) if 'GRCh38' in str(g) else np.nan,
             paper_location='Results: switch error decomposition')
    add_stat('Switch error decomposition', f'HPRC proband tSER ({g})',
             value_chm13=float(row.true_switch_error_rate) if 'CHM13' in str(g) else np.nan,
             value_grch38=float(row.true_switch_error_rate) if 'GRCh38' in str(g) else np.nan,
             paper_text=f'The measured flip error rate for HPRC probands from the 1kGP GRCh38 panel was {grch38_fer:.3f}%, while the tSER was {grch38_proband_tser:.3f}%. HPRC probands from the T2T-CHM13 panel had a flip error rate of {chm13_fer:.3f}% and a tSER of {chm13_proband_tser:.3f}% (Figure 3b).' if 'GRCh38' in str(g) else '',
             paper_location='Results: switch error decomposition')

# -- Switch-error correlation / regression values --
if isinstance(error_rate_corr, pd.DataFrame) and not error_rate_corr.empty:
    corr_export_labels = {
        'Log(Switch Error Rate)': 'SER',
        'Log(Flip Error Rate)': 'FER',
        'Log(True Switch Error Rate)': 'tSER',
    }
    corr_rows = error_rate_corr.loc[error_rate_corr['Mendelian-based error correction'].astype(str).eq('True')]
    for _, row in corr_rows.iterrows():
        label = corr_export_labels.get(row['Log error rate'], row['Log error rate'])
        add_stat(
            'Switch-error correlation',
            f'Pearson r^2 genotype-error vs {label} (HPRC probands, Mendelian-corrected)',
            value_chm13=float(row['CHM13v2.0']),
            value_grch38=float(row['GRCh38']),
            paper_location='Results: switch error decomposition',
            notes=('Separate CSV export: figures_whole_genome/correlation_error_rates_and_gt_error_rate.csv; '
                   'Pearson r^2 from log-transformed rates.'),
        )

# -- ChrX tSER --
for _, row in chrX_grouped.iterrows():
    g = row.genome
    add_stat('ChrX accuracy', f'ChrX tSER HPRC ({g})',
             value_chm13=float(row.true_switch_error_rate) if 'CHM13' in str(g) else np.nan,
             value_grch38=float(row.true_switch_error_rate) if 'GRCh38' in str(g) else np.nan,
             paper_text=f'When stratified by chromosome, we observed a {chrx_tser_fold_change:.0f}-fold reduction in tSERs on the chromosome X (GRCh38: {grch38_chrx_tser:.2f}% tSER; T2T-CHM13: {chm13_chrx_tser:.3f}% tSER).' if 'GRCh38' in str(g) else '',
             paper_location='Results: switch error decomposition')

# -- Summary / derived --
add_stat('Summary', 'Percent fewer assembly-discordant genotypes (HPRC samples; CHM13 vs GRCh38)',
         value_chm13=pct_fewer_discordant,
         paper_text=f'Among HPRC samples with published long-read-based assemblies, alignment to T2T-CHM13 resulted in {pct_fewer_discordant:.0f}% fewer assembly-discordant genotypes and {pct_fewer_switch_errors:.0f}% fewer phasing switch errors',
         paper_location='Summary; Introduction; Discussion',
         notes='Derived from HPRC-only overall concordance (per-variant TSV) and HPRC-only SER.')
add_stat('Summary', 'Percent fewer switch errors (HPRC samples; SER)',
         value_chm13=pct_fewer_switch_errors,
         paper_location='Summary; Introduction; Discussion',
         notes='Derived from HPRC-only SER using samples with published long-read-based assemblies.')
add_stat('Summary', 'Percent lower tSER (HPRC probands)',
         value_chm13=pct_fewer_tser,
         paper_text=f'Among HPRC samples with published long-read-based assemblies, true phasing error rates are {pct_fewer_tser:.0f}% lower in the T2T-CHM13 panel genome-wide and {chrx_tser_fold_change:.0f}-fold lower on chromosome X',
         paper_location='Summary; Discussion',
         notes='Derived from HPRC proband tSER and ChrX tSER.')

paper_inline_statistics_df = pd.DataFrame(paper_inline_rows)
os.makedirs('tables', exist_ok=True)
paper_inline_statistics_path = f'{tables_folder}/paper_inline_statistics_reference.tsv'
paper_inline_statistics_df.to_csv(paper_inline_statistics_path, sep='\t', index=False)

print(f'Wrote: {paper_inline_statistics_path}')
paper_inline_statistics_df


# %% [markdown]
# ### Overlap between known collapsed GRCh38 duplicates and regions of high error

# %%
import polars as pl

rolling_window_ideogram_500k = pl.read_parquet(f'{summary_statistics_folder}/rolling_stats_500k_window.parquet')
rolling_window_ideogram_500k = rolling_window_ideogram_500k.filter((pl.col('ground_truth_data_source')=='HPRC_samples') & (pl.col('method_of_phasing')=='phased_with_parents_and_pedigree'))
rolling_window_ideogram_500k = rolling_window_ideogram_500k.to_pandas(use_pyarrow_extension_array=True)

collapsed = pd.read_csv('resources/GRCh38_collapsed_duplication_FP_regions.bed.gz', sep='\t', header=None, names=['chrom','start','end'])
collapsed

# %%
from intervaltree import IntervalTree
from pyliftover import LiftOver

def get_contiguous_regions(df):
    total_seq=0
    regions = []
    for chrom, group in df.groupby('chrom', observed=True):
        start=0
        end=0
        region_start=0
        region_end=0
        for i, row in group.iterrows():
            if region_start==0:
                region_start = row['start']
                region_end = row['end']
                continue
            start = row['start']
            end = row['end']
            if start > region_end and (region_start != region_end):
                total_seq += (region_end - region_start)
                regions.append({'chrom': chrom, 'start': region_start, 'end': region_end})
                region_start = start
                region_end = end
            else:
                region_end = end
        regions.append({'chrom': chrom, 'start': region_start, 'end': region_end})
    return pd.DataFrame(regions)

def build_interval_tree(df):
    trees = {}
    for chrom in df['chrom'].unique():
        tree = IntervalTree()
        for _, row in df[df['chrom'] == chrom].iterrows():
            tree[row['start']:row['end']] = True
        trees[chrom] = tree
    return trees

def total_bases(df):
    return (df['end'] - df['start']).sum()

def overlap_trees(tree1, tree2):
    overlap = 0
    for chrom in tree1:
        if chrom in tree2:
            for iv1 in tree1[chrom]:
                for iv2 in tree2[chrom].overlap(iv1.begin, iv1.end):
                    overlap += min(iv1.end, iv2.end) - max(iv1.begin, iv2.begin)
    return overlap

chainfile = 'resources/chm13v2-grch38.chain'
lo = LiftOver(chainfile)
null_location={'lifted_chrom': None, 'lifted_start': None, 'lifted_end': None}

def liftover_coords(row, fuzzy=(10000, 500000)):
    if row['genome'] == 'CHM13v2.0':
        lifted_start = lo.convert_coordinate(row['chrom'], row['start'])
        lifted_end = lo.convert_coordinate(row['chrom'], row['end'])
        if lifted_start == [] or lifted_end == [] or lifted_start[0][0] != lifted_end[0][0]:
            return null_location
        else:
            lifted_chrom = lifted_start[0][0]
            lifted_start= lifted_start[0][1]
            lifted_end = lifted_end[0][1]
        if fuzzy is not None and type(fuzzy) == tuple and len(fuzzy) == 2:
            step, desired_len = fuzzy
            offset = lifted_start % step
            if offset < step / 2:
                rounded_start = lifted_start - offset
            else:
                rounded_start = lifted_start + (step - offset)
            rounded_end = rounded_start + desired_len
            if rounded_end-lifted_end > desired_len//2:
                return null_location
            else:
                lifted_start = rounded_start
                lifted_end = rounded_end
        if lifted_start is not None and lifted_end is not None:
            return pd.Series({
                'lifted_chrom': lifted_chrom,
                'lifted_start': lifted_start,
                'lifted_end': lifted_end
            })
    return pd.Series(null_location)


# %%
rolling_window_ideogram_500k['gt_error_rate'] = rolling_window_ideogram_500k.n_gt_errors/rolling_window_ideogram_500k.n_gt_checked * 100
rolling_window_ideogram_500k[['lifted_chrom', 'lifted_start', 'lifted_end']] = rolling_window_ideogram_500k.apply(liftover_coords, axis=1)

rolling_window_ideogram_500k_wide = (
    rolling_window_ideogram_500k.pivot(
        index=["chrom", "start", "end"],
        columns="genome",
        values="gt_error_rate"
    )
    .add_prefix("gt_error_rate_")
    .reset_index()
)

# %%
collapsed_sorted = get_contiguous_regions(collapsed.sort_values(['chrom', 'start', 'end']).reset_index(drop=True))
collapsed_trees = build_interval_tree(collapsed_sorted)
collapsed_total = total_bases(collapsed_sorted)
grch38_rolling_window_ideogram = rolling_window_ideogram_500k_wide.rename(columns={'gt_error_rate_GRCh38':'gt_error_rate'})

error_gap_overlaps = []
for error_rate_x10 in range(0, 150, 1):
    error_rate = error_rate_x10 / 10
    high_error = grch38_rolling_window_ideogram.loc[
        (grch38_rolling_window_ideogram.gt_error_rate - grch38_rolling_window_ideogram['gt_error_rate_CHM13v2.0'] > error_rate)
    ].sort_values(['chrom', 'start', 'end'])

    high_error_regions_sorted = get_contiguous_regions(high_error).sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    high_error_trees = build_interval_tree(high_error_regions_sorted)
    high_error_total = total_bases(high_error_regions_sorted)

    overlap = overlap_trees(high_error_trees, collapsed_trees)

    pct_high_error_covered_by_collapsed = 100 * overlap / high_error_total if high_error_total else 0
    pct_collapsed_covered_by_high_error = 100 * overlap / collapsed_total if collapsed_total else 0

    error_gap_overlaps.append({'error_rate_threshold': error_rate,
                     'high_error_total_bases': high_error_total,
                     'collapsed_total_bases': collapsed_total,
                     'overlap_bases': overlap,
                     'pct_high_error_covered_by_collapsed': pct_high_error_covered_by_collapsed,
                     'pct_collapsed_covered_by_high_error': pct_collapsed_covered_by_high_error})

high_error_overlaps = []
for error_rate_x10 in range(0, 150, 1):
    error_rate = error_rate_x10 / 10
    high_error = grch38_rolling_window_ideogram.loc[
        (grch38_rolling_window_ideogram.gt_error_rate > error_rate)
    ].sort_values(['chrom', 'start', 'end'])

    high_error_regions_sorted = get_contiguous_regions(high_error).sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    high_error_trees = build_interval_tree(high_error_regions_sorted)
    high_error_total = total_bases(high_error_regions_sorted)

    overlap = overlap_trees(high_error_trees, collapsed_trees)

    pct_high_error_covered_by_collapsed = 100 * overlap / high_error_total if high_error_total else 0
    pct_collapsed_covered_by_high_error = 100 * overlap / collapsed_total if collapsed_total else 0

    high_error_overlaps.append({'error_rate_threshold': error_rate,
                     'high_error_total_bases': high_error_total,
                     'collapsed_total_bases': collapsed_total,
                     'overlap_bases': overlap,
                     'pct_high_error_covered_by_collapsed': pct_high_error_covered_by_collapsed,
                     'pct_collapsed_covered_by_high_error': pct_collapsed_covered_by_high_error})

print(f"High error regions covered by collapsed regions: {pct_high_error_covered_by_collapsed:.2f}%")
print(f"Collapsed regions covered by high error regions: {pct_collapsed_covered_by_high_error:.2f}%")

# %%
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
sns.lineplot(data=pd.DataFrame(high_error_overlaps).sort_values('error_rate_threshold'),
             x='pct_collapsed_covered_by_high_error', y='pct_high_error_covered_by_collapsed', ax=ax[0])
sns.lineplot(data=pd.DataFrame(error_gap_overlaps).sort_values('error_rate_threshold'),
             x='pct_collapsed_covered_by_high_error', y='pct_high_error_covered_by_collapsed', ax=ax[1])
plt.tight_layout()

# %%
pd.DataFrame(error_gap_overlaps)

# %%
pd.DataFrame(high_error_overlaps).loc[pd.DataFrame(high_error_overlaps).error_rate_threshold==3.0]

# %% [markdown]
# ### Reviewer-response export: collapsed duplication overlap
#
# This named export supports reviewer-response text about overlap between collapsed GRCh38 duplication calls and high-error windows. The threshold sweep is exported so the exact manuscript threshold can be audited rather than inferred from a plot.
#

# %%
collapsed_duplication_overlap_df = pd.concat(
    [
        pd.DataFrame(error_gap_overlaps).assign(window_set='grch38_gt_error_minus_chm13_gt_error'),
        pd.DataFrame(high_error_overlaps).assign(window_set='grch38_absolute_gt_error'),
    ],
    ignore_index=True,
)

collapsed_duplication_overlap_df = collapsed_duplication_overlap_df[
    [
        'window_set',
        'error_rate_threshold',
        'high_error_total_bases',
        'collapsed_total_bases',
        'overlap_bases',
        'pct_high_error_covered_by_collapsed',
        'pct_collapsed_covered_by_high_error',
    ]
]

collapsed_duplication_overlap_path = f'{tables_folder}/collapsed_duplication_high_error_overlap.tsv'
collapsed_duplication_overlap_df.to_csv(collapsed_duplication_overlap_path, sep='	', index=False)

collapsed_duplication_threshold3_df = collapsed_duplication_overlap_df.loc[
    (collapsed_duplication_overlap_df['window_set'] == 'grch38_gt_error_minus_chm13_gt_error')
    & (collapsed_duplication_overlap_df['error_rate_threshold'] == 3.0)
]
collapsed_duplication_threshold3_path = f'{tables_folder}/collapsed_duplication_high_error_overlap_threshold_3pct.tsv'
collapsed_duplication_threshold3_df.to_csv(collapsed_duplication_threshold3_path, sep='	', index=False)

print(f'Wrote: {collapsed_duplication_overlap_path}')
print(f'Wrote: {collapsed_duplication_threshold3_path}')
collapsed_duplication_threshold3_df

# %% [markdown]
# ---
#
# ## Exploratory analyses (preserved; review for deletion)
#
# Cells below are intentionally preserved from exploratory work and are not required for core inline manuscript statistics. Treat them as non-published unless promoted into `PAPER_CROSS_REFERENCE.md` with a manuscript/supplement/reviewer line reference.
#

# %%
# Quick dataframe inspection
binned_maf_data

# %%
# SGDP overlap extraction for 1kGP sample list
oneKGP_meta = pd.read_csv('resources/1000_genomes_meta.tsv', sep='\t')
overlaps = oneKGP_meta.loc[
    oneKGP_meta['Sample name'].isin(all_samples.sample_id.unique())
    & oneKGP_meta['Data collections'].str.contains('Simons Genome'),
    'Sample name'
]
overlaps.to_csv('resources/sample_subsets/SGDP_samples_in_1KGP.txt', index=False, header=False)
len(overlaps)

# %%
# Additional pedigree probe
proband_or_not_trio = pedigree.loc[
    ~(pedigree.sampleID.isin(pedigree.fatherID) | pedigree.sampleID.isin(pedigree.motherID))
]
proband_or_not_trio.head()

# %%
# HPRC/HGSVC overlap check
number_of_datasets_per_sample = all_samples.loc[
    (all_samples.genome == 'GRCh38')
    & (all_samples['method_of_phasing'] == 'phased_with_parents_and_pedigree')
    & (all_samples.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples')))
].groupby(['sample_id'], observed=True).size().sort_values()

print('These samples are in both the HPRC and HGSVC datasets. We will use the HPRC genomes.')
print(number_of_datasets_per_sample.loc[number_of_datasets_per_sample > 1].index.tolist())

# %%
# Method-level quick check
all_methods.loc[
    (all_methods.method_of_phasing == 'phased_with_parents_and_pedigree')
    & (all_methods.ground_truth_data_source == 'HPRC_samples')
]

# %%
# Ancestry-level exploratory plotting
calc_ancestries = all_ancestries.loc[
    (all_ancestries.method_of_phasing == 'phased_with_parents_and_pedigree')
    & (all_ancestries.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples')))
]

calc_ancestries = calc_ancestries.groupby(['superpopulation', 'genome'], observed=True).sum(numeric_only=True)
calc_ancestries['switch_error_rate'] = calc_ancestries.n_switch_errors / calc_ancestries.n_checked * 100
calc_ancestries['flip_error_rate'] = calc_ancestries.n_flips / calc_ancestries.n_checked * 100
calc_ancestries['true_switch_error_rate'] = calc_ancestries.n_true_switch_errors / calc_ancestries.n_checked * 100
calc_ancestries = calc_ancestries.reset_index()
calc_ancestries.loc[calc_ancestries.genome == 'GRCh38', 'switch_error_ratio'] = (
    calc_ancestries.loc[calc_ancestries.genome == 'GRCh38', 'switch_error_rate']
    / calc_ancestries.loc[calc_ancestries.genome == 'GRCh38', 'switch_error_rate'].min()
)
calc_ancestries.loc[calc_ancestries.genome == 'CHM13v2.0', 'switch_error_ratio'] = (
    calc_ancestries.loc[calc_ancestries.genome == 'CHM13v2.0', 'switch_error_rate']
    / calc_ancestries.loc[calc_ancestries.genome == 'CHM13v2.0', 'switch_error_rate'].min()
)
ax = sns.barplot(data=calc_ancestries.reset_index(), x='superpopulation', hue='genome', y='switch_error_ratio')
calc_ancestries.sort_values(['genome', 'superpopulation'])[['genome', 'superpopulation', 'switch_error_ratio']]

# %%
