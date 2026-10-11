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
from pandas.api.types import CategoricalDtype
import polars as pl
import seaborn as sns
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import importlib


pedigree_file_loc = 'resources/pedigrees/gatk_1kgp.ped'
summary_statistics_folder='intermediate_data'

if os.getcwd()!='phasing_T2T_project':
    os.chdir('..') # Change to the parent directory

# %%
import glob; glob.glob(f"{summary_statistics_folder}/*")

# %%
try:
    q = pd.read_parquet('intermediate_data/per_MAF_bin.parquet')
    q = q[(q.rounded_MAF=='singleton')&(q.method_of_phasing=='phased_with_parents_and_pedigree')&(q.syntenic=='All')&(q.type=='SNPs + Indels')]
    # q.type=='SNPs + Indels'
    q
except FileNotFoundError:
    q = None
    print("WARNING: per_MAF_bin.parquet not found (imputation data unavailable)")

# %%
MAF_performance_variant_df = pl.read_parquet(f'{summary_statistics_folder}/MAF_performance_variants.parquet')
all_chroms = pd.read_parquet(f'{summary_statistics_folder}/chroms.parquet')
all_samples = pd.read_parquet(f'{summary_statistics_folder}/samples.parquet')
all_ancestries = pd.read_parquet(f'{summary_statistics_folder}/ancestries.parquet')
all_methods = pd.read_parquet(f'{summary_statistics_folder}/methods.parquet')
all_variants = pl.read_parquet(f'{summary_statistics_folder}/variants.parquet')

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

# %% [markdown]
# ### Variant Filtering Stats 
# Stats in paragraph after "Full length, telomere to telomere maps of human genetic recombination rates are largely consistent with deCODE recombination rates"

# %%
variant_filter_stats = pd.read_parquet(f"{summary_statistics_folder}/filter_summary_stats.parquet")
variant_filter_stats['gatk_filter'] = variant_filter_stats.VQSLOD_filter #| variant_filter_stats.neg_train_site_filter
chm13_variant_stats = variant_filter_stats.loc[(variant_filter_stats.genome=='CHM13v2.0')]
singletons = chm13_variant_stats.loc[chm13_variant_stats.singleton]
not_singletons = chm13_variant_stats.loc[~chm13_variant_stats.singleton]

num_affected_gatk = chm13_variant_stats.loc[~(chm13_variant_stats.GRCh38_criteria_fail) &
                                            # chm13_variant_stats.CHM13_criteria_fail &
                                            chm13_variant_stats.gatk_filter].len.sum()


singletons_gained = singletons.loc[chm13_variant_stats.GRCh38_criteria_fail & 
                                   ~(chm13_variant_stats.CHM13_criteria_fail | 
                                   variant_filter_stats.CHM13_filtered)].len.sum()
print(f"Our VQSLOD filter resulted in the loss of an additional \n\
{num_affected_gatk} \
variants compared to the NYGC’s GRCh38 1KGP variant filtering methods, but the inclusion of singleton alleles added \n\
{singletons_gained} \
unique variants to our dataset (Supplementary Figure 1).")
singletons

# %%
not_singletons.loc[(not_singletons.Syntenic)].groupby(['CHM13_filtered', 'CHM13_criteria_fail','gatk_filter']).len.sum()

# %%
variant_filter_stats.groupby(['genome','gatk_filter','GRCh38_criteria_fail'])['len'].sum().reset_index()

# %%

# %%
variant_filter_stats.loc[(~variant_filter_stats.CHM13_criteria_fail)&(variant_filter_stats.GRCh38_criteria_fail)].groupby(['GRCh38_criteria_fail','gatk_filter','genome','singleton'])['len'].sum().reset_index()    

# %% [markdown]
# Table 1 and supplemental table XX

# %%
best_phased_variants = MAF_performance_variant_df.filter((pl.col('ground_truth_data_source')!='trios')&(pl.col('method_of_phasing')=='phased_with_parents_and_pedigree'))
catVarType = pl.Enum(['SNP','Indel','SNPs + Indels'])
best_phased_variants = best_phased_variants.with_columns(type=pl.col('type').cast(catVarType))
best_phased_variants = pl.concat([best_phased_variants, best_phased_variants.with_columns(type=pl.lit('SNPs + Indels').cast(catVarType))])
best_phased_variants_by_syntenic = best_phased_variants.group_by(['type','Syntenic','method_of_phasing','ground_truth_data_source','genome']
                                                                 ).sum(
                                                                 ).with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                                                switch_error_rate = pl.col('n_switch_errors')/pl.col('n_checked')*100
                                                                 ).to_pandas()
best_phased_variants_overall = best_phased_variants.group_by(['type','method_of_phasing','ground_truth_data_source','genome']
                                                                 ).sum(
                                                                 ).with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked')*100,
                                                                                switch_error_rate = pl.col('n_switch_errors')/pl.col('n_checked')*100
                                                                 ).to_pandas()


# %%
best_phased_variants_overall.sort_values(['ground_truth_data_source','genome','type'])

# %%
best_phased_variants_by_syntenic.loc[best_phased_variants_by_syntenic.ground_truth_data_source=='HPRC_samples'].sort_values(['ground_truth_data_source','Syntenic','genome','type'])

# %% [markdown]
# ### Percent concordance in different panels

# %%
overall_subsetted = MAF_performance_variant_df.filter((pl.col('ground_truth_data_source')!='trios')&(pl.col('method_of_phasing')=='phased_with_parents_and_pedigree')).group_by(['type','Syntenic','method_of_phasing','ground_truth_data_source','genome']).sum().with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked'))
overall           = MAF_performance_variant_df.filter((pl.col('ground_truth_data_source')!='trios')&(pl.col('method_of_phasing')=='phased_with_parents_and_pedigree')).group_by(['method_of_phasing','ground_truth_data_source','genome']).sum().with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked'))
overall_syntenic  = MAF_performance_variant_df.filter((pl.col('ground_truth_data_source')!='trios')&(pl.col('Syntenic'))&(pl.col('method_of_phasing')=='phased_with_parents_and_pedigree')).group_by(['method_of_phasing','ground_truth_data_source','genome']).sum().with_columns(gt_error_rate = pl.col('n_gt_errors')/pl.col('n_gt_checked'))

# %%
overall_syntenic = overall_syntenic.to_pandas()[['Syntenic','genome','ground_truth_data_source','gt_error_rate']].sort_values(['genome','Syntenic','gt_error_rate'])
overall_syntenic['Non-reference true positive rate'] = (1-overall_syntenic.gt_error_rate) * 100
overall_syntenic['gt_error_rate'] *= 100
overall_syntenic.loc[overall_syntenic.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples','HPRC_HGSVC_all_samples'))]

# %%
overall = overall.to_pandas()[['Syntenic','genome','ground_truth_data_source','gt_error_rate']].sort_values(['genome','Syntenic','gt_error_rate'])
overall['Non-reference true positive rate'] = (1-overall.gt_error_rate) * 100
overall['gt_error_rate'] *= 100
overall.loc[overall.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples','HPRC_HGSVC_all_samples'))]

# %%
overall_subsetted = overall_subsetted.to_pandas()[['type','Syntenic','genome','ground_truth_data_source','gt_error_rate']].sort_values(['genome','type','Syntenic','gt_error_rate'])
overall_subsetted['Non-reference true positive rate'] = (1-overall_subsetted.gt_error_rate) * 100
overall_subsetted['gt_error_rate'] *= 100
overall_subsetted.loc[overall_subsetted.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples','HPRC_HGSVC_all_samples'))]

# %%
# wide_genome_ideogram['SER_diff'] = wide_genome_ideogram['SER_t2t'] - wide_genome_ideogram['SER_grch38']
# region_table = (wide_genome_ideogram.loc[(wide_genome_ideogram.ground_truth_data_source=='HPRC_samples') & 
#                          (wide_genome_ideogram.method_of_phasing=='HPRC_variation_phased_with_reference_panel') #& 
#                         & (wide_genome_ideogram.contextual_MAC_grch38 > 1e6)
#                         & (wide_genome_ideogram.contextual_MAC_t2t > 1e6)
#                         ])
# decipher_regions = ['1p36.33','1p36.32','1p36.31','1p36.23','1p36.22','1p36.21','1q21.1','1q21.2','2p21','2p16.1','2p15','2q32.3','2q33.1','2q33.2','2q33.3','2q37.3','3q29','4p16.3','5p15.33','5p15.32','5p15.31','5p15.2','5q22.3','5q23.2','5q35.2','5q35.3','7q11.22','7q21.3','8p23.1','8q21.13','9q34.3','11p13','11p11.2','12q14.2','12q14.3','12q15','15q11.2','15q12','15q13.1','15q13.2','15q13.3','15q24.1','15q24.2','15q26.3','16p13.3','16p13.11','16p12.2','16p12.1','16p11.2','17p13.3','17p12','17p11.2','17q11.2','17q12','17q21.31','21q21.3','22p13','22p12','22p11.2','22p11.1','22q11.1','22q11.21','22q11.22','22q11.23','22q13.33','Xp22.33','Xp22.31','Xp11.23','Xp11.22','Xq22.2','Xq28']
# region_table.interval_names = region_table.chrom.str.replace('chr','') + region_table.interval_names
# region_table['decipher_region'] = region_table.interval_names.isin(decipher_regions)

# region_table = (region_table.sort_values('SER_diff').iloc[:].reset_index(drop=True)
#                         [['chrom','interval_names','decipher_region','n_switch_errors_grch38','n_checked_grch38','n_switch_errors_t2t','n_checked_t2t','SER_grch38','SER_t2t','SER_diff']])

# region_table = region_table.drop(columns='chrom').rename(columns={'interval_names':'Genetic Region',
#                                                    'n_switch_errors_grch38': 'Switch Errors (GRCh38 panel)',
#                                                    'n_switch_errors_t2t': 'Switch Errors (CHM13 panel)',
#                                                    'n_checked_grch38': 'Heterozygous sites (GRCh38 panel)',
#                                                    'n_checked_t2t': 'Heterozygous sites (CHM13 panel)',
#                                                    'SER_grch38': 'SER % (GRCh38 panel)',
#                                                    'SER_t2t': 'SER % (CHM13 panel)',
#                                                    'SER_diff': 'Absolute reduction in SER (% points)'})

# pd.set_option('display.float_format', lambda x: f'{x:,.2f}%')
# region_table.set_index('Genetic Region')

# %%
rephase=all_variants.filter((pl.col('method_of_phasing') == '1kgp_variation_phased_with_reference_panel') & 
                         (pl.col('ground_truth_data_source') == 'HPRC_samples'))

singletons_rephase_sum = rephase.filter(pl.col('contextual_MAC') == 1).group_by('genome').sum().to_pandas().reset_index()
singletons_rephase_sum['SER'] = singletons_rephase_sum.n_switch_errors/singletons_rephase_sum.n_checked*100
for idx, (genome, ser) in singletons_rephase_sum[['genome','SER']].iterrows():
    print(f'Singletons in {genome}: {ser}')

rare_rephase = rephase.filter(pl.col('contextual_MAC')/pl.col('contextual_AN') < 0.01)
common_rephase = rephase.filter(pl.col('contextual_MAC')/pl.col('contextual_AN') >= 0.01)

rare_rephase_sum = rare_rephase.group_by('genome').sum().to_pandas().reset_index()
common_rephase_sum = common_rephase.group_by('genome').sum().to_pandas().reset_index()
for idx, (genome, n_switch_errors, n_checked) in rare_rephase_sum[['genome','n_switch_errors', 'n_checked']].iterrows():
    print (f"\nRare in {genome}: {n_switch_errors/n_checked*100}")
    print (f"Common in {genome}: {common_rephase_sum.iloc[idx].n_switch_errors/common_rephase_sum.iloc[idx].n_checked*100}")
syntenic = common_rephase.filter(pl.col('Syntenic') == True)
nonsyntenic = common_rephase.filter(pl.col('Syntenic') == False)

syntenic_sum = syntenic.group_by('genome').sum().to_pandas().set_index('genome')
nonsyntenic_sum = nonsyntenic.group_by('genome').sum().to_pandas().set_index('genome')

for genome, row in syntenic_sum.iterrows():
    print (f"\nCommon Syntenic in {genome}: {row.n_switch_errors/row.n_checked*100}")
    print (f"Common Nonsyntenic {genome}: {nonsyntenic_sum.loc[genome].n_switch_errors/nonsyntenic_sum.loc[genome].n_checked*100}")

# %%
synEnum = pl.Enum(categories=['All_regions', 'Syntenic', 'Nonsyntenic'])
concordance = all_variants.filter((pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree') &
                                  (pl.col('ground_truth_data_source') == 'HPRC_samples')).with_columns(Syntenic=pl.col('Syntenic').cast(str).replace('true','Syntenic').replace('false','Nonsyntenic').cast(synEnum))

# %%
concordance_summary = concordance.clone().with_columns(Syntenic = pl.lit('All_regions').cast(synEnum))
concordance_summary = pl.concat((concordance, concordance_summary)).group_by(['Syntenic','genome','type']).sum().to_pandas().reset_index()
concordance_summary = pd.concat((concordance_summary, concordance_summary.copy().assign(type = 'SNPs + Indels').groupby(['Syntenic','genome','type'], observed=True).sum(numeric_only=True).reset_index()))
concordance_summary['gt_error_rate'] = concordance_summary.n_gt_errors/concordance_summary.n_gt_checked*100
concordance_summary['gt_accuracy_rate'] = 100-concordance_summary['gt_error_rate'] 
concordance_summary['switch_error_rate'] = concordance_summary.n_switch_errors/concordance_summary.n_checked*100
concordance_summary['Syntenic'] = concordance_summary.Syntenic.replace('false', 'Nonsyntenic').replace('true','Syntenic')

# %%
pd.set_option("display.precision", 3)
concordance_summary.sort_values(['Syntenic','genome','type']).set_index(['Syntenic','genome','type'])[['switch_error_rate','gt_error_rate','n_gt_errors','n_gt_checked', 'n_switch_errors', 'n_checked']].astype(float)

# %%

# %%
# Parental-uninformative rate among proband hets (CHM13v2.0 PAR1)
from pathlib import Path

import pandas as pd

parquet_path = Path("..") / "intermediate_data" / "proband_parent_uninformative_by_MAF.parquet"
uninformative_by_maf = pd.read_parquet(parquet_path)

# Quick sanity check: show how many variants and trios were used
uninformative_by_maf[["genome", "region", "n_trios", "n_variants_seen"]].drop_duplicates().reset_index(drop=True)

# %%
# Plot breakdown of informative vs uninformative categories by MAF bin (CHM13v2.0 PAR1)
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# If the previous cell hasn't been run, load the parquet here too.
if "uninformative_by_maf" not in globals():
    parquet_path = Path("..") / "intermediate_data" / "proband_parent_uninformative_by_MAF.parquet"
    uninformative_by_maf = pd.read_parquet(parquet_path)

# Canonical rounded_MAF ordering used throughout the repo
rounded_maf_order = [
    "singleton",
    "0.00021-0.00042",
    "0.00042-0.00064",
    "0.00064-0.001",
    "0.001-0.0016",
    "0.0016-0.0022",
    "0.0022-0.003",
    "0.003-0.004",
    "0.004-0.0054",
    "0.0054-0.0072",
    "0.0072-0.0094",
    "0.0094-0.0126",
    "0.0126-0.0172",
    "0.0172-0.0244",
    "0.0244-0.037",
    "0.037-0.06",
    "0.06-0.102",
    "0.102-0.166",
    "0.166-0.255",
    "0.255-0.372",
    "0.372-0.5",
]

plot_df = uninformative_by_maf.copy()
plot_df = plot_df.loc[(plot_df.genome == "CHM13v2.0") & (plot_df.region == "chrX:1-2394410")].copy()

# Ensure all required percent columns exist (compute from counts if needed).
needed_count_cols = [
    "n_proband_hets",
    "n_uninf_both_parents_het",
    "n_uninf_both_parents_missing",
    "n_uninf_one_missing_one_het",
    "n_informative",
]
missing_counts = [c for c in needed_count_cols if c not in plot_df.columns]
if missing_counts:
    raise ValueError(f"Parquet is missing required columns: {missing_counts}")

# Fill missing MAF bins with zeros for cleaner plotting
plot_df["rounded_MAF"] = pd.Categorical(plot_df["rounded_MAF"], categories=rounded_maf_order, ordered=True)
plot_df = plot_df.set_index("rounded_MAF").reindex(rounded_maf_order).reset_index()
for col in needed_count_cols:
    plot_df[col] = plot_df[col].fillna(0).astype(float)

den = plot_df["n_proband_hets"].replace(0, np.nan)
plot_df["pct_uninf_both_parents_het"] = plot_df["n_uninf_both_parents_het"] / den * 100
plot_df["pct_uninf_both_parents_missing"] = plot_df["n_uninf_both_parents_missing"] / den * 100
plot_df["pct_uninf_one_missing_one_het"] = plot_df["n_uninf_one_missing_one_het"] / den * 100
plot_df["pct_informative"] = plot_df["n_informative"] / den * 100

x = plot_df["rounded_MAF"].astype(str).tolist()
p_both_het = plot_df["pct_uninf_both_parents_het"].fillna(0).to_numpy()
p_both_missing = plot_df["pct_uninf_both_parents_missing"].fillna(0).to_numpy()
p_one_missing_one_het = plot_df["pct_uninf_one_missing_one_het"].fillna(0).to_numpy()
p_informative = plot_df["pct_informative"].fillna(0).to_numpy()

plt.figure(figsize=(12, 4))
bottom = np.zeros_like(p_both_het)
plt.bar(x, p_both_het, bottom=bottom, label="Uninformative: both parents het")
bottom = bottom + p_both_het
plt.bar(x, p_both_missing, bottom=bottom, label="Uninformative: both parents missing")
bottom = bottom + p_both_missing
plt.bar(x, p_one_missing_one_het, bottom=bottom, label="Uninformative: one missing + one het")
bottom = bottom + p_one_missing_one_het
# plt.bar(x, p_informative, bottom=bottom, label="Informative")

plt.xlabel("rounded_MAF")
plt.ylabel("% of proband het sites")
plt.title("CHM13v2.0 PAR1: parental informativeness breakdown by MAF bin")
plt.xticks(rotation=90)
plt.legend(frameon=False, ncol=2, fontsize=9)
plt.tight_layout()
plt.show()

# %%
plot_df

# %%
