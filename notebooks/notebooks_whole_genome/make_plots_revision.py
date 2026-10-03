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
# # Whole-genome figure and table generation
#
# This notebook generates manuscript figures, supplemental figures, reviewer-response plots, and imputation summary outputs from whole-genome summary tables.
#
# ## Paper source of truth
# The publication-figure-code mapping in `../../README.md` is the source of truth for organization and published figure numbering.
#
# ## Execution order
# 1. Run the whole-genome summary dataframe creation pipeline.
# 2. Run `notebooks/notebooks_whole_genome/calc_figures_for_paper.ipynb` if paper text values need to be refreshed.
# 3. Run this notebook for figures, tables, and reviewer-response visual checks.
#
# ## Data paths
# Set `RUN_PROFILE` in the first code cell. `test` uses plain folders (`intermediate_data`, `imputation_statistics`, `figures`, `tables`) for Docker smoke/test data. `whole_genome` uses `_whole_genome` folders for publication-scale inputs and outputs.
#
# ## Published vs experimental organization
# Production manuscript and supplemental outputs are organized first. Reviewer-response additions that are now in the manuscript/supplement are labeled with their current figure numbers. Analyses that are not referenced by the manuscript, supplement, or reviewer reply are placed under the final experimental section.
#

# %%
import os
from pathlib import Path
from copy import deepcopy
import itertools
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator, PercentFormatter
import seaborn as sns
import seaborn.objects as so
from statannotations.Annotator import Annotator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from utilities import *

# Choose "test" for Docker smoke/test data or "whole_genome" for publication-scale data.
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

required_intermediate_files = (
    "binned_maf_data.parquet",
    "chroms.parquet",
    "samples.parquet",
    "ancestries.parquet",
    "methods.parquet",
    "per_variant_category_imputation_performance.parquet",
    "per_sample_imputation_performance.parquet",
    "filter_summary_stats.parquet",
)

intermediate_data_folder = require_files(
    profile_paths["intermediate_data"],
    required_intermediate_files,
    "intermediate summary data",
    RUN_PROFILE
)

imputation_statistics_folder = profile_paths["imputation_statistics"]
if not Path(imputation_statistics_folder).exists():
    raise FileNotFoundError(
        f"RUN_PROFILE={RUN_PROFILE!r} expected imputation statistics folder at "
        f"{imputation_statistics_folder}"
    )

figures_folder = profile_paths["figures"]
figdir = figures_folder  # Backward-compatible alias used throughout this notebook.
supplemental_figures_folder = f"{figures_folder}/supplemental"
tables_folder = profile_paths["tables"]
source_figures_folder = "figures"

Path(figures_folder).mkdir(parents=True, exist_ok=True)
Path(supplemental_figures_folder).mkdir(parents=True, exist_ok=True)
Path(tables_folder).mkdir(parents=True, exist_ok=True)

print(f"RUN_PROFILE: {RUN_PROFILE}")
print(f"Working directory: {repo_root}")
print(f"Intermediate data folder: {intermediate_data_folder}")
print(f"Imputation statistics folder: {imputation_statistics_folder}")
print(f"Figures folder: {figures_folder}")
print(f"Tables folder: {tables_folder}")


# %% [markdown]
# ## 1) Display settings
# 1 col fig = 3.46" wide
# 2 col fig = 7.08" wide
#
# 300dpi
#
# font: Sans serif (Helvetica, Arial)
# font size: 5-7 pt
#
# Save as AI, EPS, or PDF.
#
# Max 50MB size per fig.
#
# All text needs to be in editable vector format
#
# Legend: brief figure title, short description for each panel
#

# %%
#### Display config settings

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 6

one_col=3.34
oneandahalf_col=4.49
two_col=7.08


## Padding (in inches) around axes; defaults to 3/72 inches, i.e. 3 points.
#figure.constrained_layout.h_pad:  0.04167
#figure.constrained_layout.w_pad:  0.04167

plt.rcParams['axes.linewidth'] = 0.75
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['grid.alpha'] = 0.5
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.edgecolor'] = 'white'

plt.rcParams['figure.figsize'] = (two_col,two_col)
plt.rcParams["figure.constrained_layout.use"] = True
plt.rcParams["errorbar.capsize"] = 0.3
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.dpi"] = 300

# journal wants editable text in the vector files (no Type 3 fonts), and nothing over 7 pt except the panel letters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.titlesize'] = 7 # default is 'large', which is 7.2 pt

# Fig 3 / Fig 4c,d style for the revision - 'bars_with_points' or 'points_with_ci'
chromosome_panel_style = 'bars_with_points'
# seaborn's bootstraps and the stripplot jitter were never seeded. Seed them so the source data match the figures
seed = 42 # the one seed for every random draw in these figures
chromosome_panel_seed = seed


perc_formatter = FuncFormatter(percent_formatter)
MAF_xticks = [0.01, 0.1, 1, 10, 50]
SER_xticks = [0.1, 1, 10, 50]

genome_palette=[(0.4,0.4,0.4)] + sns.color_palette()[4:5]
type_markers = {'SNPs + Indels':',', 'SNP':'o', 'Indel':'o'}
type_dashes = {'SNPs + Indels':'', 'SNP':(1,1), 'Indel':(3,0.5)}
syntenic_markers = {'Syntenic':'o', 'Nonsyntenic':'o'}
syntenic_dashes = {'Syntenic':'', 'Nonsyntenic':(1,1)}

genome_order = ['GRCh38','CHM13v2.0']
type_order = ['SNPs + Indels', 'SNP', 'Indel', 'SNPs', 'Indels']
syntenic_order = ['Syntenic', 'Nonsyntenic']
chrom_order = ['Whole\nGenome'] + ['chr'+str(i) for i in range(1,23)] + ['PAR1','PAR2','chrX']
singletons_order = ['No Filter','no_singletons']
region_order=['simple','in_STRs','in_platinum_not_STRs','multiallelic_not_STRs_not_platinum']
HPRC_plus_samples = ['HG002', 'HG005', 'HG00733', 'HG01109', 'HG01243', 'HG02055', 'HG02080', 'HG02109', 'HG02145', 'HG02723', 'HG02818', 'HG03098', 'HG03486', 'HG03492', 'NA18906', 'NA19240', 'NA20129', 'NA21309']

region_order= [
    "All",
    "in_STRs",
    "not_in_STRs",
    "in_platinum_STRs",
    "not_in_platinum_STRs",
    "in_segdups",
    "not_in_segdups",
    "in_STRs_or_platinum",
    "not_in_STRs_or_platinum",
    "in_STRs_or_platinum_or_overlap",
    "not_in_STRs_or_platinum_or_overlap",
    "multiallelic",
    "biallelic",
]

default_style = {
    'errorbar': None,
    'palette':genome_palette,
    'markers':type_markers,
    'dashes':type_dashes,
    'markersize':3,
    'markeredgewidth':0.2,
    'markeredgecolor': 'white',
    'linewidth':1
}

formatted_ancestry_names = {'all':'All\npopulations','WestEurasia':'West\nEurasia','SouthAsia':"South\nAsia",
                            'EastAsia':'East\nAsia','Africa':'Africa', 'America':'Americas','Oceania':'Oceania',
                            'CentralAsiaSiberia':'Central Asia\n& Siberia'}


syntenic_style = deepcopy(default_style)
syntenic_style['markers'] = syntenic_markers
syntenic_style['dashes'] = syntenic_dashes

highlight_style = {'size':'highlight', 'size_order':[True, False], 'sizes':(1.5,3)}

# %% [markdown]
# ## 2) Load and fix data
#

# %%
binned_maf_data = pd.read_parquet(f"{intermediate_data_folder}/binned_maf_data.parquet")
all_chroms = pd.read_parquet(f'{intermediate_data_folder}/chroms.parquet')
all_samples = pd.read_parquet(f'{intermediate_data_folder}/samples.parquet')
all_ancestries = pd.read_parquet(f'{intermediate_data_folder}/ancestries.parquet')
all_methods = pd.read_parquet(f'{intermediate_data_folder}/methods.parquet')

per_variant_category_imputation_performance = pd.read_parquet(f'{intermediate_data_folder}/per_variant_category_imputation_performance.parquet')
per_sample_imputation_performance = pd.read_parquet(f'{intermediate_data_folder}/per_sample_imputation_performance.parquet')

sum_stats_var_filters = pd.read_parquet(f"{intermediate_data_folder}/filter_summary_stats.parquet")

# %%
per_variant_category_imputation_performance = per_variant_category_imputation_performance.replace('T2T','CHM13v2.0')
separate_panel_and_subset = per_variant_category_imputation_performance.assign(info_cutoff = per_variant_category_imputation_performance.panel
                                                                    .str.extract(r'info_cutoff_(\d+\.\d+)',expand=False).fillna('0'),
                                                   panel       = per_variant_category_imputation_performance.panel
                                                                    .str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True,)).copy()
separate_panel_and_subset = separate_panel_and_subset.assign(**{'Variant Subset':'All panel variants', 'Reference Panel':'GRCh38'})
separate_panel_and_subset.loc[separate_panel_and_subset.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
separate_panel_and_subset['panel'] = separate_panel_and_subset['panel'].str.replace('.common_variants','').str.replace('native_panel','Native').str.replace('lifted_panel','Lifted')
separate_panel_and_subset.loc[((separate_panel_and_subset.genome=='GRCh38') & (separate_panel_and_subset.panel=='Lifted')) |
                              ((separate_panel_and_subset.genome=='CHM13v2.0') & (separate_panel_and_subset.panel=='Native')), 'Reference Panel'] = 'CHM13v2.0'
separate_panel_and_subset['panel_filter'] = separate_panel_and_subset['dataset'].str.replace('SGDP', '').str.replace('pangenome', '',).str.split('_', n=1).str[1].fillna('No Filter')
separate_panel_and_subset['dataset'] = separate_panel_and_subset['dataset'].str.split('_', n=1).str[0].replace('pangenome','Pangenome')
separate_panel_and_subset = separate_panel_and_subset.rename(columns={'var_category_id':'Variant Type'})
if separate_panel_and_subset.mean_AF.max()<1:
    separate_panel_and_subset.mean_AF *= 100

# %%
# right off the bat - per-sample imputation measurement is identical when run in the syntenic vs nonsyntenic mode, so let's only handle the overall imputation bin run here
per_sample_imputation_performance = per_sample_imputation_performance.replace('T2T','CHM13v2.0').loc[~per_sample_imputation_performance.file.str.contains('yntenic')]
separate_sample_panel_and_subset = per_sample_imputation_performance.assign(info_cutoff = per_sample_imputation_performance.panel
                                                                              .str.extract(r'info_cutoff_(\d+\.\d+)',expand=False).fillna('0'),
                                                                            panel       = per_sample_imputation_performance.panel
                                                                                                .str.replace(r'\.GWAS_filtered\.info_cutoff_0\.\d+$', '', regex=True,)).copy()

separate_sample_panel_and_subset = separate_sample_panel_and_subset.assign(**{'Variant Subset':'All panel variants', 'Reference Panel':'GRCh38'})
separate_sample_panel_and_subset.loc[separate_sample_panel_and_subset.panel.str.contains('common'), 'Variant Subset'] = 'Common variants'
separate_sample_panel_and_subset['panel'] = separate_sample_panel_and_subset['panel'].str.replace('.common_variants','').str.replace('native_panel','Native').str.replace('lifted_panel','Lifted')
separate_sample_panel_and_subset.loc[((separate_sample_panel_and_subset.genome=='GRCh38') & (separate_sample_panel_and_subset.panel=='Lifted')) |
                              ((separate_sample_panel_and_subset.genome=='CHM13v2.0') & (separate_sample_panel_and_subset.panel=='Native')), 'Reference Panel'] = 'CHM13v2.0'

separate_sample_panel_and_subset['panel_filter'] = separate_sample_panel_and_subset['dataset'].str.replace('SGDP', '').str.replace('pangenome', '',).str.split('_', n=1).str[1].fillna('No Filter')
separate_sample_panel_and_subset['dataset'] = separate_sample_panel_and_subset['dataset'].str.split('_', n=1).str[0].replace('pangenome','Pangenome')
separate_sample_panel_and_subset['Synteny'] = separate_sample_panel_and_subset['file'].str.contains('syntenic').map({True:'Syntenic', False:'Nonsyntenic'})
separate_sample_panel_and_subset = separate_sample_panel_and_subset.rename(columns={'var_category_id':'Variant Type'})

# %%
#Post hoc fixes - one off events that are tricky to deal with when compiling stats

# there is one variant (I have not investigated, but likely a denovo SNP) that is a singleton variant in the HPRC dataset.
# That variant passed all filters, but makes for a very odd figure (correctly phased: 0% error rate, incorrect: 100% error rate)
# We will trim it here
binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')&
                    (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')&
                    (binned_maf_data.rounded_MAF == 'singleton'),
                    ['n_switch_errors','n_checked','switch_error_rate']] = np.nan

binned_maf_data_regions = binned_maf_data.loc[~binned_maf_data.region.str.contains('teni')]
binned_maf_data = binned_maf_data.loc[binned_maf_data.region.str.contains('teni')|(binned_maf_data.region == 'All')].rename(columns={'region':'syntenic'})

# %%
# --- Additive STR/Platinum/overlap region bins (post-hoc; all-enum approach) ---
import polars as pl

# Keep enum definitions consistent with scripts/create_summary_phasing_dataframes_polars_regional.py
methods = [
    "phased_with_parents_and_pedigree",
    "phased_without_parents_or_pedigree",
    "1kgp_variation_phased_with_reference_panel",
]
ground_truths = [
    "trios",
    "HPRC_samples",
    "HGSVC_samples",
    "HPRC_HGSVC_all_samples",
    "HPRC_HGSVC_probands",
    "HGSVC_probands",
    "HGSVC_parents",
    "HGSVC_samples_nontrios_only",
]
catMethod = pl.Enum(methods)
catGroundTruth = pl.Enum(ground_truths)
catGenomes = pl.Enum(["GRCh38", "CHM13v2.0"])
catVarType = pl.Enum(["SNP", "Indel", "SNPs + Indels"])

# Ordered r2/MAF bin labels (rounded_MAF)
r2_bin_names = (
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
)
catRoundedMAF = pl.Enum(list(r2_bin_names))

region_order=['simple','in_STRs','in_platinum_not_STRs','multiallelic_not_STRs_not_platinum']
# Region enum is defined earlier in the notebook
catRegion = pl.Enum(region_order)

maf_variants = pl.read_parquet(f"{intermediate_data_folder}/MAF_performance_variants.parquet")
print('loaded maf_variants')

# Ensure all categorical columns are cast to the same enums before any group_by/concat
maf_variants = maf_variants.with_columns(
    [
        pl.col("method_of_phasing").cast(pl.String).cast(catMethod, strict=False),
        pl.col("ground_truth_data_source").cast(pl.String).cast(catGroundTruth, strict=False),
        pl.col("genome").cast(pl.String).cast(catGenomes, strict=False),
        pl.col("type").cast(pl.String).cast(catVarType, strict=False),
        pl.col("rounded_MAF").cast(pl.String).cast(catRoundedMAF, strict=False),
    ]
 )

for col in ["in_STRs", "in_platinum_STRs", "multiallelic"]:
    maf_variants = maf_variants.with_columns(pl.col(col).fill_null(False).cast(pl.Boolean))

maf_variants = maf_variants.with_columns(pl.when(~pl.col("in_STRs") & ~pl.col("in_platinum_STRs") & ~pl.col('multiallelic'))
                                           .then(pl.lit("simple"))
                                           .when(pl.col("in_STRs"))
                                               .then(pl.lit("in_STRs"))
                                           .when(~pl.col("in_STRs") & pl.col("in_platinum_STRs"))
                                               .then(pl.lit("in_platinum_not_STRs"))
                                           .when(~pl.col("in_STRs") & ~pl.col("in_platinum_STRs") & pl.col('multiallelic'))
                                               .then(pl.lit("multiallelic_not_STRs_not_platinum"))
                                           .alias('region')
                                           .cast(catRegion)
                                           )

ascending_regions = maf_variants.select(['genome','region', 'method_of_phasing', 'ground_truth_data_source', 'type', 'rounded_MAF',
                                            'n_switch_errors', 'n_checked', 'n_gt_errors', 'n_gt_checked','MAC', 'AN']
                                        ).group_by(['genome','region','method_of_phasing','ground_truth_data_source','type','rounded_MAF']).sum(

                                        ).with_columns(
                                            [
                                                (pl.col("n_switch_errors") / pl.col("n_checked") * 100).alias("switch_error_rate"),
                                                (pl.col("n_gt_errors") / pl.col("n_gt_checked") * 100).alias("gt_error_rate"),
                                                (pl.col("MAC") / pl.col("AN") * 100).alias("MAF"),
                                            ])
maf_total_vars = maf_variants.select(['genome','method_of_phasing','ground_truth_data_source','type','rounded_MAF',
                                            'n_switch_errors', 'n_checked', 'n_gt_errors', 'n_gt_checked','MAC', 'AN']
                                        ).group_by(['genome','method_of_phasing','ground_truth_data_source','type','rounded_MAF']
                                        ).sum(
                                        ).with_columns(
                                            pl.col('n_checked').alias('total_checked_in_region'),
                                            pl.col('n_gt_checked').alias('total_gt_checked_in_region'),
                                        ).select(['genome','method_of_phasing','ground_truth_data_source','type','rounded_MAF','total_checked_in_region','total_gt_checked_in_region'])

ascending_regions = ascending_regions.join(maf_total_vars,
                                              on=['genome','method_of_phasing','ground_truth_data_source','type','rounded_MAF'],
                                              how='left',
                                             ).with_columns(cumulative_ser=pl.col('n_switch_errors')/pl.col('total_checked_in_region')*100,
                                                            cumulative_ger=pl.col('n_gt_errors')/pl.col('total_gt_checked_in_region')*100)
ascending_regions = ascending_regions.sort(['genome','region','method_of_phasing','ground_truth_data_source','type','rounded_MAF']).to_pandas()
ascending_regions = ascending_regions.loc[ (ascending_regions.method_of_phasing=='phased_with_parents_and_pedigree')
                                          &(ascending_regions.ground_truth_data_source.isin(['HPRC_samples','HGSVC_probands_no_trios']))]
del maf_variants

# %%
# I never converted syntenic representation from true/false to something named. Let's do that.
binned_maf_data = binned_maf_data.assign(syntenic=binned_maf_data.syntenic
                                                                 .astype(str)
                                                                 .fillna('All')
                                                                 .replace(True,'Syntenic')
                                                                 .replace(False,'Nonsyntenic')
                                                                 )



# Nonsyntenic GRCh38 regions are mostly GRCh38 errors, and we are not interested in the performance in these regions
binned_maf_data = binned_maf_data.loc[~((binned_maf_data.genome=='GRCh38')&(binned_maf_data.syntenic=='Nonsyntenic'))]

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

# Publication figure generation treats intermediate inputs as immutable. The
# corrected in-memory table is consumed below; regenerating pipeline parquet is
# a separate analysis step.

# %%
# Pool weighted counts within MAF bins; headline SER uses the separate summary notebook.
cat_data = binned_maf_data.loc[binned_maf_data.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_probands','HGSVC_parents', 'HGSVC_samples_nontrios_only')) &
                               (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')].copy()
cat_data.loc[cat_data.n_switch_errors == 0, 'switch_error_rate'] = np.nan
cat_data['support'] = cat_data.ground_truth_data_source.map({'HPRC_samples':'proband', 'HGSVC_probands':'proband','HGSVC_parents':'parent', 'HGSVC_samples_nontrios_only':'nontrio'})

num_probands_in_dataset = 68
num_parents_in_dataset =  13
num_nontrios_in_dataset = 19
num_probands_in_1kgp = 608
num_parents_in_1kgp  = 1195
num_nontrios_in_1kgp = 1399

# Apply the GRCh38 proband correction to both population and benchmark counts.
cat_data.loc[(cat_data.support=='proband') & (cat_data.genome=='GRCh38'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = (((num_probands_in_1kgp-1)/(num_probands_in_dataset-1)) * cat_data.loc[(cat_data.support=='proband') & (cat_data.genome=='GRCh38'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
cat_data.loc[(cat_data.support=='proband') & (cat_data.genome=='CHM13v2.0'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_probands_in_1kgp/num_probands_in_dataset) * cat_data.loc[(cat_data.support=='proband') & (cat_data.genome=='CHM13v2.0'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
cat_data.loc[cat_data.support=='parent', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_parents_in_1kgp/num_parents_in_dataset) * cat_data.loc[cat_data.support=='parent', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
cat_data.loc[cat_data.support=='nontrio', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_nontrios_in_1kgp/num_nontrios_in_dataset) * cat_data.loc[cat_data.support=='nontrio', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
cat_data.loc[(cat_data.n_switch_errors == 0), 'switch_error_rate'] = np.nan
cat_data = cat_data.groupby(['genome','type','rounded_MAF','syntenic'], observed=True).sum(numeric_only=True)
cat_data['switch_error_rate'] = cat_data.n_switch_errors/cat_data.n_checked * 100
cat_data['gt_error_rate'] = cat_data.n_gt_errors/cat_data.n_gt_checked * 100
cat_data['MAF'] = cat_data.MAC/cat_data.AN * 100
cat_data = cat_data.reset_index()
cat_data['switch_error_rate'] = cat_data['switch_error_rate'].replace(0, np.nan)
binned_maf_data = pd.concat([binned_maf_data, cat_data.assign(ground_truth_data_source='estimated', method_of_phasing='phased_with_parents_and_pedigree')])


# %%
# Pool weighted counts within regional MAF bins using the same 3,202-member panel.
regional_cat_data = binned_maf_data_regions.loc[binned_maf_data_regions.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_probands','HGSVC_parents', 'HGSVC_samples_nontrios_only')) &
                               (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')].copy()
regional_cat_data.loc[regional_cat_data.n_switch_errors == 0, 'switch_error_rate'] = np.nan
regional_cat_data['support'] = regional_cat_data.ground_truth_data_source.map({'HPRC_samples':'proband', 'HGSVC_probands':'proband','HGSVC_parents':'parent', 'HGSVC_samples_nontrios_only':'nontrio'})

num_probands_in_dataset = 68
num_parents_in_dataset =  13
num_nontrios_in_dataset = 19
num_probands_in_1kgp = 608
num_parents_in_1kgp  = 1195
num_nontrios_in_1kgp = 1399

# Keep the same GRCh38 exclusion in the regional MAF-pooled counts.
regional_cat_data.loc[(regional_cat_data.support=='proband') & (regional_cat_data.genome=='GRCh38'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = (((num_probands_in_1kgp-1)/(num_probands_in_dataset-1)) * regional_cat_data.loc[(regional_cat_data.support=='proband') & (regional_cat_data.genome=='GRCh38'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
regional_cat_data.loc[(regional_cat_data.support=='proband') & (regional_cat_data.genome=='CHM13v2.0'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_probands_in_1kgp/num_probands_in_dataset) * regional_cat_data.loc[(regional_cat_data.support=='proband') & (regional_cat_data.genome=='CHM13v2.0'), ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
regional_cat_data.loc[regional_cat_data.support=='parent', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_parents_in_1kgp/num_parents_in_dataset) * regional_cat_data.loc[regional_cat_data.support=='parent', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
regional_cat_data.loc[regional_cat_data.support=='nontrio', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']] = ((num_nontrios_in_1kgp/num_nontrios_in_dataset) * regional_cat_data.loc[regional_cat_data.support=='nontrio', ['n_switch_errors','n_checked', 'n_gt_errors','n_gt_checked','MAC','AN']]).fillna(0).astype(int)
regional_cat_data.loc[(regional_cat_data.n_switch_errors == 0), 'switch_error_rate'] = np.nan
regional_cat_data = regional_cat_data.groupby(['genome','type','rounded_MAF', 'region'], observed=True).sum(numeric_only=True)
regional_cat_data['switch_error_rate'] = regional_cat_data.n_switch_errors/regional_cat_data.n_checked * 100
regional_cat_data['gt_error_rate'] = regional_cat_data.n_gt_errors/regional_cat_data.n_gt_checked * 100
regional_cat_data['MAF'] = regional_cat_data.MAC/regional_cat_data.AN * 100
regional_cat_data = regional_cat_data.reset_index()
regional_cat_data['switch_error_rate'] = regional_cat_data['switch_error_rate'].replace(0, np.nan)
binned_maf_data_regions = pd.concat([binned_maf_data_regions, regional_cat_data.assign(ground_truth_data_source='estimated', method_of_phasing='phased_with_parents_and_pedigree')])

# %% [markdown]
# ## 3) Main manuscript figures
#
# Published figure order follows the manuscript. Paper line references below are from the extracted text in `papers/` and are summarized in `PAPER_CROSS_REFERENCE.md`.
#

# %% [markdown]
# ### Supplemental Figure 1
# Panel filtering statistics

# %%
import importlib
import upsetplot
importlib.reload(upsetplot)
from upsetplot import UpSet
print (f'Using upsetplot version {upsetplot.__file__}')
plt.rcParams["figure.dpi"] = 96        # Display rendering
plt.rcParams["savefig.dpi"] = 300      # Export rendering
class UpSetWithScaledDots(UpSet):
    """UpSet variant that scales matrix dot diameter without changing layout."""
    def __init__(self, *args, dot_scale=0.6, **kwargs):
        self._dot_scale = dot_scale
        super().__init__(*args, **kwargs)

    def plot_matrix(self, ax, fig=None):
        if hasattr(self, '_calculated_element_size'):
            self._calculated_element_size *= self._dot_scale
        # upsetplot 0.9.0 does not accept the later optional ``fig`` keyword.
        return super().plot_matrix(ax)

def move_all_passing_to_front(df):
    passing=(False,)*len(df.index[0])
    new_index = df.groupby(df.index).sum().sort_values('num_variants', ascending=False).index.drop(passing)
    new_index = new_index.insert(0, passing)
    return df.loc[new_index,:]

def reshape_data_for_upsetplot(data, genome_name, filter_columns, color_scheme):
    """Process genome-specific data for UpSet plots"""
    genome_data = data.loc[data.genome==genome_name].rename(columns={'len':'num_variants','CHM13_criteria_fail':'other_criteria_fail'} if genome_name == 'GRCh38' else {'len':'num_variants','GRCh38_criteria_fail':'other_criteria_fail'})

    if genome_name == 'CHM13v2.0':
        genome_data['grch38_fail_besides_singleton'] = genome_data[data[[x for x in ['MERR_filter', 'HWE_pop_filter','AC_filter', 'f_missing_filter', 'var_len_filter', 'alt_star_filter', 'pass_filter'] if x != 'AC_filter'] + ['MAC_filter']].columns].any(axis=1)

    genome_data['other_criteria_fail'] = genome_data['other_criteria_fail'].astype(str)

    if genome_name == 'GRCh38':
        genome_data.loc[genome_data.other_criteria_fail=='False', 'other_criteria_fail'] = 'No'
        genome_data.loc[genome_data.other_criteria_fail=='True', 'other_criteria_fail'] = 'Yes'
    else:  # CHM13v2.0
        genome_data.loc[genome_data.singleton & (genome_data.grch38_fail_besides_singleton != 'True'), 'other_criteria_fail'] = 'Yes: singleton variant'
        genome_data.loc[genome_data.other_criteria_fail=='False', 'other_criteria_fail'] = 'No'
        genome_data.loc[genome_data.other_criteria_fail=='True', 'other_criteria_fail'] = 'Yes: Other reason'

    filters = genome_data.groupby(['other_criteria_fail']+filter_columns).sum(numeric_only=True).reset_index().set_index(filter_columns)
    syntenic_annotated = genome_data.groupby(['other_criteria_fail','Syntenic']+filter_columns).sum(numeric_only=True).reset_index().set_index(filter_columns)
    syntenic_only_fail = genome_data.loc[genome_data[filter_columns].sum(axis=1)>0].groupby(['other_criteria_fail','Syntenic']+filter_columns).sum(numeric_only=True).reset_index().set_index(filter_columns)

    return {
        'filters': filters,
        'syntenic_annotated': syntenic_annotated,
        'syntenic_only_fail': syntenic_only_fail,
        'color_scheme': color_scheme
    }

def create_upset_plot(data, is_syntenic, genome_name, cutoff_perc, element_size, top_x_categories, colors, filter_column_rename, dot_scale=1):
    """Create an UpSet plot for the given data"""
    subset_data = move_all_passing_to_front(data.loc[(data.Syntenic==is_syntenic)].rename_axis(index=filter_column_rename))
    cutoff = subset_data.num_variants.sum().item() * cutoff_perc

    upset = UpSetWithScaledDots(
        subset_data[['other_criteria_fail','num_variants']],
        show_percentages=True,
        sum_over='num_variants',
        sort_by='input',
        orientation="horizontal",
        min_subset_size=cutoff,
        intersection_plot_elements=0,
        element_size=element_size,
        dot_scale=dot_scale
    )
    print(subset_data.index.names)
    upset.add_stacked_bars(
        by="other_criteria_fail",
        title="Number of variants filtered",
        sum_over='num_variants',
        elements=len(subset_data.index.names),
        colors=colors
    )

    return upset

def create_upset_figure_layout(upset_plots_config, figsize=(13,9)):
    """Create the main figure with 2x2 subplot layout and labels"""
    fig = plt.figure(figsize=figsize, layout='constrained')
    ((subA, subB), (subC, subD)) = fig.subfigures(2, 2, wspace=-0.5, hspace=-0.5)

    # upsetplot 0.9.0 accepts a Figure-like object but assumes the legacy
    # Figure size methods are present. Matplotlib SubFigure intentionally omits
    # them. Report the local subfigure size and keep its parent-managed geometry
    # fixed when upsetplot requests a resize for its 30-point elements.
    def make_upsetplot_compatible(subfigure):
        subfigure.get_figwidth = lambda: subfigure.bbox.width / subfigure.dpi
        subfigure.get_figheight = lambda: subfigure.bbox.height / subfigure.dpi
        subfigure.set_figwidth = lambda width: None
        subfigure.set_figheight = lambda height: None
        return subfigure

    font_upscale = 1.33
    # Plot each UpSet plot in its designated subfigure
    subplots = [subA, subB, subC, subD]
    configs = upset_plots_config

    for i, (subplot, config) in enumerate(zip(subplots, configs)):
        subplot = make_upsetplot_compatible(subplot)
        axes = config['upset'].plot(fig=subplot)
        axes['extra0'].legend().set_title(config['legend_title'])
        add_letter_to_subfig(subplot, chr(ord('a') + i), fontsize=12*font_upscale)

    # Add main labels
    fig.text(.01,0.25, '1KGP CHM13v2.0 panel', fontsize=9*font_upscale, weight='bold', rotation='vertical', horizontalalignment='right', verticalalignment='center')
    fig.text(.01,0.75, '1KGP NYGC GRCh38 panel', fontsize=9*font_upscale, weight='bold', rotation='vertical', horizontalalignment='right', verticalalignment='center')
    fig.text(.25,1, 'Variants in shared genomic regions', fontsize=9*font_upscale, weight='bold', rotation='horizontal', horizontalalignment='center', verticalalignment='center')
    fig.text(.75,1, 'Variants in assembly-unique regions', fontsize=9*font_upscale, weight='bold', rotation='horizontal', horizontalalignment='center', verticalalignment='center')

    return fig

# Load and prepare data
grch38_filter_columns=['AC_filter',
                       'MERR_filter', 'HWE_pop_filter', 'f_missing_filter',
                       'var_len_filter', 'alt_star_filter', 'pass_filter']
chm13_filter_columns=['VQSLOD_filter', 'MAC_filter',
                      'MERR_filter', 'HWE_pop_filter', 'f_missing_filter',
                      'var_len_filter', 'alt_star_filter', 'pass_filter']

filter_column_rename = {'HWE_pop_filter':'Low population HWE',
                        'AC_filter': 'AC <= 1',
                        'MAC_filter': 'MAC == 0',
                        'f_missing_filter':">5% missing",
                        'var_len_filter':'Variant over 50bp',
                        'alt_star_filter':'ALT allele is star (\'*\')',
                        'pass_filter': 'PASS Filter',
                        'MERR_filter': 'Mend. Error Rate > 5%',
                        'VQSLOD_filter':'VQSLOD < 0'}

sum_stats_var_filters['only_filtered_from_grch38_because_singleton'] = sum_stats_var_filters[[x for x in grch38_filter_columns if x != 'AC_filter'] + ['MAC_filter']].any(axis=1)

# Define color schemes
grch38_color_scheme={'No: Singleton variant': sns.color_palette()[6],
                     'Yes':sns.color_palette()[1],
                     'No': sns.color_palette()[0]}
chm13_color_scheme={'Yes: singleton variant': sns.color_palette()[6],
                    'Yes: Other reason':sns.color_palette()[1],
                    'No': sns.color_palette()[0]}

# Process data for both genomes
grch38_processed = reshape_data_for_upsetplot(sum_stats_var_filters, 'GRCh38', grch38_filter_columns, grch38_color_scheme)
chm13_processed = reshape_data_for_upsetplot(sum_stats_var_filters, 'CHM13v2.0', chm13_filter_columns, chm13_color_scheme)

# Configuration
element_size = 30
dot_size=1
top_x_categories = 12
cutoff_perc = 0.005

# Create UpSet plots
# is_syntenic, genome_name, cutoff_perc, element_size, top_x_categories, colors, filter_column_rename
upset_syn_grch38    = create_upset_plot(grch38_processed['syntenic_annotated'], is_syntenic=True,
                                        genome_name='GRCh38', cutoff_perc=cutoff_perc,
                                        element_size=element_size, top_x_categories=16,
                                        colors=grch38_color_scheme, filter_column_rename=filter_column_rename)
upset_nonsyn_grch38 = create_upset_plot(grch38_processed['syntenic_annotated'], is_syntenic=False,
                                        genome_name='GRCh38', cutoff_perc=cutoff_perc,
                                        element_size=element_size, top_x_categories=top_x_categories,
                                        colors=grch38_color_scheme, filter_column_rename=filter_column_rename)
upset_syn_chm13     = create_upset_plot(chm13_processed['syntenic_annotated'], is_syntenic=True,
                                        genome_name='CHM13v2.0', cutoff_perc=cutoff_perc,
                                        element_size=element_size, top_x_categories=top_x_categories,
                                        colors=chm13_color_scheme, filter_column_rename=filter_column_rename)
upset_nonsyn_chm13  = create_upset_plot(chm13_processed['syntenic_annotated'], is_syntenic=False,
                                        genome_name='CHM13v2.0', cutoff_perc=cutoff_perc,
                                        element_size=element_size, top_x_categories=top_x_categories,
                                        colors=chm13_color_scheme, filter_column_rename=filter_column_rename)

# Configure plots for figure layout
upset_plots_config = [
    {'upset': upset_syn_grch38, 'legend_title': 'Would be filtered\nin CHM13v2.0 panel?'},
    {'upset': upset_nonsyn_grch38, 'legend_title': 'Would be filtered\nin CHM13v2.0 panel?'},
    {'upset': upset_syn_chm13, 'legend_title': 'Would be filtered\nin GRCh38 panel?'},
    {'upset': upset_nonsyn_chm13, 'legend_title': 'Would be filtered\nin GRCh38 panel?'}
]

# Create and save figure
fig = create_upset_figure_layout(upset_plots_config, figsize=(13,9))

fig.savefig(f'{supplemental_figures_folder}/Supplemental 1.png', bbox_inches='tight')
fig.savefig(f'{supplemental_figures_folder}/Supplemental 1.svg', bbox_inches='tight')
fig.savefig(f'{supplemental_figures_folder}/Supplemental 1.eps', bbox_inches='tight')
fig.savefig(f'{supplemental_figures_folder}/Supplemental 1.pdf', bbox_inches='tight')

plt.rcParams["figure.dpi"] = 300        # Display rendering
plt.rcParams["savefig.dpi"] = 300      # Export rendering

# %% [markdown]
# ## Figure 1: Recombination maps
#
# Paper refs: manuscript recombination-map section; reviewer reply lines 500-507. Supplemental context: Supplemental Figures 2-6.
#
# Status: Generated in abortvin repository. The recombination-map figures are produced by the recombination workflow/R assets.
#
# ![image.png](attachment:image.png)
#

# %% [markdown]
# ## Figure 2: panel phasing accuracy by MAF

# %%
# Paper cross-reference: Figure 2; manuscript lines 257-266; reviewer reply lines 397-417 and 544-555.
# Figure 2

fig, ((axA, axB),(axC, axD)) = plt.subplots(2,2, figsize=(two_col, 5))
fig.set_layout_engine(layout='constrained', h_pad=0.15, w_pad=0.15)
sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', #hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axA, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axB, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axC, **syntenic_style)


sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axD, **syntenic_style)


for i, ax in enumerate([axA, axB, axC, axD]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    # ax.grid() # editor wants the background gridlines gone
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    if i in (2,3):
        ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    else:
        ax.set_xlabel(None)
    if i in (0,2):
        ax.set_ylabel('Switch Error Rate (%)', labelpad=0)
    else:
        ax.set_ylabel(None)

# journal max is 7 pt for everything but the panel letters (these were 10 and 9)
axA.text(-.2,0.5, 'SNPs', fontsize=7, transform=axA.transAxes, weight='bold', rotation='vertical', horizontalalignment='right', verticalalignment='center')
axC.text(-.2,0.5, 'Indels', fontsize=7, transform=axC.transAxes, weight='bold', rotation='vertical', horizontalalignment='right', verticalalignment='center')
axA.set_title('Mendelian phased HPRC samples', fontsize=7, weight='bold', pad=10)
axB.set_title('Statistically phased HGSVC samples', fontsize=7, weight='bold', pad=10)


for i, ax in enumerate([axA, axB, axC, axD]): # legend axes
    if i == 0:
        loc='upper right'
    else:
        loc='upper right'
    handles, labels = ax.get_legend_handles_labels()
    for handle in handles:
        handle.set_marker('o')
    handles = handles[1:3] + handles[5:6]
    labels = labels[1:3] + labels[5:6]
    handles[-1].set_color(handles[1].get_color())
    legend = ax.legend(ncols=1, loc=loc, handles=handles, labels=labels)
    legend.set_title(None)
    legend.get_frame().set(
        alpha=0.75, boxstyle='square', edgecolor='grey', linewidth=0)
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)

x_pos = 0.5275
y_pos = 0.5
perc_of_figure=0.95
fig.add_artist(mpl.lines.Line2D([(1-perc_of_figure)/2, 1 - (1-perc_of_figure)/2], [y_pos, y_pos], color='grey', alpha=1, linewidth=1))
fig.add_artist(mpl.lines.Line2D([x_pos, x_pos], [(1-perc_of_figure)/2, 1 - (1-perc_of_figure)/2], color='grey', alpha=1, linewidth=1))

add_letter_to_ax(axA, 'a', fontsize=8) # journal wants 8 pt bold panel letters
add_letter_to_ax(axB, 'b', fontsize=8)
add_letter_to_ax(axC, 'c', fontsize=8)
add_letter_to_ax(axD, 'd', fontsize=8)

plt.savefig(f'{figdir}/Figure 2.png', facecolor='white')
plt.savefig(f'{figdir}/Figure 2.svg', facecolor='white')
plt.savefig(f'{figdir}/Figure 2.eps', facecolor='white') # journal takes .ai/.eps/.pdf, not .svg
plt.savefig(f'{figdir}/Figure 2.pdf', facecolor='white')

# %% [markdown]
# Fig 2. Phased panel switch error rates binned by minor allele frequency. Of the 3202 participants in the 1KGP, 39 had their genomes assembled by the HPRC. An additional 61 had their genomes assembled by the HGSVC. a) The switch error rate of haplotype panel SNPs from 39 indiviudals, using their HPRC-assembled genomes as ground truth. Variants located in regions of CHM13v2.0 that are non-syntenic with GRCh38 are included in overall variant bin error rates, but are also plotted separately as CHM13v2.0 (Novel) variants. b) SER of haplotype panel SNPs from 19 individuals, using their HGSVC-assembled genomes as ground truth. Unlike the 39 HPRC samples, these 19 individuals are not part of a 1KGP trio, and therefore were phased without any mendelian-based error correction. c) Indel SER from the mendelian-phased HPRC samples. d) Indel SER from the 19 individuals phased without mendelian error correction. Switch error rates and average minor allele frequency per bin are displayed on a log scale.

# %% [markdown]
# #### Supplemental Figure 9
# Illustration describing different methods of evaluating phasing quality.
# Phasing is the process of assigning heterozygous variant calls to one of two haplotypes (in the case of diploid organisms.) One of the challenges of assessing the accuracy of statistically phased haplotypes is identifying a ground truth to compare to. In the case of trio-based phasing evaluation, the ground truth is parental genomes. While not every heterozygous site can be assigned to a haplotype of paternal or maternal origin, most can confidently be assigned to one or the other. Alternately, one can use an empirically phased reference assembly as a source of ground truth.
# When a callset’s phasing describes a variant that is from a different parent of origin than the preceeding variant, that is called a switch error. Genotyping errors, as illustrated in the fourth panel, can force statistical phasing software to assign a variant to a haplotype that is different than nearby heterozygous sites. When two switch errors occur back-to-back, that is described as a flip error. Importantly, flip errors do not affect the phasing of surrounding heterozygous sites. Discerning flip errors from true switch errors is important when evaluating phasing methods, as the majority of switch errors can frequently be found in flip errors. Counting these switches as “true” switch errors inflates switch error rates by including sites that do not affect regional phasing.
#
# ![image.png](attachment:image.png)
#

# %% [markdown]
# #### Supplemental Figure 10
# 1000 Genomes Project samples by number of relatives available to provide information for Mendelian pre-phasing.

# %% [markdown]
#

# %% [markdown]
# ### Supplemental Figure 11
#
# SNP (L column) and Indel (R column) genotyping error rates for 1KGP phased haplotype panel variation, binned by minor allele frequency. Genotyping error rates were defined as the number of inaccurate heterozygous or homozygous alternative genotypes (both alleles must be called correctly; there was no partial credit for correctly calling one allele) divided by the number of heterozygous or homozygous alternative genotypes. Correct and incorrect variant calls were determined via comparison to ground-truth assemblies for the following subsets of samples: a) Mendelian-corrected trio probands with ground truth assemblies present in the HPRC human pangenome; b) Mendelian-corrected trio probands with ground truth assemblies produced by the HGSVC; c) Mendelian-corrected trio parents with ground truth assemblies produced by the HGSVC; d) Uncorrected samples which are not part of a 1KGP trio with ground truth assemblies produced by the HGSVC. e) The genotype error rates of trio probands, trio parents, and non-trio samples were weighted by each category's prevelance in the 1000 Genomes Project to produce estimated panel-wide average genotype error rates. Genotype error rates and average minor allele frequency per bin are displayed on a log scale.

# %%
## Figure 2 with parents and estimated overall
fig, ((axA, axB), (axC, axD), (axE, axF), (axG, axH), (axI, axJ)) = plt.subplots(5,2, figsize=(oneandahalf_col,8))

# fig.set_layout_engine(layout='constrained', h_pad=0.15, w_pad=0.15)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axA, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_probands') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axC, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_parents') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axE, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axG, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'estimated') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axI, **syntenic_style)


sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axB, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_probands')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axD, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_parents')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axF, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axH, **syntenic_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'estimated') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axJ, **syntenic_style)

for i, ax in enumerate([axA, axB, axC, axD, axE, axF, axG, axH, axI, axJ]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    if i in (8,9):
        ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    else:
        ax.set_xlabel(None)
    if i in (0,2,4,6,8):
        ax.set_ylabel('Genotyping Error Rate (%)', labelpad=0)
    else:
        ax.set_ylabel(None)
    add_letter_to_ax(ax, i)


x_offset = -0.4
axA.text(x_offset, 0.5, 'Mendelian pre-phased\nHPRC samples', fontsize=10, transform=axA.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axC.text(x_offset, 0.5, 'Mendelian pre-phased\nHGSVC probands', fontsize=10, transform=axC.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axE.text(x_offset, 0.5, 'Mendelian pre-phased\nHGSVC parents', fontsize=10, transform=axE.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axG.text(x_offset, 0.5, 'Statistically phased\nHGSVC samples', fontsize=10, transform=axG.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axI.text(x_offset, 0.5, 'Estimated\npanel-wide accuracy', fontsize=10, transform=axI.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')

axA.set_title('SNPs', fontsize=9, weight='bold', pad=10)
axB.set_title('Indels', fontsize=9, weight='bold', pad=10)

for i, ax in enumerate([axA, axB, axC, axD, axE, axF, axG, axH, axI, axJ]): # legend axes
    if i in (1,3,5,7,9):
        loc='lower right'
    else:
        loc='upper right'
    handles, labels = ax.get_legend_handles_labels()
    for handle in handles:
        handle.set_marker('o')
    handles = handles[1:3] + handles[5:6]
    labels = labels[1:3] + labels[5:6]
    handles[-1].set_color(handles[1].get_color())
    legend = ax.legend(ncols=1, loc=loc, handles=handles, labels=labels)
    legend.set_title(None)
    legend.get_frame().set(
        alpha=0, boxstyle='square', edgecolor='grey', linewidth=0)
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)


plt.savefig(f'{figdir}/supplemental/Supplemental 11.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 11.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 11.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 11.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Figure 12
# Concordance of GATK Haplotypecaller variant calls with genomic assembly-derived reference variants\
#
# A) Genotype variants from the T2T-CHM13 reference haplotype panel were compared to variants derived from a joint HPRC-HGSVC pangenome released by the HGSVC. Each point is one sample-assembly pair. Two samples had assemblies generated by both the HPRC and the HGSVC: HG00733 and HG02818. The lower line connects the HG00733 datapoints, and the upper line connects the HG02818 data points. Datapoints are colored by their 1KGP assigned superpopulation group.

# %%
# Relative variant concordance with HPRC and HGSVC samples
plot_data = (
    all_samples.loc[
        (all_samples.ground_truth_data_source.isin(['HPRC_samples','HGSVC_samples'])) &
        (all_samples.method_of_phasing=='phased_with_parents_and_pedigree') &
        (all_samples.genome=='CHM13v2.0')
    ]
    .reset_index(drop=False)              # keep original row id
    .rename(columns={"index":"row_id"})
)

order = ['HPRC_samples','HGSVC_samples']
hue_order = ['AFR','AMR','EAS','SAS','EUR']

fig, ax = plt.subplots(figsize=(one_col, one_col), layout="tight")
sns.swarmplot(
    data=plot_data,
    x='ground_truth_data_source', y='gt_error_rate',
    order=order,
    hue='superpopulation', hue_order=hue_order,
    ax=ax,
    s=3,
)
ax_x = ax.get_xlim()
ax_y = ax.get_ylim()
# Add jittered x positions back onto the data
affirm_has_x = 'x_data' in plot_data.columns
if affirm_has_x:
    plot_data = add_x_pos(plot_data, ax)
else:
    plot_data = add_x_pos(plot_data, ax)  # fallback: ensure x_data is added

# Draw lines connecting datapoints for specified samples using stored x_data
samples_to_connect = ['HG00733', 'HG02818']

sample_coords = {}  # {(sample_id, source): (x, y)}
for sample_id in samples_to_connect:
    for source in order:
        sample_row = plot_data[(plot_data['sample_id'] == sample_id) &
                               (plot_data['ground_truth_data_source'] == source)]
        if len(sample_row) == 0:
            continue
        if 'x_data' not in sample_row.columns:
            continue
        x_val = sample_row['x_data'].iloc[0]
        y_val = sample_row['gt_error_rate'].iloc[0]
        sample_coords[(sample_id, source)] = (x_val, y_val)
        print(f"Using {sample_id} ({source}): x={x_val:.4f}, y={y_val:.4f}")

# Draw lines directly between recorded x_data/y values
for sample_id in samples_to_connect:
    hprc_coord = sample_coords.get((sample_id, 'HPRC_samples'))
    hgsvc_coord = sample_coords.get((sample_id, 'HGSVC_samples'))
    if hprc_coord and hgsvc_coord:
        ax.plot([hprc_coord[0], hgsvc_coord[0]], [hprc_coord[1], hgsvc_coord[1]],
                'k-', linewidth=.5, alpha=1, zorder=3)


fig,[ax]=clean_figure(fig)
ax.set_xlim(ax_x)
ax.set_ylim(ax_y)
plt.savefig(f'{figdir}/supplemental/Supplemental_12.png', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_12.pdf', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_12.svg', facecolor='white', dpi=300)


# %% [markdown]
# ### Supplemental Figure 13
# Properties of called SNPs and Indels in increasing complex genetic regions

# %%
fig, ((axA, axB),(axC, axD),(axE, axF)) = plt.subplots(3, 2, figsize=(two_col, 6), layout='tight')
fig.get_layout_engine().set(rect=(0, 0.07, 1, 1)) # leave room at the bottom for the key - it was hanging off the page and getting clipped

for ax, var_type, y_var in zip( [axA, axB, axC, axD, axE, axF],
                                ['SNP','Indel','SNP','Indel','SNP','Indel'],
                                ['n_gt_checked','n_gt_checked','cumulative_ger','cumulative_ger','cumulative_ser','cumulative_ser']):

    data = ascending_regions.loc[(ascending_regions.type == var_type) & (ascending_regions.genome == 'CHM13v2.0')
                                 & (ascending_regions.ground_truth_data_source=='HPRC_samples')]
    (
        so.Plot(data, x="rounded_MAF", y=y_var, color="region")
            .add(so.Bar(), so.Stack())
            .on(ax)
            .plot()
    )
    for label in ax.get_xticklabels():
        label.set(rotation=45, ha="right", rotation_mode="anchor")
leg = fig.legends[0]
handles = leg.legend_handles
labels  = [update_legend_values.get(t.get_text(), t.get_text()) for t in leg.get_texts()]
for leg in fig.legends:
    leg.set_visible(False)
    leg.set_in_layout(False)

fig.legend(
    handles=handles,
    labels=labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.06), # inside the strip left for it above (was 0.02, off the page)
    ncol=len(handles),
    frameon=True,
    title="Variant Region",
    columnspacing=1.4,
    handlelength=1.8,
)
# 7 pt journal max for everything but the panel letters (these were 8)
axA.text(0.5,1.2, 'SNPs', fontsize=7, transform=axA.transAxes,  rotation='horizontal', horizontalalignment='center', verticalalignment='bottom')
axB.text(0.5,1.2, 'Indels', fontsize=7, transform=axB.transAxes,  rotation='horizontal', horizontalalignment='center', verticalalignment='bottom')

axA.text(-.25,0.5, 'Number of alt variant calls', fontsize=7, transform=axA.transAxes,  rotation='vertical', horizontalalignment='right', verticalalignment='center')
axC.text(-.25,0.5, 'Genotype Error Rate (%)', fontsize=7, transform=axC.transAxes,  rotation='vertical', horizontalalignment='right', verticalalignment='center')
axE.text(-.25,0.5, 'Switch Error Rate (%)', fontsize=7, transform=axE.transAxes,  rotation='vertical', horizontalalignment='right', verticalalignment='center') # row labels at -.25 (were -.2): the 7 pt tick labels pushed the axis labels into them
fig.align_ylabels([axA, axC, axE]) # the y labels sat at different distances from the axes (each is placed from its own tick labels, and '10%' is wider than '1%')
fig.align_ylabels([axB, axD, axF])


clean_figure(fig, letter_fontsize=8) # 8 pt panel letters for the journal (this is Extended Data now)
plt.savefig(f'{figdir}/supplemental/Supplemental_13.png', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_13.pdf', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_13.svg', facecolor='white', dpi=300)


# %% [markdown]
# ### Supplemental Figure 14
#
# Genotyping and switch error rates for isolated variants in GRCh38 and T2T-CHM13 panels.
#
# Switch error rates are lower in isolated variants. Variants were stratified by whether or not they overlapped another variant in both the Bykstrom-Bishop (2022) GRCh38 reference haplotype panel and the Lalli (2025) T2T-CHM13 panel. Variants were then binned into one of 20 minor allele frequency bins, and the switch error rate per bin was calculated. As in Figure 2, we display data from four sources: a) Trio-phased SNPs from the 39 1KGP samples that are present in the HPRCv1.1 draft reference pangenome, b) Statistically-phased SNPs from the 19 1KGP non-trio individuals present in the HPRC-HGSVC pangenome, c) Trio-phased indels from the same 39 HPRC samples, and d) Statistically-phased indels from the same 19 non-trio HGSVC samples.

# %%
region_order=['biallelic','multiallelic']

biallelic_style = {k:v for k,v in default_style.items()}
biallelic_style['markers']= {'SNP': 'o', 'Indel': 'o'}
biallelic_style['dashes']=  {'SNP': '', 'Indel': (1,1)}

fig, (axA, axB) = plt.subplots(1,2, figsize=(two_col, one_col), layout='tight')

sns.lineplot(x='MAF',y='gt_error_rate', style='type', hue='genome', hue_order = genome_order, style_order = ['SNP','Indel'],
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'estimated') #trios only
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.region == 'biallelic')],
             legend=True, ax=axA, **biallelic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='type', hue='genome', hue_order = genome_order, style_order = ['SNP','Indel'],
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'estimated') #trios only
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.region == 'biallelic')],
             legend=True, ax=axB, **biallelic_style)

axA.set_yscale('log')
axB.set_yscale('log')
axA.set_xscale('log')
axB.set_xscale('log')

for i, ax in enumerate([axA, axB]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    # if i in (2,3):
    ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    # else:
    #     ax.set_xlabel(None)
    if i in (0,2):
        ax.set_ylabel('Genotype Error Rate (%)', labelpad=0)
    else:
        ax.set_ylabel('Switch Error Rate (%)', labelpad=0)


fig.suptitle('Biallelic Variant Switch Rates', fontsize=10, weight='bold', y=1.02)
axA.set_title('Genotype Error Rate', fontsize=9, weight='bold', pad=10)
axB.set_title('Switch Error Rate', fontsize=9, weight='bold', pad=10)
clean_figure(fig)
plt.savefig(f'{supplemental_figures_folder}/Supplemental_14_isolated_biallelic_estimated_error.png', facecolor='white')
plt.savefig(f'{supplemental_figures_folder}/Supplemental_14_isolated_biallelic_estimated_error.svg', facecolor='white')
plt.savefig(f'{supplemental_figures_folder}/Supplemental_14_isolated_biallelic_estimated_error.eps', facecolor='white')
plt.savefig(f'{supplemental_figures_folder}/Supplemental_14_isolated_biallelic_estimated_error.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Figure 15
# SNP (L column) and Indel (R column) switch error rates for 1KGP phased haplotype panel variation, binned by minor allele frequency.
# Switch error rates were determined by comparing estimated phasing to ground-truth assemblies for the following subsets of samples: a) Mendelian-corrected trio probands with ground truth assemblies present in the HPRC human pangenome; b) Mendelian-corrected trio probands with ground truth assemblies produced by the HGSVC; c) Mendelian-corrected trio parents with ground truth assemblies produced by the HGSVC; d) Uncorrected samples which are not part of a 1KGP trio with ground truth assemblies produced by the HGSVC. e) The switch error rates of trio probands, trio parents, and non-trio samples were weighted by each category's prevelance in the 1000 Genomes Project to produce estimated panel-wide average switch error rates. Switch error rates and average minor allele frequency per bin are displayed on a log scale.

# %%
plt.rcParams['figure.dpi'] = 300
## Figure 2 with parents and estimated overall
fig, ((axA, axB), (axC, axD), (axE, axF), (axG, axH), (axI, axJ)) = plt.subplots(5,2, figsize=(oneandahalf_col,8))

# fig.set_layout_engine(layout='constrained', h_pad=0.15, w_pad=0.15)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axA, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_probands')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axC, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_parents')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axE, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axG, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'estimated')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'SNP')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axI, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axB, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_probands')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axD, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_parents')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axF, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axH, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'estimated')
                                      & (binned_maf_data.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data.type == 'Indel')
                                      & (binned_maf_data.syntenic != 'All')],
             legend=True, ax=axJ, **syntenic_style)

for i, ax in enumerate([axA, axB, axC, axD, axE, axF, axG, axH, axI, axJ]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    if i in (8,9):
        ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    else:
        ax.set_xlabel(None)
    if i in (0,2,4,6,8):
        ax.set_ylabel('Switch Error Rate (%)', labelpad=0)
    else:
        ax.set_ylabel(None)
    add_letter_to_ax(ax, i)

x_offset = -0.4
axA.text(x_offset, 0.5, 'Mendelian pre-phased\nHPRC samples', fontsize=10, transform=axA.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axC.text(x_offset, 0.5, 'Mendelian pre-phased\nHGSVC probands', fontsize=10, transform=axC.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axE.text(x_offset, 0.5, 'Mendelian pre-phased\nHGSVC parents', fontsize=10, transform=axE.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axG.text(x_offset, 0.5, 'Statistically phased\nHGSVC samples', fontsize=10, transform=axG.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')
axI.text(x_offset, 0.5, 'Estimated\npanel-wide accuracy', fontsize=10, transform=axI.transAxes, weight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center')

axA.set_title('SNPs', fontsize=9, weight='bold', pad=10)
axB.set_title('Indels', fontsize=9, weight='bold', pad=10)

for i, ax in enumerate([axA, axB, axC, axD, axE, axF, axG, axH, axI, axJ]): # legend axes
    if i in (5,7,9):
        loc='lower right'
    else:
        loc='upper right'
    handles, labels = ax.get_legend_handles_labels()
    for handle in handles:
        handle.set_marker('o')
    handles = handles[1:3] + handles[5:6]
    labels = labels[1:3] + labels[5:6]
    handles[-1].set_color(handles[1].get_color())
    legend = ax.legend(ncols=1, loc=loc, handles=handles, labels=labels)
    legend.set_title(None)
    legend.get_frame().set(
        alpha=0, boxstyle='square', edgecolor='grey', linewidth=0)
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)


plt.savefig(f'{figdir}/supplemental/Supplemental 15.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 15.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 15.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 15.pdf', facecolor='white')

# %% [markdown]
# ## Figure 3: per-chromosome SER and chrX/PAR performance

# %%
# Paper cross-reference: Figure 3; manuscript lines 282-290.
# Figure 3: Per-chromosome error rate performance in CHM13 vs GRCh38, including chrX

all_chroms_plus_chrX = pd.concat((all_chroms_chrX, all_samples.assign(chrom = 'chr1-22+X'))).reset_index(drop=True)

full_panel_HPRC_compared = all_chroms_plus_chrX.loc[ (all_chroms_plus_chrX.method_of_phasing=='phased_with_parents_and_pedigree') &
                                                (all_chroms_plus_chrX.ground_truth_data_source.isin(('HPRC_samples','HGSVC_probands', 'HGSVC_parents', 'HGSVC_samples_nontrios_only')))].replace('chr1-22+X', 'Whole\nGenome')

fig, ax = plt.subplots(1,1,figsize=(two_col,1.5), layout='tight')

# Editor wants every sample shown as a point on the bar charts. Points go down first, then the bars + error bars go on top (zorder)
# One stripplot per genome so each genome gets its own point edge colour
for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75) # barplot desaturates its palette by .75, match the bars
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')&(full_panel_HPRC_compared.genome==genome)])

if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                # split=True, inner=None, linewidth=0.5, cut=0, density_norm='area',
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])
    # bars 75% opaque, with a solid outline the same colour as the bar
    for bar in ax.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    # same mean and 95% CI as the barplot, just drawn as lines instead of bars
    sns.pointplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])

ax.set_ylabel('Switch Error Rate')
ax.set_xlabel('Chromosome')
ax.set_ylim(0,full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')].switch_error_rate.max()*1.08)
if chromosome_panel_style == 'bars_with_points':
    legend = ax.legend()
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)

    sns.move_legend(ax, "upper left")
else:
    legend = add_points_with_ci_legend(ax, genome_palette, genome_order)
ax.yaxis.set_major_formatter(PercentFormatter(decimals=2))
# ax.grid(axis='y', alpha=0.5) # editor wants the background gridlines gone
ax.set_axisbelow(True)
if chromosome_panel_style == 'bars_with_points':
    legend.set_title('Genome')
plt.savefig(f'{figdir}/Figure 3.png', facecolor='white')
plt.savefig(f'{figdir}/Figure 3.svg', facecolor='white')
plt.savefig(f'{figdir}/Figure 3.eps', facecolor='white')
plt.savefig(f'{figdir}/Figure 3.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Figure 16
# Whole genome, pseudoautosomal region 1 (PAR1), PAR2, and non-PAR chromosome X switch error rates.
# Switch error rates for trio probands with HPRC or HGSVC assemblies from either the GRCh38 or T2T-CHM13 1KGP reference haplotype panels, with each sample's assembly used as ground truth. PAR switch error rates are calculated separately from the non-PAR region of chromosome X. Non-PAR chromosome X switch error rates are only calculated from female samples.

# %%
# Paper cross-reference: Supplemental Figure 14; supplemental lines 191-195; manuscript lines 282-290.
# fig, ax = plt.subplots(1,1,figsize=(one_col,one_col), layout='tight')
fig, ax = plt.subplots(1,1,figsize=(one_col,one_col), layout='tight')
## Supplemental: PAR1, PAR2, chrX
all_chroms_plus = pd.concat((all_chroms, all_samples.assign(chrom = 'chr1-22+X'))).reset_index(drop=True)
all_chroms_plus['chrom'] = all_chroms_plus.chrom.map(lambda x: {'chrX':'Non-PAR ChrX'}.get(x,x))

### PAR is full of outlier samples (think HG01175, with 1 switch error and 2 checked)
all_chroms_plus = all_chroms_plus.loc[all_chroms_plus.n_checked>5]

full_panel_HPRC_compared_sep_chr_X = all_chroms_plus.loc[ (all_chroms_plus.method_of_phasing=='phased_with_parents_and_pedigree') &
                                                (all_chroms_plus.ground_truth_data_source.isin(('HPRC_samples','HGSVC_probands')))].replace('chr1-22+X', 'Whole\nGenome').rename(columns={'genome':"Genome"})

ax = sns.violinplot(x='chrom', y='switch_error_rate', hue='Genome', ax=ax,
            hue_order=genome_order, palette=genome_palette, order=['Whole\nGenome','PAR1','PAR2','Non-PAR ChrX'],
            # err_kws={'linewidth':1}, capsize=0.3,
            split=True, inner=None, linewidth=0.5, cut=0, density_norm='count',
            data=full_panel_HPRC_compared_sep_chr_X.loc[full_panel_HPRC_compared_sep_chr_X.chrom.isin(('PAR1','PAR2','Non-PAR ChrX', 'Whole\nGenome'))],
            legend=True)
# ax.set_ylim(0,2)
ax.set_ylabel('Switch Error Rate')
ax.set_xlabel('Chromosome')
ax.yaxis.set_major_formatter(perc_formatter)
ax.grid(axis='y', alpha=0.5)
ax.set_axisbelow(True)
clean_figure(fig)
plt.savefig(f'{figdir}/supplemental/Supplemental 16.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 16.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 16.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 16.pdf', facecolor='white')

# %% [markdown]
# ## Figure 4: switch errors, flip errors, and true switch errors
#
# Flip error rates and true switch error rates.
#
# a) Throughout this study, we have compared our panel's phasing with ground truth reference assemblies to identify switch errors, defined as sites at which the maternal and paternal haplotypes switch strands. However, regions which are difficult to sequence can often have genotyping errors that result in flip errors, defined as two consecutive switch errors. Unlike true switch errors, flip errors do not impact the phasing accuracy of nearby variants. Eight switch errors are in the illustrated haplotypes, but only two of these are true switch errors - the other six switch errors are part of three flip errors. b) The switch error rate, flip error rate, and true switch error rate of all 100 samples with assembly-based ground truths. Data from the 19 HGSVC individuals that were phased without mendelian error correction are highlighted. Data points are plotted on top of a violin plot showing the overall distribution of error rates. c) Whole genome and per-chromosome true switch error rates from 39 HPRC-assembled samples. d) Whole genome and per-chromosome true switch error rates from 19 HGSVC-assembled samples phased without mendelian error correction.
#

# %%
### Figure 4a illustrator image

from IPython.display import Image
fig4a_image = Path(figdir) / 'Fig 4a.png'
if not fig4a_image.exists():
    fig4a_image = Path(source_figures_folder) / 'Fig 4a.png'
Image(filename=str(fig4a_image), width=500)


# %%
# Paper cross-reference: Figure 4; manuscript lines 293-321.
## Figure 4
figure_layout="""
AAABB
AAABB
AAABB
CCCCC
DDDDD
"""
fig, axes = plt.subplot_mosaic(figure_layout, figsize=(two_col, 5))#layout='tight')
axA = axes['A']
axB = axes['B']
axC = axes['C']
axD = axes['D']


## Fig 4b

error_order = ['Total SER', 'Flip rate','True SER']

compare_stats = all_samples.loc[all_samples.ground_truth_data_source.isin(('HPRC_samples','HGSVC_probands','HGSVC_parents','HGSVC_samples_nontrios_only')) &
                                (all_samples.method_of_phasing=='phased_with_parents_and_pedigree')
                               ].rename(columns={'true_switch_error_rate':'True SER', 'flip_error_rate':'Flip rate', 'flip_and_switches_SER':'Total SER', 'genome':'Genome'}
                               ).melt(value_vars=['True SER','Flip rate','Total SER'],
                                       id_vars=['sample_id', 'Genome', 'population', 'superpopulation', 'sex','trio_phased', 'method_of_phasing', 'ground_truth_data_source'],
                                      var_name='Measure of Phasing Error',
                                      value_name='Error Rate').replace('CHM13v2.0', 'T2T-CHM13')

dot_size=3
np.random.seed(chromosome_panel_seed) # stripplot jitter is random, seed it so the figure is reproducible

sns.violinplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, hue_order=['GRCh38','T2T-CHM13'], palette=genome_palette, linewidth=0.5,
               dodge=True, split=True, cut=0.5, inner=None, density_norm='count', fill=True, alpha=0.2, bw_adjust=.5,
               data = compare_stats, ax=axB, legend=True)
sns.violinplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, hue_order=['GRCh38','T2T-CHM13'], palette=genome_palette, linewidth=0.5,
               dodge=True, split=True, cut=0.5, inner=None, density_norm='count', fill=False, bw_adjust=.5,
               data = compare_stats, ax=axB, legend=False)
sns.stripplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, jitter=0.25, s=dot_size,
               dodge=True, alpha=0.6, linewidth=0.1, palette={'GRCh38':'dimgrey', 'T2T-CHM13':'rebeccapurple'},
               data = compare_stats.loc[compare_stats.trio_phased==True], ax=axB, legend=True)
sns.stripplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, jitter=0.25, s=dot_size,
               dodge=True, alpha=0.6, linewidth=0.5, palette={'GRCh38':'slategrey', 'T2T-CHM13':'plum'},
               data = compare_stats.loc[compare_stats.ground_truth_data_source=='HGSVC_samples_nontrios_only'], ax=axB, legend=True, edgecolor='black')
axB.yaxis.set_major_formatter(perc_formatter)


#Fig 4c
# points first, then bars + error bars on top (same as Fig 3)
for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75)
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')&(full_panel_HPRC_compared.genome==genome)])
if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])
    for bar in axC.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    sns.pointplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])

#Fig 4d
for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75)
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')&(full_panel_HPRC_compared.genome==genome)])
if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')])
    for bar in axD.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    sns.pointplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')])

axC.yaxis.set_major_locator(MultipleLocator(0.15))
axD.yaxis.set_major_locator(MultipleLocator(0.15))

axC.set_xlabel(None)
axD.set_xlabel('Chromosome')
axC.set_title('HPRC Mendelian Phased', pad=6)
axD.set_title('HGSVC Statistically Phased', pad=6)

handles, labels =axB.get_legend_handles_labels()
new_handles=list()
for i, (handle, label) in enumerate(zip(handles, labels)):
    handle.set_alpha(.8)
    if i == 0:
        handle.set_facecolor('#E0E0E0')
        handle.set_edgecolor('#666666')
        handle.set_alpha(1)
    elif i == 1:
        handle.set_facecolor('#E9E3F0')
        handle.set_edgecolor('#9467BD')
        handle.set_alpha(1)
    new_handles.append(handle)

labels = labels[:2] + [l + ' (Mendelian phased)' for l in labels[2:4]] + [l+ ' (Statistically phased)' for l in labels[4:6]]
axB.legend(new_handles,labels)

for ax in [axC, axD]:
    ax.set_ylabel('True Switch Error Rate')
    if chromosome_panel_style == 'bars_with_points':
        legend = ax.legend()
        for text in legend.get_texts():
            t = text.get_text()
            new_t = update_legend_values.get(t,t)
            text.set_text(new_t)

        sns.move_legend(ax, "upper left")
    else:
        legend = add_points_with_ci_legend(ax, genome_palette, genome_order)
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=2))
    # ax.grid(axis='y', alpha=0.5) # editor wants the background gridlines gone
    ax.set_axisbelow(True)
    if chromosome_panel_style == 'bars_with_points':
        legend.set_title('Genome')
axD.get_legend().remove() # with the points shown the 4d legend covers a couple of them, and it's the same as the 4c legend anyway

for i, ax in enumerate(axes.values()):
    add_letter_to_ax(ax, i, fontsize=8) # journal wants 8 pt bold panel letters


## Fig 4a - graphic
sns.despine(ax=axA, top=True, right=True, left=True, bottom=True)
axA.set_yticks([])
axA.set_xticks([])
# clean_figure(fig)
plt.savefig(f'{figdir}/Figure 4.png', facecolor='white')
plt.savefig(f'{figdir}/Figure 4.svg', facecolor='white')
plt.savefig(f'{figdir}/Figure 4.eps', facecolor='white')
plt.savefig(f'{figdir}/Figure 4.pdf', facecolor='white')

# %% [markdown]
# ### Figure 4a: drop the Illustrator graphic into Figure 4
#
# Puts the editable Illustrator version of the 4a graphic (`figures/Fig 4a.svg`) into the empty axis A, centred left-right in the space it took up in the submitted figure and centred top-to-bottom on the 4b box. Writes 'Figure 4 with graphic' as svg + pdf (text stays editable), plus a png rendered from the pdf.

# %%
import shutil
import subprocess

for ext in ['svg','pdf']:
    shutil.copy(f'{figdir}/Figure 4.{ext}', f'{figdir}/Figure 4 with graphic.{ext}')
place_panel_a(f'{figdir}/Figure 4 with graphic.svg', f'{figdir}/Figure 4 with graphic.pdf', f'{source_figures_folder}/Fig 4a.svg', 'illustrator',
              axis_top_left_points(fig, axA), centre_y=axis_vertical_centre_points(fig, axB))
subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', f'{figdir}/Figure 4 with graphic.pdf', f'{figdir}/Figure 4 with graphic'])

# %% [markdown]
# ### Extended Data: chromosome detail (Figure 3 + Figure 4c,d)
#
# Whole genome and per-chromosome switch and true switch error rates. Figure 3 and Figure 4c,d moved together into one Extended Data figure (main text is one display item over the limit), one row of chromosomes per panel.
#
# a) Switch error rates of the 39 HPRC trio probands (Mendelian phased), was Figure 3. b) True switch error rates of the same probands, was Figure 4c. c) True switch error rates of the 19 HGSVC samples that are not part of a trio (statistically phased), was Figure 4d.

# %%
# Extended Data, chromosome detail - old Figure 3 + Figure 4c,d. Same code as those two cells, just stacked into one figure
fig, axes = plt.subplot_mosaic("A\nB\nC", figsize=(two_col, 3*1.55))
ax = axes['A']
axC = axes['B']
axD = axes['C']

## a (old Figure 3)
for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75)
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')&(full_panel_HPRC_compared.genome==genome)])

if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])
    for bar in ax.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    sns.pointplot(x='chrom', y='switch_error_rate', hue='genome', ax=ax,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])

ax.set_ylabel('Switch Error Rate')
ax.set_xlabel('Chromosome')
ax.set_ylim(0,full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')].switch_error_rate.max()*1.08)
if chromosome_panel_style == 'bars_with_points':
    legend = ax.legend()
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)

    sns.move_legend(ax, "upper left")
else:
    legend = add_points_with_ci_legend(ax, genome_palette, genome_order)
ax.yaxis.set_major_formatter(PercentFormatter(decimals=2))
ax.set_axisbelow(True)
if chromosome_panel_style == 'bars_with_points':
    legend.set_title('Genome')
ax.set_xlabel(None) # only the bottom panel gets the Chromosome label
ax.set_title('HPRC Mendelian Phased', pad=6) # name the cohort like b and c do

## b, c (old Figure 4c,d)
for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75)
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')&(full_panel_HPRC_compared.genome==genome)])
if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])
    for bar in axC.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    sns.pointplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axC,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HPRC_samples')])

for genome, color in zip(genome_order, genome_palette):
    np.random.seed(chromosome_panel_seed)
    if chromosome_panel_style == 'bars_with_points':
        point_palette={g:'white' for g in genome_order}
        point_edgecolor=sns.desaturate(color, .75)
        point_linewidth=0.4
    else:
        point_palette={g:sns.desaturate(color, .75) for g in genome_order}
        point_edgecolor='white'
        point_linewidth=0.15
    sns.stripplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=point_palette, edgecolor=point_edgecolor, linewidth=point_linewidth,
                dodge=True, jitter=0.25, size=1.3, legend=False, zorder=2,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')&(full_panel_HPRC_compared.genome==genome)])
if chromosome_panel_style == 'bars_with_points':
    sns.barplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette=genome_palette,
                err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.3, seed=chromosome_panel_seed, zorder=1,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')])
    for bar in axD.patches:
        rgb = mpl.colors.to_rgb(bar.get_facecolor())
        bar.set_facecolor((*rgb, 0.75))
        bar.set_edgecolor((*rgb, 1))
        bar.set_linewidth(0.5)
else:
    sns.pointplot(x='chrom', y='true_switch_error_rate', hue='genome', ax=axD,
                order=[x for x in chrom_order if 'PAR' not in x],
                hue_order=genome_order, palette={g:'black' for g in genome_order}, dodge=0.4, linestyle='none', markers='_',
                markersize=6, markeredgewidth=0.8, err_kws={'linewidth':0.5, 'zorder':3}, capsize=0.15, seed=chromosome_panel_seed,
                legend=False, zorder=4,
                data=full_panel_HPRC_compared.loc[(full_panel_HPRC_compared.ground_truth_data_source=='HGSVC_samples_nontrios_only')])

axC.yaxis.set_major_locator(MultipleLocator(0.15))
axD.yaxis.set_major_locator(MultipleLocator(0.15))

axC.set_xlabel(None)
axD.set_xlabel('Chromosome')
axC.set_title('HPRC Mendelian Phased', pad=6)
axD.set_title('HGSVC Statistically Phased', pad=6)

for ax in [axC, axD]:
    ax.set_ylabel('True Switch Error Rate')
    if chromosome_panel_style == 'bars_with_points':
        legend = ax.legend()
        for text in legend.get_texts():
            t = text.get_text()
            new_t = update_legend_values.get(t,t)
            text.set_text(new_t)

        sns.move_legend(ax, "upper left")
    else:
        legend = add_points_with_ci_legend(ax, genome_palette, genome_order)
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=2))
    ax.set_axisbelow(True)
    if chromosome_panel_style == 'bars_with_points':
        legend.set_title('Genome')
# just one legend for the whole figure (panel a)
axC.get_legend().remove()
axD.get_legend().remove()

for i, ax in enumerate(axes.values()):
    add_letter_to_ax(ax, i, fontsize=8)

# named for what it shows, not a number - the Extended Data numbers keep moving as figures get added
plt.savefig(f'{figdir}/Extended Data - chromosome detail.png', facecolor='white')
plt.savefig(f'{figdir}/Extended Data - chromosome detail.svg', facecolor='white')
plt.savefig(f'{figdir}/Extended Data - chromosome detail.eps', facecolor='white')
plt.savefig(f'{figdir}/Extended Data - chromosome detail.pdf', facecolor='white')

# %% [markdown]
# ### Figure 4 with only panels a and b
#
# If Figure 4c,d move to the chromosome-detail Extended Data figure, Figure 4 is just the 4a graphic and the 4b violins. Same 4b code as the Figure 4 cell. Writes 'Figure 4 (a,b only)' and 'Figure 4 (a,b only) with graphic'.

# %%
# keep a and b exactly where they were in the submitted Figure 4. Without c and d under them constrained layout makes b taller and shoves a over, so place the axes by hand (positions in points, measured off the submitted figure)
fig = plt.figure(figsize=(two_col, 2.5), layout='none')
axA = fig.add_axes([36.4993/(two_col*72), 1-153.731/(2.5*72), (296.855-36.4993)/(two_col*72), (153.731-19.1315)/(2.5*72)])
axB = fig.add_axes([348.356/(two_col*72), 1-153.731/(2.5*72), (506.760-348.356)/(two_col*72), (153.731-19.1315)/(2.5*72)])
axes = {'A':axA, 'B':axB}

## Fig 4b

error_order = ['Total SER', 'Flip rate','True SER']

compare_stats = all_samples.loc[all_samples.ground_truth_data_source.isin(('HPRC_samples','HGSVC_probands','HGSVC_parents','HGSVC_samples_nontrios_only')) &
                                (all_samples.method_of_phasing=='phased_with_parents_and_pedigree')
                               ].rename(columns={'true_switch_error_rate':'True SER', 'flip_error_rate':'Flip rate', 'flip_and_switches_SER':'Total SER', 'genome':'Genome'}
                               ).melt(value_vars=['True SER','Flip rate','Total SER'],
                                       id_vars=['sample_id', 'Genome', 'population', 'superpopulation', 'sex','trio_phased', 'method_of_phasing', 'ground_truth_data_source'],
                                      var_name='Measure of Phasing Error',
                                      value_name='Error Rate').replace('CHM13v2.0', 'T2T-CHM13')

dot_size=3
np.random.seed(chromosome_panel_seed)

sns.violinplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, hue_order=['GRCh38','T2T-CHM13'], palette=genome_palette, linewidth=0.5,
               dodge=True, split=True, cut=0.5, inner=None, density_norm='count', fill=True, alpha=0.2, bw_adjust=.5,
               data = compare_stats, ax=axB, legend=True)
sns.violinplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, hue_order=['GRCh38','T2T-CHM13'], palette=genome_palette, linewidth=0.5,
               dodge=True, split=True, cut=0.5, inner=None, density_norm='count', fill=False, bw_adjust=.5,
               data = compare_stats, ax=axB, legend=False)
sns.stripplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, jitter=0.25, s=dot_size,
               dodge=True, alpha=0.6, linewidth=0.1, palette={'GRCh38':'dimgrey', 'T2T-CHM13':'rebeccapurple'},
               data = compare_stats.loc[compare_stats.trio_phased==True], ax=axB, legend=True)
sns.stripplot(x='Measure of Phasing Error', y='Error Rate', hue='Genome',
               order = error_order, jitter=0.25, s=dot_size,
               dodge=True, alpha=0.6, linewidth=0.5, palette={'GRCh38':'slategrey', 'T2T-CHM13':'plum'},
               data = compare_stats.loc[compare_stats.ground_truth_data_source=='HGSVC_samples_nontrios_only'], ax=axB, legend=True, edgecolor='black')
axB.yaxis.set_major_formatter(perc_formatter)

handles, labels =axB.get_legend_handles_labels()
new_handles=list()
for i, (handle, label) in enumerate(zip(handles, labels)):
    handle.set_alpha(.8)
    if i == 0:
        handle.set_facecolor('#E0E0E0')
        handle.set_edgecolor('#666666')
        handle.set_alpha(1)
    elif i == 1:
        handle.set_facecolor('#E9E3F0')
        handle.set_edgecolor('#9467BD')
        handle.set_alpha(1)
    new_handles.append(handle)

labels = labels[:2] + [l + ' (Mendelian phased)' for l in labels[2:4]] + [l+ ' (Statistically phased)' for l in labels[4:6]]
axB.legend(new_handles,labels)

for i, ax in enumerate(axes.values()):
    add_letter_to_ax(ax, i, fontsize=8)

## Fig 4a - graphic
sns.despine(ax=axA, top=True, right=True, left=True, bottom=True)
axA.set_yticks([])
axA.set_xticks([])
plt.savefig(f'{figdir}/Figure 4 (a,b only).png', facecolor='white')
plt.savefig(f'{figdir}/Figure 4 (a,b only).svg', facecolor='white')
plt.savefig(f'{figdir}/Figure 4 (a,b only).eps', facecolor='white')
plt.savefig(f'{figdir}/Figure 4 (a,b only).pdf', facecolor='white')

for ext in ['svg','pdf']:
    shutil.copy(f'{figdir}/Figure 4 (a,b only).{ext}', f'{figdir}/Figure 4 (a,b only) with graphic.{ext}')
place_panel_a(f'{figdir}/Figure 4 (a,b only) with graphic.svg', f'{figdir}/Figure 4 (a,b only) with graphic.pdf', f'{source_figures_folder}/Fig 4a.svg', 'illustrator',
              axis_top_left_points(fig, axA), centre_y=axis_vertical_centre_points(fig, axB))
subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', f'{figdir}/Figure 4 (a,b only) with graphic.pdf', f'{figdir}/Figure 4 (a,b only) with graphic'])

# %% [markdown]
# ### Supplemental Figure 17
#
# Association between flip error rate, switch error rate, true switch error rate, and genotyping error rate.
#
# The flip error rate, switch error rate, and true switch error rate was calculated for all chromosomes and all samples for which assembly-based ground truths existed. Genotype error rates were calculated for each chromosome in each sample. Log10 values of each error rate were taken to transform each rate to a gaussian distribution. Linear regression was sequentially performed with the formula log10(Error_Rate) ~ log10(genotype_error_rate), and r2 correlation values were calculated. a) Per chromosome per sample error rates in the GRCh38 haplotype panel for Mendelian error-corrected samples assemblies generated by the HPRC or the HGSVC. b) Per chromosome per sample error rates in the GRCh38 haplotype panel for non-trio samples with HGSVC assemblies that did not undergo mendelian error correction. c) Per chromosome per sample error rates in the T2T-CHM13 haplotype panel for Mendelian error-corrected samples assemblies generated by the HPRC or the HGSVC. d) Per chromosome per sample error rates in the T2T-CHM13 haplotype panel for non-trio samples with HGSVC assemblies that did not undergo mendelian error correction. All regression lines are plotted with a +/- 95% confidence interval.
#

# %%
plt.rcParams['figure.dpi'] = 300
hprc_evaluated_no_PAR = full_panel_HPRC_compared.loc[full_panel_HPRC_compared.chrom.str.contains('chr') & ~full_panel_HPRC_compared.chrom.str.contains('Whole')].dropna(subset=['true_switch_error_rate'])
hprc_evaluated_no_PAR[['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']] = np.log10(hprc_evaluated_no_PAR[['switch_error_rate','flip_error_rate','true_switch_error_rate','gt_error_rate']])
hprc_evaluated_no_PAR.loc[hprc_evaluated_no_PAR.log_true_switch_error_rate == -np.inf,'log_true_switch_error_rate'] = 0

def normalize_df(df):
    return ((df-df.min())/(df.max()-df.min())).values

tmp=hprc_evaluated_no_PAR.melt(id_vars=['sample_id','log_gt_error_rate','trio_phased', 'genome'], value_vars=['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate'],var_name='error_kind',value_name='log_error_rate')
tmp = tmp.rename(columns={'genome':'Panel Genome', 'trio_phased':'Mendelian Pre-phased','log_gt_error_rate':'log(Genotype Error Rate)', 'error_kind':'Measurement of panel error', 'log_error_rate':'log(Error Rate)'})
tmp = tmp.replace('log_switch_error_rate', 'log(Switch Error Rate)').replace('log_flip_error_rate', 'log(Flip Error Rate)').replace('log_true_switch_error_rate', 'log(True Switch Error Rate)').replace('CHM13v2.0','T2T-CHM13')
tmp['Mendelian Pre-phased'] = tmp['Mendelian Pre-phased'].map({True:'Mendelian pre-phased', False:'Not Mendelian pre-phased'}) # panel titles in the legend's own words
error_palette = dict(zip(['log(Switch Error Rate)', 'log(Flip Error Rate)', 'log(True Switch Error Rate)'], sns.color_palette('deep')[0:3])) # the paper's palette, as in Extended Data 1 (whose purple is the T2T-CHM13 colour)
height_width_ratio = 0.5; aspect_ratio = 1.15 #(produces a two_col by 5 inch figure, found through trial and error)
g = sns.lmplot(data=tmp, col='Mendelian Pre-phased', row='Panel Genome', col_order=['Mendelian pre-phased', 'Not Mendelian pre-phased'],
           x='log(Genotype Error Rate)',y='log(Error Rate)',hue='Measurement of panel error', palette=error_palette, scatter_kws={'alpha':0.1}, units='sample_id', height=5*height_width_ratio, aspect=aspect_ratio,
           seed=seed) # the bootstrap bands were unseeded, so they changed with every run
g.set_titles('{row_name} panel, {col_name}', size=7) # was 'Panel Genome = GRCh38 | Mendelian Pre-phased = True'; 7 pt like the other figures' titles

for i, ax in enumerate(list(itertools.chain.from_iterable(g.axes))):
    add_letter_to_ax(ax, i, points_offset=(-25, 0), fontsize=8) # 8 pt panel letters for the journal (this is Extended Data now)

g.tight_layout()

# centre the key under the middle of the two columns of axes (it sat at a hand-set .46 of the figure)
axes_middle = (g.axes[0, 0].get_position().x0 + g.axes[0, 1].get_position().x1)/2
sns.move_legend(g, loc='upper center', ncols=3, bbox_to_anchor=(axes_middle, 0))
for lh in g._legend.legend_handles: # move_legend rebuilds the key, so this comes after it
    lh.set_alpha(1)

g.savefig(f'{figdir}/supplemental/Supplemental 17.svg', facecolor='white')
g.savefig(f'{figdir}/supplemental/Supplemental 17.png', facecolor='white')
g.savefig(f'{figdir}/supplemental/Supplemental 17.eps', facecolor='white')
g.savefig(f'{figdir}/supplemental/Supplemental 17.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Table - correlation of error rates
#
# Generate correlation_error_rates_and_gt_error_rate.csv from which some text values are drawn

# %%
grch38_trio_phased_corr = (hprc_evaluated_no_PAR.loc[(hprc_evaluated_no_PAR.trio_phased)&(hprc_evaluated_no_PAR.genome=='GRCh38'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='GRCh38', mendelian_phased='True').iloc[:-1]
grch38_stat_phased_corr = (hprc_evaluated_no_PAR.loc[(~hprc_evaluated_no_PAR.trio_phased)&(hprc_evaluated_no_PAR.genome=='GRCh38'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='GRCh38', mendelian_phased='False').iloc[:-1]
grch38_all_phased_corr = (hprc_evaluated_no_PAR.loc[(hprc_evaluated_no_PAR.genome=='GRCh38'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='GRCh38', mendelian_phased='Both').iloc[:-1]
chm13_trio_phased_corr = (hprc_evaluated_no_PAR.loc[(hprc_evaluated_no_PAR.trio_phased)&(hprc_evaluated_no_PAR.genome=='CHM13v2.0'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='CHM13v2.0', mendelian_phased='True').iloc[:-1]
chm13_stat_phased_corr = (hprc_evaluated_no_PAR.loc[(~hprc_evaluated_no_PAR.trio_phased)&(hprc_evaluated_no_PAR.genome=='CHM13v2.0'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='CHM13v2.0', mendelian_phased='False').iloc[:-1]
chm13_all_phased_corr = (hprc_evaluated_no_PAR.loc[(hprc_evaluated_no_PAR.genome=='CHM13v2.0'),['log_switch_error_rate','log_flip_error_rate','log_true_switch_error_rate','log_gt_error_rate']].corr(method='pearson')**2)[['log_gt_error_rate']].assign(genome='CHM13v2.0', mendelian_phased='Both').iloc[:-1]

chm13_phased_corr = pd.concat((chm13_stat_phased_corr, chm13_trio_phased_corr, chm13_all_phased_corr)).drop(columns='genome').reset_index().rename(columns={'index':'error_rate'}).set_index(['error_rate','mendelian_phased'])
grch38_phased_corr = pd.concat((grch38_stat_phased_corr, grch38_trio_phased_corr, grch38_all_phased_corr)).drop(columns='genome').reset_index().rename(columns={'index':'error_rate'}).set_index(['error_rate','mendelian_phased'])
error_rate_corr=chm13_phased_corr.join(grch38_phased_corr, lsuffix="_chm13", rsuffix='_grch38').reset_index().rename(columns={'error_rate':'Log error rate','mendelian_phased':'Mendelian-based error correction','log_gt_error_rate_chm13':'CHM13v2.0','log_gt_error_rate_grch38':'GRCh38'})
error_rate_corr = error_rate_corr.replace('log_switch_error_rate','Log(Switch Error Rate)').replace('log_flip_error_rate','Log(Flip Error Rate)').replace('log_true_switch_error_rate','Log(True Switch Error Rate)').set_index(["Log error rate",'Mendelian-based error correction'])

error_rate_corr.reset_index().to_csv(f'{figdir}/correlation_error_rates_and_gt_error_rate.csv', index=False)

# Export the three coefficients quoted in the Results from their authoritative
# facet so they cannot be copied from a different panel or cohort.
reported_r_squared = (
    error_rate_corr
    .reset_index()
    .loc[lambda d: d['Mendelian-based error correction'] == 'True',
         ['Log error rate', 'CHM13v2.0']]
    .rename(columns={'Log error rate': 'error_rate', 'CHM13v2.0': 'r_squared'})
)
reported_r_squared['error_rate'] = reported_r_squared['error_rate'].map({
    'Log(Switch Error Rate)': 'SER',
    'Log(Flip Error Rate)': 'FER',
    'Log(True Switch Error Rate)': 'tSER',
})
reported_r_squared = reported_r_squared.assign(
    panel='T2T-CHM13',
    cohort='Mendelian-pre-phased HPRC and HGSVC probands',
    observation='sample-assembly chromosome',
    n_observations=len(hprc_evaluated_no_PAR.loc[
        hprc_evaluated_no_PAR.trio_phased
        & (hprc_evaluated_no_PAR.genome == 'CHM13v2.0')
    ]),
    statistic='squared Pearson correlation',
    transform='log10(error rate); zero tSER values represented at log10=0',
)
reported_r_squared.to_csv(f'{figdir}/reported_error_rate_r_squared.csv', index=False)


# %% [markdown]
# ### Supplemental Figure 18
#
# True switch error rate as a function of GRCh38 chromosomal completeness.
#
# Scatterplot of all chromosome-wide switch error rate in the GRCh38 and T2T-CHM13 1KGP panels for all chromosomes from HPRC-assembled samples by the percentage of each chromosome that is newly resolved in T2T-CHM13. A linear regression (true_switch_error_rate ~ %_chrom_novel), run separately for each 1KGP panel, is also plotted. Shaded areas indicate 95% ci of the regression slope. A linear regression was also performed of ΔSER% ~ %_chrom_novel. Summary statistics of all three regressions are presented in the upper left corner of the plot.

# %%
# Paper cross-reference: Supplemental Figure 1u; supplemental lines 210-216; manuscript lines 315-321.
def plot_regression(x, y, hue, data, **kws):
    fig = sns.lmplot(x=x,y=y, hue='genome', data = data, robust=False, palette=default_style['palette'], aspect=1.5, height=4,
                     facet_kws={'legend_out':False}, units='sample_id', seed=seed, **kws)#, x_jitter=0.1) # seeded - the bootstrap bands changed with every run
    ax = fig.ax
    grch = data.loc[data.genome=='GRCh38']
    t2t = data.loc[(data.genome=='CHM13v2.0')|(data.genome=='T2T-CHM13')]
    s1, i1, r1, p1, err1 = scipy.stats.linregress(grch[x], grch[y])
    s2, i2, r2, p2, err2 = scipy.stats.linregress(t2t[x], t2t[y])
    s3, i3, r3, p3, err3 = scipy.stats.linregress(grch[x], grch[y].values-t2t[y].values)
    # r2, p2 = stats.pearsonr(t2t[x], t2t[y])
    # r3, p3 = stats.pearsonr(grch[x], grch[y].values-t2t[y].values)
    # linregress gives r; report r² (square it), as the text does for Extended Data 3. And the difference is grch - t2t, so it's labeled that way round
    print ('GRCh38: ' + "r²" + f"={r1**2:#.2g}, p={p1:.2g}, slope={s1:.1e}")
    print ('CHM13v2.0: ' + "r²" + f"={r2**2:#.2g}, p={p2:.2g}, slope={s2:.1e}")
    print ('GRCh38 - CHM13v2.0: ' + "r²" + f"={r3**2:#.2g}, p={p3:.2g}, slope={s3:.1e}")
    # one block, so the three pairs of lines are evenly spaced, at the 7 pt journal max (were six separately placed 6 pt lines)
    ax.text(.05, .97, f"GRCh38\nr²={r1**2:#.2g}, p={p1:.2g}, slope={s1*100:.3f}%\n\n"
                      f"T2T-CHM13\nr²={r2**2:#.2g}, p={p2:.2g}, slope={s2*100:.3f}%\n\n"
                      f"GRCh38 - T2T-CHM13\nr²={r3**2:#.2g}, p={p3:.2g}, slope={s3*100:.3f}%",
            transform=ax.transAxes, va='top', fontsize=7)

    return fig, ax

fg, ax = plot_regression(x='percent_new',y='true_switch_error_rate', hue='genome',
                         data = all_chroms.loc[(all_chroms.method_of_phasing=='phased_with_parents_and_pedigree')
                                                &(all_chroms.ground_truth_data_source == 'HPRC_samples') &
                                                ~(all_chroms.chrom.isin(('chrX','PAR1','PAR2')))
                                                ].replace('CHM13v2.0','T2T-CHM13'),
                          scatter_kws={'linewidths':0.1})
# ax.set_ylim(0, .16)
ax.set_ylabel('True Switch Error Rate')
ax.set_xlabel('% of T2T-CHM13 chromosome\nnot present in GRCh38')
ax.xaxis.set_major_formatter(percent_formatter)
ax.yaxis.set_major_formatter(PercentFormatter(decimals=2))
# ax.grid(axis='y', alpha=.5) # editor wants the background gridlines gone
fg.legend.set_title('Genome')
ax.set_axisbelow(True)
sns.despine(ax=ax, top=False, right=False, left=False, bottom=False)
sns.move_legend(fg, "upper left", bbox_to_anchor=(.8,0.98))

fg.tight_layout()

# this is Supplemental 18 - it was saving as 'Supplemental 16' and writing over the real Supplemental 16
plt.savefig(f'{figdir}/supplemental/Supplemental 18.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 18.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 18.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 18.pdf', facecolor='white')

# %% [markdown]
# Supplemental Figure 16: True switch error rate as a function of GRCh38 chromosomal completeness. Scatterplot of all chromosome-wide switch error rate in the GRCh38 and CHM13v2.0 1KGP panels for all chromosomes from HPRC-assembled samples by the percentage of each chromosome that is newly resolved in CHM13v2.0. A linear regression (true_switch_error_rate ~ %_chrom_novel), run separately for each 1KGP panel, is also plotted. Shaded areas indicate 95% ci of the regression slope. A linear regression was also performed of ΔSER% ~ %_chrom_novel. Summary statistics of all three regressions are presented in the upper left corner of the plot.

# %% [markdown]
# ## Figure 5: out-of-panel phasing performance
#
# Switch error rates of SNPs from samples statistically phased with the 1KGP haplotype panel used as a reference.
#
# All HPRC-assembled samples and their relatives were removed from the GRCh38 1KGP and the CHM13v2.0 1KGP panels of 2504 unrelated samples, yielding 'non-HPRC' panels containing data from 2426 individuals. a) Filtered short read derived SNPs from 39 HPRC-assembled individuals were phased, using either the GRCh38 or CHM13v2.0 non-HPRC reference haplotype panel. Switch error rates were measured by mendelian concordance with parental genomes derived from short read variant calls. b) SNPs present in the HPRC human pangenome were unphased, then statistically phased with either the GRCh38 or CHM13v2.0 non-HPRC reference haplotype panel. Switch error rates were measured using the pangenome haplotypes as ground truth. Switch error rates and average minor allele frequency per bin are displayed on a log scale.

# %%
# Paper cross-reference: Figure 5; manuscript lines 325-342.
# Figure 5
fig, (axA, axB) = plt.subplots(1,2,figsize=(two_col,3), layout='tight')
sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                    & (binned_maf_data.method_of_phasing == '1kgp_variation_phased_with_reference_panel')
                                    & (binned_maf_data.syntenic!='All')
                                    & (binned_maf_data.type=='SNP')],
             legend=True, ax=axA, **syntenic_style)

sns.lineplot(x='MAF',y='switch_error_rate', style='syntenic', hue='genome', hue_order = genome_order, style_order = syntenic_order,
             data=binned_maf_data.loc[(binned_maf_data.ground_truth_data_source == 'HPRC_samples')
                                    & (binned_maf_data.method_of_phasing == '1kgp_variation_phased_with_reference_panel')
                                    & (binned_maf_data.syntenic!='All')
                                    & (binned_maf_data.type=='Indel')],
             legend=True, ax=axB, **syntenic_style)

axA.text(0.6,1.1, 'SNPs', fontsize=7, transform=axA.transAxes, weight='bold', horizontalalignment='right', verticalalignment='center') # 7 pt journal max
axB.text(0.6,1.1, 'Indels', fontsize=7, transform=axB.transAxes, weight='bold', horizontalalignment='right', verticalalignment='center')

for i, ax in enumerate([axA, axB]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    # ax.grid() # editor wants the background gridlines gone
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    if i in (2,3):
        ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    else:
        ax.set_xlabel(None)
    if i in (0,2):
        ax.set_ylabel('Switch Error Rate (%)', labelpad=0)
    else:
        ax.set_ylabel(None)
fig.supxlabel('Minor Allele Frequency (%)', fontsize=7) # default 'large' comes out at 7.2 pt


for i, ax in enumerate([axA, axB]): # legend axes
    if i == 0:
        loc='upper right'
    else:
        loc='upper right'
    handles, labels = ax.get_legend_handles_labels()
    labels=[x.replace('Syntenic','Shared genomic sequences') for x in labels]
    legend=ax.legend(handles=handles, labels=labels)
    for handle in handles:
        handle.set_marker('o')
    handles = handles[1:3] + handles[5:6]
    labels = labels[1:3] + labels[5:6]
    try:
        handles[-1].set_color(handles[1].get_color())
    except:
        print (handles)
    legend = ax.legend(ncols=1, loc=loc, handles=handles, labels=labels)
    legend.set_title(None)
    legend.get_frame().set(
        alpha=0.9, boxstyle='square', edgecolor='grey', linewidth=0)
    for text in legend.get_texts():
        t = text.get_text()
        new_t = update_legend_values.get(t,t)
        text.set_text(new_t)
    add_letter_to_ax(ax, i, points_offset=(-25, 7), fontsize=8) # 8 pt bold panel letters for the journal


# clean_figure(fig)
plt.savefig(f'{figdir}/Figure 5.svg', facecolor='white')
plt.savefig(f'{figdir}/Figure 5.png', facecolor='white')
plt.savefig(f'{figdir}/Figure 5.pdf', facecolor='white')
plt.savefig(f'{figdir}/Figure 5.eps', facecolor='white')


# %% [markdown]
# ## Figure 6: CNV-prone loci and ideograms
#
# Generated outside this notebook by the Figure 6/ideogram R script. Existing outputs are in `figures/figure6/`.
#

# %% [markdown]
# ### Supplemental Figure 19
#
# Genotype error rates of 1KGP variants located in segmental duplications (SDs) annotated in both GRCh38 and T2T-CHM13
#
# Variants from the Bykstrom-Bishop (2022) GRCh38 reference haplotype panel and this paper’s T2T-CHM13 panel were stratified by whether they overlapped segmental duplication regions as defined in Vollger (2022). They were additionally binned into 21 distinct minor allele frequency bins. Genotyping error rates of each bin were calculated as # of genotyping errors/# alt allele genotype calls in that bin. a) Genotyping error rates of SNPs from 39 trio-corrected HPRC samples, with HPRC genomes used as a source of ground truth. b) Genotyping error rates of SNPs from 19 1KGP samples that were not trio-corrected, with assembled HGSVC genomes used as a source of ground truth. c) Fold increase in the genotype error rate of SD SNPs over non-SD SNPs by MAF bin, stratified by panel. d) Genotyping error rates of indels from the same 39 HPRC samples as in figure a. e) Genotyping error rates of indels from the same 19 HGSVC samples as in panel b. f) Fold increase in indel genotyping error rate (SD/non-SD).
#

# %%
region_order=['not_in_segdups','in_segdups']

region_style = {k:v for k,v in syntenic_style.items()}
region_style['markers'] = True
del region_style['dashes']
ratio_style = {k:v for k,v in region_style.items() if k != 'markers'} # c and f are ratios, not regions: plain lines, no markers

fig, ((axA, axB, axC), (axD, axE, axF)) = plt.subplots(2,3, figsize=(two_col, 4), layout='tight')
fig.get_layout_engine().set(rect=(0, 0.15, 1, 1)) # leave room at the bottom for the keys - they were hanging off the page and getting clipped (0.13 before the keys went to 7 pt)

sns.lineplot(x='MAF',y='gt_error_rate', style='region', hue='genome', hue_order = genome_order, style_order = region_order,
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.type == 'SNP')],
             legend=True, ax=axA, **region_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='region', hue='genome', hue_order = genome_order, style_order = region_order,
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.type == 'SNP')],
             legend=False, ax=axB, **region_style)

sns.lineplot(x='MAF',y='gt_error_rate', style='region', hue='genome', hue_order = genome_order, style_order = region_order,
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'HPRC_samples')
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.type == 'Indel')],
             legend=False, ax=axD, **region_style)


sns.lineplot(x='MAF',y='gt_error_rate', style='region', hue='genome', hue_order = genome_order, style_order = region_order,
             data=binned_maf_data_regions.loc[(binned_maf_data_regions.ground_truth_data_source == 'HGSVC_samples_nontrios_only') #trios only
                                      & (binned_maf_data_regions.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (binned_maf_data_regions.type == 'Indel')],
             legend=False, ax=axE, **region_style)


id_cols=['genome','ground_truth_data_source','method_of_phasing','type','rounded_MAF']
segdups = binned_maf_data_regions.loc[binned_maf_data_regions.region.str.contains('segdup')]
segdups= segdups.pivot(index=id_cols,
                                    columns='region', values='gt_error_rate').reset_index(
               ).merge(segdups.groupby(id_cols, observed=True).first(numeric_only=True).reset_index()[id_cols + ['MAF']], on=id_cols, how='left')
segdups['diff'] = segdups.in_segdups/segdups.not_in_segdups

sns.lineplot(x='MAF',y='diff', hue='genome', hue_order = genome_order,
             data=segdups.loc[(segdups.ground_truth_data_source == 'HPRC_samples')
                                      & (segdups.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (segdups.type == 'SNP')],
                                    #   & (binned_maf_data_regions.syntenic != 'All')],
             legend=False, ax=axC, **ratio_style) # the keys at the bottom cover c and f now


sns.lineplot(x='MAF',y='diff', hue='genome', hue_order = genome_order,
             data=segdups.loc[(segdups.ground_truth_data_source == 'HPRC_samples') # d's HPRC samples, to match c (was 'HGSVC_samples_nontrios_only', e's samples)
                                      & (segdups.method_of_phasing == 'phased_with_parents_and_pedigree')
                                      & (segdups.type == 'Indel')],
                                    #   & (binned_maf_data_regions.syntenic != 'All')],
             legend=False, ax=axF, **ratio_style)



for i, ax in enumerate([axA, axB, axD, axE]):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-1, 50)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    ax.grid() # kept here: the gridlines help read these log-log panels (the editor said remove them "wherever appropriate")
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    ax.yaxis.set_major_formatter(perc_formatter)
    ax.set_xlabel('Minor Allele Frequency (%)', labelpad=3)
    ax.set_ylabel('Genotype Error Rate (%)', labelpad=0)

for i, ax in enumerate([axC, axF]):
    # ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(0, 18)
    ax.set_xlim(1e-2, 55)
    # ax.set_box_aspect(1)
    ax.grid() # kept here, as in the other four panels
    ax.set_axisbelow(True)
    ax.set_xticks(MAF_xticks)
    # ax.set_yticks(SER_xticks)
    ax.xaxis.set_major_formatter(perc_formatter)
    # ax.yaxis.set_major_formatter(perc_formatter)
    ax.set_xlabel('Minor Allele Frequency (%)')
    ax.set_ylabel('Fold increase in\ngenotype error rate') # was 'Fold increase in segmental duplication\ngenotype error rate (%)': too long (it ran into the panel letter), and a ratio has no %
    ax.axhline(1, color='grey', linestyle='--', linewidth=1, alpha=0.7)



leg = axA.get_legend()
handles = [handle for handle in leg.legend_handles]
labels  = [update_legend_values.get(t.get_text(), t.get_text()) for t in leg.get_texts()]

leg.set_visible(False)
leg.set_in_layout(False)
leg.remove()
fig.legend(
    handles=handles[1:3],
    labels=labels[1:3],
    loc="upper center",
    bbox_to_anchor=(0.14, 0.14), # inside the strip left for the keys (was 0.02, off the page); three keys now
    ncol=len(handles[1:3]),
    frameon=True,
    title=labels[0],
    columnspacing=1.0, # tighter (was 1.4) so three 7 pt keys fit across the page
    handlelength=2.5, # was 4
    markerscale=1.5, # was 2
    prop={"size":7}, # 7 pt, the journal max (was 6)
    title_fontsize=7, # 7 pt journal max
    labelspacing=1
)
fig.legend(
    handles=handles[4:6],
    labels=labels[4:6],
    loc="upper center",
    bbox_to_anchor=(0.5, 0.14),
    ncol=len(handles[4:6]),
    frameon=True,
    title=labels[3] + ' (a, b, d, e)', # the region lines are only in these four panels
    columnspacing=1.0, # tighter (was 1.4) so three 7 pt keys fit across the page
    handlelength=2.5, # was 4
    markerscale=1.5, # was 2
    prop={"size":7},
    title_fontsize=7,
    labelspacing=1
)
fig.legend( # c and f plot the ratio of the two regions, not a region - say so in their own key
    handles=[mpl.lines.Line2D([], [], color='black', **{k:v for k,v in ratio_style.items() if k not in ('palette', 'errorbar')})],
    labels=['In SD ÷ not in SD'],
    loc="upper center",
    bbox_to_anchor=(0.86, 0.14),
    frameon=True,
    title='Fold increase (c, f)',
    handlelength=2.5,
    markerscale=1.5,
    prop={"size":7},
    title_fontsize=7,
    labelspacing=1 # as in the other two keys - it also sets the gap under the title, so their rows line up
)
# the rows were never labelled: both are genotype error, a-c SNPs and d-f indels
axA.text(-.25,0.5, 'SNPs', fontsize=7, weight='bold', transform=axA.transAxes,  rotation='vertical', horizontalalignment='right', verticalalignment='center')
axD.text(-.25,0.5, 'Indels', fontsize=7, weight='bold', transform=axD.transAxes,  rotation='vertical', horizontalalignment='right', verticalalignment='center')
# axA.text(0.5,1.2, 'SNPs', fontsize=8, transform=axA.transAxes,  rotation='horizontal', horizontalalignment='center', verticalalignment='bottom')
# axB.text(0.5,1.2, 'Indels', fontsize=8, transform=axB.transAxes,  rotation='horizontal', horizontalalignment='center', verticalalignment='bottom')


clean_figure(fig, letter_fontsize=8) # 8 pt panel letters for the journal (this is Extended Data now)
plt.savefig(f'{figdir}/supplemental/Supplemental_19_segmental_duplication_genotype_error_with_ratio.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental_19_segmental_duplication_genotype_error_with_ratio.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental_19_segmental_duplication_genotype_error_with_ratio.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental_19_segmental_duplication_genotype_error_with_ratio.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Figure 20
#
# UCSC Genome Browser information on mappability of a region near the Prader-Willi breakpoint 2 (roughly chr15:24100000-24600000) that is enriched for genetic variants in the GRCh38 phased genomic panel. a) Mappability, blacklist status, location of segmental duplications, and genic regions in GRCh38. b) Mappability, location of segmental duplications, “difficult region” status from GIAB, and genic regions in CHM13v2.0.
#
# ![image.png](attachment:image.png)

# %% [markdown]
# ### Supplemental Figure 21
#
# A violin plot displaying the distribution of overall per-sample r2 values when imputing GRCh38 SGDP variation, stratified by SGDP-defined sample continental groups.
#
# Distributions are split by reference panel used during imputation: either the Byrska-Bishop 2022 GRCh38 reference haplotype panel, or the T2T-CHM13 panel lifted over to GRCh38 coordinates.
#

# %%
# Figure 8 supplements
fig, ax = plt.subplots(1,1,figsize=(two_col,3), layout='tight')
plt.rcParams["figure.dpi"] = 300

ax = sns.violinplot(x='ancestry', y='imputed_ds_rsquared', hue='Reference Panel', ax=ax,
                 hue_order=['GRCh38','T2T-CHM13\n(lifted to GRCh38)'], palette=genome_palette, #s=4,
                 data=separate_sample_panel_and_subset.loc[(separate_sample_panel_and_subset.dataset=='SGDP')
                                      & (separate_sample_panel_and_subset['Variant Subset']=='All panel variants')
                                      & (separate_sample_panel_and_subset['Variant Type']=='All variants')
                                      & (separate_sample_panel_and_subset.genome=='GRCh38')].replace({'CHM13v2.0':'T2T-CHM13\n(lifted to GRCh38)'}),
                                    #   & (separate_sample_panel_and_subset.ancestry != 'all')],
                 legend=True, alpha=0.75, split=True, cut=0, gap=0, inner='quart', saturation=1,
                 density_norm='width',common_norm=True)

formatted_ancestry_names = {'all':'All\npopulations','WestEurasia':'West\nEurasia','SouthAsia':"South\nAsia",
                            'EastAsia':'East\nAsia','Africa':'Africa', 'America':'Americas','Oceania':'Oceania',
                            'CentralAsiaSiberia':'Central Asia\n& Siberia'}

xticks = list(ax.get_xticklabels())
for tick in xticks:
    t = tick.get_text()
    tick.set_text(formatted_ancestry_names.get(t, t))
ax.set_xticklabels(xticks)

ax.set_ylim(0.8,1.01)

legend = ax.legend()
for text in legend.get_texts():
    t = text.get_text()
    new_t = update_legend_values.get(t,t)
    if new_t == 'T2T-CHM13':
        new_t = 'T2T-CHM13\n(lifted to GRCh38)'
    print (new_t)
    text.set_text(new_t)
legend.set_title('Genome') # This shouldn't be necessary, but it is and I don't really care to dig into why
sns.move_legend(ax, "lower left")#, bbox_to_anchor=(1, 0))

ax.set_ylabel('Average overall sample r²') # the mathtext superscript 2 comes out at 4.2 pt, under the journal's 5 pt minimum
ax.set_xlabel('SGDP Ancestry Grouping')
# ax.grid(axis='y', alpha=.5) # editor wants the background gridlines gone
ax.set_axisbelow(True)
clean_figure(fig, add_letters=False)
add_letter_to_ax(ax, 'b', fontsize=8) # this is panel b of the ancestry Extended Data figure now. 8 pt for the journal
# sns.despine()


plt.savefig(f'{figdir}/supplemental/Supplemental 21.png', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 21.svg', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 21.eps', facecolor='white')
plt.savefig(f'{figdir}/supplemental/Supplemental 21.pdf', facecolor='white')
supp21_fig = fig # the ancestry Extended Data figure needs this one again after Supplemental 22 is drawn

# %% [markdown]
# ## Figure 7: imputation from native and lifted reference panels
#
# Imputation of genomic variation in 256 non-1KGP Human Genome Diversity Project (HGDP) samples, using 1KGP haplotype panels as references.
#
# a) Genotyping array data was simulated by downsampling variation derived from short reads aligned to GRCh38 were downsampled to those variants present in the Infinium Omni2.5 genotyping array. We then phased this 'array' variation and imputed missing variants. At both steps, we provided a set of reference haplotypes to guide phasing and imputation - either the GRCh38 haplotype panel, or the CHM13v2.0 panel lifted to GRCh38 coordinates. R2 values were calculated using variants that were imputed by both the GRCh38 and the CHM13 panels. b) We downsampled variation by intersecting the Infinium Omni2.5 genotyping array variation (this time in CHM13v2.0 coordinates) with HGDP variation called form short reads aligned to CHM13v2.0. We then phased and imputed these samples using either the 1KGP CHM13v2.0 reference haplotype panel, or the 1KGP GRCh38 reference panel lifted to GRCh38 coordinates. All imputed variants were binned by the minor allele frequency of the variant in the reference panel used during imputation. R2 values were then calculated per bin. For both panels, R2 values were calculated using variants that were present in both the GRCh38-informed imputed variant data and the CHM13v2.0-informed imputed variant data. R2 statistics are stratified by SNPs and Indels. Average minor allele frequency per bin is displayed on a log scale.

# %%
genome_order = ['GRCh38','CHM13v2.0']
panel_order=['native_panel.common_variants','lifted_panel.common_variants']
g=sns.relplot(kind='line', data = separate_panel_and_subset.loc[(separate_panel_and_subset.Synteny == 'Syntenic')
                                                                 & (separate_panel_and_subset['Variant Subset']=='Common variants')
                                                                 & (separate_panel_and_subset.ancestry=='all')
                                                                 & (separate_panel_and_subset.dataset=='SGDP')
                                                                 & (separate_panel_and_subset.panel_filter=='No Filter')].replace('CHM13v2.0','T2T-CHM13'),
                x='mean_AF',y='imputed_ds_rsquared', hue='Reference Panel', style='Variant Type', hue_order=['GRCh38','T2T-CHM13'], style_order=['SNPs','Indels'],#style_order=['All variants','Indels','SNPs'],
                palette=default_style['palette'], col='genome', height=6, aspect=1, col_order=['GRCh38','T2T-CHM13']), #, col_wrap=4)#, , height=3, aspect=1, palette=None, row_order=None, col_order=None, hue_order=None, hue_kws=None, dropna=False, legend_out=True, despine=True, margin_titles=False, xlim=None, ylim=None, subplot_kws=None, gridspec_kws=None)

fig = (g[0].set(xscale='log', ylim=(0,1))
  .set_axis_labels("Allele Frequency", "Variant Imputation r²") # the mathtext superscript 2 came out at 4.2 pt, under the journal's 5 pt minimum
  .set_titles(template="{col_name} SGDP variation", size=7) # titles 7 pt, as in the other figures
  .tight_layout(w_pad=0)
  )

fig.figure.set_size_inches((oneandahalf_col,2))

for col_val, ax in fig.axes_dict.items():
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=1))
    # ax.grid(alpha=0.5) # editor wants the background gridlines gone
sns.move_legend(fig.figure, "lower right", bbox_to_anchor=(.92, .1))

for i, ax in enumerate(list(itertools.chain.from_iterable(fig.axes))):
    add_letter_to_ax(ax, i, points_offset=(-7, 7), fontsize=8) # 8 pt bold panel letters for the journal


# Create separate legends for hue and style
# Remove the default combined legend

# Get the first axis to extract legend information

# Create separate legend handles
axes=fig.axes[0]
handles, labels = axes[0].get_legend_handles_labels()
g[0].figure.legends[0].remove()
ax0_labels = [labels[0]] + ['GRCh38', 'T2T-CHM13\n(lifted to GRCh38)']+ labels[3:]
ax1_labels = [labels[0]] + ['GRCh38\n(lifted to T2T-CHM13)', 'T2T-CHM13'] + labels[3:]

axes[0].legend(handles=handles, labels=ax0_labels)
axes[1].legend(handles=handles, labels=ax1_labels)

plt.savefig(f'{figdir}/Figure 7.png', facecolor='white', bbox_inches='tight')
plt.savefig(f'{figdir}/Figure 7.svg', facecolor='white', bbox_inches='tight')
plt.savefig(f'{figdir}/Figure 7.eps', facecolor='white', bbox_inches='tight')
plt.savefig(f'{figdir}/Figure 7.pdf', facecolor='white', bbox_inches='tight')


# %% [markdown]
# ### Supplemental Figure 22
#
# Imputation of 256 non-1KGP HGDP GRCh38 genetic variation using 1KGP haplotype panels as references, stratified by SGDP population of origin.
#
# As previously presented, genotyping array data was simulated by downsampling variation derived from short reads aligned to GRCh38 were downsampled to those variants present in the Infinium Omni2.5 genotyping array. Variants present in our callset and in the genotyping array were then phased and imputed using either the GRCh38 1KGP haplotype panel, or the T2T-CHM13 1KGP panel lifted to GRCh38 coordinates. r2 values were calculated from the intersection of variants that were in both callsets. r2 statistics were calculated from imputed variants in samples assigned by the SGDP to the indicated ancestry group. r2 values were separately calculated in SNP variation and Indel variation. Average minor allele frequency per bin is displayed on a log scale.

# %%
per_pop_data = separate_panel_and_subset.loc[(separate_panel_and_subset.Synteny == 'All')#|(per_variant_category_imputation_performance.genome=='CHM13v2.0'))
                                                                 & (separate_panel_and_subset['Variant Subset']=='Common variants')
                                                                  &  (separate_panel_and_subset.dataset=='SGDP')
                                                                  & (separate_panel_and_subset.genome=='GRCh38')
                                                                  & (separate_panel_and_subset.panel_filter=='No Filter')].replace('CHM13v2.0','T2T-CHM13\n(lifted to GRCh38)')
per_pop_data['ancestry'] = per_pop_data.ancestry.map(lambda x: formatted_ancestry_names.get(x,x).replace('\n',' '))
panel_order=['native.common_variants','lifted.common_variants']
fig=sns.relplot(kind='line', data = per_pop_data,
                x='mean_AF',y='imputed_ds_rsquared', hue='Reference Panel', style='Variant Type', hue_order=['GRCh38','T2T-CHM13\n(lifted to GRCh38)'], style_order=['SNPs','Indels'],#style_order=['All variants','Indels','SNPs'],
                errorbar=None, palette=default_style['palette'], col='ancestry', height=4, aspect=1, col_wrap=4,facet_kws=dict(legend_out=False))
#, sharex=True, sharey=True, height=3, aspect=1, palette=None, row_order=None, col_order=None, hue_order=None, hue_kws=None, dropna=False, legend_out=True, despine=True, margin_titles=False, xlim=None, ylim=None, subplot_kws=None, gridspec_kws=None)
  # .set_ylabel('Average overall sample $\mathregular{r^2}$')
fig.figure.set_tight_layout(True)

(fig.set(xscale='log')
  .set_axis_labels("Allele Frequency", "r²") # the mathtext superscript 2 comes out at 4.2 pt, under the journal's 5 pt minimum
  .set_titles("{col_name}", size=7) # titles 7 pt, as in the other figures
  .tight_layout(w_pad=0,h_pad=5))
fig.figure.set_size_inches((two_col,4))

for col_val, ax in fig.axes_dict.items():
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=1))
    # ax.grid(alpha=0.5) # editor wants the background gridlines gone
    ax.set_ylim(0,1)

add_letter_to_ax(fig.axes[0], 'a', points_offset=(-25, 26), va='center', fontsize=8) # panel a of the ancestry Extended Data figure, on the row of the key above it
# the tight_layout above happens at the 16x8 inch relplot size, before set_size_inches, so at the real size the tick labels
# and titles hang off the page (the saved file came out 188 mm wide). Redo it at the real size
fig.tight_layout(w_pad=0,h_pad=5)

# the key goes above a now (it was under a, at a hand-set .48 of the figure), centred over the axes and on the letter's row:
# 26 pt above the top of the axes, like the letter, rather than at the figure's top edge, which moves when the layout makes room for the letter
axes_middle = (fig.axes[0].get_position().x0 + fig.axes[3].get_position().x1)/2
key_y = fig.axes[0].get_position().y1 + 26/72/fig.figure.get_figheight()
sns.move_legend(
    fig, "center",
    bbox_to_anchor=(axes_middle, key_y), ncol=6, title=None, frameon=True
)

fig.savefig(f'{figdir}/supplemental/Supplemental 22.svg', facecolor='white',bbox_inches='tight')
fig.savefig(f'{figdir}/supplemental/Supplemental 22.png', facecolor='white',bbox_inches='tight')
fig.savefig(f'{figdir}/supplemental/Supplemental 22.eps', facecolor='white',bbox_inches='tight')
fig.savefig(f'{figdir}/supplemental/Supplemental 22.pdf', facecolor='white',bbox_inches='tight')

# %% [markdown]
# ### Extended Data: ancestry imputation (Supplemental 22 over 21)
#
# Supplemental 22 (a) stacked on top of Supplemental 21 (b) on one 180 mm page, both at full size with their y axes lined up. Needs the Supplemental 21 and Supplemental 22 cells to have run first.

# %%
import subprocess

supp22_fig = fig.figure # fig is the Supplemental 22 FacetGrid from the cell above
page_width = supp21_fig.get_size_inches()[0]*72 # in points
tight22 = supp22_fig.get_tightbbox(supp22_fig.canvas.get_renderer())
pad = plt.rcParams['savefig.pad_inches']*72 # bbox_inches='tight' puts this much white space around the saved Supplemental 22

# line up the y axes of a and b, but don't let a hang off the page
x22 = axis_left_points(supp21_fig, supp21_fig.axes[0]) - axis_left_points(supp22_fig, fig.axes[0], tight22.x0)
x22 = min(max(x22, 0), page_width - tight22.width*72)

stack_figure_panels([dict(pdf=f'{figdir}/supplemental/Supplemental 22.pdf', svg=f'{figdir}/supplemental/Supplemental 22.svg',
                          clip=(pad, pad, tight22.width*72, tight22.height*72), x=x22),
                     dict(pdf=f'{figdir}/supplemental/Supplemental 21.pdf', svg=f'{figdir}/supplemental/Supplemental 21.svg',
                          clip=(0, 0, page_width, supp21_fig.get_size_inches()[1]*72), x=0)],
                    page_width, f'{figdir}/Extended Data - ancestry imputation.pdf', f'{figdir}/Extended Data - ancestry imputation.svg', gap=2) # was 8: less space between a and b
subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', f'{figdir}/Extended Data - ancestry imputation.pdf', f'{figdir}/Extended Data - ancestry imputation'])

# %% [markdown]
# ### Supplemental Figure 23
#
# Impact of including singleton variants in reference panels on imputation accuracy.
#
# A) SGDP GRCh38 variants downsampled to OMNI2.5 sites were imputed using 1KGP reference haplotypes constructed in either GRCh38 or T2T-CHM13 coordinates (the latter lifted to GRCh38). Reference panels were processed identically and included either all variants or excluded singletons (minor allele count = 1); related samples were removed, yielding 2,504 unrelated individuals per panel. Imputation accuracy was measured as genotype dosage r² against whole-genome sequencing truth and summarized across 21 bins of reference panel minor allele frequency. B) Difference in either SNP or Indel MAF bin imputation accuracy (Δ r²) between GRCh38 variants imputed using lifted T2T-CHM13 reference panels that include versus exclude singleton variants. Positive values indicate improved accuracy when singletons are included in the reference haplotype panel during imputation.
#

# %%
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(two_col, 2.5), layout='constrained')

id_cols=['genome','ancestry','dataset','panel','Reference Panel','variant_bin','Synteny','Variant Type','Variant Subset']
ax2_data = separate_panel_and_subset.loc[        (separate_panel_and_subset.dataset=='SGDP')
                                               & (separate_panel_and_subset.ancestry=='all')
                                               & (separate_panel_and_subset['Variant Subset']=='All panel variants')
                                               & (separate_panel_and_subset.Synteny=='Syntenic')
                                               & (separate_panel_and_subset.panel_filter.isin(['No Filter', 'no_singletons']))
                                               & (separate_panel_and_subset.info_cutoff=='0')
                                               & (separate_panel_and_subset.genome=='GRCh38')
                                               & (separate_panel_and_subset['Variant Type']!='All variants')].copy()

ax1_data = ax2_data.loc[ax2_data['Variant Type']=='SNPs'].copy()
mafs = ax2_data[['mean_AF'] + id_cols].drop_duplicates(id_cols).copy()
ax2_data_pivot=ax2_data.pivot(index=id_cols,
                        columns='panel_filter',
                        values='imputed_ds_rsquared').reset_index()
ax2_data_pivot = ax2_data_pivot.merge(mafs, on=id_cols, how='left')
ax2_data_pivot['rsquared_diff'] = ax2_data_pivot['No Filter'] - ax2_data_pivot['no_singletons']

### Plot data

ax1 = sns.lineplot(x='mean_AF',y='imputed_ds_rsquared', hue='Reference Panel', style='panel_filter',
             data=ax1_data, hue_order=genome_order, style_order=singletons_order, palette=genome_palette,
             legend=True, ax=ax1)

ax1.set_xscale('log')

# And ax 2
ax2=sns.lineplot(x='mean_AF', y='rsquared_diff', hue='Variant Type', #hue='Reference Panel',
             hue_order=['SNPs','Indels'],
             data=ax2_data_pivot.loc[ax2_data_pivot['Reference Panel']=='CHM13v2.0'], ax=ax2)

ax2.set_xscale('log')
ax2.axhline(0, color='grey', linestyle='--')

clean_figure(fig, letter_fontsize=8) # 8 pt panel letters for the journal (this is Extended Data now)
sns.move_legend(ax1, "upper left", bbox_to_anchor=(0, 1), frameon=True)

# Create inset axes
axins = inset_axes(ax1, width="40%", height="40%", loc='lower left',
                   bbox_to_anchor=(0.5, 0.15, .95, .95), bbox_transform=ax1.transAxes)
# Plot the same data in the inset
sns.lineplot(x='mean_AF',y='imputed_ds_rsquared', hue='Reference Panel', style='panel_filter',
                data=ax1_data, hue_order=genome_order, style_order=singletons_order, palette=genome_palette,
             legend=False, ax=axins, linewidth=1)

# Set the zoomed-in limits
axins.set_xlim(5, 50)
axins.set_ylim(0.98, 0.99)
clean_axis(axins)
axins.set_xlabel('')
axins.set_ylabel('')

# Add a box showing the zoomed region in the main plot
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.5", linestyle='--')

# no bbox_inches='tight' - its padding made the page 183 mm wide, over the journal's 180. constrained layout keeps everything on the page anyway
fig.savefig(f'{figdir}/supplemental/Supplemental 23.svg', facecolor='white')
fig.savefig(f'{figdir}/supplemental/Supplemental 23.png', facecolor='white')
fig.savefig(f'{figdir}/supplemental/Supplemental 23.eps', facecolor='white')
fig.savefig(f'{figdir}/supplemental/Supplemental 23.pdf', facecolor='white')

# %% [markdown]
# ### Supplemental Figure 24
#
# Improved liftover algorithm results in better imputation of common indel variants.
#
# a) Scatter plot of the r2 value of imputed indel variants when using a GRCh38 1KGP reference panel (x-axis) or the T2T-CHM13 panel lifted to GRCh38 coordinates using either a) GATK LiftoverVCF or b) the modified liftover algorithm implemented in this paper. Manual inspection of the high MAF variants that were poorly called by the lifted-over T2T-CHM13 panel showed that the majority of these variants are short tandem repeats that are lifted to a site that is several bases away from the same variant in the 1KGP GRCh38-native variant callset. Each dot represents one biallelic indel. The panel minor allele frequency is colored by 1KGP GRCh38 MAF.
#

# %%
import importlib
from scripts.analysis import whiskers as wsk
importlib.reload(wsk)
import subprocess



chrom='chrX'
new_genome='GRCh38'
genome='GRCh38'
old_genome='CHM13v2.0'

def get_dict(bcf):
    cmd = [
    "bcftools", "query",
    "-f", "%ID\t%INFO/SRC_CHROM\_%INFO/SRC_POS\_%INFO/SRC_REF_ALT\n",
    str(bcf)
    ]

    out=dict()
    print (' '.join(cmd))
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            key, val = line.strip().split('\t')
            out[key]=val.replace(',','_')

    ret = p.wait()
    if ret != 0:
        raise RuntimeError(f"bcftools failed with exit code {ret}")
    return out

comparative_indel_r2=list()
for liftover_method in ['LiftoverIndel','liftover_plugin']:
    basefolder=f'{imputation_statistics_folder}/imputation_results_true_false_0.05_false_1_CHM13v2.0'
    if liftover_method=='liftover_plugin':
        basefolder=basefolder+'.liftover_plugin'
    else:
        basefolder=basefolder+'.liftoverIndel_straight_omni_downsample'
    for chrom in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']:
        for genome in ['T2T','GRCh38']:
            if genome == 'GRCh38':
                old_genome='CHM13v2.0'
                bcf_genome='GRCh38'
            else:
                old_genome='GRCh38'
                bcf_genome='CHM13v2.0'

            lifted_stats=f'{basefolder}/per_chroms/{genome}.SGDP.lifted_panel.common_variants.{chrom}.glimpse2_concordance_r2_bins_r2_sites.txt.gz'
            native_stats=f'{basefolder}/per_chroms/{genome}.SGDP.native_panel.common_variants.{chrom}.glimpse2_concordance_r2_bins_r2_sites.txt.gz'
            outfile=f'{basefolder}/per_chroms/{genome}.{chrom}.combined_indel_stats.csv'

            try:
                native_stats=pl.scan_csv(native_stats, separator='\t')
            except:
                print(native_stats,'not found')
                raise
            native_indels=native_stats.filter((pl.col('allele1').str.len_chars() > 1) | (pl.col('allele2').str.len_chars() > 1)
                                              ).select(['rsid','ds_r2','maf'])
            lifted_stats=pl.scan_csv(lifted_stats, separator='\t')
            lifted_indels=lifted_stats.filter((pl.col('allele1').str.len_chars() > 1) | (pl.col('allele2').str.len_chars() > 1)
                                              ).select(['rsid','ds_r2'])
            all_stats = native_indels.join(lifted_indels, on='rsid', how='full', suffix='_lifted'
                                           ).with_columns(genome=pl.lit(genome), liftover_method=pl.lit(liftover_method)).fill_nan(0)
            comparative_indel_r2.append(all_stats)

comparative_indel_r2 = pl.concat(comparative_indel_r2).collect().to_pandas()
comparative_indel_r2 = comparative_indel_r2.rename(columns={'ds_r2':'ds_r2_native'}).drop(['rsid_lifted'], axis=1)


# %%
import importlib
from scripts.analysis import whiskers as wsk
importlib.reload(wsk)
import subprocess
import glob
from datetime import datetime
from IPython.display import Image, display

ax2_data = separate_sample_panel_and_subset.loc[ (separate_sample_panel_and_subset.dataset=='SGDP')
                                               & (separate_sample_panel_and_subset.ancestry=='all')
                                               & (separate_sample_panel_and_subset['Variant Subset']=='All panel variants')
                                               & (separate_sample_panel_and_subset.panel_filter.isin(['No Filter', 'no_singletons']))
                                               & (separate_sample_panel_and_subset.genome=='GRCh38')
                                               & (separate_sample_panel_and_subset['Reference Panel']=='CHM13v2.0')].copy()

# The submitted artwork is the lifted T2T-CHM13 panel only: 256 SGDP
# participants for every variant-type/panel-filter group. Fail if a future
# input silently pools the 256 native-GRCh38 observations back in.
supp24_group_counts = (
    ax2_data
    .groupby(['Variant Type', 'panel_filter'], observed=True)
    .agg(n_rows=('sample_name', 'size'), n_participants=('sample_name', 'nunique'))
    .reset_index()
)
if not ((supp24_group_counts.n_rows == 256) & (supp24_group_counts.n_participants == 256)).all():
    raise ValueError(f'Expected 256 lifted-panel SGDP participants per group:\n{supp24_group_counts}')

supp24_source_summary = (
    ax2_data
    .groupby(['Variant Type', 'panel_filter'], observed=True)['imputed_ds_rsquared']
    .agg(
        n='size',
        median='median',
        percentile_5=lambda values: values.quantile(0.05),
        percentile_95=lambda values: values.quantile(0.95),
    )
    .reset_index()
)
supp24_source_summary.to_csv(
    f'{figdir}/supplemental/Supplemental_24_source_summary.csv', index=False
)


fig, (ax1, ax2) = plt.subplots(1,2, figsize=(two_col, 2.5), layout='constrained')
ax1 = sns.swarmplot(x='Variant Type', y='imputed_ds_rsquared', hue='panel_filter',
                    data = ax2_data, dodge=True, s=.4,
                    legend=True, ax=ax1, hue_order=['No Filter', 'no_singletons'])

ax2 = wsk.central_whiskerplot(ax=ax2, data=ax2_data, x="Variant Type", y="imputed_ds_rsquared", hue="panel_filter",
    center="median", line_width_frac=0.95, dodge=True, palette=['black'], linewidth=1.5, caplinewidth=0.5, alpha=0.7, hue_order=['No Filter', 'no_singletons']) #,    whisker="ci",

ax2 = sns.swarmplot(x='Variant Type', y='imputed_ds_rsquared', hue='panel_filter', hue_order=['No Filter', 'no_singletons'],
                    data = ax2_data.loc[ax2_data.imputed_ds_rsquared > 0.925], dodge=True, s=0.8, legend=True, ax=ax2)#cut=0, split=True,
                    # Configure inputs and invoke the plotting script via subprocess


xlim=ax1.get_xlim()
ax1.set_xlim(xlim[0]-.25, xlim[1]+.25)

ax2.set_ylim(0.92, 1.0)
clean_figure(fig)
sns.move_legend(ax1, "lower left", frameon=True)
sns.move_legend(ax2, "lower right", frameon=True)
plt.savefig(f'{figdir}/supplemental/Supplemental_24.png', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_24.pdf', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_24.svg', facecolor='white', dpi=300)


# %% [markdown]
# ### Supplemental Figure 25
# Distribution of r2 values of imputed common indels (>10% MAF) when imputed using a reference panel lifted to GRCh38 coordinates with bcftools +liftover or LiftoverIndel.
# Each hexagonal bin represents the density of variants with a given pair of imputation r2 values (log-scaled color bar). Only bins with 3 or more variants are shown. The dashed diagonal denotes equality between methods. Points above the diagonal indicate higher imputation accuracy with LiftoverIndel. Inset: For each variant, we calculated the difference in imputed r2 value when using LiftoverIndel or Bcftools +liftover. Higher values indicate more accurate imputation when using a LiftoverIndel panel. The distribution of values for variants with an absolute difference over 0.05 are shown.
#
#

# %%
# Report common indels with MAF > 1%; do not downsample so the SVG preserves all plotted points.
s25_maf_threshold = 0.01
s25_methods = ['liftover_plugin', 'LiftoverIndel']

s25_plot_data = comparative_indel_r2.loc[
    (comparative_indel_r2.genome == 'GRCh38') &
    (comparative_indel_r2.maf > s25_maf_threshold) &
    (comparative_indel_r2.liftover_method.isin(s25_methods))
].copy()
s25_plot_counts = s25_plot_data.groupby('liftover_method').size().to_dict()
print(f"Supplemental Figure 25 plotted variants after MAF > {s25_maf_threshold:.0%}: {len(s25_plot_data):,}")

fig, axes = plt.subplots(2, 1, figsize=(two_col*0.36, two_col*0.72), layout='tight') # a over b, beside c in Extended Data 6 (were side by side, two_col by two_col/2)
sns.scatterplot(
    x='ds_r2_native',
    y='ds_r2_lifted',
    hue='maf',
    data=s25_plot_data.loc[s25_plot_data.liftover_method == 'liftover_plugin'],
    alpha=0.01,
    linewidth=0,
    s=2,
    rasterized=True, # every point is still drawn, just as an image - all vector made a 53 MB pdf, the journal limit is 30. text stays editable
    legend=True,
    ax=axes[0],
)

sns.scatterplot(
    x='ds_r2_native',
    y='ds_r2_lifted',
    hue='maf',
    data=s25_plot_data.loc[s25_plot_data.liftover_method == 'LiftoverIndel'],
    alpha=0.01,
    linewidth=0,
    s=2,
    rasterized=True,
    legend=False,
    ax=axes[1],
)

for ax in axes:
    # r² as a character - the mathtext superscript comes out at 4.2 pt, under the journal's 5 pt minimum
    ax.set_xlabel('Indel Imputation r²\nusing GRCh38 Panel')
    ax.set_ylabel('Indel Imputation r²\nusing CHM13 Panel lifted to GRCh38')

leg = axes[0].get_legend()
leg.set_title('MAF')
for h in leg.legend_handles:
    h.set_alpha(1)
    h.set_markersize(3)

axes[0].set_title(
    f"Liftover via bcftools +liftover\n"
    f"n={s25_plot_counts.get('liftover_plugin', 0):,}"
)
axes[1].set_title(
    f"Liftover via LiftoverIndel\n"
    f"n={s25_plot_counts.get('LiftoverIndel', 0):,}" # was .get('liftoverindel') - wrong case, so the title said n=0
)
clean_figure(fig, letter_fontsize=8) # 8 pt panel letters for the journal (this is Extended Data now)
plt.savefig(f'{figdir}/supplemental/Supplemental_25_liftoverindel_common_indel_scatter.png', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_25_liftoverindel_common_indel_scatter.pdf', facecolor='white', dpi=300)
plt.savefig(f'{figdir}/supplemental/Supplemental_25_liftoverindel_common_indel_scatter.svg', facecolor='white', dpi=300)


# %% [markdown]
# ### Supplemental Figure 26
#
# Distribution of r2 values of imputed common indels (>10% MAF) when imputed using a reference panel lifted to GRCh38 coordinates with bcftools +liftover or LiftoverIndel.
#
# Each hexagonal bin represents the density of variants with a given pair of imputation r2 values (log-scaled color bar). Only bins with 3 or more variants are shown. The dashed diagonal denotes equality between methods. Points above the diagonal indicate higher imputation accuracy with LiftoverIndel. Inset: For each variant, we calculated the difference in imputed r2 value when using LiftoverIndel or Bcftools +liftover. Higher values indicate more accurate imputation when using a LiftoverIndel panel. The distribution of values for variants with an absolute difference over 0.05 are shown.

# %%
wide = (
    comparative_indel_r2
      .pivot(index=["rsid", "genome",'maf'],
             columns="liftover_method",
             values=["ds_r2_native", "ds_r2_lifted"])
      .reset_index()
)

wide.columns = [
    "_".join([str(x) for x in col if x not in ("", None)])
    if isinstance(col, tuple) else str(col)
    for col in wide.columns
]


sub = wide.loc[(wide.genome == "GRCh38") & (wide.maf > 0.1)].copy()

x = sub["ds_r2_lifted_liftover_plugin"].to_numpy()
y = sub["ds_r2_lifted_LiftoverIndel"].to_numpy()
d = y - x  # Δ = improvement of LiftoverIndel over plugin

sns.set_theme(context="paper", style="white",
              rc={"axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
                  # "paper" context makes text 9-10 pt; the journal wants 5-7 pt, so the notebook's sizes: 6 pt, titles 7
                  "font.size": 6, "axes.labelsize": 6, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
                  "legend.fontsize": 6, "legend.title_fontsize": 6})

fig, ax = plt.subplots(figsize=(two_col*0.62,two_col*0.72), layout='tight') # c beside a and b (stacked) in Extended Data 6, about twice their height (was two_col*0.62 by two_col*0.52 under them; two_col by two_col as Supplemental 26)

# Background density (cleaner than 800k-point scatter)
hb = ax.hexbin(
    x, y,
    gridsize=100,
    extent=(0, 1, 0, 1),
    bins="log",
    mincnt=3,              # kills speckle; tune 5–30
    cmap="inferno",            # subdued; also good: "cividis", "viridis"
    linewidths=0,
    edgecolors="none",
    rasterized=True,
)
# Add colorbar under the plot (was on the right side)
cbar = plt.colorbar(hb, ax=ax, label='# of variants in hex bin', location='bottom', shrink=0.6, pad=0.12, aspect=30)
cbar.formatter = FuncFormatter(lambda v, pos: f'{v:,.0f}') # plain numbers - the default 10^n labels are mathtext and their exponents came out at 4.2 pt, under the journal's 5 pt minimum

# Identity line
ax.plot([0, 1], [0, 1], ls="--", lw=1.25, color="black", zorder=3)

ax.set(xlim=(0, 1), ylim=(0, 1))
ax.set_aspect("equal", adjustable="box")

# plain r² rather than mathtext: the superscript came out under the journal's 5 pt minimum
ax.set_xlabel("Imputation r² after liftover (bcftools +liftover)")
ax.set_ylabel("Imputation r² after liftover (LiftoverIndel)")
# ax.set_title("Indels (MAF > 0.1), GRCh38 truth set")

# Quantify the wedge
recovered = ((x <= 0.2) & (y >= 0.8)).mean()     # catastrophic plugin → good LiftoverIndel
harmed    = ((x >= 0.8) & (y <= 0.2)).mean()     # opposite (should be small)
win20     = (d >= 0.2).mean()
lose20    = (d <= -0.2).mean()

axins = ax.inset_axes([0.93-0.35-0.0135, 0.4-0.2-0.0135, 0.35, 0.2]) # matplotlib's own inset_axes, in the same box as the axes_grid1 call this replaces (upper right corner at (0.93, 0.4) of the axes, less its 3 pt border pad): the axes_grid1 locator re-rendered the figure during a bbox_inches='tight' save and drew the rasterized hexbin and colour bar shrunk into the bottom-left corner
sns.histplot(d[abs(d) > 0.05],kde=True, ax=axins)
axins.patch.set_alpha(1)
axins.axvline(0, ls="--", lw=.5, color="0.35")
axins.set_xlabel("Δ = r²(LI) − r²(plugin)")


sns.despine(ax=axins)
add_letter_to_ax(ax, 'c', fontsize=8) # panel c of the LiftoverIndel Extended Data figure, under Supplemental 25's a and b

fig.savefig(f'{figdir}/supplemental/Supplemental_26_liftoverindel_common_indel_hexbin.pdf', bbox_inches="tight")
fig.savefig(f'{figdir}/supplemental/Supplemental_26_liftoverindel_common_indel_hexbin.png', bbox_inches="tight")
fig.savefig(f'{figdir}/supplemental/Supplemental_26_liftoverindel_common_indel_hexbin.svg', bbox_inches="tight")


# %% [markdown]
# ### Extended Data: LiftoverIndel (Supplemental 25 beside Supplemental 26)
#
# Supplemental 25 (a over b) on the left of Supplemental 26 (c), with the letters a and c level, the pair centred on one 180 mm page. Needs the Supplemental 25 and Supplemental 26 cells to have run first.

# %%
import re
import subprocess
import fitz

page_width = two_col*72 # in points
pad = plt.rcParams['savefig.pad_inches']*72 # bbox_inches='tight' puts this much white space around the saved Supplemental 26
supp25 = f'{figdir}/supplemental/Supplemental_25_liftoverindel_common_indel_scatter'
supp26 = f'{figdir}/supplemental/Supplemental_26_liftoverindel_common_indel_hexbin'
with fitz.open(f'{supp25}.pdf') as supp25_pdf:
    w25, h25 = supp25_pdf[0].rect.width, supp25_pdf[0].rect.height # saved without bbox_inches='tight', so the whole page
    top_a = min(s['bbox'][1] for b in supp25_pdf[0].get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans'] if s['text'] == 'a')
with fitz.open(f'{supp26}.pdf') as supp26_pdf:
    pdf_size = supp26_pdf[0].rect.width, supp26_pdf[0].rect.height
    top_c = min(s['bbox'][1] for b in supp26_pdf[0].get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans'] if s['text'] == 'c') - pad
svg_size = [float(v) for v in re.search(r'<svg[^>]*width="([\d.]+)pt" height="([\d.]+)pt"', open(f'{supp26}.svg').read()).groups()]
# each save trims to the figure's tight box at that moment, and this figure's layout shifts a few points with every redraw,
# so the saved pdf and svg boxes differ - and get_tightbbox read off the figure here matched neither (it cut off the x label's descenders). Use the larger saved box
w26, h26 = max(pdf_size[0], svg_size[0]) - 2*pad, max(pdf_size[1], svg_size[1]) - 2*pad
gap = 4 # points between b's column and c
x25 = (page_width - (w25 + gap + w26))/2 # a/b and c together, centred on the page
assert x25 >= 0, f'a/b and c are {(w25 + gap + w26)/72*25.4:.1f} mm wide together, more than the page'
y26 = max(top_a - top_c, 0) # letters a and c level: c has no titles over its axes, so with the tops level its letter sat 10 pt above a's

stack_figure_panels([dict(pdf=f'{supp25}.pdf', svg=f'{supp25}.svg', clip=(0, 0, w25, h25), x=x25, y=0),
                     dict(pdf=f'{supp26}.pdf', svg=f'{supp26}.svg', clip=(pad, pad, w26, h26), x=x25 + w25 + gap, y=y26)],
                    page_width, f'{figdir}/Extended Data - LiftoverIndel.pdf', f'{figdir}/Extended Data - LiftoverIndel.svg')
subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', f'{figdir}/Extended Data - LiftoverIndel.pdf', f'{figdir}/Extended Data - LiftoverIndel'])


# %% [markdown]
# ## Table 1: HPRC error rates
#
# Paper refs: manuscript lines 207-216 and Table 1 caption lines 2037-2079.
#
# This table is restricted to HPRC samples with published long-read-based assemblies (`ground_truth_data_source == "HPRC_samples"`). It is generated from `MAF_performance_variants.parquet` and uses `n_gt_errors / n_gt_checked * 100` for genotype discordance.
#
# Output: `tables/table1_genotype_concordance.tsv`.
#

# %%
# Paper cross-reference: Table 1; manuscript lines 207-216 and caption lines 2037-2079.
# Source is restricted to HPRC samples with published long-read-based assemblies.
table1_sum_cols = ['n_gt_errors', 'n_gt_checked', 'n_switch_errors', 'n_checked']
table1_region_order = {'All_regions': 0, 'Nonsyntenic': 1, 'Syntenic': 2}
table1_genome_order = {'CHM13v2.0': 0, 'GRCh38': 1}
table1_type_order = {'SNPs + Indels': 3, 'SNP': 1, 'Indel': 2}

table1_base = (
    pl.scan_parquet(f'{intermediate_data_folder}/MAF_performance_variants.parquet')
    .select([
        'method_of_phasing',
        'ground_truth_data_source',
        'genome',
        'type',
        'Syntenic',
        *table1_sum_cols,
    ])
    .with_columns([
        pl.col('method_of_phasing').cast(pl.String),
        pl.col('ground_truth_data_source').cast(pl.String),
        pl.col('genome').cast(pl.String),
        pl.col('type').cast(pl.String),
    ])
    .filter(
        (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree')
        & (pl.col('ground_truth_data_source') == 'HPRC_samples')
        & pl.col('type').is_in(['SNP', 'Indel'])
    )
    .with_columns(
        pl.when(pl.col('Syntenic').fill_null(False))
        .then(pl.lit('Syntenic'))
        .otherwise(pl.lit('Nonsyntenic'))
        .alias('region')
    )
)

table1_per_region = table1_base.group_by(['region', 'genome', 'type']).agg([
    pl.sum(col).alias(col) for col in table1_sum_cols
])
table1_all_regions = (
    table1_base.with_columns(pl.lit('All_regions').alias('region'))
    .group_by(['region', 'genome', 'type'])
    .agg([pl.sum(col).alias(col) for col in table1_sum_cols])
)
table1_by_type = pl.concat([table1_all_regions, table1_per_region], how='vertical')
table1_combined_types = (
    table1_by_type.group_by(['region', 'genome'])
    .agg([pl.sum(col).alias(col) for col in table1_sum_cols])
    .with_columns(pl.lit('SNPs + Indels').alias('type'))
    .select(['region', 'genome', 'type', *table1_sum_cols])
)

table1_long = (
    pl.concat([table1_by_type, table1_combined_types], how='vertical')
    .with_columns([
        (pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100).alias('gt_error_rate'),
        (100 - (pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100)).alias('gt_accuracy_rate'),
        (pl.col('n_switch_errors') / pl.col('n_checked') * 100).alias('switch_error_rate'),
        pl.col('region').replace_strict(table1_region_order).alias('region_order'),
        pl.col('genome').replace_strict(table1_genome_order).alias('genome_order'),
        pl.col('type').replace_strict(table1_type_order).alias('type_order'),
    ])
    .sort(['region_order', 'genome_order', 'type_order'])
    .select([
        'region',
        'type',
        'genome',
        'n_gt_errors',
        'n_gt_checked',
        'gt_error_rate',
        'gt_accuracy_rate',
        'n_switch_errors',
        'n_checked',
        'switch_error_rate',
    ])
    .collect()
)

table1_display = (
    table1_long
    .with_columns([
        pl.col('switch_error_rate').round(3).alias('SER† (%)'),
        pl.col('gt_error_rate').round(2).alias('Discordance rate‡ (%)'),
    ])
    .select([
        pl.col('region').alias('Genomic Region'),
        pl.col('genome').alias('Reference Genome'),
        pl.col('type').alias('Variant type'),
        'SER† (%)',
        'Discordance rate‡ (%)',
        pl.col('n_switch_errors').alias('# of haplotype switches'),
        pl.col('n_checked').alias('# of heterozygous sites'),
        pl.col('n_gt_errors').alias('# of discordant variant calls'),
        pl.col('n_gt_checked').alias('# of evaluated genotypes'),
    ])
    .to_pandas()
)

table1_path = Path(tables_folder) / 'table1.tsv'
table1_display.to_csv(table1_path, sep='\t', index=False)
print(f'Wrote: {table1_path}')

table1_display

# %%
# Percent drop from GRCh38 -> CHM13v2.0 using table1_display
# Limited to all regions; reported for SNPs, Indels, and Both (SNPs + Indels)

all_regions = table1_display.loc[
    table1_display["Genomic Region"] == "All_regions",
    ["Reference Genome", "Variant type", "SER† (%)", "Discordance rate‡ (%)"],
].copy()

rows = []
for vt in ["SNP", "Indel", "SNPs + Indels"]:
    sub = all_regions.loc[all_regions["Variant type"] == vt]

    grch38_ser = float(sub.loc[sub["Reference Genome"] == "GRCh38", "SER† (%)"].iloc[0])
    chm13_ser = float(sub.loc[sub["Reference Genome"] == "CHM13v2.0", "SER† (%)"].iloc[0])

    grch38_dis = float(
        sub.loc[sub["Reference Genome"] == "GRCh38", "Discordance rate‡ (%)"].iloc[0]
    )
    chm13_dis = float(
        sub.loc[sub["Reference Genome"] == "CHM13v2.0", "Discordance rate‡ (%)"].iloc[0]
    )

    rows.append(
        {
            "Variant type": "Both" if vt == "SNPs + Indels" else vt,
            "GRCh38 SER† (%)": grch38_ser,
            "CHM13v2.0 SER† (%)": chm13_ser,
            "SER % drop": (grch38_ser - chm13_ser) / grch38_ser * 100,
            "GRCh38 Discordance rate‡ (%)": grch38_dis,
            "CHM13v2.0 Discordance rate‡ (%)": chm13_dis,
            "Discordance % drop": (grch38_dis - chm13_dis) / grch38_dis * 100,
        }
    )

percent_drop_all_regions = pd.DataFrame(rows)
percent_drop_all_regions[["SER % drop", "Discordance % drop"]] = (
    percent_drop_all_regions[["SER % drop", "Discordance % drop"]].round(2)
)

percent_drop_all_regions

# %% [markdown]
# ## Table 2: cytobands improved by T2T-CHM13 reference-panel phasing
#
# Paper refs: manuscript lines 343-354 and Table 2 caption lines 2149-2156.
#
# Outputs: `tables/1kgp_variation_HPRC_ground_truth_summarized.tsv`, `tables/1kgp_variation_HPRC_ground_truth.tsv`, `tables/1kgp_variation_parental_haplotypes_ground_truth.tsv`, `tables/decipher_top10_cytoband_summary.tsv`, `tables/decipher_cytoband_improvement_concentration.tsv`, `tables/decipher_cytoband_delta_ser_permutation.tsv`, and DECIPHER CNV overlap tables.
#

# %%
cyto = pd.read_parquet(f"{intermediate_data_folder}/per_cytoband_variant_data.parquet")
cnv = pd.read_parquet(f"{intermediate_data_folder}/per_genome_cnv_regions.parquet")
cnvslop = pd.read_parquet(f"{intermediate_data_folder}/per_genome_cnv_regions_with_slop.parquet")

def concat_overlapping_sydromes(df):
    df['Syndrome'] = df['Syndrome'] + '/'
    df = df.groupby(['n_switch_errors', 'n_checked', 'n_gt_errors', 'n_gt_checked', 'method_of_phasing', 'ground_truth_data_source', 'genome', 'MAC', 'AN', 'start', 'end', 'switch_error_rate']).sum()
    df['Syndrome'] = df['Syndrome'].str[:-1]
    return df.reset_index()

cnv = concat_overlapping_sydromes(cnv)
cnvslop = concat_overlapping_sydromes(cnvslop)

decipher_cnvs = pd.read_csv('resources/decipher_syndromes.txt', sep='\t')
decipher_cnvs.loc[decipher_cnvs.end_grch38 < decipher_cnvs.start_grch38, ['start_grch38','end_grch38']] = decipher_cnvs.loc[decipher_cnvs.end_grch38 < decipher_cnvs.start_grch38, ['end_grch38','start_grch38']].values
decipher_cnvs.loc[decipher_cnvs.end_chm13 < decipher_cnvs.start_chm13, ['start_chm13','end_chm13']] = decipher_cnvs.loc[decipher_cnvs.end_chm13 < decipher_cnvs.start_chm13, ['end_chm13','start_chm13']].values
decipher_band_annotation = decipher_cnvs.drop(columns=[x for x in decipher_cnvs.columns if 'start' in x or 'end' in x])
decipher_band_annotation['cytoband'] = decipher_band_annotation.cytoband.str.split(',')
decipher_band_annotation = decipher_band_annotation.explode('cytoband').reset_index(drop=True)

cyto['interval'] = cyto.chrom.str.replace('chr','') + cyto['interval'].astype(str)
cyto['decipher'] = cyto.interval.isin(decipher_band_annotation.cytoband)

cyto['switch_error_rate'] = cyto.n_switch_errors/cyto.n_checked * 100
cyto['genotype_error_rate'] = cyto.n_gt_errors/cyto.n_gt_checked * 100

cnv['switch_error_rate'] = cnv.n_switch_errors/cnv.n_checked * 100
cnv['genotype_error_rate'] = cnv.n_gt_errors/cnv.n_gt_checked* 100

cnvslop['switch_error_rate'] = cnvslop.n_switch_errors/cnvslop.n_checked * 100
cnvslop['genotype_error_rate'] = cnvslop.n_gt_errors/cnvslop.n_gt_checked * 100

merge_columns=['method_of_phasing','ground_truth_data_source','chrom','interval']

cyto_grch38 = cyto.loc[cyto.genome=='GRCh38']
cyto_chm13 = cyto.loc[cyto.genome=='CHM13v2.0']
cnv_grch38 = cnv.loc[cnv.genome=='GRCh38']
cnv_chm13 = cnv.loc[cnv.genome=='CHM13v2.0']
cnvslop_grch38 = cnvslop.loc[cnvslop.genome=='GRCh38']
cnvslop_chm13 = cnvslop.loc[cnvslop.genome=='CHM13v2.0']

cyto = cyto_grch38.merge(cyto_chm13, on=['method_of_phasing','ground_truth_data_source','chrom','interval','decipher'], suffixes=['_grch38','_chm13'])
cnv = cnv_grch38.merge(cnv_chm13, on=['method_of_phasing','ground_truth_data_source','Syndrome'], how='outer', suffixes=('_grch38', '_chm13'))
cnvslop = cnvslop_grch38.merge(cnvslop_chm13, on=['method_of_phasing','ground_truth_data_source','Syndrome'], how='outer', suffixes=('_grch38', '_chm13'))

cyto['abs_diff_SER'] = cyto.switch_error_rate_grch38.fillna(0) - cyto.switch_error_rate_chm13.fillna(0)
cnv['abs_diff_SER'] = cnv.switch_error_rate_grch38.fillna(0) - cnv.switch_error_rate_chm13.fillna(0)
cnvslop['abs_diff_SER'] = cnvslop.switch_error_rate_grch38.fillna(0) - cnvslop.switch_error_rate_chm13.fillna(0)

cyto['rel_drop_SER'] = (1- (cyto.switch_error_rate_chm13.fillna(0) / cyto.switch_error_rate_grch38.fillna(0))) * 100
cnv['rel_drop_SER'] = (1- (cnv.switch_error_rate_chm13.fillna(0) / cnv.switch_error_rate_grch38.fillna(0))) * 100
cnvslop['rel_drop_SER'] = (1- (cnvslop.switch_error_rate_chm13.fillna(0) / cnvslop.switch_error_rate_grch38.fillna(0))) * 100

cyto['abs_diff_GER'] = cyto.genotype_error_rate_grch38.fillna(0) - cyto.genotype_error_rate_chm13.fillna(0)
cnv['abs_diff_GER'] = cnv.genotype_error_rate_grch38.fillna(0) - cnv.genotype_error_rate_chm13.fillna(0)
cnvslop['abs_diff_GER'] = cnvslop.genotype_error_rate_grch38.fillna(0) - cnvslop.genotype_error_rate_chm13.fillna(0)

cyto['rel_drop_GER'] = (1- (cyto.genotype_error_rate_chm13.fillna(0) / cyto.genotype_error_rate_grch38.fillna(0))) * 100
cnv['rel_drop_GER'] = (1- (cnv.genotype_error_rate_chm13.fillna(0) / cnv.genotype_error_rate_grch38.fillna(0))) * 100
cnvslop['rel_drop_GER'] = (1- (cnvslop.genotype_error_rate_chm13.fillna(0) / cnvslop.genotype_error_rate_grch38.fillna(0))) * 100

# %%
# Paper cross-reference: Table 2; manuscript lines 343-354 and caption lines 2149-2156.
filtered_cyto = cyto.loc[(cyto.method_of_phasing=='1kgp_variation_phased_with_reference_panel')&(cyto.ground_truth_data_source=='HPRC_samples')&(cyto.n_checked_grch38 > 50)&(cyto.n_checked_chm13 > 50),
        ['interval','n_switch_errors_grch38','n_switch_errors_chm13','n_checked_grch38','n_checked_chm13','switch_error_rate_grch38','switch_error_rate_chm13',
         'n_gt_checked_grch38','n_gt_checked_chm13','genotype_error_rate_grch38','genotype_error_rate_chm13',
         'abs_diff_GER','rel_drop_GER','abs_diff_SER','rel_drop_SER', 'decipher']
         ].sort_values('rel_drop_SER',ascending=False).reset_index(drop=True)
filtered_cyto.to_csv(f'{tables_folder}/1kgp_variation_HPRC_ground_truth_summarized.tsv', sep='\t', index=False)

cyto.loc[(cyto.method_of_phasing=='1kgp_variation_phased_with_reference_panel')&(cyto.ground_truth_data_source=='HPRC_samples')
         ].sort_values('rel_drop_SER',ascending=False).reset_index(drop=True).to_csv(f'{tables_folder}/1kgp_variation_HPRC_ground_truth.tsv', sep='\t', index=False)
cyto.loc[(cyto.method_of_phasing=='1kgp_variation_phased_with_reference_panel')&(cyto.ground_truth_data_source=='trios')
         ].sort_values('rel_drop_SER',ascending=False).reset_index(drop=True).to_csv(f'{tables_folder}/1kgp_variation_parental_haplotypes_ground_truth.tsv', sep='\t', index=False)

table2_columns = {
    'interval': 'interval',
    'decipher': 'decipher',
    'n_gt_checked_grch38': 'GRCh38 Num. of alt genotypes',
    'n_checked_grch38': 'GRCh38 Num. of het sites',
    'genotype_error_rate_grch38': 'GRCh38 GT error rate (%)',
    'switch_error_rate_grch38': 'GRCh38 Switch error rate (%)',
    'n_gt_checked_chm13': 'CHM13 Num. of alt genotypes',
    'n_checked_chm13': 'CHM13 Num. of het sites',
    'genotype_error_rate_chm13': 'CHM13 GT error rate (%)',
    'switch_error_rate_chm13': 'CHM13 Switch error rate (%)',
    'rel_drop_SER': 'Relative reduction in SER (%)',
}

table2_display = filtered_cyto[list(table2_columns)].rename(columns=table2_columns).copy()
for col in [
    'GRCh38 GT error rate (%)',
    'GRCh38 Switch error rate (%)',
    'CHM13 GT error rate (%)',
    'CHM13 Switch error rate (%)',
    'Relative reduction in SER (%)',
]:
    table2_display[col] = table2_display[col].map(lambda value: f"{value:.2f}%" if pd.notna(value) else np.nan)

filtered_cyto.loc[filtered_cyto.n_switch_errors_chm13 > 10]
filtered_cyto.loc[filtered_cyto.interval=='15q22.1']

# %%
# Descriptive Table 2 summary: how many of the top cytobands overlap DECIPHER CNV syndromes.
top10_cytobands = (
    filtered_cyto
    .sort_values('rel_drop_SER', ascending=False)
    .head(10)
    .copy()
)
top10_cytobands['rank'] = np.arange(1, len(top10_cytobands) + 1)

decipher_top10_cytoband_summary = pd.DataFrame([{
    'rank_metric': 'rel_drop_SER',
    'n_top_bands': 10,
    'n_decipher_top_bands': int(top10_cytobands['decipher'].sum()),
    'top_decipher_bands': ','.join(top10_cytobands.loc[top10_cytobands['decipher'], 'interval']),
    'top10_bands': ','.join(top10_cytobands['interval']),
}])

decipher_top10_cytoband_summary.to_csv(
    f'{tables_folder}/decipher_top10_cytoband_summary.tsv',
    sep='	',
    index=False,
)
decipher_top10_cytoband_summary

# %%
# Primary Table 2 test: base-pair-weighted concentration of positive CHM13 improvement in DECIPHER cytobands.
IMPROVEMENT_CONCENTRATION_N_PERMUTATIONS = 100_000
IMPROVEMENT_CONCENTRATION_RANDOM_SEED = 42
IMPROVEMENT_CONCENTRATION_BATCH_SIZE = 5_000


def _permuted_decipher_mass_shares(mass, decipher_mask,
                                   n_permutations=IMPROVEMENT_CONCENTRATION_N_PERMUTATIONS,
                                   seed=IMPROVEMENT_CONCENTRATION_RANDOM_SEED,
                                   batch_size=IMPROVEMENT_CONCENTRATION_BATCH_SIZE):
    mass = np.asarray(mass, dtype=float)
    decipher_mask = np.asarray(decipher_mask, dtype=bool)
    n_regions = mass.size
    n_decipher = int(decipher_mask.sum())
    total_mass = float(mass.sum())
    if total_mass <= 0:
        raise ValueError('Improvement concentration requires positive total mass')

    rng = np.random.default_rng(seed)
    shares = []
    for start in range(0, n_permutations, batch_size):
        batch_n = min(batch_size, n_permutations - start)
        selected = np.argpartition(rng.random((batch_n, n_regions)), n_decipher - 1, axis=1)[:, :n_decipher]
        permuted_masks = np.zeros((batch_n, n_regions), dtype=bool)
        permuted_masks[np.arange(batch_n)[:, None], selected] = True
        shares.append((permuted_masks @ mass) / total_mass)

    return np.concatenate(shares)


improvement_concentration_cyto = cyto.loc[
    (cyto.method_of_phasing == '1kgp_variation_phased_with_reference_panel')
    & (cyto.ground_truth_data_source == 'HPRC_samples')
    & (cyto.n_checked_grch38 > 50)
    & (cyto.n_checked_chm13 > 50),
    [
        'interval', 'decipher', 'start_grch38', 'end_grch38', 'start_chm13', 'end_chm13',
        'switch_error_rate_grch38', 'switch_error_rate_chm13', 'rel_drop_SER'
    ],
].sort_values('rel_drop_SER', ascending=False).reset_index(drop=True)

improvement_concentration_cyto['delta_SER'] = (
    improvement_concentration_cyto['switch_error_rate_grch38']
    - improvement_concentration_cyto['switch_error_rate_chm13']
)
improvement_concentration_cyto['mean_region_bp'] = (
    (improvement_concentration_cyto['end_grch38'] - improvement_concentration_cyto['start_grch38'])
    + (improvement_concentration_cyto['end_chm13'] - improvement_concentration_cyto['start_chm13'])
) / 2

positive_improvement_mass = (
    np.maximum(improvement_concentration_cyto['delta_SER'].to_numpy(dtype=float), 0)
    / 100
    * improvement_concentration_cyto['mean_region_bp'].to_numpy(dtype=float)
)
decipher_mask = improvement_concentration_cyto['decipher'].to_numpy(dtype=bool)
observed_share = float(positive_improvement_mass[decipher_mask].sum() / positive_improvement_mass.sum())
permuted_shares = _permuted_decipher_mass_shares(positive_improvement_mass, decipher_mask)
n_greater_equal = int((permuted_shares >= observed_share).sum())

decipher_cytoband_improvement_concentration = pd.DataFrame([{
    'analysis': 'base_pair_weighted_positive_improvement_concentration',
    'question': 'share of total positive CHM13-associated SER improvement in DECIPHER cytobands',
    'delta_SER_definition': 'switch_error_rate_grch38 - switch_error_rate_chm13',
    'delta_SER_transform': 'max(delta_SER, 0)',
    'exposure_method': 'mean_region_bp',
    'exposure_definition': 'mean of GRCh38 and CHM13 cytoband span in bp',
    'n_regions': int(len(improvement_concentration_cyto)),
    'n_decipher_regions': int(decipher_mask.sum()),
    'n_non_decipher_regions': int((~decipher_mask).sum()),
    'decipher_region_share': float(decipher_mask.mean()),
    'decipher_exposure_share': float(
        improvement_concentration_cyto.loc[decipher_mask, 'mean_region_bp'].sum()
        / improvement_concentration_cyto['mean_region_bp'].sum()
    ),
    'n_permutations': IMPROVEMENT_CONCENTRATION_N_PERMUTATIONS,
    'random_seed': IMPROVEMENT_CONCENTRATION_RANDOM_SEED,
    'total_improvement_mass': float(positive_improvement_mass.sum()),
    'decipher_improvement_mass': float(positive_improvement_mass[decipher_mask].sum()),
    'observed_decipher_improvement_share': observed_share,
    'permuted_share_mean': float(permuted_shares.mean()),
    'permuted_share_median': float(np.median(permuted_shares)),
    'permuted_share_95pct': float(np.quantile(permuted_shares, 0.95)),
    'enrichment_ratio_vs_permuted_median': observed_share / float(np.median(permuted_shares)),
    'n_permuted_shares_greater_equal': n_greater_equal,
    'empirical_p_value_greater_equal_plus_one': (n_greater_equal + 1) / (IMPROVEMENT_CONCENTRATION_N_PERMUTATIONS + 1),
}])

decipher_cytoband_improvement_concentration.to_csv(
    f'{tables_folder}/decipher_cytoband_improvement_concentration.tsv',
    sep='	',
    index=False,
)
decipher_cytoband_improvement_concentration

# %%
# Corroborating Table 2 test: assessed-heterozygous-site-weighted mean delta-SER difference.
DELTA_SER_N_PERMUTATIONS = 100_000
DELTA_SER_RANDOM_SEED = 42
DELTA_SER_BATCH_SIZE = 5_000


def _weighted_mean_difference(values, group_mask, weights):
    return float(
        np.average(values[group_mask], weights=weights[group_mask])
        - np.average(values[~group_mask], weights=weights[~group_mask])
    )


def _permuted_weighted_mean_differences(values, n_decipher, weights,
                                         n_permutations=DELTA_SER_N_PERMUTATIONS,
                                         seed=DELTA_SER_RANDOM_SEED,
                                         batch_size=DELTA_SER_BATCH_SIZE):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n_regions = values.size
    rng = np.random.default_rng(seed)
    weighted_values = values * weights
    total_weight = float(weights.sum())
    total_weighted_value = float(weighted_values.sum())

    differences = []
    for start in range(0, n_permutations, batch_size):
        batch_n = min(batch_size, n_permutations - start)
        selected = np.argpartition(rng.random((batch_n, n_regions)), n_decipher - 1, axis=1)[:, :n_decipher]
        permuted_masks = np.zeros((batch_n, n_regions), dtype=bool)
        permuted_masks[np.arange(batch_n)[:, None], selected] = True

        decipher_weights = permuted_masks @ weights
        decipher_weighted_sums = permuted_masks @ weighted_values
        decipher_means = decipher_weighted_sums / decipher_weights
        non_decipher_means = (total_weighted_value - decipher_weighted_sums) / (total_weight - decipher_weights)
        differences.append(decipher_means - non_decipher_means)

    return np.concatenate(differences)


delta_ser_cyto = filtered_cyto.copy()
delta_ser_cyto['delta_SER'] = delta_ser_cyto['switch_error_rate_grch38'] - delta_ser_cyto['switch_error_rate_chm13']
delta_ser_cyto['harmonic_mean_n_checked'] = (
    2
    / ((1 / delta_ser_cyto['n_checked_grch38'].astype(float))
       + (1 / delta_ser_cyto['n_checked_chm13'].astype(float)))
)

values = delta_ser_cyto['delta_SER'].to_numpy(dtype=float)
weights = delta_ser_cyto['harmonic_mean_n_checked'].to_numpy(dtype=float)
decipher_mask = delta_ser_cyto['decipher'].to_numpy(dtype=bool)
n_decipher_regions = int(decipher_mask.sum())
observed_difference = _weighted_mean_difference(values, decipher_mask, weights)
permuted_differences = _permuted_weighted_mean_differences(values, n_decipher_regions, weights)
n_greater_equal = int((permuted_differences >= observed_difference).sum())

decipher_cytoband_delta_ser_permutation = pd.DataFrame([{
    'analysis': 'harmonic_n_checked_weighted_delta_SER_mean_difference',
    'effect_definition': 'weighted_mean(delta_SER_DECIPHER) - weighted_mean(delta_SER_non_DECIPHER)',
    'delta_SER_definition': 'switch_error_rate_grch38 - switch_error_rate_chm13',
    'weight_method': 'harmonic_mean_n_checked',
    'weight_definition': 'harmonic mean of GRCh38 and CHM13 assessed heterozygous sites',
    'n_regions': int(values.size),
    'n_decipher_regions': n_decipher_regions,
    'n_non_decipher_regions': int(values.size) - n_decipher_regions,
    'n_permutations': DELTA_SER_N_PERMUTATIONS,
    'random_seed': DELTA_SER_RANDOM_SEED,
    'mean_delta_SER_decipher': float(np.average(values[decipher_mask], weights=weights[decipher_mask])),
    'mean_delta_SER_non_decipher': float(np.average(values[~decipher_mask], weights=weights[~decipher_mask])),
    'delta_SER_difference': observed_difference,
    'n_permuted_differences_greater_equal': n_greater_equal,
    'empirical_p_value_greater_equal_plus_one': (n_greater_equal + 1) / (DELTA_SER_N_PERMUTATIONS + 1),
    'permuted_difference_mean': float(permuted_differences.mean()),
    'permuted_difference_median': float(np.median(permuted_differences)),
    'permuted_difference_95pct': float(np.quantile(permuted_differences, 0.95)),
}])

decipher_cytoband_delta_ser_permutation.to_csv(
    f'{tables_folder}/decipher_cytoband_delta_ser_permutation.tsv',
    sep='	',
    index=False,
)
decipher_cytoband_delta_ser_permutation

# %%
# Paper cross-reference: Table 2 DECIPHER CNV overlap; manuscript lines 343-354.
cnv_clean = cnv.loc[(cnv.method_of_phasing=='1kgp_variation_phased_with_reference_panel')&(cnv.ground_truth_data_source=='HPRC_samples'),
                            ['Syndrome','n_switch_errors_grch38','n_switch_errors_chm13','n_checked_grch38','n_checked_chm13','switch_error_rate_grch38','switch_error_rate_chm13',
                            'n_gt_checked_grch38','n_gt_checked_chm13','genotype_error_rate_grch38','genotype_error_rate_chm13',
                            'abs_diff_GER','rel_drop_GER','abs_diff_SER','rel_drop_SER']
                    ].sort_values('rel_drop_SER',ascending=False).reset_index(drop=True)
cnv_clean.to_csv(f'{tables_folder}/decipher_cnv_region_performances_rephased_1KGP_variants_HPRC_ground_truth.tsv', sep='\t', index=False)
cnv_clean.sort_values('rel_drop_SER',ascending=False)

# %%
# Paper cross-reference: Table 2 DECIPHER CNV +/- 1 Mb overlap; manuscript lines 343-354.
cnv_slop_clean = cnvslop.loc[(cnvslop.method_of_phasing=='1kgp_variation_phased_with_reference_panel')&(cnvslop.ground_truth_data_source=='HPRC_samples')&(cnvslop.n_checked_grch38 > 50)&(cnvslop.n_checked_chm13 > 50),
                                ['Syndrome','n_switch_errors_grch38','n_switch_errors_chm13','n_checked_grch38','n_checked_chm13','switch_error_rate_grch38','switch_error_rate_chm13',
                                            'n_gt_checked_grch38','n_gt_checked_chm13','genotype_error_rate_grch38','genotype_error_rate_chm13',
                                            'abs_diff_GER','rel_drop_GER', 'abs_diff_SER','rel_drop_SER']
                            ].sort_values('rel_drop_SER',ascending=False).reset_index(drop=True)
cnv_slop_clean.to_csv(f'{tables_folder}/decipher_cnv_region_performances_rephased_1KGP_variants_HPRC_ground_truth_plusminus_1mb.tsv', sep='\t', index=False)
cnv_slop_clean.sort_values('rel_drop_SER',ascending=False)
