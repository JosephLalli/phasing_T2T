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
# # Paper inline statistics (per-variant source notebook)
#
# This notebook computes **per-variant and filter-derived statistics used directly in manuscript text and tables**.
# Each section prints plug-and-play f-string sentences that use the exact manuscript sentence structure with
# computed values inserted. These sentences are the primary deliverable -- copy-paste them into the paper and
# use the surrounding context to ctrl-F to the correct section.
#
# ## Paper source of truth
# The manuscript PDFs in `papers/` are the source of truth for organization and published figure numbering.
# See `notebooks/notebooks_whole_genome/PAPER_CROSS_REFERENCE.md` for the full figure/data-point map.
#
# ## Execution order
# 1. Run the whole-genome summary dataframe creation pipeline.
# 2. Run this notebook first (heavy per-variant calculations).
# 3. Run `notebooks/notebooks_whole_genome/make_plots.ipynb`.
# 4. Run `notebooks/notebooks_whole_genome/calc_figures_for_paper.ipynb` (summary fast-path and consolidated exports).
#
# ## Data paths
# Set `RUN_PROFILE` in the first code cell. `test` uses plain folders (`intermediate_data`, `tables`);
# `whole_genome` uses `_whole_genome` folders for publication-scale inputs and outputs.
#
# ## Organization
# Production manuscript statistics are organized first. Exploratory checks are at the end under a dedicated heading.
# The central statistics reference TSV is exported to `tables/paper_per_variant_inline_statistics_reference.tsv`.

# %%
import os
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import polars as pl

from utilities import find_repo_root, require_files

# Choose "test" for Docker test data or "whole_genome" for publication-scale data.
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

summary_statistics_folder = profile_paths["intermediate_data"]
tables_folder = profile_paths["tables"]

required_summary_files = (
    "MAF_performance_variants.parquet",
    "samples.parquet",
    "variants.parquet",
)

summary_statistics_folder = require_files(
    profile_paths["intermediate_data"],
    required_summary_files,
    "per-variant summary data",
    RUN_PROFILE,
)

Path(tables_folder).mkdir(parents=True, exist_ok=True)

print(f"RUN_PROFILE: {RUN_PROFILE}")
print(f"Working directory: {repo_root}")
print(f"Summary statistics folder: {summary_statistics_folder}")
print(f"Tables folder: {tables_folder}")

# %% [markdown]
# ## 1) Load data

# %%
MAF_performance_variant_df = pl.read_parquet(f'{summary_statistics_folder}/MAF_performance_variants.parquet')
all_samples = pd.read_parquet(f'{summary_statistics_folder}/samples.parquet')
all_variants = pl.read_parquet(f'{summary_statistics_folder}/variants.parquet')

# %% [markdown]
# ## 2) Variant filtering statistics
#
# Reconstructs CHM13 variant filter flags from bcftools query data to compute the proportion
# of variants removed by each filter criterion in T2T-CHM13-unique regions.

# %%
# Build CHM13 variant ID set for filter reconstruction
chm13_ids = pl.scan_parquet(f"{summary_statistics_folder}/MAF_performance_variants.parquet") \
    .select(['variant_id','genome']) \
    .filter(pl.col('genome')=='CHM13v2.0') \
    .select('variant_id') \
    .unique() \
    .collect() \
    .to_series()
chm13_ids = set(chm13_ids)

# Load bcftools query data and compute filter flags
variant_filter_stats = (
    pl.scan_parquet(f"{summary_statistics_folder}/bcftools_query_variant_data.parquet")
    .with_columns(
        MAC_original=pl.min_horizontal(
            pl.col('AN_original_panel') - pl.col('AC_original_panel'),
            pl.col('AC_original_panel'),
        )
    )
    .filter(pl.col('genome') == 'CHM13v2.0')
    .with_columns(
        VQSLOD_filter=(pl.col('VQSLOD') < 0),
        MERR_filter=~(pl.col('MERR') > (602 * 0.05)),
        in_panel=pl.col('ID').is_in(chm13_ids),
        HWE_pop_filter=(
            (pl.col('HWE_EUR') < 1e-10) & (pl.col('HWE_AFR') < 1e-10) &
            (pl.col('HWE_EAS') < 1e-10) & (pl.col('HWE_AMR') < 1e-10) &
            (pl.col('HWE_SAS') < 1e-10)
        ).cast(pl.Boolean),
        MAC_filter=(pl.col('MAC_original') == 0).cast(pl.Boolean),
        AC_filter=(pl.col('AC_original_panel') <= 1),
        f_missing_filter=(pl.col('F_MISSING') > 0.05).cast(pl.Boolean),
        var_len_filter=(
            (pl.col('ID').str.split('_').list[2].str.len_chars().cast(pl.Int64)
             - pl.col('ID').str.split('_').list[3].str.len_chars().cast(pl.Int64)
            ).abs() >= 50 - 1  # minus 1 because bcftools ILEN does not count the anchor base
        ).cast(pl.Boolean),
        alt_star_filter=(pl.col('ALT') == '*').cast(pl.Boolean),
        singleton=(pl.col('MAC_original') == 1).cast(pl.Boolean),
        pass_filter=(pl.col('FILTER') != 'PASS').cast(pl.Boolean),
    )
    .collect()
)

variant_filter_stats = variant_filter_stats.with_columns(
    CHM13_criteria_fail=(
        pl.col('MAC_filter') | pl.col('VQSLOD_filter') |
        pl.col('f_missing_filter') | pl.col('pass_filter') |
        pl.col('HWE_pop_filter') | pl.col('MERR_filter') |
        pl.col('var_len_filter') | pl.col('alt_star_filter')
    ).cast(pl.Boolean)
)

# %% [markdown]
# ## 3) Table 1: Genotype concordance
#
# Computes genotype discordance/concordance and switch error summaries by variant class
# (SNP, Indel, combined), genome (GRCh38 vs CHM13v2.0), syntenic/non-syntenic regions,
# and ground truth source (HPRC/HGSVC/combined).

# %%
best_phased_variants = MAF_performance_variant_df.filter(
    (pl.col('ground_truth_data_source') != 'trios') &
    (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree')
)
catVarType = pl.Enum(['SNP', 'Indel', 'SNPs + Indels'])
best_phased_variants = best_phased_variants.with_columns(type=pl.col('type').cast(catVarType))
best_phased_variants = pl.concat([
    best_phased_variants,
    best_phased_variants.with_columns(type=pl.lit('SNPs + Indels').cast(catVarType)),
])

best_phased_variants_by_syntenic = (
    best_phased_variants
    .group_by(['type', 'Syntenic', 'method_of_phasing', 'ground_truth_data_source', 'genome'])
    .sum()
    .with_columns(
        gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100,
        switch_error_rate=pl.col('n_switch_errors') / pl.col('n_checked') * 100,
    )
    .to_pandas()
)

best_phased_variants_overall = (
    best_phased_variants
    .group_by(['type', 'method_of_phasing', 'ground_truth_data_source', 'genome'])
    .sum()
    .with_columns(
        gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100,
        switch_error_rate=pl.col('n_switch_errors') / pl.col('n_checked') * 100,
    )
    .to_pandas()
)

best_phased_variants_overall.sort_values(['ground_truth_data_source', 'genome', 'type'])

# %%
best_phased_variants_by_syntenic.loc[
    best_phased_variants_by_syntenic.ground_truth_data_source == 'HPRC_samples'
].sort_values(['ground_truth_data_source', 'Syntenic', 'genome', 'type'])

# %% [markdown]
# ## 4) Percent concordance summaries

# %%
overall_subsetted = (
    MAF_performance_variant_df
    .filter(
        (pl.col('ground_truth_data_source') != 'trios') &
        (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree')
    )
    .group_by(['type', 'Syntenic', 'method_of_phasing', 'ground_truth_data_source', 'genome'])
    .sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked'))
)

overall = (
    MAF_performance_variant_df
    .filter(
        (pl.col('ground_truth_data_source') != 'trios') &
        (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree')
    )
    .group_by(['method_of_phasing', 'ground_truth_data_source', 'genome'])
    .sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked'))
)

overall_syntenic = (
    MAF_performance_variant_df
    .filter(
        (pl.col('ground_truth_data_source') != 'trios') &
        (pl.col('Syntenic')) &
        (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree')
    )
    .group_by(['method_of_phasing', 'ground_truth_data_source', 'genome'])
    .sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked'))
)

# %%
overall_syntenic = overall_syntenic.to_pandas()[['Syntenic', 'genome', 'ground_truth_data_source', 'gt_error_rate']].sort_values(['genome', 'Syntenic', 'gt_error_rate'])
overall_syntenic['Non-reference true positive rate'] = (1 - overall_syntenic.gt_error_rate) * 100
overall_syntenic['gt_error_rate'] *= 100
overall_syntenic.loc[overall_syntenic.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples', 'HPRC_HGSVC_all_samples'))]

# %%
overall = overall.to_pandas()[['Syntenic', 'genome', 'ground_truth_data_source', 'gt_error_rate']].sort_values(['genome', 'Syntenic', 'gt_error_rate'])
overall['Non-reference true positive rate'] = (1 - overall.gt_error_rate) * 100
overall['gt_error_rate'] *= 100
overall.loc[overall.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples', 'HPRC_HGSVC_all_samples'))]

# %%
overall_subsetted = overall_subsetted.to_pandas()[['type', 'Syntenic', 'genome', 'ground_truth_data_source', 'gt_error_rate']].sort_values(['genome', 'type', 'Syntenic', 'gt_error_rate'])
overall_subsetted['Non-reference true positive rate'] = (1 - overall_subsetted.gt_error_rate) * 100
overall_subsetted['gt_error_rate'] *= 100
overall_subsetted.loc[overall_subsetted.ground_truth_data_source.isin(('HPRC_samples', 'HGSVC_samples', 'HPRC_HGSVC_all_samples'))]

# %% [markdown]
# ## 5) Singleton and doubleton accuracy

# %%
# Paper cross-reference: manuscript lines 210-214; Figure 2b,d;
# reviewer reply lines 348-354 (singleton SER ~30%).

final_dataset_singletons = MAF_performance_variant_df.filter(
    (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree') &
    (pl.col('ground_truth_data_source') == 'HGSVC_samples_nontrios_only') &
    (pl.col('rounded_MAF').cast(pl.String) == 'singleton')
)

final_dataset_singleton_ser = (
    final_dataset_singletons
    .group_by(['genome', 'type'])
    .sum()
    .with_columns(SER=pl.col('n_switch_errors') / pl.col('n_checked') * 100)
    .sort(['genome', 'type'])
    .to_pandas()
    .reset_index(drop=True)
)

# Genotype-discordance rows for the manuscript text about singleton and doubleton SNP concordance.
# The doubleton values use all HPRC/HGSVC assembly samples; the singleton value uses HGSVC non-trio samples.
first_non_singleton_maf_bin = '0.00021-0.00042'
paper_maf_concordance = (
    MAF_performance_variant_df
    .filter(
        (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree') &
        (pl.col('type') == 'SNP') &
        (
            ((pl.col('ground_truth_data_source') == 'HPRC_HGSVC_all_samples') &
             (pl.col('rounded_MAF').cast(pl.String) == first_non_singleton_maf_bin)) |
            ((pl.col('ground_truth_data_source') == 'HGSVC_samples_nontrios_only') &
             (pl.col('rounded_MAF').cast(pl.String) == 'singleton'))
        )
    )
    .group_by(['ground_truth_data_source', 'rounded_MAF', 'genome', 'type'])
    .sum()
    .with_columns(
        gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100,
        gt_accuracy_rate=100 - (pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100),
    )
    .sort(['ground_truth_data_source', 'rounded_MAF', 'genome', 'type'])
    .to_pandas()
    .reset_index(drop=True)
)

# %% [markdown]
# ## 6) Concordance summary (Table 1 source)

# %%
synEnum = pl.Enum(categories=['All_regions', 'Syntenic', 'Nonsyntenic'])
concordance = all_variants.filter(
    (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree') &
    (pl.col('ground_truth_data_source') == 'HPRC_samples')
).with_columns(
    Syntenic=pl.col('Syntenic').cast(str).replace('true', 'Syntenic').replace('false', 'Nonsyntenic').cast(synEnum)
)

concordance_summary = concordance.clone().with_columns(Syntenic=pl.lit('All_regions').cast(synEnum))
concordance_summary = pl.concat((concordance, concordance_summary)).group_by(['Syntenic', 'genome', 'type']).sum().to_pandas().reset_index()
concordance_summary = pd.concat((
    concordance_summary,
    concordance_summary.copy().assign(type='SNPs + Indels').groupby(['Syntenic', 'genome', 'type'], observed=True).sum(numeric_only=True).reset_index(),
))
concordance_summary['gt_error_rate'] = concordance_summary.n_gt_errors / concordance_summary.n_gt_checked * 100
concordance_summary['gt_accuracy_rate'] = 100 - concordance_summary['gt_error_rate']
concordance_summary['switch_error_rate'] = concordance_summary.n_switch_errors / concordance_summary.n_checked * 100
concordance_summary['Syntenic'] = concordance_summary.Syntenic.replace('false', 'Nonsyntenic').replace('true', 'Syntenic')

pd.set_option("display.precision", 3)
concordance_summary.sort_values(['Syntenic', 'genome', 'type']).set_index(['Syntenic', 'genome', 'type'])[
    ['switch_error_rate', 'gt_error_rate', 'n_gt_errors', 'n_gt_checked', 'n_switch_errors', 'n_checked']
].astype(float)

# %% [markdown]
# ## 7) STR and isolated variant stratification

# %%
# Compute indel discordance within STRs (GIAB + Platinum) and for isolated (biallelic,
# non-STR, non-overlapping) indels, grouped by genome.

indels_only = MAF_performance_variant_df.filter(
    (pl.col('type') == 'Indel') &
    (pl.col('method_of_phasing') == 'phased_with_parents_and_pedigree') &
    (pl.col('ground_truth_data_source') == 'HPRC_samples')
)

# STR indels: in GIAB or Platinum STR regions
str_indels = (
    indels_only
    .filter(pl.col('in_STRs') | pl.col('in_platinum_STRs'))
    .group_by('genome').sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100)
    .to_pandas().set_index('genome')
)

# multiallelic indels: in GIAB, or Platinum STR regions, or overlapping another variant
multiallelic_indels = (
    indels_only
    .filter(pl.col('in_STRs') | pl.col('in_platinum_STRs') | pl.col('multiallelic'))
    .group_by('genome').sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100)
    .to_pandas().set_index('genome')
)


# Isolated indels: not in STRs, not overlapping other variants (biallelic)
isolated_indels = (
    indels_only
    .filter(~pl.col('in_STRs') & ~pl.col('in_platinum_STRs') & ~pl.col('multiallelic'))
    .group_by('genome').sum()
    .with_columns(gt_error_rate=pl.col('n_gt_errors') / pl.col('n_gt_checked') * 100)
    .to_pandas().set_index('genome')
)

chm13_str_indel_disc = str_indels.loc['CHM13v2.0', 'gt_error_rate']
grch38_str_indel_disc = str_indels.loc['GRCh38', 'gt_error_rate']
chm13_multiallelic_indel_disc = multiallelic_indels.loc['CHM13v2.0', 'gt_error_rate']
grch38_multiallelic_indel_disc = multiallelic_indels.loc['GRCh38', 'gt_error_rate']
chm13_isolated_indel_disc = isolated_indels.loc['CHM13v2.0', 'gt_error_rate']
grch38_isolated_indel_disc = isolated_indels.loc['GRCh38', 'gt_error_rate']

print(f'STR indel discordance  -- GRCh38: {grch38_str_indel_disc:.2f}%, CHM13: {chm13_str_indel_disc:.2f}%')
print(f'Multiallelic indel discordance -- GRCh38: {grch38_multiallelic_indel_disc:.2f}%, CHM13: {chm13_multiallelic_indel_disc:.2f}%')
print(f'Isolated indel discordance -- GRCh38: {grch38_isolated_indel_disc:.2f}%, CHM13: {chm13_isolated_indel_disc:.2f}%')


# %% [markdown]
# ## 8) Plug-and-play manuscript text
#
# The sentences below use the exact manuscript sentence structure with computed values inserted.
# Copy-paste these into the paper and use the surrounding context to ctrl-F to the correct section.

# %%
# Helper: extract a value from concordance_summary
def cs(genome, vtype, syntenic, col='gt_error_rate'):
    mask = (
        (concordance_summary.genome == genome) &
        (concordance_summary.type == vtype) &
        (concordance_summary.Syntenic == syntenic)
    )
    return float(concordance_summary.loc[mask, col].iloc[0])

# Helper: extract singleton SER
def sing_ser(genome, vtype):
    mask = (final_dataset_singleton_ser.genome == genome) & (final_dataset_singleton_ser.type == vtype)
    return float(final_dataset_singleton_ser.loc[mask, 'SER'].iloc[0])

# Helper: extract doubleton/singleton discordance from paper_maf_concordance
def maf_disc(gt_source, maf_bin, genome):
    mask = (
        (paper_maf_concordance.ground_truth_data_source == gt_source) &
        (paper_maf_concordance.rounded_MAF == maf_bin) &
        (paper_maf_concordance.genome == genome)
    )
    return float(paper_maf_concordance.loc[mask, 'gt_error_rate'].iloc[0])


# ---- Overall concordance (Table 1) ----
chm13_all_disc = cs('CHM13v2.0', 'SNPs + Indels', 'All_regions')
grch38_all_disc = cs('GRCh38', 'SNPs + Indels', 'All_regions')

print('--- Results: 1kGP short-read variant calls (Table 1 / concordance) ---')
print()
print(f'Overall, T2T-CHM13 1kGP variants were {100 - chm13_all_disc:.2f}% concordant '
      f'with assembly-derived variation, while GRCh38 1kGP variants were '
      f'{100 - grch38_all_disc:.2f}% concordant (Table 1).')
print()

# ---- Shared-region concordance ----
chm13_snp_syn = cs('CHM13v2.0', 'SNP', 'Syntenic')
grch38_snp_syn = cs('GRCh38', 'SNP', 'Syntenic')
chm13_indel_syn = cs('CHM13v2.0', 'Indel', 'Syntenic')
grch38_indel_syn = cs('GRCh38', 'Indel', 'Syntenic')

print(f'SNPs from the T2T-CHM13 panel were less discordant with assembled genomes '
      f'({chm13_snp_syn:.2f}%) than those from the GRCh38 panel ({grch38_snp_syn:.2f}%).')
print()
print(f'Indels in shared regions followed the same pattern '
      f'(T2T-CHM13: {chm13_indel_syn:.2f}%, GRCh38: {grch38_indel_syn:.2f}%, Table 1).')
print()

# ---- Doubleton/singleton discordance ----
grch38_doubleton = maf_disc('HPRC_HGSVC_all_samples', first_non_singleton_maf_bin, 'GRCh38')
chm13_doubleton = maf_disc('HPRC_HGSVC_all_samples', first_non_singleton_maf_bin, 'CHM13v2.0')
chm13_singleton = maf_disc('HGSVC_samples_nontrios_only', 'singleton', 'CHM13v2.0')

print(f'Genotyping accuracy did not vary substantially by minor allele frequency (0.02%-50%), '
      f'consistent with the VQSR method\'s design to maintain a uniform false positive rate '
      f'regardless of allele frequency. SNPs with only two minor alleles present in the panel '
      f'(doubletons) were called almost as accurately as SNPs in the overall panel '
      f'(GRCh38: {grch38_doubleton:.2f}% discordant, T2T-CHM13: {chm13_doubleton:.2f}% discordant), '
      f'and singleton SNPs in the T2T-CHM13 panel were only {chm13_singleton:.2f}% discordant '
      f'with assembly genotypes.')
print()

# ---- STR / isolated indel discordance ----
print('--- Results: STR and isolated variant stratification ---')
print()
print(f'Manual inspection of pangenome-discordant indels showed that they were frequently '
      f'found within multiallelic short tandem repeats (STRs) -- a trend that we confirmed '
      f'with genome-wide analysis (GRCh38: {grch38_multiallelic_indel_disc:.2f}%, T2T-CHM13: '
      f'{chm13_multiallelic_indel_disc:.2f}% indel discordance within STRs, Extended Data Figure 1).')
print()
print(f'In contrast, isolated indels that did not overlap with other variant classes in the '
      f'panel showed substantially lower discordance rates (GRCh38: {grch38_isolated_indel_disc:.2f}%, '
      f'T2T-CHM13: {chm13_isolated_indel_disc:.2f}%, Supplementary Figure 13).')
print()

# ---- Singleton SER ----
print('--- Results: Singleton phasing accuracy ---')
print()
snp_sing_ser = sing_ser('CHM13v2.0', 'SNP')  # manuscript uses CHM13 values
indel_sing_ser = sing_ser('CHM13v2.0', 'Indel')
print(f'Singleton variants were statistically phased with a switch error rate (SER) of ~30% '
      f'({snp_sing_ser:.2f}% for SNPs; {indel_sing_ser:.2f}% for indels) (Figure 2b,d).')
print()

# ---- Summary/derived stats ----
pct_fewer_discordant = (1 - chm13_all_disc / grch38_all_disc) * 100
print('--- Summary / Introduction ---')
print()
print(f'Among HPRC samples with published long-read-based assemblies, alignment to T2T-CHM13 resulted in '
      f'{pct_fewer_discordant:.0f}% fewer assembly-discordant genotypes')

# %% [markdown]
# ## 9) Export: paper_per_variant_inline_statistics_reference.tsv
#
# Consolidated export with one row per unique statistic (both genome values on the same row)
# and a `Paper_Text` column containing the plug-and-play manuscript sentence.

# %%
type_abbrev_map = {'SNP': 'SNP', 'Indel': 'Ind', 'SNPs + Indels': 'All'}
syn_abbrev_map = {'All_regions': 'All', 'Syntenic': 'Syn', 'Nonsyntenic': 'Nsyn'}

rows = []

# --- Table 1 / concordance: one row per (type, syntenic) combination ---
for (vtype, syn), grp in concordance_summary.groupby(['type', 'Syntenic'], observed=True):
    type_abbrev = type_abbrev_map.get(vtype, vtype)
    syn_abbrev = syn_abbrev_map.get(syn, syn)
    stat_code = f'C_{type_abbrev}_{syn_abbrev}'
    by_genome = grp.set_index('genome')
    val_chm13 = float(by_genome.loc['CHM13v2.0', 'gt_error_rate']) if 'CHM13v2.0' in by_genome.index else np.nan
    val_grch38 = float(by_genome.loc['GRCh38', 'gt_error_rate']) if 'GRCh38' in by_genome.index else np.nan
    rows.append({
        'Stat_Code': stat_code,
        'Section': 'Per-variant concordance',
        'Statistic_Description': f'{syn} | {vtype} | gt_error_rate',
        'Value_CHM13': val_chm13,
        'Value_GRCh38': val_grch38,
        'Percent_Change': (val_chm13 - val_grch38) / val_grch38 * 100 if (val_chm13 and val_grch38) else np.nan,
        'Paper_Text': '',
        'Paper_Location': 'Results: Table 1 / concordance paragraphs',
        'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
        'Notes': 'ground_truth_source=HPRC_samples (phased_with_parents_and_pedigree)',
    })

# Attach Paper_Text to the overall concordance row
for r in rows:
    if r['Stat_Code'] == 'C_All_All':
        r['Paper_Text'] = (
            f'Overall, T2T-CHM13 1kGP variants were {100 - r["Value_CHM13"]:.2f}% concordant '
            f'with assembly-derived variation, while GRCh38 1kGP variants were '
            f'{100 - r["Value_GRCh38"]:.2f}% concordant (Table 1).'
        )
    elif r['Stat_Code'] == 'C_SNP_Syn':
        r['Paper_Text'] = (
            f'SNPs from the T2T-CHM13 panel were less discordant with assembled genomes '
            f'({r["Value_CHM13"]:.2f}%) than those from the GRCh38 panel ({r["Value_GRCh38"]:.2f}%).'
        )
    elif r['Stat_Code'] == 'C_Ind_Syn':
        r['Paper_Text'] = (
            f'Indels in shared regions followed the same pattern '
            f'(T2T-CHM13: {r["Value_CHM13"]:.2f}%, GRCh38: {r["Value_GRCh38"]:.2f}%, Table 1).'
        )

# --- Final-dataset singleton SER: one row per type ---
for vtype in ['SNP', 'Indel']:
    type_abbrev = type_abbrev_map.get(vtype, vtype)
    ser_by_genome = final_dataset_singleton_ser.loc[final_dataset_singleton_ser.type == vtype].set_index('genome')
    val_chm13 = float(ser_by_genome.loc['CHM13v2.0', 'SER']) if 'CHM13v2.0' in ser_by_genome.index else np.nan
    val_grch38 = float(ser_by_genome.loc['GRCh38', 'SER']) if 'GRCh38' in ser_by_genome.index else np.nan
    rows.append({
        'Stat_Code': f'SING_SER_{type_abbrev}',
        'Section': 'Singleton SER',
        'Statistic_Description': f'Final-dataset singleton SER | {vtype}',
        'Value_CHM13': val_chm13,
        'Value_GRCh38': val_grch38,
        'Percent_Change': np.nan,
        'Paper_Text': '',
        'Paper_Location': 'Results: statistically phased singleton text; Figure 2b,d',
        'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
        'Notes': 'method_of_phasing=phased_with_parents_and_pedigree; ground_truth_source=HGSVC_samples_nontrios_only; rounded_MAF=singleton',
    })

# Attach Paper_Text to singleton SNP row
for r in rows:
    if r['Stat_Code'] == 'SING_SER_SNP':
        indel_row = next(x for x in rows if x['Stat_Code'] == 'SING_SER_Ind')
        r['Paper_Text'] = (
            f'Singleton variants were statistically phased with a switch error rate (SER) of ~30% '
            f'({r["Value_CHM13"]:.2f}% for SNPs; {indel_row["Value_CHM13"]:.2f}% for indels) (Figure 2b,d).'
        )

# --- Singleton/doubleton genotype-discordance ---
maf_concordance_defs = {
    ('HPRC_HGSVC_all_samples', first_non_singleton_maf_bin, 'SNP'): {
        'Stat_Code': 'GT_DISC_SNP_DOUBLETON',
        'Statistic_Description': 'Doubleton SNP genotype discordance (first non-singleton MAF bin)',
        'Paper_Location': 'Results: genotype concordance paragraph, manuscript lines 215-218',
        'Notes': (
            f'method_of_phasing=phased_with_parents_and_pedigree; '
            f'ground_truth_source=HPRC_HGSVC_all_samples; rounded_MAF={first_non_singleton_maf_bin}; '
            f'concordance=100-discordance'
        ),
    },
    ('HGSVC_samples_nontrios_only', 'singleton', 'SNP'): {
        'Stat_Code': 'GT_DISC_SNP_SINGLETON_FINAL',
        'Statistic_Description': 'Singleton SNP genotype discordance (HGSVC non-trio final dataset)',
        'Paper_Location': 'Results: genotype concordance paragraph, manuscript lines 215-218',
        'Notes': (
            'method_of_phasing=phased_with_parents_and_pedigree; '
            'ground_truth_source=HGSVC_samples_nontrios_only; rounded_MAF=singleton; '
            'concordance=100-discordance'
        ),
    },
}

for key, group in paper_maf_concordance.groupby(['ground_truth_data_source', 'rounded_MAF', 'type'], observed=True):
    definition = maf_concordance_defs.get(tuple(key))
    if definition is None:
        continue
    by_genome = group.set_index('genome')
    val_chm13 = float(by_genome.loc['CHM13v2.0', 'gt_error_rate']) if 'CHM13v2.0' in by_genome.index else np.nan
    val_grch38 = float(by_genome.loc['GRCh38', 'gt_error_rate']) if 'GRCh38' in by_genome.index else np.nan
    rows.append({
        'Stat_Code': definition['Stat_Code'],
        'Section': 'MAF-bin genotype concordance',
        'Statistic_Description': definition['Statistic_Description'],
        'Value_CHM13': val_chm13,
        'Value_GRCh38': val_grch38,
        'Percent_Change': np.nan,
        'Paper_Text': '',
        'Paper_Location': definition['Paper_Location'],
        'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
        'Notes': definition['Notes'],
    })

# Attach Paper_Text to the doubleton row (which also covers singleton in same sentence)
for r in rows:
    if r['Stat_Code'] == 'GT_DISC_SNP_DOUBLETON':
        sing_row = next(x for x in rows if x['Stat_Code'] == 'GT_DISC_SNP_SINGLETON_FINAL')
        r['Paper_Text'] = (
            f'SNPs with only two minor alleles present in the panel (doubletons) were called '
            f'almost as accurately as SNPs in the overall panel (GRCh38: {r["Value_GRCh38"]:.2f}% discordant, '
            f'T2T-CHM13: {r["Value_CHM13"]:.2f}% discordant), and singleton SNPs in the T2T-CHM13 panel '
            f'were only {sing_row["Value_CHM13"]:.2f}% discordant with assembly genotypes.'
        )

# --- STR / isolated indel discordance ---
rows.append({
    'Stat_Code': 'STR_INDEL_DISC',
    'Section': 'STR indel discordance',
    'Statistic_Description': 'Indel discordance within GIAB/Platinum STR regions',
    'Value_CHM13': chm13_str_indel_disc,
    'Value_GRCh38': grch38_str_indel_disc,
    'Percent_Change': np.nan,
    'Paper_Text': (
        f'Manual inspection of pangenome-discordant indels showed that they were frequently '
        f'found within multiallelic short tandem repeats (STRs) -- a trend that we confirmed '
        f'with genome-wide analysis (GRCh38: {grch38_str_indel_disc:.2f}%, T2T-CHM13: '
        f'{chm13_str_indel_disc:.2f}% indel discordance within STRs, Extended Data Figure 1).'
    ),
    'Paper_Location': 'Results: STR stratification; manuscript lines 222-229',
    'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
    'Notes': 'ground_truth_source=HPRC_samples; type=Indel; in_STRs or in_platinum_STRs',
})

rows.append({
    'Stat_Code': 'ISOLATED_INDEL_DISC',
    'Section': 'Isolated indel discordance',
    'Statistic_Description': 'Isolated (non-STR, biallelic) indel discordance',
    'Value_CHM13': chm13_isolated_indel_disc,
    'Value_GRCh38': grch38_isolated_indel_disc,
    'Percent_Change': np.nan,
    'Paper_Text': (
        f'In contrast, isolated indels that did not overlap with other variant classes in the '
        f'panel showed substantially lower discordance rates (GRCh38: {grch38_isolated_indel_disc:.2f}%, '
        f'T2T-CHM13: {chm13_isolated_indel_disc:.2f}%, Supplementary Figure 13).'
    ),
    'Paper_Location': 'Results: isolated variant stratification; manuscript lines 222-229',
    'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
    'Notes': 'ground_truth_source=HPRC_samples; type=Indel; ~in_STRs & ~in_platinum_STRs & ~multiallelic',
})

# --- Summary derived stat ---
rows.append({
    'Stat_Code': 'PCT_FEWER_DISCORDANT',
    'Section': 'Summary derived',
    'Statistic_Description': 'Percent fewer assembly-discordant genotypes (HPRC samples; CHM13 vs GRCh38)',
    'Value_CHM13': chm13_all_disc,
    'Value_GRCh38': grch38_all_disc,
    'Percent_Change': pct_fewer_discordant,
    'Paper_Text': (
        f'Among HPRC samples with published long-read-based assemblies, alignment to T2T-CHM13 resulted in '
        f'{pct_fewer_discordant:.0f}% fewer assembly-discordant genotypes'
    ),
    'Paper_Location': 'Summary; Introduction; Discussion',
    'Source_Notebook': 'calc_per_variant_figures_for_paper.ipynb',
    'Notes': 'Derived from HPRC samples only: (1 - chm13_disc/grch38_disc) * 100; All_regions SNPs+Indels gt_error_rate',
})

# --- Build and export ---
paper_inline_statistics_per_variant = pd.DataFrame(rows)
out_path = f'{tables_folder}/paper_per_variant_inline_statistics_reference.tsv'
paper_inline_statistics_per_variant.to_csv(out_path, sep='\t', index=False)
print(f'Wrote {len(paper_inline_statistics_per_variant)} rows to {out_path}')
paper_inline_statistics_per_variant
