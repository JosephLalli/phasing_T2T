"""
Contains functions for loading and processing data for the whole genome phasing paper, as well as some plotting utilities.
"""
import os
from pathlib import Path
from copy import deepcopy
from string import ascii_lowercase as lowercase

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter, MultipleLocator, PercentFormatter
from matplotlib.transforms import ScaledTranslation


## Import functions
def find_repo_root(start=None):
    """Locate the repository root from a notebook subdirectory."""
    start = Path.cwd().resolve() if start is None else Path(start).resolve()
    for candidate in (start, *start.parents):
        if (candidate / "notebooks").exists() and (candidate / "resources").exists():
            return candidate
    raise RuntimeError(f"Could not locate repository root from {start}")


def require_files(folder, required_files, label, RUN_PROFILE):
    folder_path = Path(folder)
    print(f"Checking for required files in {folder_path}")
    print (Path(os.path.join(folder_path, required_files[0])))
    missing = [required_file for required_file in required_files if not Path(os.path.join(folder_path, required_file)).exists()]
    if missing:
        missing_display = "\n  - ".join(missing)
        raise FileNotFoundError(
            f"RUN_PROFILE={RUN_PROFILE!r} expected {label} at {folder_path}, "
            f"but missing required files:\n  - {missing_display}"
        )
    return str(folder_path)


### Plotting utilities

#### Initial constants for plotting:

update_legend_values={'genome':'Genome',
                      'SNPs + Indels': 'All Variants',
                      'type': 'Variant Type',
                      'syntenic':'Genomic Region',
                      'Syntenic':'Shared genomic regions',
                      'Nonsyntenic':'T2T-CHM13\n(Previously unresolved)',
                      'CHM13v2.0':'T2T-CHM13',
                      'GRCh38':'GRCh38',
                      'not_in_STRs':'Not in GIAB STR region',
                      'in_STRs':'In GIAB STR region',
                      'in_segdups':'In Segmental Duplication',
                      'not_in_segdups':'Not in Segmental Duplication',
                      'biallelic':'Isolated, biallelic variant',
                      'multiallelic':'Variant overlaps another variant',
                      'not_in_platinum_STRs':'Not in Platinum Genomes STR region',
                      'in_platinum_STRs':'In Platinum Genomes STR region',
                      'no_singletons':'No',
                      'No Filter':'Yes', 
                      'panel_filter':'Imputed with singletons?',
                      'imputed_ds_rsquared': "Variant Imputation $\mathregular{r^2}$",
                      'rsquared_diff': "Difference in Variant Imputation $\mathregular{r^2}$",
                      'mean_AF': "Minor Allele Frequency",
                      "in_STRs_or_platinum": "In GIAB or Platinum STRs",
                      "not_in_STRs_or_platinum": "Not in GIAB/Platinum STRs",
                      "in_STRs_or_platinum_or_overlap": "In GIAB/Platinum STRs or overlapping variants",
                      "not_in_STRs_or_platinum_or_overlap": "Not in STR/overlap regions",
                      "simple": "Isolated variants",
                      'in_platinum_not_STRs': "Non-GIAB Platinum STR region",
                      'multiallelic_not_STRs_not_platinum': "Non-STR variants that overlap another variant",
                      "region":"Genomic region",
                      "n_gt_checked": "# alt alleles called",
                      "cumulative_ser": "Contribution to\nSwitch error rate (%)",
                      "cumulative_ger": "Contribution to\nGenotype error rate (%)",
                      'rounded_MAF': "MAF bin",
                      'superpopulation':'Superpopulation',
                      'gt_error_rate':'% of Haplotypecaller variant calls\nnot present in genome assembly',
                      'HPRC_samples':'HPRC samples',
                      'HGSVC_samples':'HGSVC samples',
                      'ground_truth_data_source':'Assembly Source'
                    }

### plotting functions

def add_letter_to_ax(ax: mpl.axes, label: str | int, points_offset=(-25, 7),
                         va='bottom', fontsize=10, weight='bold', fontfamily='sans', alphabet=lowercase, **kwargs) -> mpl.axes:
    # Use ScaledTranslation to put the label
    # - at the top left corner (axes fraction (0, 1)),
    # - offset 20 pixels left and 7 pixels up (offset points (-20, +7)),
    # i.e. just outside the axes.
    if type(label) == int:
        label=alphabet[label]
    ax.text(
        0.0, 1.0, label, transform=(
            ax.transAxes + ScaledTranslation(points_offset[0]/72, points_offset[1]/72, ax.get_figure().dpi_scale_trans)),
            va=va, fontsize=fontsize, weight=weight, fontfamily=fontfamily, **kwargs)
    return ax

def add_letter_to_subfig(subfig, label, points_offset=(-25, 7),
                         va='bottom', fontsize=10, weight='bold', fontfamily='sans', alphabet=lowercase, **kwargs):
    # Use ScaledTranslation to put the label
    # - at the top left corner (axes fraction (0, 1)),
    # - offset 20 pixels left and 7 pixels up (offset points (-20, +7)),
    # i.e. just outside the axes.
    if type(label) == int:
        label=alphabet[label]
    subfig.text(
        0.0, 1.0, label, transform=(
            subfig.transSubfigure + ScaledTranslation(points_offset[0]/72, points_offset[1]/72, subfig.get_figure().dpi_scale_trans)),
            va=va, fontsize=fontsize, weight=weight, fontfamily=fontfamily, **kwargs)
    return subfig


def percent_formatter(x, pos, max_decimals=4):
    x = str(x)
    if '.' in x:
        x = '.'.join([x.split('.')[0], x.split('.')[1][:max_decimals]])
    while x[-1] == '0':
        x = x[:-1]
    if x[-1] == '.':
        x = x[:-1]
    return x+'%'

perc_formatter = FuncFormatter(percent_formatter)

def axis_has_legend(ax) -> bool:
    return ax.get_legend() is not None  # or: getattr(ax, "legend_", None) is not None

def figure_has_legend(ax) -> bool:
    return axis_has_legend(ax) or bool(ax.figure.legends)

def add_grid(ax, **kwargs):
    return ax.grid(alpha=0.5, **kwargs)

def highlight_value(df, var, val):
    df['highlight'] = False
    df.loc[df[var] == val, 'highlight'] = True
    return df

def clean_legend(
    ax,
    update_values=update_legend_values,
    drop_labels=('highlight', 'True', 'False'),
    dedupe=True,
    title=None,
    normalize_symbol_markers=True,
    handle_alpha=1.0,
    symbol_markersize=3,
    **legend_kwargs,
):
    """Clean a seaborn/matplotlib legend by dropping junk labels and renaming.

    Also (optionally) normalizes scatter-style legend handles (PathCollection)
    to a consistent size and alpha.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis whose legend should be cleaned.
    update_values : dict
        Mapping from raw legend labels -> cleaned display labels.
    drop_labels : iterable[str]
        Legend labels to remove entirely (e.g. seaborn's size/hue helper entries).
    dedupe : bool
        If True, drop duplicate labels while preserving first occurrence.
    title : str | None
        If provided, force legend title to this value. If None, preserve existing
        legend title (optionally remapped via `update_values`).
    normalize_symbol_markers : bool
        If True, detect PathCollection legend handles and set marker size/alpha.
    symbol_markersize : float
        Default marker size (points) for PathCollection legend entries.
    symbol_alpha : float
        Default alpha for PathCollection legend entries.
    **legend_kwargs
        Passed to `ax.legend(...)`.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return ax

    existing_legend = ax.get_legend()
    existing_title = (
        existing_legend.get_title().get_text() if existing_legend is not None else ""
    )

    # Preserve title by default; only remap if a non-empty title exists.
    if title is None:
        if existing_title and update_values is not None:
            new_title = update_values.get(existing_title, existing_title)
        else:
            new_title = existing_title
    else:
        new_title = title

    drop_set = set(drop_labels) if drop_labels is not None else set()
    new_handles = []
    new_labels = []
    seen = set()
    for handle, label in zip(handles, labels):
        handle.set_alpha(handle_alpha)
        try:
            handle.set_markersize(symbol_markersize)
        except AttributeError:
            pass
        
        if label in drop_set:
            continue

        clean_label = update_values.get(label, label) if update_values is not None else label

        if dedupe:
            if clean_label in seen:
                continue
            seen.add(clean_label)

        new_handles.append(handle)
        new_labels.append(clean_label)

    if new_handles:
        ax.legend(handles=new_handles, labels=new_labels, title=new_title, **legend_kwargs)

    return ax


def add_letter_to_ax(
    ax,
    label,
    points_offset=(-25, 7),
    va='bottom',
    fontsize=10,
    weight='bold',
    fontfamily='sans',
    alphabet=lowercase,
    **kwargs,
):
    """Add a panel letter to a specific Axes (not a SubFigure)."""
    if type(label) == int:
        label = alphabet[label]

    ax.text(
        0.0,
        1.0,
        label,
        transform=(
            ax.transAxes
            + ScaledTranslation(
                points_offset[0] / 72,
                points_offset[1] / 72,
                ax.figure.dpi_scale_trans,
            )
        ),
        va=va,
        fontsize=fontsize,
        weight=weight,
        fontfamily=fontfamily,
        **kwargs,
    )
    return ax


def clean_axis(
    ax,
    update_values=update_legend_values,
    clean_labels=True,
    clean_legend=True,
    legend_drop_labels=('highlight', 'True', 'False'),
    legend_dedupe=True,
    percent_formatter_cols=(
        'Minor Allele Frequency (%)',
        'Switch Error Rate (%)',
        'Genotype Error Rate (%)',
        'Genotyping Error Rate (%)',
        'Switch Error Rate',
        'Genotype Error Rate',
        'Genotyping Error Rate',
        'True Switch Error Rate',
        '% of CHM13v2.0 chromosome\nnonsyntenic with GRCh38',
        '% of T2T-CHM13 chromosome\nnonsyntenic with GRCh38',
        'switch_error_rate',
        'true_switch_error_rate',
        'genotyping_error_rate',
        'maf',
        'MAF',
        'mean_AF',
        'cumulative_ser',
        'cumulative_ger',
        'gt_error_rate',
    ),
    x_formatter=None,
    y_formatter=None,
    **legend_kwargs,
):
    """Clean a plot axis (labels + legend + formatters).

    - Title is intentionally NOT modified (set titles per-panel).
    - Axis labels are replaced via `update_values` when they match keys.
    - `perc_formatter` is applied to axes whose label matches `percent_formatter_cols`.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    update_values : dict
        Mapping from raw/internal labels to cleaned labels.
    clean_labels : bool
        If True, remap `ax.get_xlabel()` and `ax.get_ylabel()` via `update_values`.
    clean_legend : bool
        If True and the axis has a legend, rebuild it with cleaned labels.
    percent_formatter_cols : iterable[str] | None
        Labels (or column names used as labels) that should use `perc_formatter`.
        If None, no automatic percent formatting is applied.
    x_formatter, y_formatter : matplotlib.ticker.Formatter or None
        If provided, explicitly set major formatter for x or y.
    **legend_kwargs
        Passed through to `ax.legend(...)` when cleaning legends.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    raw_xlab = ax.get_xlabel()
    raw_ylab = ax.get_ylabel()
    raw_x_tick_labels = ax.get_xticklabels()
    new_x_tick_labels = []
    if clean_labels and update_values is not None:
        if raw_xlab:
            ax.set_xlabel(update_values.get(raw_xlab, raw_xlab))
        if raw_ylab:
            ax.set_ylabel(update_values.get(raw_ylab, raw_ylab))
        if raw_x_tick_labels:
            for tick in raw_x_tick_labels:
                tick.set_text(update_values.get(tick.get_text(), tick.get_text()))
                new_x_tick_labels.append(tick)
            ax.set_xticklabels(new_x_tick_labels)
    if clean_legend and axis_has_legend(ax):
        _clean_legend_fn = globals().get('clean_legend', None)
        if not callable(_clean_legend_fn):
            raise TypeError(
                "Expected global 'clean_legend' to be callable; got "
                f"{type(_clean_legend_fn)}"
            )
        _clean_legend_fn(
            ax,
            update_values=update_values,
            drop_labels=legend_drop_labels,
            dedupe=legend_dedupe,
            **legend_kwargs,
        )

    # Apply explicit formatters last (after label cleaning)
    if x_formatter is not None:
        ax.xaxis.set_major_formatter(x_formatter)
    if y_formatter is not None:
        ax.yaxis.set_major_formatter(y_formatter)

    # Apply percent formatter only for chosen labels/cols
    if percent_formatter_cols is not None:
        targets = set(percent_formatter_cols)
        cur_xlab = ax.get_xlabel()
        cur_ylab = ax.get_ylabel()

        if x_formatter is None and (raw_xlab in targets or cur_xlab in targets):
            ax.xaxis.set_major_formatter(perc_formatter)
        if y_formatter is None and (raw_ylab in targets or cur_ylab in targets):
            ax.yaxis.set_major_formatter(perc_formatter)

    return ax


def clean_figure(
    fig,
    axes=None,
    update_values=update_legend_values,
    add_letters=True,
    alphabet=lowercase,
    letter_points_offset=(-25, 7),
    **clean_axis_kwargs
):
    """Clean all axes in a figure and optionally add panel letters.

    Axes are lettered in row-major order (top-to-bottom, left-to-right).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axes : list[matplotlib.axes.Axes] or None
        Axes to clean; defaults to all non-colorbar axes in the figure.
    update_values : dict
        Passed to `clean_axis`.
    add_letters : bool
        If True, add letters to each axis.
    alphabet : str
        Sequence of letters (e.g. `string.ascii_lowercase`).
    letter_points_offset : (float, float)
        Offset for letter placement.
    clean_axis_kwargs : dict or None
        Extra kwargs passed to `clean_axis`.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : list[matplotlib.axes.Axes]
    """
    if axes is None:
        axes = [ax for ax in fig.axes if ax.get_label() != '<colorbar>']

    clean_axis_kwargs = clean_axis_kwargs or {}

    for ax in axes:
        clean_axis(ax, update_values=update_values, **clean_axis_kwargs)

    if add_letters:
        positioned = [(ax, ax.get_position()) for ax in axes]
        # Row-first order: top row -> bottom row, left -> right
        positioned.sort(key=lambda x: (-x[1].y0, x[1].x0))

        for i, (ax, _bbox) in enumerate(positioned):
            if i >= len(alphabet):
                raise ValueError(f"Not enough letters in alphabet for {len(positioned)} axes")
            add_letter_to_ax(ax, alphabet[i], points_offset=letter_points_offset)

    return fig, axes

def add_highlights_legend(ax, var, val, sizes=(1.5,3), update_legend_values=update_legend_values):

    def is_variable_title(handle):
        return handle.get_color()=='w'

    new_handles = list()
    new_labels = list()
    handles, labels = ax.get_legend_handles_labels()
    in_higlight_var_section=False
    for handle, label in zip(handles, labels):
        if is_variable_title(handle):
            in_higlight_var_section = label == var

        if label not in ['highlight','True','False']:
            if in_higlight_var_section:
                if label == val:
                    handle.set_linewidth(sizes[1])
                else:
                    handle.set_linewidth(sizes[0])
            new_labels.append(update_legend_values.get(label,label))
            new_handles.append(handle)
        else:
            print(label)
    ax.legend(handles=new_handles, labels=new_labels)
    return ax


def add_x_pos(df, ax):
    all_points = []
    for collection in ax.collections:
        if hasattr(collection, 'get_offsets'):
            offsets = collection.get_offsets()
            all_points.extend([(x, y) for x, y in offsets])
    data_points=pd.DataFrame(np.array(sorted(all_points, key=lambda x: x[1])), columns=['x_data','gt_error_rate'])
    df=df.merge(data_points,on='gt_error_rate')
    return df