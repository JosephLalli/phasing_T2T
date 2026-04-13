import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def _categorical_order(vec, order=None):
    """Like seaborn: preserve order if provided; else use pandas categorical or sorted uniques."""
    if order is not None:
        return list(order)
    if pd.api.types.is_categorical_dtype(vec):
        return list(vec.cat.categories)
    # seaborn uses categorical_order which preserves original appearance order for object sometimes;
    # here we keep stable unique-in-appearance order:
    return list(pd.Index(vec).dropna().unique())

def _bootstrap_ci(values, center_fn, level=95, n_boot=2000, seed=None):
    rng = np.random.default_rng(seed)
    v = np.asarray(values)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan
    boots = rng.choice(v, size=(n_boot, v.size), replace=True)
    stat = np.apply_along_axis(center_fn, 1, boots)
    alpha = (100 - level) / 2
    lo = np.percentile(stat, alpha)
    hi = np.percentile(stat, 100 - alpha)
    return lo, hi

def infer_hue_geometry_from_ax(ax):
    """
    If you've already drawn a dodged seaborn categorical plot on `ax`,
    this tries to infer per-category hue centers and per-hue width.
    Returns dict with keys: 'centers' (sorted unique x centers),
    'delta' (typical spacing between hue centers within a category).
    """
    xs = []

    # Boxes/patches
    for p in ax.patches:
        if hasattr(p, "get_x") and hasattr(p, "get_width"):
            xs.append(p.get_x() + p.get_width() / 2)

    # Lines (caps/medians etc)
    for ln in ax.lines:
        xdata = ln.get_xdata()
        if xdata is None or len(xdata) == 0:
            continue
        xs.extend(list(np.asarray(xdata)))

    xs = np.asarray(xs, dtype=float)
    xs = xs[np.isfinite(xs)]
    if xs.size == 0:
        return None

    # Round a bit to collapse near-identical values
    xs_r = np.round(xs, 6)
    uniq = np.unique(xs_r)

    # Try to infer typical small spacing (hue dodge) as the smallest nonzero diff
    diffs = np.diff(np.sort(uniq))
    diffs = diffs[diffs > 1e-6]
    delta = np.nan if diffs.size == 0 else np.percentile(diffs, 10)

    return {"centers": np.sort(uniq), "delta": delta}

def central_whiskerplot(
    *,
    data,
    x,
    y,
    hue=None,
    ax=None,
    order=None,
    hue_order=None,
    palette=None,
    dodge=True,
    width=0.8,
    center="median",          # "median" or "mean"
    line_width_frac=0.6,      # horizontal line length relative to per-hue width
    whisker="percentile",     # "percentile" | "ci" | "std"
    percentiles=(5, 95),     # used if whisker="percentile"
    ci_level=95,              # used if whisker="ci"
    n_boot=2000,              # used if whisker="ci"
    std_n=1.0,                # used if whisker="std"
    cap_frac=0.3,            # cap width relative to per-hue width (seaborn-ish)
    linewidth=1,
    caplinewidth=None,
    marker_size=0,            # set >0 if you want a center marker
    center_marker="o",
    zorder=10,
    err_kws=None,             # dict passed to Line2D (e.g., alpha)
    seed=None,
    **kwargs
):
    """
    Draws a seaborn-like categorical plot with a central line (median/mean)
    and whiskers extending from that line to a chosen statistic.
    """
    if ax is None:
        ax = plt.gca()
    plot_kws = dict(kwargs)
    if err_kws is not None:
        plot_kws.update(err_kws)   # err_kws wins

    if caplinewidth is None:
        caplinewidth = linewidth

    df = data.copy() if isinstance(data, pd.DataFrame) else pd.DataFrame(data)

    x_levels = _categorical_order(df[x], order)
    if hue is not None:
        hue_levels = _categorical_order(df[hue], hue_order)
    else:
        hue_levels = [None]

    n_hue = len(hue_levels)
    # per-hue "slot" width when dodging
    if (hue is not None) and dodge:
        w = width / n_hue
        offsets = (np.arange(n_hue) - (n_hue - 1) / 2.0) * w
    else:
        w = width
        offsets = np.array([0.0])

    # palette handling
    if hue is not None:
        colors = sns.color_palette(palette, n_hue) if palette is not None else sns.color_palette(n_colors=n_hue)
        color_map = dict(zip(hue_levels, colors))
    else:
        single = sns.color_palette(palette, 1)[0] if palette is not None else sns.color_palette(n_colors=1)[0]
        color_map = {None: single}

    # center fn
    if center == "median":
        center_fn = np.nanmedian
    elif center == "mean":
        center_fn = np.nanmean
    else:
        raise ValueError("center must be 'median' or 'mean'")

    if whisker == "percentile":
        if len(percentiles) != 2:
            raise ValueError("percentiles must contain exactly two values")
        p_lo, p_hi = np.sort(np.asarray(percentiles, dtype=float))
        if p_lo < 0 or p_hi > 100:
            raise ValueError("percentiles values must be between 0 and 100")
    else:
        p_lo = p_hi = np.nan

    def whisker_bounds(vals):
        v = np.asarray(vals)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return np.nan, np.nan, np.nan

        c = center_fn(v)

        if whisker == "percentile":
            lo, hi = np.percentile(v, [p_lo, p_hi])
        elif whisker == "ci":
            lo, hi = _bootstrap_ci(v, center_fn, level=ci_level, n_boot=n_boot, seed=seed)
        elif whisker == "std":
            m = np.nanmean(v)
            s = np.nanstd(v, ddof=1) if v.size > 1 else np.nan
            lo, hi = m - std_n * s, m + std_n * s
            # center line is still median/mean per `center`
        else:
            raise ValueError("whisker must be 'percentile', 'ci', or 'std'")

        if np.isfinite(lo) and np.isfinite(hi) and lo > hi:
            lo, hi = hi, lo

        return c, lo, hi

    # draw
    for i, xlvl in enumerate(x_levels):
        x0 = float(i)

        if hue is None:
            groups = [(None, df.loc[df[x] == xlvl, y])]
        else:
            groups = []
            for h in hue_levels:
                mask = (df[x] == xlvl) & (df[hue] == h)
                groups.append((h, df.loc[mask, y]))

        for j, (h, series) in enumerate(groups):
            xpos = x0 + offsets[j]
            col = color_map[h]

            c, lo, hi = whisker_bounds(series.values)

            # central horizontal line
            requested_half = (w * line_width_frac) / 2.0
            category_cap = 0.5 - abs(offsets[j])
            half = min(requested_half, category_cap)

            ax.plot([xpos - half, xpos + half], [c, c],
                    color=col, linewidth=linewidth, solid_capstyle="butt", zorder=zorder, **plot_kws)

            # vertical whiskers from center to bounds
            ax.plot([xpos, xpos], [lo, c], color=col, linewidth=caplinewidth, zorder=zorder, **plot_kws)
            ax.plot([xpos, xpos], [c, hi], color=col, linewidth=caplinewidth, zorder=zorder, **plot_kws)
            # caps
            cap_half = (w * cap_frac) / 2.0
            ax.plot([xpos - cap_half, xpos + cap_half], [lo, lo], color=col, linewidth=caplinewidth, zorder=zorder, **plot_kws)
            ax.plot([xpos - cap_half, xpos + cap_half], [hi, hi], color=col, linewidth=caplinewidth, zorder=zorder, **plot_kws)
            # optional center marker
            if marker_size and marker_size > 0:
                ax.plot([xpos], [c], marker=center_marker, markersize=marker_size, color=col, zorder=zorder)

    # categorical ticks/labels
    ax.set_xticks(np.arange(len(x_levels)))
    ax.set_xticklabels(x_levels)

    return ax
