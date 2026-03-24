#!/usr/bin/env python3
"""
stitch_svgs.py - Vertically stitch two SVG files (top + bottom) into one.

Single-pair mode:
    python stitch_svgs.py top.svg bottom.svg output.svg [--spacing -8] [--pdf]

Batch mode (pairs grch38_*.svg with t2t_*.svg by matching suffix):
    python stitch_svgs.py --batch figures/figure6/ [--spacing -8] [--pdf]

Figure 6 block mode (build the pre-colorbar chr15 + chr22 regional block):
    python stitch_svgs.py --grid figures/figure6/ [--pdf]
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

try:
    from fontTools.ttLib import TTFont
except ImportError:  # pragma: no cover - optional runtime dependency
    TTFont = None

from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
DEFAULT_REGIONAL_SPACING = -8.0
DEFAULT_WHOLE_GENOME_SPACING = -4.0
DEFAULT_GRID_GAP = 3.0
BASE_GRID_TOP_Y = 20.0
GRID_LEFT_PAD = 0.0
LABEL_X = 48.0
DEFAULT_FONT_FAMILY = "Arial"
LABEL_FONT_SIZE_PX = 8.0 * (96.0 / 72.0)
TITLE_FONT_SIZE_PX = 6.4
TITLE_TOP_MARGIN_PX = 5.333333333333333
TITLE_GAP_ABOVE_PANEL_PX = 2.2
TITLE_LINE_STEP_PX = 7.2
TITLE_Y_SHIFT_PX = TITLE_LINE_STEP_PX
PANEL_TAG_FONT_SIZE_PX = 10.0
PANEL_TAG_FONT_WEIGHT = "bold"
PANEL_TAG_METRICS_FAMILY = "Arial:style=Bold"
PANEL_TAG_LOCAL_X = -0.68 - (10.0 * (72.0 / 96.0))
# Tag positions are stored in the figure's pt-like user units.
PANEL_TAG_LOCAL_Y = -1.621719 + (10.0 * (72.0 / 96.0))
FIGURE6_BLOCK_LEFT_PAD = 12.0
GRID_BOTTOM_PADDING = 10.72
FIGURE6_BOTTOM_PADDING = 0.9
BASE_COLORBAR_Y = 254.18
GENOME_STYLES = {
    "grch38": ("GRCh38", "#8499B1"),
    "t2t": ("T2T-CHM13", "#593F62"),
    "chm13": ("T2T-CHM13", "#593F62"),
}
GRID_LAYOUT = {
    "width": 495.0 + GRID_LEFT_PAD,
    "block_width": 228.0,
    "left_x": 48.0 + GRID_LEFT_PAD,
    "right_x": 267.0 + GRID_LEFT_PAD,
    "top_y": BASE_GRID_TOP_Y,
}
GRID_ROOT_LABELS = [
    ("T2T-CHM13", "#593F62", 61.0 + GRID_LEFT_PAD, 66.656),
    ("GRCh38", "#8499B1", 61.0 + GRID_LEFT_PAD, 183.296),
]
GRID_TITLE_SPECS = [
    {"lines": ["Prader-Willi/Angelman", "syndromes"]},
    {
        "font_size_px": 6.4,
        "lines": ["22q11", "deletion/duplication", "syndromes"],
    },
    {
        "font_size_px": 6.4,
        "lines": ["22q11.2", "distal deletion", "syndrome"],
    },
]
FIGURE6_LAYOUT = {
    "width": 495.0 + GRID_LEFT_PAD,
}
FIGURE6_PNG_DPI = 300
FIGURE6_EXTRA_PNG_DPIS = (600,)


def _float_attr(element, attr):
    """Return a float attribute from an SVG element, stripping 'pt' or 'px' suffixes."""
    raw = element.get(attr, "")
    if not raw:
        raise ValueError(f"SVG is missing required '{attr}' attribute")
    return float(raw.rstrip("ptpxcmin"))


def _nsmap_for_root(*roots):
    """Merge namespace maps from one or more roots."""
    ns = {}
    for root in roots:
        if root.nsmap:
            ns.update(root.nsmap)
    return ns


def _prefix_ids(root, prefix):
    """Prefix all element ids in-place and update references to keep the merged SVG valid."""
    id_map = {}
    for element in root.iter():
        elem_id = element.get("id")
        if elem_id:
            new_id = f"{prefix}-{elem_id}"
            id_map[elem_id] = new_id
            element.set("id", new_id)

    if not id_map:
        return

    url_pattern = re.compile(r"url\(#([^)]+)\)")

    for element in root.iter():
        for attr_name, attr_value in list(element.attrib.items()):
            updated = url_pattern.sub(
                lambda match: f"url(#{id_map.get(match.group(1), match.group(1))})",
                attr_value,
            )
            if updated.startswith("#"):
                ref_id = updated[1:]
                updated = f"#{id_map.get(ref_id, ref_id)}"
            if updated != attr_value:
                element.set(attr_name, updated)


def _spacing_for_paths(top_svg, bottom_svg, override=None):
    """Use separate default spacing for regional and whole-genome panels."""
    if override is not None:
        return override

    names = f"{top_svg.name.lower()} {bottom_svg.name.lower()}"
    if "whole_genome" in names:
        return DEFAULT_WHOLE_GENOME_SPACING
    return DEFAULT_REGIONAL_SPACING


def _infer_panel_style(svg_path):
    """Infer label text and color from the input filename prefix."""
    stem = svg_path.stem.lower()
    for prefix, style in GENOME_STYLES.items():
        if stem.startswith(prefix):
            return style
    return (svg_path.stem, "#000000")


def _svglite_text_style(font_size_px, font_family=DEFAULT_FONT_FAMILY, font_weight=None):
    """Mirror svglite's CSS text serialization for stitched text."""
    style = [f"font-size: {font_size_px:.2f}px", f'font-family: "{font_family}"']
    if font_weight:
        style.append(f"font-weight: {font_weight}")
    return "; ".join(style) + ";"


@lru_cache(maxsize=None)
def _font_metrics(font_family):
    if TTFont is None:
        return None

    try:
        result = subprocess.run(
            ["fc-match", "--format", "%{file}\n", font_family],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None

    font_path = result.stdout.strip()
    if not font_path:
        return None

    font = TTFont(font_path)
    cmap = font.getBestCmap() or {}
    hmtx = font["hmtx"].metrics
    units_per_em = font["head"].unitsPerEm
    return cmap, hmtx, units_per_em


def _text_length_px(text, font_size_px, font_family=DEFAULT_FONT_FAMILY):
    metrics = _font_metrics(font_family)
    if metrics is None:
        return None

    cmap, hmtx, units_per_em = metrics
    advance_units = 0
    fallback_glyph = ".notdef" if ".notdef" in hmtx else None

    for char in text:
        glyph_name = cmap.get(ord(char), fallback_glyph)
        if glyph_name is None or glyph_name not in hmtx:
            return None
        advance_units += hmtx[glyph_name][0]

    return (advance_units / units_per_em) * font_size_px


def _apply_svglite_text_metrics(
    element,
    text,
    font_size_px,
    font_family=DEFAULT_FONT_FAMILY,
):
    if not text:
        return

    text_length = _text_length_px(text, font_size_px, font_family=font_family)
    if text_length is None:
        return

    element.set("textLength", f"{text_length:.2f}px")
    element.set("lengthAdjust", "spacingAndGlyphs")


def _add_panel_label(root, label, color, y, x=LABEL_X):
    text = etree.SubElement(
        root,
        f"{{{SVG_NS}}}text",
        attrib={
            "x": f"{x:g}",
            "y": f"{y:.3f}",
            "text-anchor": "end",
            # Emit the same CSS-style text metadata svglite uses upstream so the
            # stitched labels share the same font metrics as the panel text.
            "style": _svglite_text_style(LABEL_FONT_SIZE_PX),
            "fill": color,
        },
    )
    text.text = label
    _apply_svglite_text_metrics(text, label, LABEL_FONT_SIZE_PX)


def _add_panel_tag(root, label, x, y):
    text = etree.SubElement(
        root,
        f"{{{SVG_NS}}}text",
        attrib={
            "x": f"{x:.3f}",
            "y": f"{y:.3f}",
            "text-anchor": "start",
            "style": _svglite_text_style(
                PANEL_TAG_FONT_SIZE_PX,
                font_weight=PANEL_TAG_FONT_WEIGHT,
            ),
            "fill": "black",
        },
    )
    text.text = label
    _apply_svglite_text_metrics(
        text,
        label,
        PANEL_TAG_FONT_SIZE_PX,
        font_family=PANEL_TAG_METRICS_FAMILY,
    )
    return text


def _add_multiline_text(
    root,
    lines,
    x,
    y_start,
    font_size_px,
    fill="black",
    line_step=TITLE_LINE_STEP_PX,
):
    text = etree.SubElement(
        root,
        f"{{{SVG_NS}}}text",
        attrib={
            "x": f"{x:.3f}",
            "text-anchor": "middle",
            "style": _svglite_text_style(font_size_px),
            "fill": fill,
        },
    )
    for i, line in enumerate(lines):
        tspan = etree.SubElement(
            text,
            f"{{{SVG_NS}}}tspan",
            attrib={"x": f"{x:.3f}", "y": f"{y_start + i * line_step:.3f}"},
        )
        tspan.text = line
        _apply_svglite_text_metrics(tspan, line, font_size_px)
    return text


def _viewbox_tuple(element):
    viewbox = element.get("viewBox")
    if viewbox:
        parts = viewbox.split()
        if len(parts) != 4:
            raise ValueError(f"Unexpected viewBox on element {element.tag}: {viewbox}")
        return tuple(float(part) for part in parts)
    return (0.0, 0.0, _float_attr(element, "width"), _float_attr(element, "height"))


def _find_svg_by_id_suffix(root, suffix):
    for element in root.findall(f".//{{{SVG_NS}}}svg"):
        element_id = element.get("id", "")
        if element_id.endswith(suffix):
            return element
    raise ValueError(f"Could not find nested SVG with id suffix '{suffix}'")


def _svg_local_x_to_parent_x(svg_element, local_x):
    viewbox_x, _, viewbox_width, _ = _viewbox_tuple(svg_element)
    rendered_width = _float_attr(svg_element, "width")
    parent_x = float(svg_element.get("x", "0"))
    scale = rendered_width / viewbox_width
    return parent_x + (local_x - viewbox_x) * scale


def _svg_local_y_to_parent_y(svg_element, local_y):
    _, viewbox_y, _, viewbox_height = _viewbox_tuple(svg_element)
    rendered_height = _float_attr(svg_element, "height")
    parent_y = float(svg_element.get("y", "0"))
    scale = rendered_height / viewbox_height
    return parent_y + (local_y - viewbox_y) * scale


def _title_text_height(spec):
    line_step = spec.get("line_step", TITLE_LINE_STEP_PX)
    font_size_px = spec.get("font_size_px", TITLE_FONT_SIZE_PX)
    return font_size_px + (len(spec["lines"]) - 1) * line_step


def _title_band_height(specs):
    return (
        TITLE_TOP_MARGIN_PX
        + max(_title_text_height(spec) for spec in specs)
        + TITLE_GAP_ABOVE_PANEL_PX
    )


def _figure6_title_specs_from_grid(grid_root):
    chr15_block = _find_svg_by_id_suffix(grid_root, "chr15-block")
    chr22_block = _find_svg_by_id_suffix(grid_root, "chr22-block")

    chr15_center = _svg_local_x_to_parent_x(chr15_block, _extract_cnv_centers(chr15_block, 1)[0])
    chr22_centers = [
        _svg_local_x_to_parent_x(chr22_block, center)
        for center in _extract_cnv_centers(chr22_block, 2)
    ]

    return [
        {
            "x": chr15_center + GRID_TITLE_SPECS[0].get("x_offset", 0.0),
            "font_size_px": GRID_TITLE_SPECS[0].get("font_size_px", TITLE_FONT_SIZE_PX),
            "line_step": GRID_TITLE_SPECS[0].get("line_step", TITLE_LINE_STEP_PX),
            "lines": GRID_TITLE_SPECS[0]["lines"],
        },
        {
            "x": chr22_centers[0] + GRID_TITLE_SPECS[1].get("x_offset", 0.0),
            "font_size_px": GRID_TITLE_SPECS[1].get("font_size_px", TITLE_FONT_SIZE_PX),
            "line_step": GRID_TITLE_SPECS[1].get("line_step", TITLE_LINE_STEP_PX),
            "lines": GRID_TITLE_SPECS[1]["lines"],
        },
        {
            "x": chr22_centers[1] + GRID_TITLE_SPECS[2].get("x_offset", 0.0),
            "font_size_px": GRID_TITLE_SPECS[2].get("font_size_px", TITLE_FONT_SIZE_PX),
            "line_step": GRID_TITLE_SPECS[2].get("line_step", TITLE_LINE_STEP_PX),
            "lines": GRID_TITLE_SPECS[2]["lines"],
        },
    ]


def _figure6_panel_tag_specs_from_grid(grid_root):
    chr15_block = _find_svg_by_id_suffix(grid_root, "chr15-block")
    chr22_block = _find_svg_by_id_suffix(grid_root, "chr22-block")

    specs = []
    for block, label in (
        (chr15_block, "a"),
        (chr22_block, "b"),
    ):
        x = _svg_local_x_to_parent_x(block, PANEL_TAG_LOCAL_X)
        y = _svg_local_y_to_parent_y(block, PANEL_TAG_LOCAL_Y)
        specs.append({"label": label, "x": x, "y": y})

    return specs


def _figure6_colorbar_x(grid_root, colorbar_root):
    chr15_block = _find_svg_by_id_suffix(grid_root, "chr15-block")
    chr22_block = _find_svg_by_id_suffix(grid_root, "chr22-block")
    left_panel_center = _svg_local_x_to_parent_x(chr15_block, 108.0)
    right_panel_center = _svg_local_x_to_parent_x(chr22_block, 108.0)
    axis_width = right_panel_center - left_panel_center
    return left_panel_center - ((_float_attr(colorbar_root, "width") - axis_width) / 2.0)


def _new_svg_root(width, height, nsmap):
    if None not in nsmap:
        nsmap[None] = SVG_NS

    root = etree.Element(
        f"{{{SVG_NS}}}svg",
        nsmap=nsmap,
        attrib={
            "version": "1.1",
            "width": f"{width}pt",
            "height": f"{height}pt",
            "viewBox": f"0 0 {width} {height}",
        },
    )
    etree.SubElement(
        root,
        f"{{{SVG_NS}}}rect",
        attrib={"width": "100%", "height": "100%", "fill": "white"},
    )
    return root


def _merge_defs(new_root, roots):
    combined_defs = etree.SubElement(new_root, f"{{{SVG_NS}}}defs")
    for root in roots:
        for defs in root.findall(f"{{{SVG_NS}}}defs"):
            for child in defs:
                combined_defs.append(deepcopy(child))
    return combined_defs


def _append_svg_children(target_group, source_root):
    for child in source_root:
        tag = etree.QName(child.tag).localname if child.tag != etree.Comment else None
        if tag == "defs":
            continue
        if _is_background_rect(child):
            continue
        child_copy = deepcopy(child)
        _strip_background_rects(child_copy)
        target_group.append(child_copy)


def _strip_background_rects(element):
    for child in list(element):
        tag = etree.QName(child.tag).localname if child.tag != etree.Comment else None
        if tag == "rect" and _is_background_rect(child):
            element.remove(child)
            continue
        _strip_background_rects(child)


def _is_background_rect(element):
    return (
        etree.QName(element.tag).localname == "rect"
        and element.get("height") == "100%"
        and (_extract_fill_color(element) or "").upper() in {"#FFFFFF", "WHITE"}
    )


def _extract_fill_color(element):
    fill = element.get("fill")
    if fill:
        return fill.strip()
    style = element.get("style", "")
    match = re.search(r"fill:\s*([^;]+)", style)
    if match:
        return match.group(1).strip()
    return None


def _merge_intervals(intervals, tolerance=1.0):
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + tolerance:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _extract_cnv_centers(root, expected_count):
    intervals = []
    for rect in root.findall(f".//{{{SVG_NS}}}rect"):
        fill = _extract_fill_color(rect)
        if fill != "#FFDDDD":
            continue
        x = rect.get("x")
        width = rect.get("width")
        if not x or not width:
            continue
        try:
            start = float(x)
            rect_width = float(width)
        except ValueError:
            continue
        if rect_width < 5:
            continue
        intervals.append((start, start + rect_width))
    merged = _merge_intervals(sorted(intervals))
    if len(merged) < expected_count:
        raise ValueError(
            f"Expected at least {expected_count} CNV intervals, found {len(merged)}"
        )
    merged = sorted(merged, key=lambda item: item[0])
    selected = merged[:expected_count]
    return [((start + end) / 2.0) for start, end in selected]


def _normalize_figure6_panel(root):
    """Fix legacy short-panel SVG base-number text that lands just below the canvas."""
    panel_height = _float_attr(root, "height")
    for text in root.findall(f".//{{{SVG_NS}}}text"):
        text_value = "".join(text.itertext()).strip()
        y = text.get("y")
        if not y:
            continue
        try:
            y_value = float(y)
        except ValueError:
            continue
        if text_value.isdigit() and y_value > panel_height:
            text.set("y", "115.52")


def _normalize_figure6_block(root):
    """Normalize text that falls just below the block canvas after export."""
    block_height = _float_attr(root, "height")
    for text in root.findall(f".//{{{SVG_NS}}}text"):
        text_value = "".join(text.itertext()).strip()
        y = text.get("y")
        if not y:
            continue
        try:
            y_value = float(y)
        except ValueError:
            continue
        if text_value.isdigit() and y_value > block_height:
            text.set("y", f"{block_height - 0.83:.2f}")


def stitch(top_svg: Path, bottom_svg: Path, out_svg: Path, spacing=None):
    """Parse two SVG files and combine them vertically into out_svg."""

    top_tree = etree.parse(str(top_svg))
    bot_tree = etree.parse(str(bottom_svg))

    top_root = top_tree.getroot()
    bot_root = bot_tree.getroot()

    _prefix_ids(top_root, "top")
    _prefix_ids(bot_root, "bottom")

    top_w = _float_attr(top_root, "width")
    top_h = _float_attr(top_root, "height")
    bot_w = _float_attr(bot_root, "width")
    bot_h = _float_attr(bot_root, "height")

    spacing = _spacing_for_paths(top_svg, bottom_svg, spacing)
    combined_w = max(top_w, bot_w)
    combined_h = top_h + spacing + bot_h
    if combined_h <= 0:
        raise ValueError(
            f"Computed non-positive combined height {combined_h} for {top_svg.name} and {bottom_svg.name}"
        )

    nsmap = _nsmap_for_root(top_root, bot_root)
    new_root = _new_svg_root(combined_w, combined_h, nsmap)

    _merge_defs(new_root, (top_root, bot_root))

    # Wrap top content (everything except <defs>) in a <g id="top">
    top_g = etree.SubElement(new_root, f"{{{SVG_NS}}}g", attrib={"id": "top"})
    _append_svg_children(top_g, top_root)

    # Wrap bottom content in <g id="bottom" transform="translate(0, offset)">
    offset = top_h + spacing
    bot_g = etree.SubElement(
        new_root,
        f"{{{SVG_NS}}}g",
        attrib={"id": "bottom", "transform": f"translate(0,{offset})"},
    )
    _append_svg_children(bot_g, bot_root)

    top_label, top_color = _infer_panel_style(top_svg)
    bot_label, bot_color = _infer_panel_style(bottom_svg)
    _add_panel_label(new_root, top_label, top_color, top_h / 2)
    _add_panel_label(new_root, bot_label, bot_color, offset + (bot_h / 2))

    out_tree = etree.ElementTree(new_root)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write(
        str(out_svg),
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=False,
    )
    print(f"  Written: {out_svg}")
    return out_svg


def make_figure6_region_block(grch38_svg: Path, t2t_svg: Path, out_svg: Path):
    """Stack short Figure 6 source panels into one canonical region block."""

    grch38_root = etree.parse(str(grch38_svg)).getroot()
    t2t_root = etree.parse(str(t2t_svg)).getroot()

    _normalize_figure6_panel(grch38_root)
    _normalize_figure6_panel(t2t_root)
    _prefix_ids(grch38_root, "grch38")
    _prefix_ids(t2t_root, "t2t")

    width = _float_attr(grch38_root, "width")
    height = _float_attr(grch38_root, "height")
    t2t_width = _float_attr(t2t_root, "width")
    t2t_height = _float_attr(t2t_root, "height")
    if width != t2t_width or height != t2t_height:
        raise ValueError(
            f"Figure 6 source size mismatch: {grch38_svg.name} is {width}x{height}, "
            f"{t2t_svg.name} is {t2t_width}x{t2t_height}"
        )

    combined_height = height + t2t_height
    root = _new_svg_root(width, combined_height, _nsmap_for_root(grch38_root, t2t_root))
    source_viewbox = grch38_root.get("viewBox", f"0 0 {width} {height}").split()
    if len(source_viewbox) == 4:
        root.set(
            "viewBox",
            f"{source_viewbox[0]} {source_viewbox[1]} {source_viewbox[2]} {combined_height}",
        )
    background = root.find(f"{{{SVG_NS}}}rect")
    if background is not None and FIGURE6_BLOCK_LEFT_PAD:
        background.set("x", f"{-FIGURE6_BLOCK_LEFT_PAD}")
        background.set("width", f"{width}")
    _merge_defs(root, (grch38_root, t2t_root))

    top_group = etree.SubElement(root, f"{{{SVG_NS}}}g", attrib={"id": "t2t"})
    _append_svg_children(top_group, t2t_root)

    bottom_group = etree.SubElement(
        root,
        f"{{{SVG_NS}}}g",
        attrib={
            "id": "grch38",
            "transform": f"translate(0,{height})",
        },
    )
    _append_svg_children(bottom_group, grch38_root)

    out_tree = etree.ElementTree(root)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write(
        str(out_svg),
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=False,
    )
    print(f"  Written: {out_svg}")
    return out_svg


def make_chr15_chr22_grid(chr15_block_svg: Path, chr22_block_svg: Path, out_svg: Path):
    """Create the pre-colorbar Figure 6 chr15 + chr22 block from two assembled region blocks."""

    chr15_root = etree.parse(str(chr15_block_svg)).getroot()
    chr22_root = etree.parse(str(chr22_block_svg)).getroot()

    _prefix_ids(chr15_root, "chr15")
    _prefix_ids(chr22_root, "chr22")

    chr15_w = _float_attr(chr15_root, "width")
    chr15_h = _float_attr(chr15_root, "height")
    chr22_w = _float_attr(chr22_root, "width")
    chr22_h = _float_attr(chr22_root, "height")

    if chr15_w != GRID_LAYOUT["block_width"] or chr22_w != GRID_LAYOUT["block_width"]:
        raise ValueError(
            f"Figure 6 block width mismatch. Expected {GRID_LAYOUT['block_width']}pt, "
            f"found {chr15_w}pt and {chr22_w}pt"
        )
    if chr15_h != chr22_h:
        raise ValueError(f"Figure 6 block height mismatch. Found {chr15_h}pt and {chr22_h}pt")

    root = _new_svg_root(
        GRID_LAYOUT["width"],
        GRID_LAYOUT["top_y"] + chr15_h + GRID_BOTTOM_PADDING,
        _nsmap_for_root(chr15_root, chr22_root),
    )
    _merge_defs(root, (chr15_root, chr22_root))

    block_specs = [
        ("chr15-block", chr15_root, GRID_LAYOUT["left_x"], GRID_LAYOUT["top_y"]),
        ("chr22-block", chr22_root, GRID_LAYOUT["right_x"], GRID_LAYOUT["top_y"]),
    ]
    for block_id, block_root, x, y in block_specs:
        block_viewbox = block_root.get(
            "viewBox",
            f"0 0 {GRID_LAYOUT['block_width']} {chr15_h}",
        )
        block_group = etree.SubElement(
            root,
            f"{{{SVG_NS}}}svg",
            attrib={
                "id": block_id,
                "x": f"{x}",
                "y": f"{y}",
                "width": f"{GRID_LAYOUT['block_width']}",
                "height": f"{chr15_h}",
                "viewBox": block_viewbox,
                "overflow": "visible",
            },
        )
        _append_svg_children(block_group, block_root)

    for label, color, x, y in GRID_ROOT_LABELS:
        _add_panel_label(root, label, color, y, x=x)

    out_tree = etree.ElementTree(root)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write(
        str(out_svg),
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=False,
    )
    print(f"  Written: {out_svg}")
    return out_svg


def assemble_figure6(grid_svg: Path, colorbar_svg: Path, out_svg: Path):
    """Assemble the final Figure 6 SVG from the chr15/chr22 grid and colorbar."""

    grid_root = etree.parse(str(grid_svg)).getroot()
    colorbar_root = etree.parse(str(colorbar_svg)).getroot()

    _prefix_ids(grid_root, "grid")
    _prefix_ids(colorbar_root, "colorbar")

    title_specs = _figure6_title_specs_from_grid(grid_root)
    panel_tag_specs = _figure6_panel_tag_specs_from_grid(grid_root)
    colorbar_x = _figure6_colorbar_x(grid_root, colorbar_root)
    content_height = max(
        _float_attr(grid_root, "height"),
        BASE_COLORBAR_Y + _float_attr(colorbar_root, "height") + FIGURE6_BOTTOM_PADDING,
    )

    root = _new_svg_root(
        FIGURE6_LAYOUT["width"],
        content_height,
        _nsmap_for_root(grid_root, colorbar_root),
    )
    _merge_defs(root, (grid_root, colorbar_root))

    grid_group = etree.SubElement(root, f"{{{SVG_NS}}}g", attrib={"id": "figure6-grid"})
    _append_svg_children(grid_group, grid_root)

    colorbar_group = etree.SubElement(
        root,
        f"{{{SVG_NS}}}g",
        attrib={
            "id": "figure6-colorbar",
            "transform": f"translate({colorbar_x},{BASE_COLORBAR_Y})",
        },
    )
    _append_svg_children(colorbar_group, colorbar_root)

    for spec in panel_tag_specs:
        _add_panel_tag(root, spec["label"], spec["x"], spec["y"])

    panel_top_y = BASE_GRID_TOP_Y
    bottom_baseline_y = panel_top_y - TITLE_GAP_ABOVE_PANEL_PX
    for spec in title_specs:
        line_step = spec.get("line_step", TITLE_LINE_STEP_PX)
        y_start = (
            bottom_baseline_y
            - (len(spec["lines"]) - 1) * line_step
            + TITLE_Y_SHIFT_PX
        )
        _add_multiline_text(
            root,
            spec["lines"],
            spec["x"],
            y_start,
            spec["font_size_px"],
            line_step=line_step,
        )

    out_tree = etree.ElementTree(root)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write(
        str(out_svg),
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=False,
    )
    print(f"  Written: {out_svg}")
    return out_svg


def to_pdf(svg_path: Path):
    """Convert an SVG to PDF using cairosvg if available, else print rsvg-convert hint."""
    pdf_path = svg_path.with_suffix(".pdf")
    try:
        import cairosvg  # type: ignore

        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
        print(f"  PDF:     {pdf_path}")
    except ImportError:
        print(
            f"  cairosvg not found — to convert to PDF run:\n"
            f"    rsvg-convert -f pdf -o {pdf_path} {svg_path}"
        )


def to_png(svg_path: Path, png_path: Path, dpi: int = FIGURE6_PNG_DPI):
    """Convert an SVG to a PNG with explicit DPI metadata."""
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_handle:
            tmp_path = Path(tmp_handle.name)

        try:
            import cairosvg  # type: ignore

            cairosvg.svg2png(url=str(svg_path), write_to=str(tmp_path), dpi=dpi)
        except ImportError:
            rsvg_convert = shutil.which("rsvg-convert")
            if rsvg_convert is None:
                raise RuntimeError(
                    "Neither cairosvg nor rsvg-convert is available for PNG export"
                )
            subprocess.run(
                [
                    rsvg_convert,
                    "-d",
                    str(dpi),
                    "-p",
                    str(dpi),
                    "-o",
                    str(tmp_path),
                    str(svg_path),
                ],
                check=True,
            )

        try:
            from PIL import Image
        except ImportError as exc:
            raise RuntimeError("Pillow is required to tag the PNG with 300 dpi metadata") from exc

        png_path.parent.mkdir(parents=True, exist_ok=True)
        with Image.open(tmp_path) as image:
            image.save(png_path, dpi=(dpi, dpi))
        print(f"  PNG:     {png_path} ({dpi} dpi)")
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def run_single(top: Path, bottom: Path, output: Path, spacing: float, pdf: bool):
    stitch(top, bottom, output, spacing)
    if pdf:
        to_pdf(output)


def run_batch(directory: Path, spacing: float, pdf: bool):
    grch38_files = sorted(directory.glob("grch38_*.svg"))
    if not grch38_files:
        print(f"No grch38_*.svg files found in {directory}", file=sys.stderr)
        sys.exit(1)

    for grch38_svg in grch38_files:
        suffix = grch38_svg.name[len("grch38_"):]  # everything after "grch38_"
        t2t_svg = directory / f"t2t_{suffix}"
        if not t2t_svg.exists():
            print(f"  Warning: no matching t2t_{suffix}, skipping {grch38_svg.name}")
            continue

        stem = Path(suffix).stem  # strip .svg
        out_svg = directory / f"combined_{stem}.svg"
        print(f"Stitching {grch38_svg.name} + {t2t_svg.name} -> {out_svg.name}")
        stitch(grch38_svg, t2t_svg, out_svg, spacing)
        if pdf:
            to_pdf(out_svg)

    run_grid(directory, spacing, pdf)


def run_grid(directory: Path, spacing: float, pdf: bool):
    source_paths = {
        "grch38_chr15": directory / "figure6_grch38_chr15_anglemans.svg",
        "t2t_chr15": directory / "figure6_t2t_chr15_anglemans.svg",
        "grch38_chr22": directory / "figure6_grch38_chr22_22q11.svg",
        "t2t_chr22": directory / "figure6_t2t_chr22_22q11.svg",
    }
    missing = [path.name for path in source_paths.values() if not path.exists()]
    if missing:
        print(
            "Missing Figure 6 source SVGs: " + ", ".join(sorted(missing)),
            file=sys.stderr,
        )
        sys.exit(1)

    chr15_block = directory / "figure6_block_chr15_anglemans.svg"
    print(
        "Overlaying "
        f"{source_paths['grch38_chr15'].name} + {source_paths['t2t_chr15'].name} "
        f"-> {chr15_block.name}"
    )
    make_figure6_region_block(
        source_paths["grch38_chr15"],
        source_paths["t2t_chr15"],
        chr15_block,
    )

    chr22_block = directory / "figure6_block_chr22_22q11.svg"
    print(
        "Overlaying "
        f"{source_paths['grch38_chr22'].name} + {source_paths['t2t_chr22'].name} "
        f"-> {chr22_block.name}"
    )
    make_figure6_region_block(
        source_paths["grch38_chr22"],
        source_paths["t2t_chr22"],
        chr22_block,
    )

    out_svg = directory / "grid_chr15_chr22.svg"
    print(f"Placing {chr15_block.name} + {chr22_block.name} -> {out_svg.name}")
    make_chr15_chr22_grid(chr15_block, chr22_block, out_svg)

    colorbar_svg = directory / "colorbar.svg"
    if colorbar_svg.exists():
        figure6_svg = directory / "figure6.svg"
        print(f"Assembling {out_svg.name} + {colorbar_svg.name} -> {figure6_svg.name}")
        assemble_figure6(out_svg, colorbar_svg, figure6_svg)

        figure6_png = directory / "figure6.png"
        print(f"Rasterizing {figure6_svg.name} -> {figure6_png.name} ({FIGURE6_PNG_DPI} dpi)")
        to_png(figure6_svg, figure6_png, dpi=FIGURE6_PNG_DPI)

        figure6_300dpi_png = directory / f"figure6_{FIGURE6_PNG_DPI}dpi.png"
        shutil.copy2(figure6_png, figure6_300dpi_png)
        print(f"  PNG:     {figure6_300dpi_png} ({FIGURE6_PNG_DPI} dpi)")

        for dpi in FIGURE6_EXTRA_PNG_DPIS:
            dpi_png = directory / f"figure6_{dpi}dpi.png"
            print(f"Rasterizing {figure6_svg.name} -> {dpi_png.name} ({dpi} dpi)")
            to_png(figure6_svg, dpi_png, dpi=dpi)

    if pdf:
        to_pdf(out_svg)


def main():
    parser = argparse.ArgumentParser(
        description="Stitch regional SVGs and rebuild the pre-colorbar Figure 6 block."
    )
    parser.add_argument("--batch", metavar="DIR", type=Path,
                        help="Batch mode: process all grch38_*/t2t_* pairs in DIR")
    parser.add_argument("--grid", metavar="DIR", type=Path,
                        help="Create the pre-colorbar Figure 6 chr15 + chr22 regional block in DIR")
    parser.add_argument("top",    nargs="?", type=Path, help="Top SVG (single-pair mode)")
    parser.add_argument("bottom", nargs="?", type=Path, help="Bottom SVG (single-pair mode)")
    parser.add_argument("output", nargs="?", type=Path, help="Output SVG (single-pair mode)")
    parser.add_argument(
        "--spacing",
        type=float,
        default=None,
        help=(
            "Vertical gap between the two halves in SVG points "
            f"(default: {DEFAULT_REGIONAL_SPACING:g} regional, "
            f"{DEFAULT_WHOLE_GENOME_SPACING:g} whole-genome)"
        ),
    )
    parser.add_argument("--pdf", action="store_true",
                        help="Also export combined PDF via cairosvg (or print rsvg-convert hint)")

    args = parser.parse_args()

    if args.grid:
        if not args.grid.is_dir():
            parser.error(f"--grid path is not a directory: {args.grid}")
        run_grid(args.grid, args.spacing, args.pdf)
    elif args.batch:
        if not args.batch.is_dir():
            parser.error(f"--batch path is not a directory: {args.batch}")
        run_batch(args.batch, args.spacing, args.pdf)
    else:
        if not (args.top and args.bottom and args.output):
            parser.error("Provide top.svg bottom.svg output.svg, or use --batch DIR")
        run_single(args.top, args.bottom, args.output, args.spacing, args.pdf)


if __name__ == "__main__":
    main()
