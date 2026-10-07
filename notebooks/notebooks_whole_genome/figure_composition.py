"""Compose saved vector figures while preserving their PDF and SVG geometry."""
import os
from string import ascii_lowercase as lowercase

import numpy as np


# Original Figure 4 geometry maps to current Figure 3; coordinates are points from the figure top-left, measured from
# the Word-embedded SVG: axis a's top-left corner and the illustration's bounding box.
ORIGINAL_FIG4_AXIS_A_TOP_LEFT = (36.4993, 19.1315)
ORIGINAL_FIG4_ILLUSTRATION_BOX = (45.9, 4.2, 292.4, 173.5)
FIG4_COMPOSITE_GROUP_ID = '4a-4'


def lay_out_as_in_pdf(fig):
    """Lay the figure out with the PDF renderer, leaving its axes where they are in a PDF saved from it.

    Constrained layout places axes up to 0.3 pt differently for each renderer, and plt.savefig and the notebook's
    inline display both redraw the figure for the screen after the last save, so positions read from a figure are
    those of whichever renderer drew it last. A throwaway PDF save re-runs the layout with the PDF renderer; for
    Figure 4 the positions it gives do not depend on which renderer drew the figure before it.
    """
    import io
    fig.savefig(io.BytesIO(), format='pdf')


def axis_top_left_points(fig, ax):
    """Axis top-left corner in points from the figure's top-left corner (SVG and PDF page coordinates), as in the PDF."""
    lay_out_as_in_pdf(fig)
    pos = ax.get_position()
    width_pt, height_pt = fig.get_size_inches() * 72
    return pos.x0 * width_pt, (1 - pos.y1) * height_pt


def axis_vertical_centre_points(fig, ax):
    """Vertical centre of an axes box in points from the figure's top edge, as in the PDF."""
    lay_out_as_in_pdf(fig)
    pos = ax.get_position()
    return (1 - (pos.y0 + pos.y1) / 2) * fig.get_size_inches()[1] * 72


def axis_left_points(fig, ax, origin_x_in=0.0):
    """x of an axes' left edge in points, as in the PDF, measured from a saved page whose left edge is at origin_x_in (inches)."""
    lay_out_as_in_pdf(fig)
    return (ax.get_position().x0 * fig.get_size_inches()[0] - origin_x_in) * 72


def ink_bbox(page, pad=0.0, zoom=4):
    """Bounding box (points) of drawn content on a PDF page, from a transparent rendering."""
    import fitz
    pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=True)
    alpha = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[..., -1]
    ys, xs = np.nonzero(alpha)
    return fitz.Rect(xs.min() / zoom - pad, ys.min() / zoom - pad, (xs.max() + 1) / zoom + pad, (ys.max() + 1) / zoom + pad)


def _prefixed(text, prefix):
    """Prefix ids, id references and CSS classes so inserted artwork cannot collide with matplotlib's."""
    import re
    for i in sorted(set(re.findall(r'\bid="([^"]+)"', text)), key=len, reverse=True):
        text = re.sub(rf'\bid="{re.escape(i)}"', f'id="{prefix}{i}"', text)
        text = text.replace(f'#{i}"', f'#{prefix}{i}"').replace(f'#{i})', f'#{prefix}{i})')
    return re.sub(r'\bcls-(\d+)', rf'{prefix}cls-\1', text)


def register_svg_namespaces():
    from xml.etree import ElementTree as ET
    for prefix, uri in [('', 'http://www.w3.org/2000/svg'), ('xlink', 'http://www.w3.org/1999/xlink'),
                        ('dc', 'http://purl.org/dc/elements/1.1/'), ('cc', 'http://creativecommons.org/ns#'),
                        ('rdf', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#')]:
        ET.register_namespace(prefix, uri)


def panel_a_as_standalone_svg(source_svg, kind, workdir):
    """Write the panel-a artwork as an SVG whose user units are points; return its path."""
    from xml.etree import ElementTree as ET
    register_svg_namespaces()
    root = ET.parse(source_svg).getroot()
    if kind == 'composite':
        # Keep only the illustration group of the submitted composite
        for child in list(root):
            if child.tag.split('}')[1] != 'defs' and child.attrib.get('id') != FIG4_COMPOSITE_GROUP_ID:
                root.remove(child)
    vb = root.attrib['viewBox'].split()
    root.attrib['width'] = vb[2] + 'pt'
    root.attrib['height'] = vb[3] + 'pt'
    out = os.path.join(workdir, f'panel_a_{kind}.svg')
    ET.ElementTree(root).write(out, xml_declaration=True, encoding='utf-8')
    return out


def place_panel_a(svg_path, pdf_path, source_svg, kind, axis_top_left, centre_y=None):
    """Fit the panel-a artwork into the region it occupies in the submitted composite, in both SVG and PDF.

    The artwork is scaled to fit that region and centred in it horizontally. Vertically it is centred in the region,
    or, with centre_y (points from the page top), its drawn content is centred on that height (panel b's axes box).
    """
    import subprocess
    import tempfile
    from xml.etree import ElementTree as ET
    import fitz
    dx = axis_top_left[0] - ORIGINAL_FIG4_AXIS_A_TOP_LEFT[0]
    dy = axis_top_left[1] - ORIGINAL_FIG4_AXIS_A_TOP_LEFT[1]
    box = fitz.Rect(ORIGINAL_FIG4_ILLUSTRATION_BOX) + (dx, dy, dx, dy)
    with tempfile.TemporaryDirectory() as tmp:
        lone_svg = panel_a_as_standalone_svg(source_svg, kind, tmp)
        lone_pdf = os.path.join(tmp, 'panel_a.pdf')
        subprocess.run(['rsvg-convert', '-f', 'pdf', '-o', lone_pdf, lone_svg], check=True)
        art = fitz.open(lone_pdf)
        clip = ink_bbox(art[0], pad=0.5)
        scale = min(box.width / clip.width, box.height / clip.height)
        target = fitz.Rect(box.x0, box.y0, box.x0 + clip.width * scale, box.y0 + clip.height * scale)
        shift_x = (box.width - target.width) / 2
        shift_y = (box.height - target.height) / 2 if centre_y is None else centre_y - (box.y0 + target.height / 2)
        target += (shift_x, shift_y, shift_x, shift_y)

        # PDF: place the artwork as a vector form object (text stays text)
        doc = fitz.open(pdf_path)
        doc[0].show_pdf_page(target, art, 0, clip=clip, keep_proportion=True, overlay=True)
        doc.save(pdf_path + '.tmp', garbage=3, deflate=True)
        doc.close()
        os.replace(pdf_path + '.tmp', pdf_path)

        # SVG: nest the artwork with a viewBox equal to the same clip, at the same target rectangle
        lone_text = _prefixed(open(lone_svg, encoding='utf-8').read(), 'panelA_')
        lone_root = ET.fromstring(lone_text.split('?>', 1)[-1])
        lone_root.attrib.update({'x': f'{target.x0:.3f}', 'y': f'{target.y0:.3f}', 'width': f'{target.width:.3f}',
                                 'height': f'{target.height:.3f}',
                                 'viewBox': f'{clip.x0:.3f} {clip.y0:.3f} {clip.width:.3f} {clip.height:.3f}',
                                 'preserveAspectRatio': 'xMidYMid meet', 'id': 'panel_a_artwork'})
        tree = ET.parse(svg_path)
        tree.getroot().append(lone_root)
        tree.write(svg_path, xml_declaration=True, encoding='utf-8')
    return scale


def stack_figure_panels(panels, page_width, out_pdf, out_svg, gap=8):
    """Stack saved figures top to bottom on one page, unscaled and `gap` points apart, in PDF and SVG.

    Each panel is a dict with its saved 'pdf' and 'svg', 'clip' = (x0, y0, width, height) of the part of its saved
    page to show (points), and 'x', its left edge on the new page. A panel with a 'y' (its top edge, points) goes there
    instead of under the panel before it, so panels can also sit side by side. Returns the page height in points.
    """
    from pathlib import Path
    from xml.etree import ElementTree as ET
    import fitz
    y, placed = 0.0, []
    for p in panels:
        top = p.get('y', y)
        placed.append((p, top))
        y = top + p['clip'][3] + gap
    height = max(top + p['clip'][3] for p, top in placed)

    out = fitz.open()
    page = out.new_page(width=page_width, height=height)
    for p, top in placed:
        cx, cy, cw, ch = p['clip']
        page.show_pdf_page(fitz.Rect(p['x'], top, p['x'] + cw, top + ch), fitz.open(p['pdf']), 0,
                           clip=fitz.Rect(cx, cy, cx + cw, cy + ch))
    out.save(out_pdf, garbage=3, deflate=True)

    register_svg_namespaces()
    root = ET.Element('{http://www.w3.org/2000/svg}svg', {
        'width': f'{page_width:.3f}pt', 'height': f'{height:.3f}pt', 'viewBox': f'0 0 {page_width:.3f} {height:.3f}',
        'version': '1.1'})
    for (p, top), letter in zip(placed, lowercase):
        text = _prefixed(Path(p['svg']).read_text(encoding='utf-8'), f'panel{letter.upper()}_')
        part = ET.fromstring(text.split('?>', 1)[-1])
        cx, cy, cw, ch = p['clip']
        part.attrib.update({'x': f"{p['x']:.3f}", 'y': f'{top:.3f}', 'width': f'{cw:.3f}', 'height': f'{ch:.3f}',
                            'viewBox': f'{cx:.3f} {cy:.3f} {cw:.3f} {ch:.3f}', 'id': f'panel_{letter}'})
        root.append(part)
    ET.ElementTree(root).write(out_svg, xml_declaration=True, encoding='utf-8')
    return height
